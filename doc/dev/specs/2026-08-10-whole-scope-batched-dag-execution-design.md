# Whole-scope batched DAG execution: a scope-tree scheduler and pure-realizer executor

_2026-08-10. Branch context: `evaleev/feature/multimode-batched-eval` (SeQuant),
paired MPQC batched-CCk work. Motivated by the water-20 PNO-CCSD 2.2x slowdown
investigation, which isolated the root architectural limit: recursive forest
descent cannot share a K-slice across trees._

## Motivation

The current batched runtime is **recursive forest descent**: it evaluates the
eval forest one tree at a time. Cross-tree common-subexpression reuse survives
only at the **top scope** (the one cache that persists across trees); a value
shared by two trees whose lifetime is BELOW the top scope is re-formed per tree,
because a K-slice of it cannot persist in a per-tree scope across trees.

The water-20 investigation quantified the consequence (job 658937 vs 631196,
identical energies, 2.2x wall):
- The DP batches K-carrying gC composites `Κ:con` across the board (36 vs the old
  ~12 distinct shapes), and each is rebuilt per consumer batch group.
- One concrete case: `g(μ̃,i₃,Κ)·C(i₃,i₂,μ̃;a₃<i₂,i₃>)` = `g{i₃,a₃<i₂,i₃>;K}` was
  built ONCE and shared (`store=10`, `access=30`) by the old code and rebuilt 15x
  (`store=150`) by the new code -- a genuine placement/sharing regression.

The analysis showed a real performance impact from RECOMPUTATION of shared
intermediates that forest descent forces (a K-slice cannot persist across trees)
and that a whole-scope DAG schedule would eliminate. That recompute is what this
project targets.

The fix is the execution model the DAG-scope machinery (home scopes,
`home_modes`/`home_depth`, the router, remat, the peak oracle) has always assumed
but no runtime has honored: **whole-scope descent** -- descend the fused forest
scope by scope, do all of every tree's work homed at each scope for each batch
block, so a shared sub-result is built once at its home scope and reused by its
whole subtree.

## Two execution models (from the water-20 analysis)

- **Recursive forest descent (current).** One tree at a time; cross-tree CSE only
  at the top scope. The static peak oracle (`peak_profile_sweep` over
  `home_modes` co-residency) does NOT model this runtime -- which is why the
  oracle and the replay never matched.
- **Whole-scope / DAG descent (this design).** One fused loop-nest walk for the
  whole forest; a value homed at a scope node is built once per block there and
  reused by its subtree across all trees. The co-residency peak oracle models
  THIS runtime. This design makes the runtime and the oracle agree.

## The scope tree

Loops/scopes form a **scope tree** (a.k.a. loop tree): the root is the top scope
(outside all loops); each child of a node is a batch loop entered at that scope;
every root->node path is one particular loop nest. A value is homed at some node
in the tree, built when that node's loop block is active, and reused by its whole
subtree.

**Placement = the LCA of a value's consumers in the scope tree.** A value consumed
only within one branch homes deep in that branch; a value whose consumers live in
different branches homes at their lowest common ancestor, possibly the top scope.

**Narrow (current) vs general (future) trees.** The existing placement machinery
emits a NARROW scope tree: one loop per index type in a single canonical order --
effectively one linear chain the whole forest shares. For a single global order,
the current set-based `home_modes` (intersection of enclosing batch mode-SETS) IS
the correct LCA. The GENERAL tree has per-branch loop orders, where set-meet
breaks:

> Consumer A in path `[occ, K]`, consumer B in path `[K, occ]`: tree-LCA = root
> (they diverge at the first loop), but set-meet = `{occ,K} ∩ {K,occ}` = `{occ,K}`
> -- non-empty, WRONG (it would home the value in a nest neither consumer shares).

The general branching tree, with tree-LCA placement replacing set-meet, is a
later placement upgrade. The executor is built GENERAL so this is a placement
change, not an executor rewrite.

## Design

### Pipeline -- three stages, two of them reuse existing machinery

1. **Fuse + place (reuse).** The forest fuses into one DAG by hash CSE, which
   `compute_dag_boulevard` already does (occurrences grouped into `ValueCell`s),
   and each value gets its home scope from the existing DAG-scope placement. Home
   is a mode-set (narrow-chain level) today; a tree-path later.
2. **Schedule (new).** Build the scope tree (union of home paths -- one chain
   today) and assign each value to its home node. The cross-scope VALUE ORDERING
   need not be pre-computed at first: the executor can materialize an above-homed
   value lazily on first reference (see Open Questions). So the initial scheduler
   deliverable is the scope tree + per-value home; an explicit ahead-of-time topo
   order is a later optimization.
3. **Execute (new driver, reused primitives).** Walk the scope tree.

### The load-bearing insight: replace the DRIVER, not the primitives

The batched-scratch mechanics already exist in `eval.hpp` and the current
opportunistic co-evaluation uses them: per-block scratch, ACCUMULATE (contracted
mode), SCATTER (external mode), and SLICE-ON-USE (a consumer slices an
ancestor-homed value to the current block). Whole-scope descent swaps the DRIVER
-- from per-tree recursion + the co-evaluation trigger, to a scope-tree walk over
the whole fused forest. The per-op eval (`prod`/`sum`) and the memory primitives
are untouched. This contains the surface area and the correctness risk.

### The executor walk (scope-tree DFS)

```
walk(scope_node, block_context):
  for each block b of scope_node.mode:
    bind context (scope_node.mode = b)
    build every value HOMED at scope_node for block b
        -- operands sliced-on-use from ancestors + already-built siblings,
           shared across ALL trees (the CSE forest descent could not do)
    for child in scope_node.children:
        walk(child, block_context + b)
    on block exit:
        contracted-mode results ACCUMULATE into the parent scope's result;
        external-mode results SCATTER into their pre-sized result slice
  free scope_node's scratch
```

A value homed at a node is built once per block and reused by its whole subtree.
The root (top scope) is walked once (a single "block"); values homed at the root
are the cross-tree shared, K-independent intermediates built once for the run.

### The executor is a pure realizer

It honors the schedule's placement exactly -- no demote-on-pressure, no ad hoc
spill, no runtime placement decisions. All peak/recompute policy lives upstream in
placement (remat / the budget / the cost model). The executor's only job is to
walk the given scope tree and run whatever placement it is handed.

### The paradox resolved, and the oracle bound

The water-20/C60 work found placement was INERT on peak under forest descent --
because forest descent never co-resides cross-tree sets, so where a value is homed
does not change the realized peak. Under whole-scope descent, the co-resident set
at a scope block is exactly "everything homed at or below that block, across all
trees, for the current block" -- so PLACEMENT IS THE DOMINANT PEAK LEVER, precisely
as the DAG-scope model always assumed. The co-residency peak oracle (the Task-1b
`home_modes` model, which models whole-scope descent) becomes the accurate
predictor of the realized peak. This design is what makes the runtime and this
session's cost-model / remat / router machinery finally describe the same thing.

## The peak/recompute tradeoff is RELOCATED, not eliminated (honest risk)

Whole-scope descent RAISES co-residency -- it holds all trees' work at a scope
block. For C60 that peak is exactly why batching exists. So this design does NOT
give free CSE at no peak cost. It RELOCATES the peak/recompute tradeoff into
principled placement: home a shared value high (great CSE, high co-residency) vs
deep (low co-residency, rebuilt per block = recompute). If the budget forces deep
homes, recompute returns -- but as a PRICED placement choice, not a forest-descent
artifact. The wins are (1) cross-tree sharing is finally POSSIBLE, and (2) the cost
model predicts what actually runs. The win is not "no peak cost."

Crucially, whether a given placement is feasible -- e.g. whether C60's top scope
holds too much data to fit -- is a PLACEMENT / OPTIMIZER question, NOT an executor
question. The executor evaluates a fixed schedule; it does not produce or change
it. If the top scope is over budget, the placement / budget / optimizer homes
values deeper (accepting priced recompute) until it fits; the executor faithfully
runs whatever it is handed. C60 feasibility therefore belongs to the future
scope-tree-aware optimizer (out of scope, below), not to this project.

## Scope boundary

**In scope:**
- The scope-major SCHEDULER: fuse + read placement + build scope tree + emit a
  walkable schedule (cross-scope value ordering lazy at first; explicit topo
  later).
- The general scope-tree EXECUTOR (pure realizer): the scope-tree DFS walk reusing
  the existing accumulate/scatter/slice-on-use/scratch primitives.
- Validation on the NARROW canonical-chain placement the current machinery emits.

**Out of scope (explicit):**
- General branching-tree PLACEMENT -- per-branch loop orders, tree-LCA homes
  replacing the set-meet. The executor is built general to accept it later.
- A more complex OPTIMIZER that shapes the scope-tree schedule for a batched DAG
  and prices INTRALOOP CSE correctly. The current DP prices per-tree and cannot
  see cross-tree sharing; a scope-tree-aware optimizer is a separate, larger
  project. This design realizes the EXISTING placement; it does not optimize it.
- Retiring the current forest-descent path -- the two coexist until whole-scope is
  validated.

## Validation

- **Correctness.** The whole-scope executor must produce results that agree with
  forest descent to within SMALL NUMERICAL NOISE -- NOT byte-identical: the
  schedule changes the evaluation/contraction order, so floating-point
  non-associativity yields a tiny diff. Validate on the water-20 and C60 CSV-CCSD
  doubles residuals (converged energies agree to a tight tolerance, not to the
  bit).
- **Cost model.** The dry-run cost model must predict the whole-scope realized
  peak/recompute. Under whole-scope descent the realized peak IS the placement
  co-residency, so the `home_modes` co-residency oracle becomes the faithful
  predictor -- the direct test that the runtime now honors the DAG-scope model.
- **The regression it targets.** On water-20, the shared K-carrying composites
  (e.g. `g{i₃,a₃<i₂,i₃>;K}`) should be built once at their home scope and reused,
  not rebuilt per group -- recovering the old code's sharing without its
  per-tree-descent accident.

## Open questions (resolve during implementation)

- **Cross-scope ordering: lazy first, explicit topo later.** A value homed at
  scope S is needed by consumers in deeper scopes. The FIRST implementation need
  not build an ahead-of-time topological order across scopes: it can reuse the
  existing LAZY mechanism -- when the walk hits a node inside a batch loop whose
  home is ABOVE, compute it in the form its home needs (e.g. full extent), place
  it in the home's cache, and access it by slicing thereafter. On-demand
  home-materialization + slice-on-access handles the ordering implicitly, exactly
  as the current evaluator already does for a mid-loop reference to an
  above-homed value. An explicit ahead-of-time topo sort (DFS of the scope tree +
  per-node value ordering by data dependence) is a later optimization, not a
  prerequisite.
- **Fused-DAG construction.** `compute_dag_boulevard` already groups occurrences
  by hash into `ValueCell`s; confirm this is a sufficient fused-DAG representation
  for the scheduler, or what it must add (explicit inter-value edges, per-value
  home node).
- **Accumulate/scatter across the fused nest.** The primitives exist per-tree;
  confirm they compose correctly when a single scope block's result aggregates
  contributions from MULTIPLE trees (the accumulate target and scatter slice must
  be shared/pre-sized across trees).
- **Coexistence with forest descent.** A flag selecting whole-scope vs forest
  descent, and how the dry-run cost model selects the matching peak model.

## Reproduce / evidence (from the water-20 investigation)

- Forest-descent limit + fragmentation: `[.][dryrun-water-frag]` (36 Κcon members,
  gC priced `rf==1`), `[.][dryrun-water-remat-sweep]` (remat demote-only, seed
  homes 17/21 gC in-loop, no peak_threshold hoists). Both committed on this branch.
- The sharing regression: cache `store`/`access` for `g{i₃,a₃<i₂,i₃>;K}` = 10/30
  (old, hoisted) vs 150/30 (new, rebuilt 15x).
- Two-execution-models framing + the oracle-vs-replay gap:
  `doc/dev/specs/2026-08-09-remat-into-cost-profile-design.md`, Task 1b.
