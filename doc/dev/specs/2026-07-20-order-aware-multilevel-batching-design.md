# Order-aware multilevel batching: cost model + lifetime management (design spec)

**Status:** designed 2026-07-20 (brainstorming). Supersedes the runtime-only
framing of `mpqc4/doc/dev/specs/2026-07-17-nested-batch-group-join-design.md`,
which it folds in and generalizes. This is the foundational cost+runtime piece:
without it, the factorizer cannot produce a correct batching structure at all.

> **RETRACTED / CORRECTED (2026-07-27) -- SECOND, separate confound (role-predicate).**
> The 2026-07-22 correction in Section 1 fixed the heuristic-fallback confound. A
> DISTINCT, additional confound remains: the D5 table's "external-occ" arm ran with
> occupied indices in the CONTRACTED-role batchability predicate (only
> `is_batchable_index` set; `is_batchable_external_index` never supplied), so the DP
> ALSO batched contracted occ (cost cells sliced the contracted occ pair
> `i_3 i_4`). Every "external-occ" number here was "external-occ ON TOP OF
> contracted-occ"; the IDENTICAL peak across the two rows is the tell. Post
> role-split (SeQuant `79540c831`..`7f5240014`, MPQC `1685c7ac1e`..`5a3d94dc63`)
> honest re-measure of `[.][dryrun-occ-veto]` (nterms=55): peak ~= **5860.9 GB**
> (NOT 2302), avoidable_time ~= **44.6%**, `contracted_occ_stamps == 0`,
> `external_occ_stamps == 244`. The ~2302 GB floor cited in the success criterion
> (Section 7) and out-of-scope (Section 10) is therefore mis-attributed. As the
> 2026-07-22 note already states, the spec's CORE argument (order-blindness is a
> representability defect) is independent of these magnitudes and is NOT retracted;
> only the C60 numbers and their attribution are. Root cause:
> `.superpowers/sdd/contamination-role-predicate.md`.

**Repos:** SeQuant (`core/optimize/cost_model.hpp`, `core/eval/eval.hpp`,
`core/eval/cache_manager.hpp`, `core/eval/eval_expr.{hpp,cpp}`). Consumed by MPQC
CCk (CSV/PNO-CCSD); any behavior change needs an MPQC repin.

**Instrument / reproducer:** `tests/unit/test_eval_dryrun.cpp -> [.][dryrun-occ-veto]`
(SeQuant `5118c7738`). Replays the real C60 residual forest in seconds and reports
`avoidable_time` per policy. The perf-first `peak_threshold` ceiling made the
dry-run faithfully mirror the runtime, so this is a trustworthy gate.

---

## 1. Problem and motivation

The perf-first batched objective bounds the dominant particle-particle-ladder giant
`W` by **external-mode batching** (see `2026-07-20-external-mode-batching-design.md`,
D1-D3): it batches `W`'s external occ `i_1,i_2` work-neutrally. The D5 gate proved
that is **necessary but not sufficient** on the C60 residual forest:

| config | replay ops | avoidable_time | scatter | modeled peak |
|---|---|---|---|---|
| contracted-occ | 92848 | 75.1% | 0 | 2302.30 GB |
| external-occ | 66039 | 44.7% | 46 | 2302.30 GB |

> **[CORRECTION 2026-07-22 — this table is CONFOUNDED. See
> `.superpowers/sdd/oamb-a0-note.md` sections 7, 12, 14, 16.]** Its two arms
> flip `batch_spectator_indices` **and** `suppress_heuristic_fallback`
> together, so the contracted arm ran the legacy runtime heuristic and the
> external arm did not. Holding the flag fixed:
>
> | config (single variable) | replay ops | avoidable_time |
> |---|---|---|
> | contracted-occ, fallback ON | 120424 | 74.53% |
> | **contracted-occ, fallback OFF (the missing cell)** | **61275** | **43.72%** |
> | external-occ, fallback OFF | 83573 | 43.96% |
>
> So 44.7% is not "what external batching left behind" — it is what you get
> with the fallback off, external or not. **External-mode batching contributes
> ~0 here**, and is measured neutral-to-slightly-worse. The heuristic fallback
> has since been removed outright.
>
> Three further corrections bear on this section's framing:
> - **External batching is not under-priced, and needs no cost-model work.**
>   `seeded_forest_peak` declines any seed that is not work-neutral, so a seed
>   whose invariant subtrees would be recomputed per block is never adopted and
>   no `External` stamp is emitted. Verified: on a minimal network where seeding
>   would cost 1200 vs 720, External stamps = **0**.
> - **The 76% is retracted** (witness misconfiguration; MPQC suppresses that
>   fallback), and job 631467 never exhibited it — its log has zero heuristic
>   firings. Its actual pathology was the broad caching veto, fixed by
>   `2a52e063c` two days after the run.
> - **The ~2302 GB peak is not addressable by any caching or cost change.** It
>   is one transient 2231 GB operand formed inside an already-sliced `Κ_2` loop,
>   full in the free PAO index and the proto-occ pair, never cached.
>
> **None of this weakens the spec's core argument.** Order-blindness is a
> *representability* defect — a loop **set** cannot encode a loop **tree**, so
> the free-hoist vs charged-hoist distinction is inexpressible on any witness.
> That is independent of C60 magnitudes and is unaffected by every correction
> above. What is retracted is the prioritisation rhetoric attached to it.
>
> **[ADDENDUM 2026-07-27 -- a SEPARATE, additional confound.]** Beyond the
> heuristic-fallback confound corrected above, the table also had occ in the
> CONTRACTED-role predicate (role-predicate confound), so both arms batched
> contracted occ and the ~2302 GB peak is mis-attributed. Honest post-role-split
> re-measure: peak ~= 5860.9 GB, avoidable ~= 44.6%. See the top-of-doc note.

External batching engages and strictly improves both metrics but neither eliminates
the recompute (44.7% remains) nor bounds the forest peak. The surviving offender is a
**contracted middle-gap** node: an intermediate that sits inside a batch loop over a
mode it does NOT carry, so it is rebuilt per batch of that mode. On the C60
occ+aux run, **76% of modeled execution time is avoidable recomputation** of this
class -- the worst offender is `g(μ̃,μ̃,Κ_2) * C -> I(i_2,i_1,μ̃,Κ_2;a_1<i_1,i_2>)`
(the half-transformed DF factor `gC`, carrying the aux `Κ_2`, rebuilt across an occ
loop it does not carry).

### Metric: avoidable-recompute TIME

The only meaningful waste measure is the share of modeled execution TIME spent
recomputing a value already built. Op counts do not work: batching a mode into `n`
batches legitimately turns one op into `n` work-neutral slice-ops, and a cheap node
rebuilt 1000x and an expensive one rebuilt twice look opposite by count. Replay the
forest through the real eval loop; key each op by `(expression, {slice of each
batched mode that OCCURS in the op})`; the first build of a key is necessary, any
later build with the same key is avoidable.

    avoidable_time = Sigma exec_cost over duplicate-key builds / Sigma exec_cost over all builds

`avoidable_time == 0` is a perfect-sharing evaluator under the same schedule.

---

## 2. Root cause: the DP is order-blind

`dense_time_space` is a single joint subset DP over the contraction tree
(`PeakBatchedModel::relax`). It charges a batch-recompute term into the flops
(`cost_model.hpp:910`, gated by `charge_batch_recompute`, default on):

```cpp
std::size_t const esc = B & ~ctx.open_modes[n];   // batched modes n does NOT carry
for k: if (esc >> k & 1) rf *= ctx.nbatches[k];    // charged nbatches times each
```

Today `B` is a **set** (a bitmask of batched modes) and `esc = B & ~carried` is a
set; the DP enumerates every `B`. For `gC` -- which carries the aux `Κ_2` but not the
contracted occ `i_3` -- the min-cost choice is `B = {Κ_2}`: dropping `i_3` gives
`esc = {}` and **`rf = 1`**, treating `gC` as a CSE built once and shared. The runtime
cannot honor that: `gC`'s unbatched footprint is 4.7 TB (it carries `Κ_2`, which must
stay batched), so it is rebuilt once per `i_3` batch -- **DP `rf = 1` vs runtime
`rf ~ nBatch(i_3) = 15`**.

Whether `gC` can be hoisted -- and therefore whether `rf = 1` is real -- depends on
the **nesting order**, which a set cannot represent:

> **Place an intermediate in the cache of the outermost batching context that
> contains all of its carried batched modes. Its recompute factor is the product of
> `nBatch` over the batched modes nested OUTER to that placement that it does not
> carry.**

- `[Κ_2 outer, i_3 inner]`: `gC` placed at the top `Κ_2` level, above the `i_3`
  loop, reused across `i_3` -> **`rf = 1`** (legitimate: rebuilt only per `Κ_2` batch,
  which it carries).
- `[i_3 outer, Κ_2 inner]` (what we see): `gC` placed at the inner `Κ_2` level,
  rebuilt per `i_3` batch -> **`rf = nBatch(i_3)`**.

`esc = B & ~carried` is the set-collapsed version of this rule -- correct only when
every carried mode is batched at or below the node. Otherwise it both over-charges the
hoistable order and lets the DP drop an outer mode from `B` to escape the charge (the
`rf = 1` phantom). The fix (Section 3) charges each node from its placement in the
subtree-bound loop tree instead of the collapsed set.

The fix has **two coupled parts**, joined by one annotation contract:
**(A) an order-aware cost model** and **(B) multilevel lifetime management**. This
spec designs both; the annotation contract (the DP emits, the runtime obeys) is the
interface between them.

---

## 3. The order-aware cost model (DP)

The batching structure is a **loop tree** mirroring the contraction tree -- **not** a
single global nest over the term. Each batchable index's loop wraps the **subtree**
where that index lives; a node's recompute factor is set by the loops that enclose
*it*, not by a term-wide order. (An earlier draft assumed one linear nest per term; the
`R{i,m}` analysis below shows why that is wrong.)

### Subtree-binding: where each batch loop goes

A batchable index `x`'s loop must wrap the **minimal subtree containing every
occurrence of `x`** -- wrapping any larger subtree only forces needless recompute of
the extra, non-`x` branches. Two regimes:

- **Contracted index** (summed inside the term: `j`, `l`, `k` below): its minimal
  subtree is exactly the subtree rooted at its contraction node, so its loop placement
  is **determined by the contraction tree**. There is no order to search and no
  interchange legality to check -- both are automatic from "wrap the contraction
  subtree." Choosing the tree chooses the contracted-index loop structure, so this
  folds into the existing tree DP; the DP's job is only to **charge each tree
  correctly**.
- **External index** (in the term's result: `i`, `m` below): it occurs at leaves and
  at the root, so its minimal covering subtree can still contain branches it does not
  carry. **Which subtree its loop wraps is a genuine optimization** -- the open, harder
  case, the same problem as the forest-peak follow-on (Gap 1) and the external-mode
  scatter (D1-D3). This spec models its cost; choosing it optimally is deferred there.

### Multiple contracted indices at one node

A contraction node sums over one *or more* indices (every index shared by its two
children and closed at the parent). Batching several of them puts several loops at that
one subtree -- a **sub-nest**. Its **enclosing** order is fixed by the tree (an index
contracted *above* is outer); the order **among the indices co-contracted at that node**
is free. That free order is **cost-neutral at the roofline level** -- both operands
carry every co-contracted index, so an `(K, L, ...)` slice block has the same
operand/result footprints and the same flop count whichever loop is outer; peak and
flops are symmetric in it. So the DP does **not** search it. (The only shape that could
break neutrality -- a one-sided/unary reduction -- is not representable in the current
tensor IR: unary reductions are ARRAY operations, always performed eagerly, and never
fall out of tensor products. So neutrality is unconditional here.) But it **must be pinned by
a canonical, deterministic rule** (e.g. ascending batched-slot ordinal), never left to
hash/iteration accident: an unspecified order is both a scheduling nondeterminism and a
**bit-reproducibility** hazard (loop order changes the batch accumulation order, hence
rounding). The annotation (Section 4) therefore emits the co-contracted loop order, not
just the set of batched modes. (Finer locality/communication effects of that order are
below this roofline model's resolution.)

### The charge: per-node hoist-vs-recompute, against the budget

For a batch loop over `x` enclosing a node `n`:

- if `n` carries `x`, `n` is sliced per batch (`x1`, no recompute);
- if `n` does **not** carry `x` (`n` is invariant to the loop), the DP chooses per node:
  - **hoist** `n` above the loop -- built once, but **resident** across it
    (peak += `n`'s footprint at its placement); or
  - **recompute** `n` per batch -- `flops(n) x= nBatch(x)`, no residency.

The choice minimizes time subject to `peak <= budget` (Section 7). A node's recompute
factor is `Prod nBatch(x)` over the enclosing loops `x` it does not carry **and is not
hoisted above** -- a tree-structural quantity read off the loop tree, replacing the
order-blind `esc = B & ~open_modes[n]` (which, being a set, cannot say which loops
actually enclose the node, so it both misplaces and mis-charges).

**Materialization cost is gated by volatility.** A (re)built leaf/subtree costs only if
it is actually recomputed: a *volatile* (lazily recomputed) leaf is charged
`volatile_weight x` per replay, a fully-materialized one is free; batched, that
multiplier becomes the recompute factor above. The model already carries this hook
(`cost_model.hpp` `volatile_mask` / `volatile_weight`) -- it needs the loop-tree
structure to supply the multiplier.

### Peak model: a resident scan across every root-leaf path (DIRECTIVE 2026-07-23)

The charge above prices the *flops* side of hoist-vs-recompute; the **peak** side needs
a model the DP does not have today. `PeakBatchedModel`'s peak is a purely *local* max at
each node (`szlp + szrp + szn + contrib`, crossed with child peaks) -- it never accounts
for what is still **alive above** the node. That is why A3a's `Carr == 0` exemption is an
under-charge: it sets `rf = 1` and adds **no** peak term, so it trades a flops
over-charge for a peak under-charge. Hoisting is never free; the model simply cannot see
its cost.

The correct single-term peak is a **max-peak search over every root-leaf path**. At each
node along a path:

    peak(node) = local co-residency here
               + Sum over enclosing loop levels L of (footprint of the intermediate
                 resident at L)

An intermediate hoisted/accumulating at an enclosing batch-loop level is resident across
that loop's whole body, so it co-exists with everything evaluated inside. **Each
enclosing batch loop contributes one resident charge**; a node that contracts more than
one mode has more than one loop level there, hence more than one charge. The term peak is
the max of `peak(node)` over all nodes; hoisting a node above a loop adds its footprint
to the resident set for that loop's subtree, which the path-scan now sees -- so the
`(peak, flops)` hoist-vs-recompute trade is finally *representable*, for BOTH the
inner-escaped (small, still-sliced) and outer-escaped (un-sliced, the 8.7 GB -> 2231 GB
case) cases. This makes peak context-dependent (it depends on the enclosing loop tree,
i.e. the ordered `B` + placements) rather than a pure bottom-up local max, and
`reconstructed_batched_peak`'s `sim` oracle must mirror the same scan.

Scope: this models **within-term** peak correctly. Cross-term residency (a persistent
intermediate alive across terms) is under-charged **by definition** -- the peak model is
single-term. That limitation is accepted here (Section 9), not solved.

### External loop placement is in scope (DIRECTIVE 2026-07-23)

Earlier drafts treated external-mode placement as a follow-on and left external modes out
of `B` entirely (they are contracted nowhere, so `C = B | aprime` can never admit them),
which silently assumes an external loop wraps the **whole term** -- emit stamps `External`
on every carrying node and the runtime scatters over the entire tree. That is wrong in
general: in `(AB)(CD)` an external mode's loop can legitimately wrap only a subtree, and
its position relative to the contracted loops is exactly what the loop tree expresses. So
an external mode **must be able to enter `B` and be ordered/placed** like any other. Two
consequences: (a) the ordered-key cap is set on all `B`-entering modes including
externals (C60 `m_B` rises from 3 to 5), and (b) `seeded_forest_peak`'s work-neutrality
guard, which today declines a seed that is non-work-neutral *at root level*, must be
evaluated at the seed's actual placement subtree -- a seed non-work-neutral over the whole
term may be work-neutral over a subtree it carries everywhere, so proper placement makes
currently-declined seeds adoptable (a capability gain, not a tidy-up).

### No interchange-legality question

Because contracted-index loops are subtree-bound (previous section), there is no
separate "which orders are legal" predicate for them: a nest that would read a partial
accumulator is simply not expressible -- the accumulated mode's loop *is* its
(inner) contraction subtree, so it cannot float outside. Order and legality both fall
out of the tree. (For a fixed tree only some batch structures are realizable; the tree
DP already ranges over the trees, so this is covered by searching trees, not by a
bolt-on order search.)

### Worked examples

**`R{i,m} = (A{i,j} B{j,k}) (C{k,l} D{l,m})`**, tree `I1{i,k} = Σ_j A B`,
`I2{k,m} = Σ_l C D`, `R = Σ_k I1 I2`:

- **batch `k`** (carried by `I1` and `I2`, contracted at `R`): the loop wraps the whole
  tree; each nest does its per-`k`-slice work; **flop-neutral** -- nothing invariant is
  enclosed.
- **batch `k,j,l`**: `j`'s loop wraps only `I1`'s subtree, `l`'s only `I2`'s -- each at
  its contraction subtree, both inside the `k` loop. Flop-neutral *if leaves are
  materialized*; if `A` is volatile, its slice is rebuilt `#k` times (it sits inside the
  `k` loop but does not carry `k`) -- the volatility penalty above.
- **batch `i,k` (the middle gap)**: `i` occurs in `I1` and `R` but **not** `I2`.
  Naively wrapping the `i` loop around the whole tree recomputes `I2` `#i` times. The
  model's fix is to **hoist `I2` out of the `i` loop** (place it at the `k` level, above
  `i`): built once per `k`, read `nBatch(i)` times from inside the `i` loop. That
  placement + count is exactly the lifetime-scope / effective-use-count of Sections
  4-5.

**C60 `gC`** (`g(μ̃,μ̃,Κ_2)·C -> I(...,Κ_2)`, batched aux `Κ_2` and contracted occ
`i_3`): `gC` carries `Κ_2`, not `i_3`. Unlike `I2`, hoisting `gC` above the `i_3` loop
needs it whole over `Κ_2` (4.7 TB > budget), so it is **forced to recompute** --
`rf = nBatch(i_3) ~ 15`. Today's order-blind `esc` drops `i_3` and prices `rf = 1` (a
hoist the peak cannot afford) -- the phantom this model removes.

---

## 4. The annotation contract (DP emits, runtime obeys)

**Decision 0: the DP annotates, the runtime obeys verbatim -- never re-derives.**
The DP authored `(B, order, tree, nBatch-per-loop)`, so both annotations below are
static functions of the schedule; that is the only way realized scope cannot drift
from priced scope. Per node the DP emits:

1. **lifetime-scope (storage level)** = the node's placement in the loop tree: the
   outermost enclosing loop level that contains all its carried batched modes. Selects
   which cache holds the node and whose per-batch reset re-arms it. The loops enclosing
   any one node form a chain, so its scopes nest monotonically
   **run ⊃ term ⊃ L0 ⊃ L1 ⊃ ...**; the
   existing **P (run)** and **NP (term)** are the two outermost, batching adds inner
   ones. The context is keyed by **nest-position -> canonical slot**, not index label
   (only plain indices are batched; a proto-composite `a_1<i_1,i_2>` never is), so
   equivalent nests **coincide across terms** -- preserving cross-term sharing.

2. **effective use count** = `Sigma over consumers c of Prod nBatch(L) over loops L
   enclosing c but not n`. The dynamic read count. Drives **tight release**
   (decrement on every read -- fall-through reads DO decay -- release at zero), with
   **scope-reset as the always-correct backstop** (if runtime reads are fewer than the
   dense count, e.g. screening, the node clears at scope close: a peak penalty, never
   a correctness bug; reads can never exceed the dense count).

**P is the sole no-decay case.** The one quantity the DP cannot supply is
#CC-iterations-to-convergence (data-dependent), which is exactly P's cross-iteration
reuse, so P alone releases only at teardown. `persistent_` stays **one bit**.

**Emit mechanism.** The per-node stamp already carries `batched_here()` via
`reconstruct_batched_modes` / `set_batched_here`. The contract **extends that same
stamp** with the scope level, the effective count, and -- where a node contracts
several batched indices -- the **pinned co-contracted loop order** (Section 3), so the
runtime never picks it by accident. No new channel. One invariant:
an accumulator's internal read-modify-write is **not** a counted consumer read (the
evaluator owns accumulation and must not route it through the counted `access`), else
`I` releases before its post-loop consumer reads it.

---

## 5. The runtime: visibility + lifetime

Two orthogonal mechanisms; keeping them separate dissolves the down-propagation
blocker that killed the earlier attempt.

### Visibility -- the scope chain (walk up on read)

Each scratch/cache gets a `parent_` pointer following the loop tree (a node's enclosing
loops form a chain). `access(key)`
checks local, then walks **up** (inner scratch -> outer scratch -> real cache),
returning the first hit. A hoisted node lives in **exactly one** cache (its
lifetime-scope's); every deeper body finds it by walking up. **This is the committed
inert Task-1 primitive** (`cache_manager.hpp` `parent_` / `set_parent` /
fall-through `access`, `75d240079`) finally given its consumer.

- **Reads walk up; writes aim.** `store` targets the one cache at the node's
  DP-emitted scope level (walk up to that depth -- one target), never pushed down
  into descendants (the many-target blocker).
- **Eviction stays per-level and blunt.** A loop's per-batch reset clears only its
  own scratch; ancestors untouched. Each level holds only its own-scope nodes, so a
  fall-through read can never surface a stale slice of a batch-variant node (those
  live-and-die in the loop cache). The "slow-not-wrong" failure mode is preserved.

### Lifetime -- effective count + reset

Two release triggers: **count-zero** (tight early release, decrement on each read) and
**scope reset** (always-correct backstop; over-holding is a peak penalty, never
wrong). **P is the sole no-decay case.**

### Where the evaluator stores invariants

At a batching node, instead of building each loop-invariant descendant on the per-batch
scratch that resets under it, the evaluator `store`s it at its DP-emitted context
level with the DP's effective count -- loop-invariant code motion on the batched
replay.

### Falls out for free / must be respected

- **Accumulators need no new machinery.** The standard rule puts an accumulator's
  scope at the loop's *parent* (consumed after the loop), i.e. above the reset line,
  so per-batch resets never disturb accumulation. Nested (`I += ...` inside an `i_5`
  loop) -> `I`'s scope is the `i_5` cache: survives inner resets, cleared by its own.
- **The free-batchable veto is reinterpreted, not removed.** It encodes
  "batched => not run-scope," which the lifetime-scope rule states directly; it must
  not evict a hoisted (batch-invariant, non-run-scope) node from its context level.
- **Two gotchas from Attempt 1.** Build a whole invariant node on a **fresh** cache
  via `evaluate(n, le)` (variadic -> fresh empty cache; `evaluate(n, le, cache)`
  re-enters the batched evaluator and SIGSEGVs). Invariance is a **subtree** property
  (does the subtree touch the batched mode -- Contracted `batched_here` at this level
  AND below -- not the node's result indices).

---

## 6. Constraints any fix must respect

- **Exactness.** Sharing must be slice-exact: identical canonical batch-mode position
  AND element range (and, for an External mode, identical external slice -- a node
  carrying an outer External mode is not batch-invariant under it and must not be
  seeded whole). Batching relies on `sum_K = sum_{batches} sum_{K in batch}` per
  mode; nesting composes these.
- **Peak.** The point of batching is to never materialize a batched mode's full
  extent at once. Restoring sharing by caching a shared node *whole* (across its
  batched mode) defeats the memory saving. Sharing must be of *sliced* values, per
  outer batch -- caught by the run's reported peak.
- **No regression on working paths.** Unbatched (~0%) and aux-only (~0.48%) both work
  today and must not regress. The `peak_threshold = +inf` / no-annotation legacy
  heuristic path must stay **byte-identical**.
- **Cross-repo.** Validate on the local `[.][dryrun-occ-veto]` repro before any Owl
  run; behavior changes need an MPQC repin.

---

## 7. Success criterion

The committed `[.][dryrun-occ-veto]` witness asserts

    both.peak_gb          < 100.0    // PRIMARY gate; FAILS today at ~2302 GB
    aux.avoidable_time()  < 0.10     // aux-only is clean
    both.avoidable_time() < 0.10     // FAILS today at 0.437

with **peak still sliced, not materialized whole**.

> **[CORRECTION 2026-07-22.]** The earlier criterion read `< 0.05` against a 76%
> baseline and claimed "roughly a 4x per-iteration speedup on C60 -- the primary
> lever". Both are withdrawn: the 76% was a runtime misconfiguration (43.72% is
> the production-representative figure) and the 4x was inferred from it plus the
> assumption that C60's 4-hour wall was this recompute, which it was not.
>
> Peak is now the PRIMARY gate: this is a memory problem, and a time metric
> alone cannot distinguish a real fix from a schedule that shares less and holds
> more. Note that peak is a *slicing* property, not a caching one, so Phase B
> cannot be judged by it -- see the note's section 14.2.
>
> The replay aggregates are **diagnostics, not targets** (the witness header
> states the trust boundary and its five open defects). Phase A should be gated
> on DP-side structural facts, which do not pass through the replay: the
> committed RED gate `[.][loop-tree]` asserts on the DP cell `st[n][B]` directly
> and fails 4000 vs 1000 on the free-hoist shape.
>
> **[ADDENDUM 2026-07-27.]** The "~2302 GB" figure in this criterion is the
> contaminated "external-occ ON TOP OF contracted-occ" peak (role-predicate
> confound). Honest post-role-split value is ~5860.9 GB. The peak-gate criterion
> stands as a criterion; only its cited magnitude is mis-attributed. See the
> top-of-doc note.

---

## 8. Implementation phasing (cost -> runtime, one spec)

- **Phase A -- order-aware cost model.** Charge each contraction tree by its
  subtree-bound loop-tree recompute (per-node hoist-vs-recompute, placement-aware `rf`
  read off the loop tree) in place of the order-blind `esc`; emit
  `{lifetime-scope, effective-count}` on the per-node stamp. *Checkpoint
  (dry-run-measurable, before any runtime):* the DP forces `gC` to recompute
  (`rf = nBatch(i_3)`, no `rf = 1` phantom), hoists genuinely-hoistable invariants
  (`I2`-class), and its modeled avoidable reflects the true cost.
- **Phase B -- runtime.** Wire the inert Task-1 scope chain as consumer + store-at-
  context + effective-count release + veto reinterpretation. *Checkpoint:*
  `both.avoidable_time() < 0.05`, peak sliced.
- **Phase C -- cross-repo.** MPQC repin + local `[dryrun-occ-veto]`; Owl C60
  (separate session): per-iteration time down, energy unchanged.

---

## 9. Known limitations (carry forward)

- **External-index loop placement is modeled, not optimized here.** For a *contracted*
  index the loop is subtree-bound (tree-determined); for an *external* index the
  wrapping subtree is a free choice this spec prices but does not search -- optimal
  external-index placement is the forest-peak follow-on (Gap 1) / external-mode scatter.
  Contracted-index batching (the `gC` middle gap) is fully in scope.
- **Per-node hoist-vs-recompute is decided locally.** Each (node, enclosing loop) picks
  hoist or recompute against the budget; a globally-joint optimum could differ. The
  gate's residual localizes any case where the local choice is wrong.
- **Cross-term coupling.** Loop trees are per term; an intermediate shared across terms
  (persistent CSE) needs compatible placement in both -- not solved here.
- **Grouping generalization** -- co-batching a *carried*-mode shared descendant across
  its consumers (built once per batch by co-evaluation) is the complementary regime;
  the current single-partition/first-axis join is not sequence-aware. Deferred.
- **Peak-model verification** that resident accounting equals held-at-context
  footprint -- a verification task.

## 10. Out of scope / adjacent

- **Forest-level / multi-mode co-batching** to drive the forest **peak** under budget
  (Gap 1 from the D5 gate: the 2302 GB peak set by a term that batching the proto-occ
  pair does not reach). This is the explicit **follow-on** that stands on this
  foundation -- correct multilevel placement pricing is its prerequisite.
  > **[RETRACTED 2026-07-27 -- see the top-of-doc note.]** The "2302 GB" here is the
  > contaminated "external-occ ON TOP OF contracted-occ" figure; honest
  > post-role-split peak is ~5860.9 GB, and the driver is the cross-pair
  > two-PNO-leg giants external slicing cannot reduce. Re-scope this follow-on
  > against a role-separated config.
- The External / occ-forest emission path (`2026-07-20-external-mode-batching-design.md`)
  is done; it is what bounds `W`, and it is the schedule this spec's contracted-occ
  middle-gap sits inside.

## 11. To consult

- **CoNST (Raje et al.)** -- believed to be recent work on generating optimized
  loop-nest / communication schedules for tensor-contraction trees; plausibly directly
  relevant to the external-index subtree-placement calculus (Section 3) and the
  materialize-vs-recompute trade. Pull the actual paper and mine it before relying on
  it; this pointer is from memory and unverified.
