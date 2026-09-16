# Batched array DAG evaluation -- as-built design (phase 1)

## 1. Purpose and scope

This document describes, as built, the machinery SeQuant uses to evaluate a whole
evaluation FOREST as one fused DAG under batch (memory-blocking) loops: the batched
dynamic program that decides where loops are opened, the loop identity and occurrence
facts that name those loops, the ordered schedule that realizes them as a block tree,
the explicit cell table that states every value form and every read, the table-driven
executor that runs it, and the dry-run backend that meters it without touching a real
tensor. It supersedes the incremental design drafts that preceded it; where a draft and
the code disagreed, the code is what is written here.

**What phase 1 does.** A batched cost model (`PeakBatchedModel`,
`SeQuant/core/optimize/cost_model.hpp`): ORDERED nests of batched modes per node,
contracted AND external opens priced together, recompute and bound-persistence charged,
selection on a three-objective `(peak, flops, nsl)` frontier. A loop-identity pass
(`compute_dag_boulevard`, `SeQuant/core/eval/peak_profile.hpp`): physical loops as
connected components over OCCURRENCES, value identity defined from them. A schedule
builder (`build_ordered_schedule`, `SeQuant/core/eval/ordered_schedule.hpp`): those
loops realized as a nested block tree, with per-nest pass splits and multi-level escape
chains. An explicit cell table (`cell_table.hpp`, `cell_table_builder.hpp`) plus a
static validator, so a schedule is checked in milliseconds before it runs. A
table-driven executor (`ordered_executor.hpp`) with no inference of its own: every read,
slice, life, persistence and assemble is read off the table. And a dry-run backend and
metered walk (`backends/dryrun/`) driving the SAME driver entry a wet run does.

**What phase 1 does NOT do.** It does not choose fusion assignments optimally (where
producer/consumer connectivity conflicts, the union-find records A valid choice); does
not reorder loop groups across spaces for cost (the realized nesting order is the one
the DP's occurrences witness); does not do cost-aware materialize-vs-recompute across a
pass split (builder rule 4 materializes unconditionally); does not share transients
across reads (a value the table holds no cell for is recomputed inside each production
tree that needs it); does not support an Assemble scattering more than one position of
one value at one loop instance (a joint sub-block write -- the executor THROWS; see
section 8.5); and does not change backends' tile geometry or the factorizer.

## 2. Terminology

The canonical vocabulary -- mode, index, the TAPP index roles, external vs contracted,
batch vs slice, `BatchModeType`, node-local vs forest-level batching, and the
new-vs-legacy identifier map -- lives in `doc/dev/batching-mode-terminology.md`. Read it
first; this document uses those terms without redefining them.

Three terms are specific to this document. A **value** is an identity in the fused DAG
-- phase 1 defines it as `EvalExpr::value_key()` (section 5), not the bare node id. An
**occurrence** is one use-site of a value, the triple (value, consumer value, leg),
recorded as `OccurrenceRec` (`peak_profile.hpp`). A **cell** is one FORM of a value
resident at one scope -- which positions are sliced, which loops it is a partial sum
over, where it lives, how it is produced (`TableCell`, `cell_table.hpp`).

## 3. Array modes and canonical layout

An eval node's ARRAY MODE LIST is its `canon_indices()`
(`SeQuant/core/eval/eval_expr.hpp`, `eval_expr.cpp`), produced by the core canonicalizer
-- `TensorNetwork::canonicalize_slots` for a proto-carrying tensor,
`TensorBlockCanonicalizer` for a plain one. For layout the list is split into PROTO-FREE
OUTER modes first and PROTO-CARRYING INNER composites second: exactly what
`EvalExpr::indices_annot()` emits (`outer + ";" + inner`), the annotation every backend
receives.

Three rules follow, and all three are load-bearing:

1. **Labels bind modes for composition only.** An index label is meaningful inside
   one occurrence's frame: it says which mode of the left operand pairs with which
   mode of the right at one contraction, and nothing across occurrences or trees.
   Every cross-occurrence or cross-frame decision here is by POSITION in
   `canon_indices()`, never by label -- see `detail::home_modes_in_cell_frame`
   (`cell_table_builder.hpp`), which translates a node's own home into the rich
   cell's frame positionally.
2. **A leaf is laid out per the REQUESTING node's `canon_indices()`.** The leaf
   evaluator is handed that list and must return the array in that layout. The
   fixtures' `rand_tensor_yield` (`tests/unit/test_eval_ta.cpp:205`) generates one
   random array per tensor in a fixed slot layout and serves each request as a
   permutation of it, so every occurrence sees the same tensor in the layout it
   asked for.
3. **Tensor SLOT ordinals are never used for array construction.** Rebuilding a
   node's mode list from the tensor's bra/ket/aux slots (explicit indices plus
   appended composite protos) can repeat an index and hand the backend a
   degenerate contraction. `canon_indices()` is the only source.

For a composite (CSV/PNO) index `a<i,j>` the array carries mode `a` over an `<i,j>`-tied
range. A batch loop is always over a PLAIN mode; slicing an `i` loop does nothing to
mode `a`. `eval::detail::slot_modes_of` (`lifetime_mask.hpp`) encodes this: a composite
slot contributes `a`, never its proto pair.

## 4. The batched DP

Entry points: `PeakBatchedModel` (`SeQuant/core/optimize/cost_model.hpp`) --
`build_context`, `Context::build_cells`, `relax`, `select_root`,
`reconstruct_batched_modes`, `is_external_mode` -- driven by `run_single_term_opt`
(`optimize/single_term.hpp`).

### 4.1 Cells are ordered nests over batchable modes

`Context::batchable_modes` is the ordered, deduplicated list of batchable indices; bit
`k` is mode `k`, `m` their count. A DP CELL is not a subset but an ORDERED SEQUENCE of
batched modes, outermost first, identified by `id`. `Context::cell_union(id)` recovers
the order-independent bitmask the footprint tables are indexed by (slicing is a
footprint change, so size depends only on WHICH modes are in the union);
`Context::descend(id, Ap)` appends a set as inner positions; `Context::escaped_outer(id,
carried)` gives the enclosing modes a node is actually re-executed across -- only those
OUTER to its innermost carried placement, none at all when the node carries no enclosing
mode and hoists above the whole nest.

`Context::build_cells` enumerates every sequence up to `Context::cap` = `min(m, 3)`,
grown by `popcount(external_mask)` under spectator batching so a contracted nest of the
previous depth still fits under an external one. On an enumeration blowup (estimate
above 100000 cells) `Context::ordered` is cleared, `id` degenerates to the plain
bitmask, and every helper reduces to set arithmetic: a bounded-nest schedule is still
CORRECT, only optimality is lost.

### 4.2 Per-node opens

At node `n` in enclosing cell `B`, `relax` forms

    contracted_here = (open_modes[lp] | open_modes[rp]) & ~open_modes[n]
    ext_here        = external_mask & open_modes[n] & ~cell_union(B)

and enumerates every subset `S` of their union, splitting it into `Ap = S &
~external_mask` (modes contracted here) and `E = S & external_mask` (external loops
OPENED here). Children are read at `descend_pt(B, E, Ap)` -- the external loops sit
OUTSIDE the contracted ones opened at the same node, exactly the nesting order the
schedule builder realizes. Both gates carry the same unlimited-budget revert: an
infinite `peak_threshold` forces `contracted_here == 0` and `ext_here == 0`, so the
batched model produces the same schedule as the unbatched one. `contracted_here` is
additionally forced to 0 when `batch_persistent_only` is set and the subset holds a
volatile leaf (`volatile_mask & n`); `ext_here` additionally requires
`batch_spectator_indices`. `is_external_mode` decides membership of `external_mask`: a
mode OPEN on the root subset and open in every subset where a carrying leaf is present
-- contracted at no node, so slicing it is work-neutral.

### 4.3 Costs

- **Footprint.** `both = szlp + szrp + szn + contrib`. An external open produces the
  node per batch (sized on `cell_union(B) | E`) and scatters into the full, pre-sized
  result held across the loop; with a contracted mode also sliced here the per-batch
  block accumulates and is weighted by `accumulation_factor`. A purely contracted open
  charges `accumulation_factor * szn` for the in-flight contribution co-residing with
  the accumulator. Any non-empty `S` also puts `szn` into the resident-scan term, so
  each ancestor batching node's accumulator stacks onto every descendant's peak.
- **Recompute** (`charge_batch_recompute`, default on): `rf` multiplies `nbatches[k]`
  for every `k` in `escaped_outer(B, open_modes[n])`.
- **Bound persistence** (`charge_bound_persistence`, default on): at a VOLATILE node
  with a non-empty `S`, a NON-volatile child carrying a mode of `S` is scaled by
  `volatile_weight` -- the loop runs once per production of the node that CLOSES it, so
  such a child is rebuilt every replay and cannot persist. The scale keys on the whole
  `S`, contracted and external alike (an external loop the DP opens here is closed here
  too, the full result being assembled outside it). A nest closed by a NON-volatile node
  runs once and is charged once.

### 4.4 Objective and frontier

A frontier point (`BFrontPoint`) carries `peak`, `flops`, the two children's frontier
indices, `aprime` (contracted opens), `eopen` (external opens), and `nsl`, the
cumulative count of batchable modes sliced anywhere in this realization;
`pareto_insert_ceiling` keeps `nsl` as a third objective only when `perf_first &&
isfinite(peak_threshold)`.

`select_root` on the perf-first arm treats `peak_threshold` as a CEILING: among points
whose byte peak fits, fewest flops, ties to fewer `nsl`, then lower peak -- so a term is
batched only when its cheapest schedule would exceed the budget, with no free slicing
below the ceiling. If nothing fits, the tiebreak INVERTS to min flops then min PEAK (the
most-sliced realization); keeping the least-sliced schedule there would blow the budget
by the most. `reconstruct_batched_modes` then walks the chosen realization and stamps
each node: `node_slice_mask` gets `External` for every external mode of the node's
enclosing cell or opened at it that the node carries plus `Contracted` for `aprime`;
`opened_here` gets `eopen` as `External` (only under `batch_spectator_indices`) AND
`aprime` as `Contracted`, a contracted loop opening at its unique reduction node. There
is no post-DP external placement pass.

### 4.5 `BatchPolicy` fields

`SeQuant/core/batch_policy.hpp`:

| field | meaning |
|---|---|
| `is_batchable_contracted_index` | spaces batchable where SUMMED; feeds the DP's contracted candidates |
| `is_batchable_external_index` | spaces batchable where open on the term root; feeds `external_mask` |
| `is_batchable_index()` | DERIVED union, never settable; the runtime's accept predicate |
| `batch_target_size` | per-index per-batch slice size, an UPPER BOUND rounded down to a tile multiple (one-tile floor) |
| `is_volatile_leaf` | which leaves are amplitude-dependent |
| `batch_spectator_indices` | admits external modes as nestable cell modes; without it the DP opens none |
| `persistent_only` | decline to batch any subtree holding a volatile leaf (default off: batch across the board) |
| `accumulation_factor` | in-flight batch-contribution footprint multiplier |
| `scheduler` | `forest_descent` (default, legacy) or `ordered` (this pipeline) |
| `peak_threshold` | budget in BYTES; a finite value is what ENABLES batching at all |

## 5. Loop identity and occurrence facts

### 5.1 Identity versus layout

In `SeQuant/core/eval/dag_scope.hpp`, `LoopKey{depth, loop_slot}` is the STABLE IDENTITY
of a physical batch loop: `depth` says which loop GROUP (separating even two groups of
one space -- an external occupied group beside a contracted occupied one -- which
`space` cannot), `loop_slot` which MEMBER of that group. `LoopKey::color()` packs both
into one opaque `size_t` where a single integer must distinguish loops. `DagScopeLevel`
adds LAYOUT: `altitude_ordinal` (nesting rank of the slot within its group) and
`latitude_ordinal` (the PASS index within a forced-split nest). `space` is a color for
group matching and is NEVER identity.

The cell table, the per-occurrence atlas and `ordered_n_batches_by_loop` (one batch
count per LOOP, not per loop group) key on IDENTITY. `CellScope::encloses` deliberately
compares full `(LoopKey, latitude)` entries, two passes being genuinely different
scopes; residency and read multiplicity compare identity alone through
`detail::same_key`.

### 5.2 `compute_dag_boulevard`: loops first, values second

`compute_dag_boulevard` (`peak_profile.hpp`) stamps `stamp_lifetime_masks` and
`stamp_occurrence_homes` (`lifetime_mask.hpp`), walks every tree in post-order assigning
monotone static points, and records one `NodeRec`/`OccurrenceRec` per occurrence. It
then numbers loops BEFORE defining values -- the ordering that breaks the catch-22,
since a value keyed by node id alone forces two physically distinct loops through one
shared node.

A loop instance is a connected component of a CONFLICT-AWARE union-find whose nodes are
`(occurrence, position)` and, for a mode contracted in batches, `(occurrence, reduced
mode)`. Edges join a home-sliced position to the same physical mode at the parent
occurrence within a tree; across trees, occurrences sharing a group key (node id plus
home-sliced positions) FOLD position by position, each fold a `try_unite` attempt. Four
constraints reject a fold:

1. **Position conflict** -- a component never holds two positions of ONE
   occurrence (one batch loop slices one mode of one array). This is the only
   PHYSICAL constraint; one loop may slice different positions of one NODE in
   different occurrences.
2. **Forbidden partner** (`forbid`) -- a value's reduction loop never unites with
   a loop enclosing, or opened by, an occurrence that reads that value COMPLETE;
   otherwise the reader would sit inside the loop it waits for.
3. **Kind** (`comp_kind`) -- External and Contracted components never unite: one
   scatters into the result, the other sums, and a realized block is one kind.
4. **Nesting** (`inner_of`) -- two components one of which must (transitively)
   enclose the other never unite; this keeps the builder's nesting constraints
   acyclic.

Components are numbered per space in first-seen order; the number is the `loop_slot`.
Each occurrence then carries `loop_slot[p]` per carried position (-1 where not
loop-sliced) and `reduced_slot` per mode reduced in batches; `RichSchedule::loop_kind`
records each instance's kind and `RichSchedule::loop_order` the `(outer, inner)` pairs
read off the occurrences' enclosing chains.

### 5.3 Value identity

With every occurrence's positions and reductions stamped, the value key is computed
bottom-up over the post-order records:

    value key = hash(node id,
                     (position, slot) of every home-sliced position,
                     (index, slot) of every mode reduced in batches,
                     value key of the left and right operands)

or the node id alone when nothing in the subtree is sliced. It is stamped on the forest
node (`EvalExpr::value_key`, read via `value_key_of`) so the boulevard,
`analyze_legality`, the builder and the executor all agree; occurrences group into
`ValueCell`s by it (`ValueCell::key` beside `ValueCell::hash`). Nothing changes for a
schedule with no home-sliced value: every key equals the hash.

### 5.4 Per-occurrence home, and the three mask fields

`home_scope(n)` returns `n->occurrence_home()` -- the loops opened at or above this node
filtered to its own result slots, stamped per occurrence by `stamp_occurrence_homes`. It
is NOT a cross-occurrence meet; `[lifetime_mask][seed]` pins both halves.

Three distinct fields on `EvalExpr` are routinely confused:

- `node_slice_mask()` -- the batchable indices the DP chose to slice AT this node, each
  tagged `Contracted`/`External`, stamped on EVERY carrying node. A BUILD-SITE role
  source, never a per-occurrence use-induced fact (it is empty on a top-homed leaf's
  deep fetch).
- `batch_loops_opened_here()` -- the subset of `node_slice_mask()` for which this node
  is the LOOP-OPEN site, each physical loop appearing exactly once. Anything
  reconstructing a loop NEST must read this or it multi-counts a loop once per carrying
  node.
- `sliced_modes()` -- the cross-occurrence MEET: modes slicing EVERY occurrence of the
  canonical node, accumulated top-down from opens (`stamp_lifetime_masks`) and filtered
  to the node's own result slots. Empty = all-full / loop-invariant. The forest path's
  residency reads this.

## 6. Ordered schedule builder

`build_ordered_schedule` (`ordered_schedule.hpp`) lowers `RichSchedule` plus
`LegalitySchedule` (`legality.hpp`: `build_site_of`, `classify_axis`,
`analyze_legality`, per-value `CellLegality::per_axis` of `LoopLocal` / `Reduction` /
`LoopCarried`) into an `OrderedSchedule`: a tree of `ScopeBlock`s holding `Step`s (a
`BuildStep` or a nested block), each block carrying `axis`, `level` (its
`DagScopeLevel`), `kind` (`Contracted`/`External`) and `outputs` (`value_id ->
OutputKind`).

### 6.1 Per-instance chain, then clusters

The realized chain is PER LOOP INSTANCE, not per space: for each space, the distinct
FUSION `loop_slot`s appearing on it across all cells each become one realized depth
(`slots_of_space`, `types`, `type_slot`, resolved through `fusion_slot`). A doubles
amplitude carrying two same-space externals therefore gets TWO loops, not one -- the
collapse fix. `fusion_slot` returning -1 makes a `LoopCarried` caller fall back to slot
0, but a `Reduction` mode with no slot is a HARD ERROR and throws. The chain is nested
to satisfy every `RichSchedule::loop_order` pair, the space-major slot-ascending order
being only a stable tie-break; a cycle is a builder error (loop identity fused two loops
that nest in opposite orders).

**Co-occurrence clusters (un-fuse).** Two loop members co-occur iff some value is
home-sliced on both; members that never co-occur go into DISJOINT nests. Two residual
sub-DAGs that both batch the same space but are connected only through a whole symmetric
intermediate are genuinely separate loop groups, and realizing them as one over-deep
chain deadlocks. A union-find over member depths (`type_cluster`, `cluster_min`) cuts
the single chain into independent nests, concatenated at root.

Each block's steps are ordered by a REAL topological sort
(`ordered_schedule_topo_sort_steps`) over dependency edges recovered from `RichSchedule`
alone (`ordered_schedule_dep_graph`, from every occurrence's `consumer_point`); a nested
child block's `requires_` is its whole subtree's bubbled external need, so a need
satisfiable only several levels out still surfaces where it can be met. Ties break on
`first_use`.

### 6.2 Placement, escapes, and escape chains

A value whose `per_axis` is entirely `LoopLocal` is a plain `BuildStep` at the deepest
depth any of its LoopLocal modes resolves to (`local_home_depth`) and gets NO `outputs`
entry anywhere: `Transient` is the ABSENCE of an escape, never a recorded `outputs`
entry (a recorded one would double-produce it under `well_formed`).

A value with at least one `Reduction` or `LoopCarried` entry ESCAPES: it has no
`BuildStep`, and is an `outputs` entry of each escaped instance's block --
`AccumulateSum` for a Reduction, `AccumulateScatter` for a LoopCarried. With more than
one escape this forms a multi-level ESCAPE CHAIN: raw production at the DEEPEST escape
site, pure forwarding at every shallower one. Such a chain may SKIP a level the value is
invariant on: the fused chain nests a term inside loops of OTHER terms, and a value
crosses such an interposed loop not by an escape of its own but by RESIDENCY plus
`produce_if_absent` (section 7).

### 6.3 Forced producer/consumer split, per nest

`forced_split_levels` assigns every value an integer PASS, global over every batched
space. Two edge kinds bump a reader's pass, both keyed on the OPERAND: (1) the operand
is `LoopCarried` on some axis -- its full array exists only after its own loop closes,
so every direct reader is bumped; (2) the operand is a `Reduction` on some instance AND
the reader is produced INSIDE that instance (the caller-supplied `inside(reader,
source)` predicate resolves this exactly as escape placement does), so it would
otherwise read the current batch's partial sum.

A forward sweep computes each value's base pass; a reverse sweep LIFTS a non-pinned
value with consumers to the minimum of its consumers' passes, so a value whose readers
all sit later is built with them. Sources of a bumping edge are PINNED and skipped by
the lift. A nest holding members of more than one pass emits ONE SIBLING BLOCK PER PASS,
run in schedule order, disambiguated by `latitude_ordinal`. Sibling pass blocks sit at
the same nesting level; each pass's cells are distinct because scopes carry their
latitude, and a later pass's read names the earlier pass's assembled cell explicitly.
`fork_subchain` produces the predicate-false / predicate-true copies of an already-built
inner sub-chain when a nest must be split below.

Rule 4: a value with a LATER-PASS reader produced in the SAME nest is materialized --
every instance of that nest it is loop-local on and not already scatter-escaped is
escaped as a Scatter (a Reduction escape on such an instance is UPGRADED to Scatter,
since a loop-local mode's batches are disjoint and summing them would combine values
that must stay separate). A later-pass reader produced in a CONFIDENTLY DIFFERENT nest
throws, as does a reader of a reduction inside the reduction's own loop in the SAME pass
-- the pass levels should have bumped that one, so it stands as a tripwire on the
levels, not an expected outcome.

## 7. Cell table

`cell_table.hpp` defines the model; `cell_table_builder.hpp` derives it from the
finished schedule (`build_cell_table` over `CellTableInputs`).

### 7.1 Cells

A `TableCell` is one FORM of one value at one scope: `value_id`; `scope` (a `CellScope`
-- the enclosing `(LoopKey, latitude)` path, outermost first, empty = root); `sliced`
(carried POSITION -> the loop instance slicing it, empty = whole); `partial_over` (the
loop instances it is a partial sum over -- a reduced mode has no carried position, so it
cannot appear in `sliced`); `production` (`Build`, `Leaf`, or `Assemble` with `assemble`
= `Sum` | `Scatter`, a `source` cell, and for a Scatter a `scatter_map` of target
position -> the loop instance whose batches are scattered into it); plus
`produce_if_absent`, `persistent`, `life`.

Three helpers are the single source of three decisions. `detail::bound_instances(c)` =
`sliced`'s instances plus `partial_over` -- the one place "is this cell bound to
instance k" is answered. `detail::residency_scope(s)` is the residency CEILING, by
production kind: `Leaf` -> root; `Build` -> its own scope always (a step's value lives
in its own block and dies with it, whole or not); `Assemble` -> the prefix of its scope
path ending at the DEEPEST bound instance, or ROOT when bound to none -- that last case
is what lets a term's complete partial cross an interposed foreign loop.
`detail::deepest_visible_form` is the "which form is visible here" rule: among a value's
cells, the deepest-scoped one whose residency scope encloses the query, earlier
candidate winning a tie. The builder's read-source selection and `CellRegistry::cell_of`
both use it, so table and runtime cannot disagree.

`detail::read_multiplicity(source, consumer_scope, n_batches_of)` weights one read by
`max(1, n_batches(k))` for every loop on the consumer's scope path NOT on the SOURCE'S
RESIDENCY scope path; the builder's life bookkeeping and validator rule 4 both call it,
so they cannot drift.

### 7.2 Reads

A `Read` is one per LEG of a consumer's production tree, not one per distinct operand
value: a consumer contracting one value with itself carries TWO reads of it. It names
`consumer`, `operand_value_id`, `source`, the `slice` (operand position -> loop
instance), and `invariant_on` (loop instances on which the seam explicitly recorded that
no slicing decision bound this read). The seam's facts are per leg --
`SlicedModeAssignment::occ_facts` / `occ_invariant` carry the consumer's leg as their
last element -- and the builder pairs each leg's read with its own.

### 7.3 `produce_if_absent` and persistence

`detail::set_residency_flags` (`cell_table_builder.hpp`) sets both from the cell's bound
instances against its enclosing path:

    produce_if_absent = !path.empty() && !bound_on_innermost
    persistent (candidacy) = !volatile_value && !bound_on_any_enclosing

`produce_if_absent` means the cell is produced on the FIRST iteration of a loop it is
not bound to and reused on the rest -- the mechanism by which a value crosses an
interposed foreign loop. `CellRegistry::clear_bound_to` clears such a cell only on the
instance it IS bound to, so it survives the batches it is invariant across.
`detail::apply_persistence_frontier` then narrows persistence to the FRONTIER, demoting
a candidate that has at least one consumer but NO volatile consumer: such a cell is
never read again, because those consumers are themselves skipped by cache-halt on every
later evaluation. As built:

    persistent = !volatile
                 && bound_instances(cell) is empty
                 && (some consumer's value is volatile || no consumer at all)

This is NOT the legacy runtime's "survives the reset of its home scope" flag, which
tests only the INNERMOST home loop: a cell bound to an OUTER loop passes that test and
fails this one.

### 7.4 Validator

`validate_cell_table(table, root, n_batches_of)` returns `CellViolation{rule, what}`
entries; `assert_valid_cell_table` throws on any. Five rules:

1. **visibility** -- every read's source must have been produced earlier (or be a
   Leaf) and its RESIDENCY scope must enclose the consumer's scope; an Assemble's
   source must be produced.
2. **form** -- per Build cell and per loop instance it is BOUND to, every operand
   bound to that GROUP (keyed on `depth` alone, since group members share a
   depth) must be bound to the SAME instance; an UNDECIDED whole operand on such a
   group is a mismatch. A read recorded `invariant_on` that instance is neither
   bound nor whole and never contributes. NECESSARY, NOT SUFFICIENT: the dry-run
   range check remains ground truth.
3. **chain** -- an Assemble must strictly enclose its source's scope and share its
   value id; a Scatter's every `scatter_map` entry must name a position the source
   is sliced on by that instance, and an empty map is a violation. A Sum's CLOSING
   instance is the loop just below the Assemble's own scope ALONG THE SOURCE'S
   PATH (not the source's innermost loop -- foreign loops may be interposed below
   it); the source must be a partial over it and bound to NONE of the deeper
   interposed loops. **(3b)** A cell with a non-empty `partial_over` may be
   consumed only by the Assemble that closes it.
4. **life** -- reads weighted by `read_multiplicity`, plus one per Assemble
   consuming the cell, must equal `life`; Leaf cells are skipped. A produced cell
   nobody reads is dead work UNLESS its RESIDENCY scope is root -- those are the
   schedule's results.
5. **uniqueness** -- at most one Build per (value, scope).

`CellTable::unresolved` is diagnostics, not model: (cell, position) pairs in the value's
sliced modes that matched no enclosing instance of the same space. Such a position is
recorded WHOLE, so a non-empty list on a schedule expected to be fully resolved is to be
treated like a violation.

## 8. Executor

`ordered_executor.hpp`. `sequant::evaluate(forest, policy, layout, leaf_evaluator,
cache, mode_order, make_scope_guard)` routes on `policy.scheduler`;
`BatchScheduler::ordered` goes to `eval::evaluate_ordered`, which builds schedule and
table, validates the table (`assert_valid_cell_table` THROWS on any violation), and
walks the block tree.

### 8.1 Reads go through the resolver

`CellReadResolver::fetch(value_key, ctx)` names the source cell for the current
consumer's next declared Read, spends one of its declared lives, applies the Read's
declared slices against the batch context, and returns the value in the registry's
CANONICAL orientation; `compute_cell` converts it to the operand node's own orientation
with `apply_canon_phase`. `CellRegistry::read` moves the value out on the last
non-persistent read rather than copying, and `operand_drained(key)` is the table-side
answer to the in-place-accumulation provenance question: true only when the most recent
fetch drained a non-persistent source.

### 8.2 `compute_cell` and `apply_one_op`

`compute_cell(node, cell, resolver, leaf_evaluator, cache, ctx)` evaluates ONE
production; every operand leg is a table Read. Two cases do not resolve to a held cell.
On **a leaf's FIRST touch** `fetch` defers (leaving the Read unconsumed), the leaf
evaluator runs on the WHOLE leaf, the result is recorded as that leaf's cell in
canonical orientation, and the Read is then served -- which is what applies the declared
slice; serving the whole leaf directly would hand the consumer a whole operand where the
schedule says a batch slice. **A node the table holds no cell for** is a TRANSIENT of
this production tree, computed in place from its own operands, recursively, by the same
rules; the executor never hands a production back to the tree-walking engine.

The op itself is always `apply_one_op_traced` (shaped-product hook, recompute tally and
trace event included), or `sum_in_place_traced` for an accumulating Sum whose left
operand is an internal node the table says is drained. The production's own node and
every transient of its tree go through the SAME gate, so they cannot answer the
provenance question differently.

### 8.3 Blocks, assembles, reset, and block close

`run_ordered_contracted_block` realizes one `ScopeBlock`. `Contracted` and `External`
blocks are realized UNIFORMLY; their difference is carried entirely by each Assemble
cell's own kind. Any other `BatchModeType` is refused loudly.

Per batch, in order: (1) `registry.clear_bound_to(block.level.key())` -- the per-batch
reset expressed on cells, dropping every cell bound to THIS block's own loop instance so
a stale prior-batch cell is never read as this batch's; (2) push the batch range onto
the batch context; (3) run the block's steps -- a `BuildStep` resolves its own Build
cell AT THIS EXACT SCOPE (a miss means table and schedule tree disagree), calls
`compute_cell`, and `registry.set`s the canonical orientation, while a child block
recurses; (4) run each output as an Assemble step, which reads the source's per-batch
cell, spending its own declared read, then accumulates. `Sum` takes the value by move on
the first batch when that read drained its source, else clones it; `Scatter` allocates
its destination from the Assemble cell's OWN form -- the value's canonical descriptor
narrowed by the Assemble's own sliced bindings against the current batch context, with
no separate destination-sizing inference -- and writes at the scatter map's positions.
At BLOCK CLOSE the Assemble CELL is recorded (`registry.set(out_cells[k], ...)`), not a
raw result, so the next link of the escape chain reads a cell like any other.

The `CacheManager` this route wires uses `min_repeats = 1` (2 on master): the table
derives every life from its own exact per-cell read count, which a cache declining
single-use nodes would falsify.

Cache-halt is a skip set over cells computed once per evaluation from the persistent
store, closed under "a cell all of whose consumers are skipped is skipped", with a
per-visit seed from held, unbound `produce_if_absent` cells
(`ordered_visit_skip_seedable`). A skipped cell FORGOES its declared reads with the
multiplicity of the loops the skip collapses, so its sources are released exactly as if
the reads had happened.

### 8.4 The scope guard, and why it exists

`run_ordered_contracted_block` takes a `ScopeGuardFactory` and instantiates
`make_scope_guard(batches.size())` AT BLOCK ENTRY; nested blocks' guards stack by RAII.

This is not optional bookkeeping. TiledArray screens result tiles on a Cauchy-Schwarz
bound; for a partial over 1/n of a reduced index that bound is roughly 1/n of the
full-sum bound, so a tile that survives the full contraction can be screened away in
EVERY batch and vanish from the blocked sum. MPQC's CSV path runs the CC equations under
a RAISED sparse threshold and compensates by handing the evaluator a scope-guard factory
that DIVIDES the threshold by the batch count while a batched partial runs. The executor
accepted the factory but never called it; the omission stayed invisible until the DP
placed a contracted loop, inside an occupied nest, around a node with many
near-threshold tiles, and the water-8 wet gate lost 7.8e-6 Eh from the second iteration
on. Calling the factory at block entry is the fix: guarded partials keep marginally more
tiles than one full product, agreeing with the unguarded reference to 5e-10. The unit
fixtures cannot reproduce this -- they run at TiledArray's default threshold -- so it is
a WET gate and stays one.

### 8.5 Refusals

`SEQUANT_ASSERT` compiles to `do {} while(0)` unless the build sets
`SEQUANT_ASSERT_BEHAVIOR` to ABORT or THROW, which `CMakeLists.txt` does only for Debug.
Every refusal of an input the pipeline cannot express is therefore a `throw
sequant::Exception`, holding in the Release / RelWithDebInfo builds that run production
work: an unsupported `ScopeBlock` batch-mode kind or `Step` variant or escape
`OutputKind`; an Assemble whose `scatter_map` is empty or names more than one position;
an Assemble step that realized zero batches; a missing `BackendArrayOps`; a value id
absent from the forest's value-node map; a scheduled value never produced; and
`build_ordered_schedule`'s `well_formed` self-check. Genuine internal invariants --
range and size checks that cannot fail unless the code is wrong -- stay
`SEQUANT_ASSERT`. `[ordered-executor][refuse]` pins one refusal in whatever build the
suite is compiled with.

## 9. Dry-run backend, metered walk and report

`SeQuant/core/eval/backends/dryrun/` provides a sizing-only backend (`eval_expr.hpp`,
`result.hpp`, `cost_model_object.hpp`, `size_regime.hpp`): no real tensor is allocated,
every op reports flops and footprint through the `CostModel`, and the strict range
checks (lobound and extent on shared labels in prod / sum / accumulate / scatter) plus
the fill-once tripwire run on every step -- the FIRST gate on any schedule or builder
change, milliseconds, no execution.

`meter()` (`backends/dryrun/meter.hpp`) runs the POLICY-SELECTED executor over a dry-run
forest with a fresh, `PeakMonitor`-wired, build-tallying cache and returns a
`MeterReport`: `peak_bytes` and the op hash where the peak occurred, persistent/volatile
FLOPs and `CostModel` exec-time split (classified by `compute_volatility`, the same
bottom-up rule the gated `cache_manager` factory uses), `builds_total`, the per-value
`HomeFidelity` list (builds versus where the value is homed and used), and the `scheduler`
the report describes. It drives the SAME driver entry (`sequant::evaluate(Nodes const&,
BatchPolicy const&, ...)`) a real solve uses, selected by the same `policy.scheduler`, so
the metered replay is the run the policy describes, not a proxy; the forest-descent
evaluator is installed at `Trace::On`, since `note_working_set` (which feeds `peak_bytes`)
is compile-time gated on it. The `[dryrun-2iter-report]` fixture runs two iterations, cold
and warm, forest versus ordered, reporting builds / FLOPs / peak at water-20 residual
scale; it is how persistence and cache-halt effects are measured.

## 10. Diagnostics

Every environment-gated dump the ordered evaluator has lives in
`SeQuant/core/eval/ordered_dump.hpp`, whose file doc block carries the COMPLETE table of
variables and what each prints, plus a pointer to the cache / runtime / optimize knobs
documented with their own facility; it is not duplicated here. Two design properties:
every dump is OFF unless its variable is set (an unset variable costs one `std::getenv`
test and the result is byte-identical), and the header takes its operands as template
parameters rather than including the schedule / table / executor headers, so it sits
BELOW all of them with no include cycle. The comparison hooks the scope-guard
investigation used -- `Result::norm2`, `layout_desc`, `tile_diff` -- remain on `Result`
(`result.hpp:353-367`); the environment-gated cross-check probes that drove them do not.

## 11. Test fixtures and gates

`tests/unit/test_cell_table.cpp` (20 cases) pins the model and the validator on
hand-built tables: scope prefix and equality; residency (a whole Assemble is resident at
root, a plain Build inside a loop is confined to its own scope regardless of
`persistent`); each of the five rules' violation and its non-violating twin; that life
weighs a read against the SOURCE'S RESIDENCY scope, not its raw scope; that an
explicitly invariant read is not a form mismatch while the same read without the marker
is; that two legs of one consumer reading one value are two reads; and that an escape
chain may skip a level the value is invariant to. `tests/unit/test_legality.cpp` (8
cases) pins `build_site_of` and `classify_axis`: loop-local versus reduction on a value
carrying one occupied index and contracting another in batches; the two same-space
cases; and the `forced_split_axes` fixpoint.

`tests/unit/test_ordered_executor.cpp` -- (H) = a Catch2 `[.]`-HIDDEN case, run only by
naming its tag, never by the default suite:

| tag | what it proves |
|---|---|
| `[b-full]` | cells are computed through `apply_one_op` ONLY: on an unbatched scalar forest, one build per `BuildStep` and zero cache probes |
| `[w20-auxocc-walk]` | the water-20 aux+occ residual dry-run walk completes with no vanished home value -- the multi-space, multi-instance schedule |
| `[cell_table][ordered]` | the w20 default and INPUT-MIRRORED configurations each derive a VALID table; the mirrored one RUNS the strict walk to completion |
| `[cache-halt]` | a resident persistent composite is not re-formed, and its dead batch-block prerequisites are skipped |
| `[block-skip]` | a loop-invariant escape is not re-formed on later batches |
| `[pia-rebind]` | a `produce_if_absent` cell bound to an enclosing loop IS re-produced on that loop's batches; only an UNBOUND one may seed a per-visit skip set |
| `[cell_registry]` | a leaf's first touch inside a batch loop is served through its Read, with the declared slice |
| `[assemble-dest]` | an Assemble's scatter destination is sized from its own cell form |
| `[ordered-executor-witness-water20]` (H) | ordered and forest descent model the same result shape |
| `[meter]` (H, `b3`) | `meter()`'s ordered replay peak equals a direct `evaluate_ordered_schedule` replay EXACTLY (one-byte margin, deliberately not a tolerance) |
| `[dryrun-2iter-report]` (H) | section 9 |

`tests/unit/test_eval_ta.cpp`, real TiledArray arrays, none hidden, ground truth = the
unbatched evaluation of the same forest:

| tag | what it proves |
|---|---|
| `[mixed-open]` | two External occupied loops plus one Contracted occupied loop opened at ONE node; the contracted loop nests between the externals, the value is built inside all three, partial over the contracted loop, escaping scatter -> sum -> scatter. Dense |
| `[mixed-open-pia]` | the same shape over an aux-contracted inner product INVARIANT on the inner external loop, so the `produce_if_absent` path runs |
| `[mixed-open-pia-sparse]` | the same, under a SPARSE policy |
| `[mixed-open-tot]` | the nested-array (ToT) replica of the culprit term: pair composites, the amplitude sliced on both outer modes under external > contracted > external |

The WET gates live on the MPQC side and are named here only: the water-8 gate
(RelWithDebInfo, asserts live; lossless energy AND per-iteration residual norm, for
aux-c, aux-c/occ-e and aux-c/occ-e/occ-c), and the water-20 run in Release on the Linux
box.

## 12. Legacy route kept, and open items

### 12.1 The batched forest descent

`BatchScheduler::forest_descent` remains the DEFAULT and is kept for now: an
unconditional forward to the pre-existing `sequant::evaluate(forest, layout,
leaf_evaluator, cache)` overload (`eval.hpp`), with no schedule built and nothing from
sections 5-8 running. It evaluates one tree at a time, hoists a value to an enclosing
scope through `ensure_hoist_slot`, and slices a hoisted value to the current batch on
fetch through `slice_to_use` (the "hops" model), placement driven by the
cross-occurrence meet `sliced_modes()` rather than a table. With a finite
`peak_threshold` the DP still batches, so this is a genuine BATCHED route, not an
unbatched fallback -- simply the non-DAG one: a value shared across trees is rebuilt per
tree. The `stamp_lifetime_masks` meet, the node-keyed cache, the shared-buffer in-place
test and the `PeakSink` metering (hence the evaluator's `Trace` parameter) are kept for it.

### 12.2 Known limitations and open items

- **`Context::ordered == false` is untested.** The bounded-nest degradation
  (`build_cells`' blowup guard, above 100000 estimated cells) is believed
  correct-but-suboptimal, but no fixture forces it.
- **The Adjoint trace label differs between paths.** `apply_one_op_traced` (the ordered
  op path) appends ` | L.<trange> O.<trange>`, ` R.<trange>` interposed for a binary op;
  `evaluate_impl`'s `Stage::NeedLeftAdj` arm emits the bare `log::label`. Same value.
- **An Assemble scattering more than one position is unsupported** -- it would need a
  joint sub-block write; the executor THROWS on `scatter_map.size() != 1` (see section
  8.5).
- **`validate_cell_table`'s block-tree walk is not implemented.** `root` is reserved for
  the intended visibility rule; as built, visibility comes from cell scopes and emission
  order (7.4 rule 1).
- **The form rule is necessary, not sufficient** -- it accepts a read as whole on
  `Read::invariant_on`; the dry-run range check is ground truth.
- **`CellTableInputs::operands_of` has a lossy fallback.** Without per-leg operands the
  builder falls back to the DE-DUPLICATED `depends_on`, so a self-contracting consumer
  gets one read, not two.
- **`residency_scope` is a CEILING.** The legacy close-store walk stops EARLIER (at a
  consumer-internal loop, and at the cache chain's end), so a deeper hold is no defect.
- **Fusion choices are not optimal.** Where connectivity conflicts `try_unite` records A
  valid assignment (5.2), the loser read transposed; minimizing transpositions and
  reordering loop groups across spaces are future work.
- **A reader of a reduction inside its own loop in the same pass is unsupported;**
  `build_ordered_schedule` throws, a tripwire on the pass levels. **Cost-aware
  materialize-vs-recompute** is likewise absent; rule 4 always materializes.
- **The ordered-vs-forest FLOPs gap at tight budgets is open.** The fused chain nests
  terms under other terms' loops and rebuilds what those loops do not skip
  (`SEQUANT_UT_TOP_BUILDS` is the tool); the batched cost models also run with NO
  subnetwork CSE, costing terms standalone.
- **`subnet_cse` is asserted off in Debug builds ONLY** (`SEQUANT_ASSERT`,
  `single_term.hpp:104,129`): a Release caller must not set it with a batched objective.
- **Thirteen `[blocked-layers-1-2]` fixtures are hidden and NOT retargeted**
  (`test_ordered_schedule.cpp` 7, `test_eval_ta.cpp` 5, `test_eval_dryrun.cpp` 1). Each
  encodes the older loop-identity / layout model the identity rework superseded, not a
  known-good test switched off, and their "blocked on Layers 1-2" reason has EXPIRED
  (6.2 / 6.3 / 7 are built). At least two FAIL when run: `test_eval_ta.cpp`'s "batched
  ToT External occ loop" and `test_ordered_schedule.cpp`'s "forced-split occ axis
  realizes TWO ordered sibling blocks" -- the latter sets only a `BatchPolicy` role
  predicate and never stamps `node_slice_mask`, so the builder realizes NO loop and
  finds zero occ blocks: its INPUT contract is stale and 6.3 stands as written. Each
  file carries a note saying so; none should be un-hidden without rewriting its
  expectation first.
