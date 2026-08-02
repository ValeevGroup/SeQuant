# Placement as register allocation: values, instances, cells

_Status: design. Author discussion: 2026-08-02._

## 1. Problem

The batched evaluator maintains a **cache hierarchy**: one persistent
*run-scope* cache (survives every iteration; holds cross-iteration
intermediates) plus a nest of transient *batch-scratch* caches, one per batch
loop, re-created on each iteration of the loops enclosing it. A value's correct
placement is whichever level it is invariant across — genuinely shared
intermediates high (run scope), batch-variant ones deep (the scratch of the loop
that varies them).

Today that placement is produced by **three mechanisms that do not agree with
each other**, and it does not respect the memory budget:

1. **Run-scope population** (`cache_manager.hpp`) is a static use-count CSE walk:
   count each node's consumer edges, cache those `>= min_repeats`, then **veto**
   the batch-variant ones (non-empty cross-occurrence lifetime mask) so they do
   not sit in the persistent cache. The veto is the *negative* half of a
   placement — "not here" with no statement of where.

2. **Runtime hoist** (`eval.hpp`, `place_at_this_level`/`collect`) decides, per
   batch level and lazily during evaluation, which nodes to hoist into that
   level's scratch, refusing a node with `has_demoted_external`.

   These two read *different* signals and diverge on the demoted cross-pair
   giant: its cross-occurrence meet empties the external slot, so
   `mask_all_full()` is true and the run-scope veto **admits** it, while its own
   `batched_here` still carries the External stamp, so `has_demoted_external` is
   true and the runtime **refuses** to hoist it. A value that wants to live in
   two places at once, with no way to name the second.

3. The **use count is batching-blind**: it counts static consumer edges, so a
   value re-read every iteration of a loop it is invariant to is counted `1`,
   falls below `min_repeats`, is not cached, and is **recomputed per iteration**
   — a direct component of the avoidable recompute the dry-run cost model
   reports, invisible to the count.

And the deeper problem: hoisting a shared value to a common scope to reuse it
**lengthens its live range**, so it is co-resident with the inner loops' working
sets for that whole span — which can violate `peak_threshold`. The factorizer
cannot see this, because it is **per-term** and CSE (the sharing that triggers
the hoist) is **cross-term / cross-occurrence**, hence invisible to it. Its
per-term peak is a lower bound; the true peak appears only once CSE + placement
are resolved over the whole forest.

## 2. Goal

Stop asking "which single home does this node get" (ill-posed for a shared,
budget-constrained value) and instead solve the actual problem:

> given the CSE'd value-DAG and a memory budget, choose how many copies of each
> value to materialize, where each lives, and which consumers it serves, to
> **minimize recompute subject to the peak staying under `peak_threshold`**.

This is **register allocation with loop-invariant code motion and
rematerialization** (§6). The rest of this note builds the vocabulary (§3-4),
the reuse measure (§5), the decision and its shape (§6), the pass (§7), the
correctness cases (§8), and the mapping onto SeQuant (§9).

## 3. Three identities: value, instance, cell

- **Value** — what SeQuant hashes (`TreeNode`). The mathematical object: a
  canonical (sub)expression. Perfect CSE is exactly "one value = one node," so
  the perfectly-reused DAG is a DAG of *values*.
- **Instance / occurrence** — a single *use-site* of a value: a position in a
  specific term's tree. One value has many instances (what CSE folds together).
- **Cell** — a *chosen physical copy* of a value serving some **subset of its
  instances**: one build, one live range, one cache entry. A cell is a
  **partition-cell of the value's instances**.

```
value ──1:N──► instances ──partitioned into──► cells ──cover──► instances
```

The decision is a **partition of each value's instances into cells**. Its
extremes are the two DAGs:

- **perfect reuse** → one cell per value (cells ≡ values) → the DAG of unique
  nodes;
- **no reuse** → one cell per instance (cells ≡ instances) → the original
  per-term forest.

Everything between is a **materialization DAG**: the value-DAG with the split
values *duplicated* — a *partial un-CSE*. So the materialization DAG may hold
several cells carrying the **same value** (several copies of one node), each
wired to its own consumer subset. That is the crux: once peak forces un-sharing,
"the node" is no longer one thing — its instances split into cells, and
placement is a property of the **cell**, not the value.

## 4. Cell identity

A cell needs a cache key that keeps two cells of the same value from collapsing
back into one copy. Its identity has up to three coordinates:

```
cell = ( value , home-scope , split-index )
```

- **value** — the hash. *What* it computes.
- **home-scope** — *which loop level of the batch nest the cell lives at*. This
  is the coordinate batching adds; the flat value-DAG has none. Two cells of one
  value at different scopes are different copies (one held full at run scope,
  one sliced per block in a scratch).
- **split-index** — a discriminator for two cells of one value at the **same**
  scope, split purely for peak. This is the register allocator's move: splitting
  a live range gives each piece a **fresh name** (an SSA subscript); the
  split-index is that name.

Perfect CSE fixes `split-index ≡ 0`, `home-scope ≡ meet`, so the key degenerates
to the value — today's hash-keyed cache. At runtime a cell homed at a batch
scope is *re-instantiated per iteration* of the loops above it (the scratch
rebuilds each pass), but that is the **same static cell** — one materialization
decision, many temporal realizations.

**Realization: `(home-scope, split-index)` selects a store, not a key.** The
`home-scope` and `split-index` coordinates are *not* folded into the cache's
key. Instead there is **one value-keyed store per `(home-scope, split-index)`**
(see §7a), so the cache map stays keyed by the value (`TreeNode`) exactly as
today, and a cell is identified by *which store* it lives in. `home-scope`
already corresponds to a store (the run-scope cache vs. a per-loop scratch);
`split-index` adds a second store at the *same* scope only where a peak split
occurs. Perfect CSE is byte-identical because it uses one store per scope with
`split-index ≡ 0`.

**Terminology.** *Home scope* = the scope where a cell is built and lives; every
read is at-or-below it. It is the same thing the code calls the *lifetime scope*
and one might call the *store scope* — this note standardizes on **home scope**
(and flags in code that "lifetime scope" is the synonym). *Use scope* is the
consumer's scope. The pair **home-scope vs. use-scope** is the read's two ends.

The three split reasons (why a value gets more than one cell):

- **temporal** — free from the scope hierarchy: a cell homed in a scratch is
  rebuilt per outer-iteration, so instances in different iterations are already
  separate cells;
- **structural** — from the cross-occurrence meet: instances binding a mode to
  incompatible blocks are demoted, so they *cannot* share a full value (the
  demoted giant);
- **peak** — the new one: instances that *could* share, split anyway because the
  shared cell's live range would exceed the budget. This is the only reason that
  needs new mechanism (a `split-index` and a group-scoped cache key); (a) and
  (b) already fall out of the scope hierarchy and the meet.

## 5. `W(cell)`: the batching-aware reuse measure

Cache-worthiness of a cell is "would building it once and reusing it save work?"
The honest measure is **how much of the value's data is actually read within the
cell**, in units of one full value. An access reads a fraction `f`; the count is
`Σ f` over the cell's runtime accesses.

For a consumer `C` enclosed in batch loops with extents `B_L`, partition the
loops between the cell's home and `C` into

- **S** = loops whose mode the value **carries** (the access is a slice), and
- **I** = loops whose mode the value does **not** carry (the access reads full).

`C` runs `∏_all B_L` times, each reading `1/∏_S B_L` of the value, so contributes
`∏_all / ∏_S = ∏_I B_L`. Summing over the cell's consumers:

```
W(cell) = Σ over the cell's consumers C of  ∏_{loops L between home(cell) and C : value does not carry L.mode}  B_L
```

= the number of full-value-worths of data read within the cell. Properties:

- **batching-decoupled** — a `B=1` loop contributes 1; a slicing loop (value
  carries the mode) contributes 1 (its `B` disjoint slices tile one full read).
  Batching on carried modes leaves `W` identical to the unbatched schedule.
- **redundancy counted** — a loop of extent `B` the value is invariant to reads
  full `B` times → contributes `B`.
- **strict generalization** — no batching → `W` is the integer consumer count,
  so `min_repeats >= 2` keeps its meaning: `W = 1` (read once, no reuse) → do not
  materialize a shared cell; `W >= 2` → worth it.

`W` is **rational** (`f = 1/∏_S B_L`); use Boost rational (already a SeQuant
dependency) so the products/sums stay exact.

Crucially, `W(cell) - 1` times the build cost is the **avoidable recompute** the
cost model reports for that cell: the cache-worthiness count, the split decision,
and the avoidable-recompute measure are one quantity.

## 6. The per-value decision is 2-D — and it is register allocation

For each value the pass chooses, jointly:

1. **how far to hoist** each cell (its `home-scope`), and
2. **how to partition** the instances into cells.

The two trade **peak against recompute**:

- coarser partition / higher hoist → more reuse (less recompute) but longer live
  ranges (higher peak);
- finer partition / lower hoist → shorter ranges (lower peak) but more recompute.

This is **register allocation** with the batching extensions:

- **scope = loop nest.** A live range lives at a *level*; hoisting a value out of
  a loop to reuse it is **loop-invariant code motion** — LICM under register
  pressure is exactly the hard coupling flat RA avoids.
- **pressure is per-scope-per-iteration.** The program point whose live set is
  measured is `(static location × loop iteration)`; the peak is the max over the
  whole nest, not a linear timeline.
- **partial hoist (slicing).** A value carrying a batched mode is a family of
  slices; a "register" can hold 1/B of a tensor, and rematerialization can rebuild
  a *slice* (cheap) rather than the whole value — a granularity flat RA lacks.

Splitting a shared cell = **live-range splitting with rematerialization** (spill
by recompute, not reload). Minimizing recompute under a hard memory bound on a
DAG is NP-hard in general (the pebble game / `revolve` / gradient-checkpointing
family), so the design targets a **heuristic** with a **detection safety net**,
not a closed form.

## 7. The pass, the objective, and the honest boundary

**Peak is a placement constraint, not a factorizer constraint.** The factorizer
owns *arithmetic* (per-term, CSE-blind, unchanged); a whole-forest,
**post-CSE** placement pass owns *memory*. `peak_threshold` is the placement
pass's constraint. Only that pass sees the value-DAG (CSE) and the global live
set.

Objective:

```
minimize  Σ over cells of (W(cell) - 1) * build_cost(value)      # total recompute
subject to  peak-profile(materialization DAG) <= peak_threshold  # whole-forest peak
```

`W` is the objective; the peak profile is the constraint; the `meet` partition
(perfect CSE) is the constraint-free optimum. The pass:

1. **seed** with perfect CSE: one cell per value at its `meet` home (max reuse,
   max pressure);
2. while the peak profile exceeds threshold, **split the worst live range**:
   pick the cell whose materialization contributes the most peak per unit
   recompute cost and split it (partition its instances / lower its home /
   slice it further), inserting the recompute — the register allocator's
   spill-the-highest-pressure move, adapted to rematerialization;
3. **stop** when the profile fits, or when no split reduces peak.

**Detection safety net & the boundary.** The whole-forest peak is exactly what
`cost_profile()`'s replay already computes, so a violation is *detectable*. If
even full splitting (recompute everything locally, no shared cells) still exceeds
threshold, the peak is **inherent to the factorization** — a single large
intermediate placement cannot fix. That is the boundary of a two-phase design:
detect it via the replay peak and report it back (or feed a hint to
re-factorization); you cannot repair a factorization choice from the placement
pass.

## 7a. Runtime realization: cell stores + an explicit router

The placement decision is made concrete by two pieces, decoupling *placement*
(the router) from *mechanism* (dumb stores):

- **Cell stores.** One value-keyed store per `(home-scope, split-index)` — the
  existing per-scope `CacheManager`, keyed by `TreeNode` unchanged; a peak split
  adds a second store at the same scope. Scope *lifecycle* still nests along the
  loop nest (a scratch store is created per outer-iteration and reset), so a
  cell's lifetime is its store's.
- **Router (new).** An explicit map `{value, use-site} -> (home-scope,
  split-index)` — the placement pass's output, and the materialization decision
  of §3 made first-class. The default entry is `{value} -> (home, 0)` covering
  *all* use-sites; peak splits are per-use-site *overrides*, so the finer
  `use-site` granularity costs nothing where nothing is split. The `use-site` is
  the consumer's own position in the eval traversal (occurrence-id), finer than
  the tree-id because a value can occur more than once within one tree.

This **replaces the implicit `parent_` fall-through search**. Today a read walks
the scope chain (`access_at`, counting `hops`) to find a value stored at an
outer scope. With the router the read is direct:

```
read N at use-site U in use-scope:
  (home-scope, split) = router[{N, U}]
  full_or_partial      = store(home-scope, split).access(N)
  result               = slice full_or_partial by  (use-scope - home-scope) INTERSECT carried(N)
```

The slice is **the existing Enter-stage formula** (`eval.hpp:595`, "(use scope
MINUS lifetime scope) INTERSECT carried"): only the loops between home and use
*that N carries* need narrowing; the rest leave N full. The single change is
that the store scope is taken **directly from the router** instead of being
derived from `hops` (`d - hops`). So the slicer logic is unchanged; routing is
explicit; and the `parent_` walk and `hops` are subsumed.

## 7b. The greedy split pass (O2)

The pass is a **spill loop** in the register-allocation sense: start at the
recompute optimum and relax it just enough to fit the memory budget.

**Scope — a second phase on a FIXED schedule.** O2 runs *after* the min-time
factorizer, which is per-term and picks both the factorization *and* the
batch-loop assignments (`batched_here` — which modes are batched at which
nodes), targeting the **per-term** peak. O2 takes that factorization and loop
nest as **fixed** and decides only the whole-forest **eval/placement** strategy:
where each cell lives (home-scope) and the cell partition (router). It never
adds, removes, or re-assigns a batch loop. The reason the second phase exists is
that the per-term, CSE-blind factorizer only bounds the *per-term* peak;
CSE-hoisting can blow the *whole-forest* peak, which only a whole-forest,
post-CSE pass sees.

**State.** A placement = the router `{value, use-site} -> (home-scope,
split-index)` (§7a). **Seed** = perfect CSE: one cell per value at its `meet`
home. This is the **recompute-minimal** placement (maximum sharing, maximum
hoisting) and therefore the **peak-maximal** one.

**Objective / readout.** Minimize recompute subject to `peak <= threshold`. The
recompute of a placement is exactly the avoidable-recompute the cost model
already reports: `Σ (W(cell) - 1) * build_cost`, i.e. `cost_profile()`'s
`avoidable_flops`. So `cost_profile()` returns *both* the objective
(`avoidable_flops`) and the constraint (`peak_bytes`) for any placement — the
greedy needs no new measurement, only the ability to attribute the peak to the
cells alive at the peak point.

**The moves.** The peak at a program point is the sum of the footprints of the
cells **alive** there. To reduce a cell's contribution, either shrink it or
evict it from that point:

- **Shrink — home into an existing carried loop.** A cell homed *above* a batch
  loop `m` whose mode it *carries* holds `m` full; re-homing it *inside* that
  (already existing) loop lets the loop slice it, dropping its footprint by `B_m`
  at ~no recompute (`W` counts carried-mode slicing as free, §5). This is a
  **placement** choice on the fixed nest -- it uses a loop the factorizer already
  placed; it does NOT batch a new mode. (Deciding to batch a mode that is *not*
  batched -- adding a loop -- is the factorizer's lever, not O2's; see the
  boundary below.) Nearly free; try first.
- **Evict — shorten the live range.** For the peak a slice cannot reach:
  - **delay / un-hoist** a cell *invariant* to a loop it is currently held
    across: lower its home into (or toward) that loop so it is built lazily at
    its use instead of held idle from the loop's start. Footprint unchanged;
    lifetime shortened, so it is no longer co-resident at the earlier iterations.
    Recompute = rebuild once per outer block (the `∏_I` factor of §5).
  - **split instances** of a cell whose consumers span a long range: partition
    its use-sites into groups (a new `split-index`), each group its own short-
    lived cell at its own `meet`. The long shared live range breaks into short
    ones. Recompute = one extra build per new group.

  Shrink applies to a loop the cell *carries* (re-home into it, the loop slices
  it); delay/un-hoist is for a loop the cell is *invariant* to (re-home into it,
  rebuilt per outer block, same size). Both are placement on the fixed nest --
  **all O2 moves are CSE-aware placement the per-term factorizer cannot make**
  (they act on shared cells / cross-occurrence live ranges it never sees). What
  O2 does *not* touch is `batched_here` (which modes are batched) -- that is the
  factorizer's lever, reached only via the re-batch feedback below.

**The loop.**

```
seed = perfect CSE (router: {value} -> (meet, 0) for all use-sites)
while peak(placement) > threshold:
    p     = the binding peak point (argmax live-set)
    cands = cells alive at p
    move  = best over cands of:  prefer any zero-cost shrink;
                                 else max  ΔPeak / ΔRecompute  (the spill metric)
    if no move reduces peak(p): break        # infeasible: see Termination below
    apply move to the router; incrementally update the peak profile
report (avoidable_flops, peak_bytes)          # objective + constraint
```

Focusing candidates on the **binding peak point** (the max-pressure point) is the
register allocator's "spill from where the pressure is." `ΔRecompute` is the
added `avoidable_flops`; `ΔPeak` is the drop in the profile's max. Apply,
**incrementally** re-cost the affected lifetime interval (not the whole forest),
repeat.

**Termination.** `peak <= threshold` (success), or no placement move reduces the
binding point — two distinct failure modes, both detected via `peak_bytes` and
neither repairable *within* O2:

- **factorization-inherent** — a single intermediate is larger than the budget
  even fully placed; the factorization itself must change;
- **re-batch-needed** — the fixed batching left placement too little room (e.g. a
  *shared* cell would need slicing on a mode the per-term factorizer never
  batched, because no single term saw the sharing). This is not an O2 move; it
  feeds back to the factorizer to add batch loops (§7, O6).

So the honest structure is `factorize + batch (per-term peak) -> O2 place
(whole-forest peak, fixed batching) -> if infeasible: re-batch`.

**Sub-items.**

- **O2a** the incremental peak-profile update after one move (which lifetime
  intervals change), so the loop is not `O(moves × forest)`.
- **O2b** the exact `ΔPeak`/`ΔRecompute` estimator per move, and whether the
  greedy needs lookahead (a shrink that enables a later cheaper evict) or a
  single pass suffices in practice.
- **O2c** relation to the existing DP external-slice pass: subsume it as the
  shrink move, or run it first (DP-local shrink) then this pass (whole-forest
  evict) — the latter is the smaller change.

## 7c. Sizing cells and the peak profile (O3)

O2's `ΔPeak` needs two things: the footprint of a cell at a placement, and the
peak of a whole placement, attributable to the cells that set it.

**Cell footprint is home-relative.** A cell homed at scope `S` holds each index
it carries at the extent seen *at* `S`. For a carried mode `m` batched at the
(fixed) loop `L(m)`:

```
extent_at_S(m) = block_extent(m)   if L(m) ENCLOSES S   (cell built inside the loop -> sliced)
                 full_extent(m)     otherwise            (cell built above L(m) or m unbatched -> full)
footprint(cell) = product over carried modes m of  extent_at_S(m)
```

This is exactly the existing moment-aware `memsize` (CSV/composite moments) with
home-relative extent overrides — a cell held *full* above a loop it carries (the
demoted-giant option) is sized full; re-homing it inside that loop (the shrink
move, §7b) sizes it by `block_extent` — so O2's `ΔPeak` for shrink is just the
footprint delta. Non-carried modes do not contribute.

**The peak profile is max weighted-interval overlap.** Linearize the *static*
schedule tree in execution order (one representative iteration per loop; same-loop
different-iteration cells are never co-resident, and the single-iteration
linearization captures that, while a cell homed above a loop spans that loop's
whole subtree). Each cell is an interval `[first-use, last-use]` in that order
(from the router's use-sites + the schedule order) weighted by its `footprint`. A
cell contributes its footprint to *every* point in its live range (it is stored
at its home granularity; deeper consumers slice on access, they do not shrink the
stored cell). Then

```
peak = max over static points p of  Σ footprint(cell) over cells live at p
```

a sweep line over weighted intervals — and the argmax point is exactly the
**binding peak point** O2 spills from, with its live set the spill candidates.

**This corrects today's under-count.** `cost_profile()`'s `peak_bytes` is
`max(batched-scratch hwmark, cache working_set_hwmark)` — a `max`, so it
under-counts when a persistent cross-term cell co-resides with a batched-inner
transient (§1 flagged it a lower bound). The sweep **sums** all live cells at a
point, so it is the *true* peak; O3 both models the peak for O2 and fixes the
measurement. The replay remains the ground-truth oracle: the sweep must equal a
replay that sums (not `max`es) co-resident residency — that equality is O3's
validation.

**Incrementality (feeds O2a).** The profile is weighted-interval overlap, so an
O2 move changes one cell's interval endpoints (a home change moves its live
range) and its weight (a footprint change), and the sweep updates over just the
affected span — not the whole forest.

**Sub-items.** O3a the sweep/interval structure and its incremental update; O3b
the replay-with-summed-co-residency oracle to validate the sweep against; O3c
composite (proto-indexed) sizing under home-relative overrides (reuse
`memsize_counter`'s moment path).

## 7d. Home scope and the demotion fold (O5)

`home_scope` is the seed of everything (the perfect-CSE placement O2 starts
from), and it must be *correct* before any peak move — it is a structural, not a
budget, decision. Definition:

```
home_scope(value) = deepest scope enclosing the loops of  ( sliced_modes  ∪  demoted_external_modes )
```

- **`sliced_modes`** — the cross-occurrence *meet* (`stamp_lifetime_masks`): the
  modes on which the value is *consistently* sliced across all its occurrences.
  Homing inside their loops is the highest (max-reuse) legal placement — the
  perfect-CSE upper bound.
- **`demoted_external_modes`** — External `batched_here` stamps the meet
  *demoted* out of `sliced_modes` (the `has_demoted_external` condition):
  occurrences bind such a mode to *incompatible* blocks, so the value's
  per-occurrence value is block-specific and **cannot be shared as one full value
  above that loop**. Home must therefore also be *inside* those loops.

**The fold is adding `demoted_external_modes` to the meet set** — and it is
exactly what removes the split-authority bug (§1). Without it, a demoted value
has an all-full meet mask, so the cache veto reads "run-scope cacheable" while
the runtime's `has_demoted_external` refuses to hoist it — two answers. With it,
`home_scope` is inside the external loop, the cache never gets a run-scope cell,
the runtime slices it, and **both sides read the one `home_scope`.**

**Realization is temporal, not a split-index.** A demoted value homed at the
external loop is *one* static cell, re-instantiated per block by the loop (the
free temporal split of §4); each occurrence's read lands on its own block via the
router + the `(use - home) INTERSECT carried` slice (§7a). No `split-index` — that
coordinate is reserved for O2's *peak*-driven same-scope splits. Demotion is
structural.

**Order.** `home_scope` (with the fold) is computed from the meet *before* O2:
it fixes the perfect-CSE seed's homes correctly (structural correctness); O2 then
only lowers homes *further* for peak (§7b). Consistency with `W`: a demoted value
is sliced per block on its external mode (its loop encloses the home), so `W`
counts that mode as free tiling (§5) — the demotion adds no avoidable recompute,
as it should.

**Sub-items.** O5a confirm `has_demoted_external` is the exact demoted-external
signal and the edge cases (a value with *both* consistently-sliced and demoted
modes; nested/independent external loops); O5b the seed router construction
(meet home + demotion fold) as the input to O2.

## 8. Worked cases (correctness tests)

Each must be a unit test.

**(a) Invariant value under one looped consumer.** Invariant to loop `m`
(extent `B`), one consumer inside it. Today: edge count `1` → not cached →
recomputed `B` times. Here: `W = B >= 2`, one cell hoisted to its home if it
fits the budget → built once, sliced/fetched per iteration. The concrete
avoidable-recompute win.

**(b) Demoted cross-pair giant.** External slot demoted by the meet. The value's
instances **partition into per-block cells** (structural split); each cell homed
at the external loop with `W = 1` → sliced per block in scratch, never a full
run-scope copy. The mask-vs-`has_demoted_external` disagreement disappears
because the value simply has more than one cell, each with a name.

**(c) Genuinely shared gC.** `home == run`, several run-scope consumers,
`W >= 2`, fits budget → one cell at run scope. Must be unaffected.

**(d) One value, two consumers at different depths.** If the shared cell fits:
one cell at `meet`, the deeper consumer's accesses weight against home (`W`
inflated by its extra invariant loops), slice-on-use serves its fetch. If it
does not fit: the pass splits into two cells (per consumer), trading the reuse
for peak. Confirms both the cross-depth `W` and the peak-driven split.

## 9. Mapping onto SeQuant / what changes

- **value** = `TreeNode` (hash); **perfect CSE** = the hash-keyed
  `cache_manager` (its use-count walk, `cache_manager.hpp:618`).
- the **home-scope axis** already exists as the cache *hierarchy* (persistent
  run-scope cache vs per-loop scratch); it is simply not treated as a coordinate
  of a cell.
- the **cross-occurrence meet** (`stamp_lifetime_masks`, `lifetime_mask.hpp`)
  already performs the *structural* split (b); `place_at_this_level`/`collect`
  already perform *temporal* placement.
- **the demoted giant** is a value that already wants two cells; the current
  veto-vs-hoist split fumbles it because there is no way to *name* the second
  cell.

Changes required:

1. **Cell stores + router, not a wider cache key** (O1, resolved — see §7a). The
   cache stays `TreeNode`-keyed; a cell is a value-keyed store per
   `(home-scope, split-index)`, and an explicit router `{value, use-site} ->
   (home-scope, split-index)` (the placement pass's output) replaces the
   implicit `parent_` fall-through search. Reads route via the map then slice by
   `(use-scope - home-scope) INTERSECT carried(N)` — the existing Enter-stage
   slicer, fed the home scope directly instead of via `hops`. The default
   `{value} -> (home, 0)` for every use-site is byte-identical to today.
2. **Batching-aware `W`.** Replace the integer edge count with the rational
   `W(cell)` (§5), computed in the static walk (needs `home` and the loop nest,
   both static).
3. **A whole-forest, post-CSE placement pass** (§7) that seeds with perfect CSE
   and greedily splits live ranges against the peak profile, replacing the veto
   and unifying with `place_at_this_level`.
4. `BatchPolicy::is_batchable_contracted_index` (the batching *decision*) and
   `BatchPolicy::is_batchable_index()` (the eval accept union) are untouched.

Precursor already landed: disjunct (a) of the batch-variant veto and its
`is_batchable_index` parameter were removed as structurally dead (commit
`eval: remove the dead free-batchable-mode caching veto (disjunct a)`), leaving
only the mask-based disjunct (b) — the seam this design absorbs into cells.

## 10. Prior art / relation to the literature

We are **not** inventing the core problem; it is well studied in several
communities, and the design should borrow formulations rather than reinvent. The
core -- *store an intermediate vs. recompute it, under a hard memory budget, to
minimize compute* -- is the **rematerialization / recompute-vs-memory
scheduling** problem:

- **Register allocation.** Rematerialization instead of spilling (Briggs,
  Cooper, Torczon, PLDI 1992); graph-coloring allocation (Chaitin). The flat
  prototype; our `split-index` is its live-range-split rename.
- **Reverse-mode AD / checkpointing.** Griewank & Walther's `revolve` (ACM TOMS
  2000) is optimal for the *chain* case; the ML rematerialization line -- Chen
  et al., "Training Deep Nets with Sublinear Memory Cost" (2016); Jain et al.,
  **Checkmate**, "Breaking the Memory Wall with Optimal Tensor
  Rematerialization" (MLSys 2020) -- formulates recompute-under-a-memory-budget
  as an **ILP** over a tensor DAG. That is essentially our objective (§7), for a
  *fixed* DAG.
- **Tensor Contraction Engine (our exact domain).** The Sadayappan group's work
  on quantum-chemistry tensor contractions is the closest match: Cociorva,
  Baumgartner, Sadayappan et al., **"Space-Time Trade-Off Optimization for a
  Class of Electronic Structure Calculations"** (PLDI 2002) does the
  memory-vs-recompute tradeoff with loop fusion *for electronic-structure tensor
  contractions*; and "Memory-Constrained Data Locality Optimization for Tensor
  Contractions" (LCPC 2003). This is the literature to mine first -- same domain,
  same tiling/memory tradeoff.
- **Pebble games.** The space/space-time lower bounds and hardness (Sethi 1975;
  Hopcroft-Paul-Valiant) -- why the exact problem is NP-hard and heuristics are
  the norm.

So the *pieces* are off-the-shelf: adopt the **rematerialization framing**
(Checkmate-style ILP as an exact oracle for small forests; a greedy spill for
production), and reuse the **TCE space-time-tradeoff** algorithms for the
tiling/memory part.

What is **genuinely ours** (not directly in those works) is the *combination*:

- **cross-term / cross-occurrence CSE partition** (choosing cells) -- the
  checkpointing/AD work assumes a fixed DAG with no sharing-partition choice;
- **partial (sliced) rematerialization** -- textbook rematerialization is
  all-or-nothing (store the whole value or rebuild the whole value); batching
  lets a cell hold `1/B` and rebuild a *slice*, a granularity the general
  formulations lack;
- **placement across a batch loop nest** driven by the SeQuant cost model.

So: build on the established formulations; the contribution is the
tensor+batching+CSE cost model feeding a known memory-bounded-scheduling engine,
not a bespoke solver. Caveat: the exact ILP/DP are NP-hard / expensive for large
DAGs, so production ships a heuristic (greedy live-range spill), with the ILP a
reference oracle for validation on small cases.

## 11. Open items

- **O1 — RESOLVED (§7a).** Not a wider cache key: value-keyed stores per
  `(home-scope, split-index)` + an explicit `{value, use-site} -> (home-scope,
  split-index)` router; reads route then reuse the existing `(use - home)
  INTERSECT carried` slicer. Residual sub-items: (i) a cheap **use-site /
  occurrence id** for the router key (likely the eval-tree node pointer /
  position); (ii) **audit** the `parent_` / `access_at` / `hops` /
  `batch_context` users and re-express slicing via home-vs-use scope difference;
  (iii) **standardize naming** on "home scope" (code's "lifetime scope").
- **O2 — DESIGNED (§7b).** Spill loop from the perfect-CSE seed; moves are
  *shrink* (slice a carried mode -- the existing external-slice) and *evict*
  (delay/un-hoist an invariant cell, or split its instances -- the new
  CSE-aware move); candidates are cells alive at the binding peak point, ranked
  by ΔPeak/ΔRecompute; objective+constraint are `cost_profile()`'s
  `avoidable_flops`/`peak_bytes`. Residual: O2a incremental peak-profile update,
  O2b the per-move ΔPeak/ΔRecompute estimator and lookahead question, O2c
  subsume-vs-run-after the DP external-slice pass.
- **O3 — DESIGNED (§7c).** Cell footprint = home-relative `memsize` (a carried
  mode is sliced iff its fixed loop encloses the home, else full). Peak profile =
  max weighted-interval overlap (each cell an `[first-use, last-use]` interval
  weighted by footprint; sweep line), whose argmax is O2's binding peak point.
  This *sums* co-resident live cells, correcting today's `max(scratch, cache)`
  under-count (§1). Residual: O3a the incremental sweep structure, O3b a
  replay-with-summed-co-residency oracle to validate against, O3c composite
  (proto-indexed) sizing under home-relative overrides.
- **O4 — `W`'s computation order.** `W(cell)` depends on the partition, which the
  pass is choosing — a fixed point / two-pass (seed at `meet`, then refine).
- **O5 — DESIGNED (§7d).** `home_scope(value)` = deepest scope enclosing the
  loops of `(sliced_modes ∪ demoted_external_modes)`; the demotion fold (adding
  the demoted externals) is what unifies the current cache-veto-vs-
  `has_demoted_external` split into one authority. Structural, computed from the
  meet before O2 (which only lowers homes further for peak); per-block is
  temporal (no split-index). Residual O5a-b: confirm the exact signal / edge
  cases, and the seed router construction.
- **O6 — factorization feedback.** What to do when detection shows the peak is
  factorization-inherent vs. re-batch-needed (§7b): report vs. hint back to the
  DP (add batch loops / re-factorize).
