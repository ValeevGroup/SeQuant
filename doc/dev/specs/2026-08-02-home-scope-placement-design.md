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

1. **Group-scoped cache keys.** The cache map keys by value hash; give it a
   `(hash, home-scope, split-index)` key so two cells of one value do not
   collapse. This is the implementation friction point (O1): the hasher/
   comparator (`TreeNodeHasher`/`TreeNodeEqualityComparator`) and the map must
   admit the extra coordinates without breaking the perfect-CSE default
   (`split-index ≡ 0`, `home-scope ≡ meet` reproduces today's keying byte-for-
   byte).
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

- **O1 — group-scoped keying.** How `(hash, home-scope, split-index)` threads
  through `cache_manager`'s hasher/comparator and the scratch hierarchy without
  regressing the perfect-CSE path. The part with real teeth.
- **O2 — the greedy split move.** What "split the worst live range" does
  concretely: partition instances vs. lower home vs. slice further, and the
  peak-per-recompute ranking that drives it.
- **O3 — per-placement footprint / peak profile.** How a partially-sliced cell
  (some carried modes full, some sliced) is sized, and how the profile is
  evaluated, so it matches the replay peak `cost_profile()` already reports.
- **O4 — `W`'s computation order.** `W(cell)` depends on the partition, which the
  pass is choosing — a fixed point / two-pass (seed at `meet`, then refine).
- **O5 — canonical `home_scope`(cell)** as `meet` of the cell's instances (upper
  bound), refined downward by the peak constraint, with the demotion fold.
- **O6 — factorization feedback.** What to do when detection shows the peak is
  factorization-inherent (report vs. hint back to the DP).
