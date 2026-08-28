# Frame-correct use-induced sliced-mode map for the ordered DAG evaluator

**Status:** design draft — 2026-08-27
**Area:** `SeQuant/core/eval` — `compute_sliced_mode_assignment`, `build_ordered_schedule`, `placement_remat`, the runtime slice seam.
**Motivates:** aux+occ (external-occ) CCSD wet evaluation, which crashes today; aux-only works.

## 1. Problem

The ordered (DAG) evaluator slices values two different ways, and the runtime
seam (`occ_facts` / `LoopColoredSliceSeam`) that drives `slice_to_use` must
report **both**, per occurrence, in the occurrence's own index-frame:

1. **Home-slicing** — a value homed *inside* a batch loop is stored block by
   block (`A[_]=0; for i: A[i]=…`). Its sliced modes are recorded at home
   placement as `OrderedSchedule::home_mode_depth[vid]` (mode → DAG-scope
   *depth*, in the value's own frame) and consumed by the value-id coloring
   (`value_id_coloring`). This is the correct *concept*, but its current
   implementation is broken **the same two ways** as the use-induced seam and
   simply has not been exercised yet (aux-only, or cases with no sibling-ordinal
   split):
   - **Ordinal-blind, but NOT a cache-collision bug.** `home_mode_depth` stores
     `int depth`, discarding `ordinal` (the walk does `inner.push_back(
     {child->axis, child->level.depth})`); `value_id_coloring` colors by that
     `depth`. This does *not* cause collisions: each batched scope has its OWN
     scratch (`make_batched_scratch`), so values homed in distinct loops
     (distinct `ordinal`s) live in distinct caches regardless of color — a remat
     split into two slices homed at sibling loops is safe by construction.
     Within one scratch the level is constant, so what the color must separate is
     *which mode is sliced* (`ctx_modes`); the `depth` is at best a proxy.
     Whether the home coloring needs the level at all — versus just the sliced
     *slot* — is an OPEN question (§9), not a settled fix.
   - **Space-matched.** It resolves mode→carried by matching each level's
     *space* to the value's first-unconsumed carried mode of that space
     (lines 1410-1419) — the same `uniq`/space anti-pattern (§3). If two carried
     modes of one space are homed at one scope it picks arbitrarily; it should
     resolve the mode via `per_axis` (own frame, per-instance) instead. This is a
     real (if latent) correctness weakness, independent of the ordinal question.

   So the mode-resolution fix is **symmetric** across home coloring and
   use-induced seam (both should use `per_axis`, not space-matching); the
   *level*-keying question is settled only for the use-induced seam (needs full
   `{space,depth,ordinal}`) and OPEN for the home coloring (§9).

2. **Use-induced slicing** — a value homed at (or above) some scope, then
   **pulled into a deeper enclosing loop nest** at a use, is sliced *on fetch*
   to the current block for the loops crossed between its home and the use
   (`build_scope[consumer] \ home_scope[value]`, the "hops" loops). A top-homed
   leaf `g` fetched two loops deep (`home_d=0, hops=2`) is the archetype. This
   is what `occ_facts` is *for*, and it is where every current and historical
   construction is wrong for the aux+occ case.

Aux-only works because the aux (contracted) batch is decided by the single-term
optimizer and lands in `batch_loops_opened_here` → `OccurrenceRec::ectx`, which
`compute_sliced_mode_assignment` reads. External-occ batching is a **forest /
DAG-scheduler** decision (not single-term), so it is absent from `ectx` and
`node_slice_mask` and lives only in the ordered block tree.

## 2. The governing principle

**A DAG scope loop is identified by `DagScopeLevel = {space, depth, ordinal}`,
never by an index label.** Index labels are meaningful only inside a particular
value/occurrence's frame. `ScopeBlock::axis` is a per-*space* canonical
representative shared across a type's loops (verified: a forced split emits two
sibling blocks with the *same* `axis` and *same* `depth`, differing only in
`ordinal`), so it is **not** a loop identity and must never be matched against a
value's modes. `LoopId` (the index into `SlicedModeAssignment::levels`) *is* the
identity, because it carries `{space, depth, ordinal}`.

Corollary — `depth` alone is not a loop identity: sibling loops over one space
(e.g. a forced split's producer/consumer passes) share `space` **and** `depth`
and differ only by `ordinal`.

## 3. Why every prior attempt failed (do not revisit)

All measured on w8 with forced occ batching (`csv-cck-w8.json`); each crashed.

- **`OccurrenceRec::ectx` positional pairing** (the shipped baseline): `ectx` is
  a *pre-remat forest* artifact built from single-term `batch_loops_opened_here`.
  Remat moves values between scopes and the DAG adds external-occ loops the
  forest never opened, so `ectx` disagrees with the realized `build_scope`
  (occ loop simply missing for a top-homed fetch). Not the enclosure authority.

- **`node_slice_mask`**: a single-term *node* annotation (`eval_expr.cpp:612`,
  from `binarize`). It "does not deal with values" — it is populated on some
  forest nodes of a value and empty on others (empty on the `vid6` fetch
  occurrences, present on `i@o1` ones). Not a per-occurrence value fact.

- **`home_mode_depth`**: correct, but records **home**-placement slicing only —
  a value homed *inside* a loop. Top-homed leaves and hoisted intermediates
  (`home_d=0`) that are sliced *on fetch* have empty `home_mode_depth`. It is
  the wrong half (home, not use-induced).

- **Matching by index label ("exact" — `L.axis` in `occ.carried`)**: a
  cross-frame label match against the shared canonical `axis`; frame-dependent,
  the very thing to avoid. Under-slices a divergent occurrence whose external
  is relabeled.

- **Matching by space + uniqueness ("uniq")**: a hack that is only ever right
  when a value has exactly one mode per space. It **over-slices** a *contracted*
  index that merely shares the occ space (verified: it sliced `I(μ̃,i_2,Κ)`'s
  `i_2`, contracted in `g(i_2,i_1,Κ)·I(μ̃,i_2,Κ)→I(i_1,μ̃)`, while `g`'s `i_2`
  stayed full → BLAS overran → `corrupted size vs. prev_size`).

The lesson: the seam must **read** a frame-correct binding, never reconstruct it
by matching labels or spaces across frames.

## 4. The frame-correct data that already exists

- **`LegalitySchedule::cells[vid].per_axis`** — `AxisClass{axis, role}`, one per
  the value's own build-site axis, **in the value's own frame**, **per
  Index-instance** (an outer product `A{;i_3}·A{;i_4}` lists *both* `i_3` and
  `i_4`), with `role ∈ {LoopLocal, Reduction, LoopCarried}`:
  - `LoopLocal` — consumed within the value's build → home-sliced (this is
    `home_floor` / `home_mode_depth`).
  - `Reduction` — summed at/below the value → contracted; never sliced on fetch.
  - `LoopCarried` — survives into the value's own result → an external the value
    carries; sliced on fetch when pulled into the loop batching it.
  This is exactly the external-vs-contracted, per-instance, frame-correct
  classification the heuristics were failing to reconstruct.

  **Caveat (audit finding):** `per_axis` is the right *shape* (per the value's
  own build-site axis, per-instance) but its **producer**, `classify_axis`
  (`legality.hpp:234-278`), is **not** frame-pure: it derives the role by
  `base_key` (space, line 239), by walking `occ.ectx` (252-255), and by an
  across-occurrence label compare (`c == ectx_it->first`, 272). Within one
  non-divergent cell whose occurrences are position-aligned this happens to
  yield the correct role, but it is exactly the space/ectx/label machinery this
  design forbids, and it can misclassify when occurrences diverge. So the role
  is a *value* property that is well-defined in the value's own frame, but the
  current mechanism computing it is suspect. **The plan therefore does not build
  on `per_axis` blindly: Task 1 verifies (by measurement on the w8 repro) that
  the roles `classify_axis` emits are correct per-instance before anything
  downstream consumes them, and repairs the producer to a frame-pure derivation
  if they are not.**

- **`build_scope[occurrence's consumer]`** — the ordered block tree walk
  (`populate_build_scope_walk`) already records, per value, the enclosing blocks
  as `{axis, level}` with the **full `DagScopeLevel` including `ordinal`**. This
  is the per-occurrence enclosure, correct and available.

## 5. The gap

`build_ordered_schedule`'s per-value placement (~lines 1053-1068) maps each
`per_axis` mode to a **depth by space** (`depth_of_type(ac.axis)`) and
**collapses same-space instances**: two carried occ externals `{i_1, i_2}`
(a doubles amplitude) fold into **one** escape at the occ depth ("Same-depth
axis-classes … collapse to ONE escape at that depth, with LoopCarried
dominating"). The per-value **mode→(which loop)** binding is destroyed here, and
nothing downstream can recover it. This same collapse is why the original w20
external-occ scatter was single-mode (`dm[k]` = one axis).

## 6. Design: distinct per-instance loops (no collapse)

The factorizer (single-term DP) already batches per **index instance**, and this
is the model the runtime must preserve end-to-end. Verified in
`sliced_footprints`: a node's sliced-set is a **subset** (bitmask) of `aux_list`
— the *distinct* batchable indices — and each sliced index independently drops
its own mode to `min(extent, batch_target_size(ix))`, so slicing `{i_1, i_2}`
reduces `i_1`'s extent **and** `i_2`'s (multiplicative). `{i_1, i_2}` is therefore
**two distinct, independent loops**, 1-to-1 with the indices — not a group
collapsed to one, and neither is a spectator.

The *only* freedom among same-space loops is their **nesting order** (`i_1`
outer / `i_2` inner commutes with the reverse): a group of loops interchangeable
in **order**, but distinct in **identity**. Scheduling/remat breaks that order
equivalence by assigning each loop a definite `DagScopeLevel = {space, depth,
ordinal}`; the loop's *identity* (which index instance) is fixed throughout.
(A separate, orthogonal source of `ordinal` is the producer/consumer forced
split of a *single* index — §7/§9 — not the interchangeable-instance order.)

With loops kept distinct, `mode → loop` is **direct and per-instance**, no slot
or permutation indirection:

- **mode** = the value's own-frame index instance, read from `per_axis` (§4) —
  frame-correct, per-instance, never a label/space match.
- **loop** = that instance's realized `DagScopeLevel`/`LoopId` from the schedule
  (where remat fixed the order). Per occurrence, *which enclosing copy* of that
  loop (build vs consumer pass) is read from the occurrence's placement
  `build_scope[consumer]`, which already carries the full level with `ordinal`.

Genuinely divergent occurrences (different physical pairs) are un-folded by remat
`apply_split` into distinct non-divergent cells, so each cell's occurrences are
position-aligned and no per-occurrence mode table is needed; the per-occurrence
dimension is only "which enclosing loop copy," already in `build_scope`.

The current bug is precisely the **collapse** (§5): `build_ordered_schedule`
merges the DP's distinct per-instance loops into one-per-space and folds two
carried same-space instances into a single escape, destroying the 1-to-1 the
map (and the external-occ scatter) needs.

## 7. Code changes

1. **`build_ordered_schedule`: realize the DP's per-instance loops; stop
   collapsing.** This is the core change and the biggest one, because the
   *level/depth structure itself* assumes one loop per space:
   - `types` (line 992-999) keeps one representative Index **per space**, so the
     nest has exactly one depth per space. It must instead host **one loop per
     batched index instance** — two carried occ externals `{i_1, i_2}` become
     two distinct occ loops at **distinct depths** (nested; the DP's model is a
     multiplicative grid), interchangeable only in nesting *order*.
   - The escape collapse (lines 1049-1068: "two carried occ indices collapse to
     ONE escape") must be removed: each per-instance carried mode escapes to its
     own loop.
   - `depth_of_type` (a space→depth lookup) must become an instance→level
     resolution, since a space now spans several depths.
   - **Exactly how the depth structure grows to host per-instance same-space
     loops (nesting order, ordinal interaction with the producer/consumer split,
     `DagScopeLevel` enumeration) is the immediate open work — see §9.**

   *Home side (~1395-1447):* `home_mode_depth` must resolve mode→carried via
   `per_axis` rather than space-matching (§1). Whether it must key on the full
   `DagScopeLevel`/`LoopId` instead of `int depth` is OPEN (§9): per-scope
   scratches already prevent cross-ordinal collision, so the level-color may be
   redundant. Decide before touching the coloring.

2. **`compute_sliced_mode_assignment`: read, do not match.** For each
   occurrence, for each of its `LoopCarried` modes (from `per_axis`, own frame,
   per-instance), the loop is that instance's realized `LoopId` — read from the
   schedule, composed per-occurrence with the enclosing loop copy from
   `build_scope[consumer]` (which carries `ordinal`). No `axis`, space, or label
   match anywhere. Home-sliced (`LoopLocal`) modes stay with the value-id
   coloring and are excluded from the use-induced seam.

3. **External-occ scatter (`run_ordered_contracted_block`, `write_array_into_mode`):
   per-instance.** With per-instance loops restored, a scatter output carrying
   two batched occ externals must restrict *both* modes in its destination block
   (today `dm[k]`/`member_axis` handle exactly one). This is the runtime face of
   the same collapse and is fixed by the same 1-to-1 data.

The runtime seam and `slice_to_use` are unchanged in shape; they consult the
now-correct, per-instance `occ_facts`.

## 8. Invariants and validation

- **Frame purity:** no code path outside a single value/occurrence frame may
  compare index labels or space keys to assign a mode to a loop. Loops are
  addressed only by `LoopId` / `{space, depth, ordinal}`.
- **Completeness guard:** a value whose mode is *not* batched by any enclosing
  loop of a fetch (invariant there) is correctly unsliced; the guard must not
  fire on it merely because another occurrence is sliced. Re-express the guard
  against "this occurrence carries a mode bound to the crossed loop", not a
  hash-level participates().
- **Validation:** w8 with forced occ batching (`csv-cck-w8.json`,
  `occ_target_size` small) must run crash-free with zero unresolved facts, and
  its batched CCSD energy must match the unbatched reference (`peak_threshold`
  ≈ ∞) to full precision. Then w20 aux+occ (the original failure) lossless.
  Existing `[ordered-schedule]` unit tests must stay green.

## 9. Open items to resolve during planning

- **Depth structure for per-instance same-space loops (the core work):** `types`
  (one per space) and `depth_of_type` (space→depth) assume one loop per space. To
  host two independent occ loops (`i_1, i_2`) at distinct nested depths, resolve:
  how the depth chain is enumerated when a space contributes several loops;
  what fixes their nesting *order* (the interchangeable-order equivalence — is it
  arbitrary/canonical, or a placement/remat decision); how that interacts with the
  producer/consumer `ordinal` split (which is a *different* use of `ordinal`, over
  one index); and how `DagScopeLevel`/the `LoopId` enumeration change. This is the
  next investigation.

  **RESOLVED (hinge read, 2026-08-27):** the factorizer does **not** emit a nest
  for same-space external loops. `NodeBatchAnnotation::opened_here`
  (`cost_model.hpp:1999-2016`) is a **flat list** of `{Index, BatchModeType}`:
  external opens on the default root-seed path are *all* pushed at the term ROOT
  (`opened_ext_mask = (n==root) ? (chosen_seed_mask & open_modes[n]) : 0`), so
  `{i_1, i_2}` are co-opened at one node with **no inter-order** — a flat
  per-instance set, not a nest. (Contracted opens land at each mode's `aprime`
  node and *can* differ in depth; the node-level external path can too, but it is
  a default-off flag and not the failing case.) Therefore the nest is
  **synthesized at Stage 4** (`build_ordered_schedule`), not preserved from the
  factorizer: the factorizer's contract is *which* per-instance indices are
  batched, with roles; Stage 4 owns turning that flat set into distinct nested
  `DagScopeLevel`s and *chooses* the (interchangeable) order. This confirms §6
  ("only freedom is nesting order … Scheduling/remat breaks that order
  equivalence") and makes the collapse in §5 the single defect to remove: Stage 4
  currently collapses the flat per-instance set to one-loop-per-space instead of
  expanding it to one-loop-per-instance.
- **Instance identity across values:** the loops are shared across the values in
  a scope (CSE); confirm how a value's own-frame instance (`per_axis`) is bound to
  the *shared* physical loop it names — the residual's externals give a common
  frame within one residual, but pin the cross-value binding explicitly.
- **Home coloring level-keying:** does `ValueIdColoring` need the full
  `{space, depth, ordinal}` per sliced mode, or only the sliced *slot*
  (`ctx_modes`)? Per-scope scratches already separate cross-ordinal homes, so
  the `depth` color may be redundant. Establish what distinguishes two values of
  one node that must coexist in a *single* scratch — that, and only that, is
  what the color must encode. (The mode-resolution fix — `per_axis` over
  space-matching — is needed regardless of this answer.)

## 10. Non-goals

- No change to the single-term optimizer, `binarize`, `node_slice_mask`, or the
  aux (contracted) batch path, which already works via `ectx`. (Whether the aux
  path should also migrate to this scheme is a follow-up, not this change.)
