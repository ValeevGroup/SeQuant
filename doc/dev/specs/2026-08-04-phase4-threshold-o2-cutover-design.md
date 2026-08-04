# Phase 4: threshold-driven O2 spill + the CAT-1 cutover

Status: design, agreed via brainstorming 2026-08-04. Consumes Phase 1 (peak oracle),
Phase 2 (router read+store seam), Phase 3a (`home_scope` seed predictor), Phase 3b
(static `PeakProfile`), all complete on `evaleev/feature/multimode-batched-eval`
(`...f31bda794`). Realizes O2 (spec `2026-08-02-home-scope-placement-design.md` section
7b) and the CUTOVER. **Supersedes** the Phase-3 design doc's roadmap item "Phase 4 --
O2 greedy + CUTOVER ... lands COUPLED behind a flag": there is NO flag; `peak_threshold`
is the sole control (see Section 1).

**Naming.** The pass this document calls **O2** (the parent spec's ordinal for the
whole-forest greedy placement-relaxation pass) is named **`remat`**
(rematerialization -- https://en.wikipedia.org/wiki/Rematerialization) in the code
(`SeQuant/core/eval/placement_remat.hpp`), because its mechanism is rematerialization,
not register-spilling: it reduces peak by REBUILDING values in smaller/shorter-lived
pieces (trading recompute for peak), never by moving them to slower storage. "O2" and
"remat" are used interchangeably below; O2 is the spec ordinal, `remat` the code name.

## 1. End-state model: threshold-driven O2, no flag

After Phase 4 the runtime has ONE placement model:

- The **unified all-modes meet** is the single `residency` (drop the `External` filter
  in `ext_modes_of`; DELETE `contracted_modes`; `sliced_modes` becomes the meet of ALL
  batched modes on a node's result slots -- this is exactly `seed_residency` from 3a,
  promoted to the runtime path). `home_scope` = deepest enclosing loop over that residency.
- **O2 always runs**; `peak_threshold` is the sole control. `threshold = inf` => O2 makes
  no move => the raw perfect-CSE seed (peak-maximal); a finite threshold => O2 spills the
  demoted giants until `peak <= threshold`. The "off" state falls out of `threshold = inf`;
  there is no separate on/off flag and no parallel heuristic path.
- `has_demoted_external` (the demotion veto) is DELETED; its job -- keep the scattered
  giant from materializing full -- is now an O2 spill move (shrink), chosen cost-based.

This is a genuine CUTOVER, not an additive change: it replaces `place_at_this_level`'s
residency computation and its veto for EVERY placement. Correctness (CCk energies) is the
load-bearing gate; the peak is expected to CHANGE (drop toward the documented-RED witness
targets). This is the deliberate retirement of the heuristic the whole effort targets
(memory LESSON: "stop treating the current place_at_this_level heuristic as ground truth").

## 2. Structure: 4a (additive O2 pass) then 4b (cutover)

### Phase 4a -- O2 as a standalone static pass (additive, NON-regressing)

A pure function over the Phase-3b analysis, with ZERO runtime wiring (like 3b):

```
O2(PeakProfile, threshold) -> router overrides        // {value,use-site} -> (home, split)
seed = perfect CSE (home_scope for every value; the 3b Schedule)
while peak(placement) > threshold:
    p     = binding_point                              // argmax live-set (from the sweep)
    cands = live_at_binding                            // the cells alive at p
    move  = prefer any zero-cost SHRINK among cands;
            else argmax over cands of  DeltaPeak / DeltaRecompute
    if no move reduces peak(p): report infeasible (Section 4); break
    apply move to the router; re-sweep (Section 3d)
report (avoidable_flops, peak_bytes)                   // objective + constraint
```

O2 lives in a new component (consuming `PeakProfile`, the router, and the sizing
primitives from 3b) and is validated STANDALONE: run it on the witness forests' static
profiles and assert the resulting modeled peak drops under threshold (or reports an honest
infeasibility). No runtime path changes in 4a.

### Phase 4b -- the cutover (behavior-changing)

Wire seed+O2 into `place_at_this_level` and do the CAT-1 deletions (design doc
`2026-08-03-...phase3-design.md` section 8 lists file:line):

- Drop the `External` filter in `ext_modes_of` (`lifetime_mask.hpp:73-83`) so
  `stamp_lifetime_masks` meets ALL batched modes; `sliced_modes` is now the unified meet.
- DELETE `contracted_modes` end-to-end (field/accessors `eval_expr.hpp`;
  `NodeBatchAnnotation`; cost-model emission `cost_model.hpp:1885-1911`; the `in_union`
  union half and `residency_all_outer` contracted conjunct in `eval.hpp`;
  `schedule_dump.hpp`; tests). `residency_all_outer` stays as the structural
  "home = deepest residency loop" walk, now over the unified `sliced_modes` alone.
- DELETE `has_demoted_external` (`eval.hpp:1720-1743,1748,1784`) and the veto coupling;
  its expectations invert in the affected `test_eval_ta.cpp` cases (they now assert the
  hoisted+spilled placement, not the vetoed-local one).
- `place_at_this_level` consults O2's router overrides for the home (the Phase-2 store
  seam), falling to the plain `home_scope` (deepest residency loop) when O2 made no move.

## 3. O2 internals (YAGNI-first)

### 3a. Moves (spec 7b)
- **Shrink -- the workhorse.** A demoted cell homed ABOVE a batch loop `m` it CARRIES
  holds `m` full; re-home it INSIDE that (already-placed) loop so the loop slices it,
  dropping the footprint by `B_m` at ~zero recompute (carried-mode slicing is free per the
  `W` model). This directly undoes the seed's peak-maximal hoisting and is the C60 /
  water-20 lever. Build FIRST; prefer it (zero-cost) before any evict.
- **Un-hoist** (evict, invariant loop): lower a cell's home toward a loop it is invariant
  to, so it is rebuilt lazily per outer block instead of held idle. Footprint unchanged,
  lifetime shortened. Recompute = one rebuild per outer block.
- **Split** (evict, long live range): partition a cell's use-sites into groups (a new
  `split_index`), each its own short-lived cell at its own meet. Recompute = one extra
  build per group.

Build shrink first; if it alone brings the witnesses under their targets, un-hoist/split
are lower-priority (specced, validated on constructed forests, but not on the critical
path). None of these touch `batched_here` (which modes are batched) -- that is the
factorizer's lever, reached only via the re-batch feedback (Section 4, Phase 5).

### 3b. The greedy (O2b)
Single pass: at the binding point, prefer any zero-cost shrink, else take the move
maximizing `DeltaPeak / DeltaRecompute`. `DeltaRecompute` = the added `avoidable_flops`
(the `W` model / `cost_profile().avoidable_flops`); `DeltaPeak` = the drop in the sweep's
max. Lookahead (a shrink that unlocks a cheaper later evict) is DEFERRED -- add only if a
single pass demonstrably underperforms on the witnesses.

### 3c. Relation to the existing DP external-slice pass (O2c)
Run the existing per-term DP external-slice pass FIRST (unchanged), then O2 as the
whole-forest evict pass on TOP -- the smaller change. O2 sees the post-DP placement and
only makes the CSE-aware, cross-occurrence moves the per-term factorizer cannot.

### 3d. Incremental re-sweep (O2a) -- DEFERRED
Start with a FULL re-sweep of the 3b `Schedule` after each move (correctness-first; the
sweep already exists; moves are bounded). Add the incremental update (a move changes one
cell's interval endpoints + weight; re-sweep only the affected span) ONLY if
`moves x forest` is measurably too slow on the witnesses.

## 4. Termination and the two honest infeasibilities

O2 stops on `peak <= threshold` (success) or when no move reduces the binding point. The
latter is NOT repairable within O2 and splits into two reported outcomes (spec 7b), both
feeding Phase 5:
- **factorization-inherent** -- a single intermediate exceeds the budget even fully
  placed; the factorization itself must change.
- **re-batch-needed** -- a shared cell would need slicing on a mode the per-term
  factorizer never batched (no single term saw the sharing). Feeds back to add batch loops
  (Phase 5 / spec O6). O2 reports it; it does not add loops.

## 5. The default-threshold policy (cutover risk)

Deleting the veto means an UNTUNED run (no explicit `peak_threshold`) would get
`threshold = inf` => seed-alone => peak-maximal => a DEFAULT REGRESSION into the
C60-OOM class the veto used to prevent. So 4b MUST define a default so O2 spills by
default: `peak_threshold` defaults to a finite budget (e.g. a fraction of addressable
memory, or the existing DP/perf-first default if one exists) rather than infinity. The
exact default is a 4b decision; the REQUIREMENT is that "no explicit threshold" does NOT
mean "peak-maximal." Document the chosen default and its rationale in the cutover commit.

## 6. Consequences

- **Behavior-changing (not byte-identical).** Every placement changes. Correctness (CCk
  energies) is the gate; peak is expected to drop. Validate CCk energies UNCHANGED, then
  peak drops toward the documented-RED witness targets.
- **No fallback.** The heuristic is gone; there is no non-O2 production path. O2 must be
  robust; infeasibility is reported (Section 4), not silently fallen back.
- **The unified meet corrects a latent bug** (the per-occurrence `contracted_modes`
  under-meet of aux modes that differ across occurrences, `A[c,_]` vs `A[_,c]`), so some
  placements legitimately differ from today beyond the demotion change.

## 7. Validation

- **4a (standalone, additive):** O2 on the witness static profiles lowers the modeled peak
  under threshold; `threshold = inf` is a no-op (== seed-alone == the raw 3b profile);
  each move's `DeltaPeak`/`DeltaRecompute` matches a hand-computed expectation on
  constructed forests; infeasible forests report the right mode.
- **4b (cutover):**
  1. **Correctness first** -- the CCk validation suite energies UNCHANGED (the load-bearing
     gate; placement must not change results).
  2. `threshold = inf` reproduces EXACTLY the raw perfect-CSE seed placement (O2 no-op).
  3. Finite threshold: the witness peaks DROP toward their documented-RED targets
     (`[.][dryrun-occ-veto]` 525.9 GB -> under the <100 GB target), with energies still
     correct.
  4. The deleted-veto `test_eval_ta.cpp` cases re-baselined to the hoisted+spilled
     placement (expectations inverted deliberately, justified in the commit).

## 8. Roadmap after Phase 4

- **Phase 5 -- feedback (spec O6):** consume O2's reported infeasibilities (re-batch-needed)
  to drive the factorizer to add batch loops; reconcile the static `PeakProfile` estimate
  against measured runtime peak; tune the default threshold and move ordering on real
  workloads.

## 9. Phase 4b decomposition, and 4b-1 (the unified meet) detailed design (2026-08-04)

Section 2's Phase 4b (the cutover) is too large and too risky for one change, so it is
decomposed into three slices, ordered safest-first:

- **4b-1 -- the unified meet (this section):** drop the `External` filter so `sliced_modes`
  is the all-batched-modes meet; DELETE `contracted_modes` end-to-end; retire the now-
  duplicate `seed_residency`. KEEP the veto (`has_demoted_external`). Correctness-gated
  (CCk energies unchanged); no remat wiring. This is the first phase whose placement is NOT
  byte-identical to today's.
- **4b-2 -- hash threading + router emission (additive):** thread the value hash /
  occurrence-key through `linearize_rich` -> `ValueCell` (currently computed in
  `hash_to_cell` and discarded), and add `remat_to_router(placement) -> PlacementRouter`
  overrides. Produces overrides; not yet driving runtime.
- **4b-3 -- the coupled behavior change (its own design pass):** wire a dry-run remat
  pre-pass -> populate the router -> `place_at_this_level` consults it -> DELETE
  `has_demoted_external`; set the default-threshold policy (Section 5). The runtime-
  integration architecture (how the dry-run remat placement reaches real eval via the
  occurrence-keyed router) is designed when 4b-3 is reached, informed by 4b-1/4b-2.

### 9.1 4b-1: what changes

1. **Drop the `External` filter** in `stamp_lifetime_masks` (`lifetime_mask.hpp:148-153`):
   its `ext_modes_of` selector becomes the all-batched-modes selector (identical to
   `stamp_seed_residency`'s `all_batched_modes_of`). `sliced_modes` is now the cross-
   occurrence meet of ALL batched modes on a node's result slots.
2. **Delete `contracted_modes` end-to-end** (the CAT-1 list with file:line in
   `2026-08-03-meet-based-home-scope-phase3-design.md` section 8): field/accessors
   (`eval_expr.hpp`), `NodeBatchAnnotation`, the cost-model emission (`cost_model.hpp`),
   the `in_union` union-half + `residency_all_outer` contracted conjunct + doc
   (`eval.hpp`), `schedule_dump.hpp`, and tests. `residency_all_outer` stays as the
   structural "home = deepest residency loop" walk, now over `sliced_modes` alone. Because
   the unified `sliced_modes` already contains every batched mode on the slots, the old
   `sliced UNION contracted` is subsumed -- deleting `contracted_modes` is net-corrective
   (it also fixes the latent per-occurrence aux under-meet, e.g. `A[c,_]` vs `A[_,c]`).
3. **Retire the duplicate `seed_residency`:** delete `stamp_seed_residency` +
   `EvalExpr::seed_residency_` (+ accessors); point `home_scope` at `sliced_modes()`;
   `peak_profile.hpp`'s `linearize_rich` calls `stamp_lifetime_masks` instead of
   `stamp_seed_residency`; the 3a/4a tests that read `seed_residency` move to
   `sliced_modes`. After step 1, `sliced_modes == ` the old `seed_residency` exactly, so
   this retirement is BYTE-IDENTICAL for `peak_profile`/remat -- a behavior-neutral de-dup.

### 9.2 4b-1: what stays

`has_demoted_external` (the veto) -- giants are still un-hoisted to local, so 4b-1 does NOT
peak-regress and stays correctness-preserving (the veto's deletion is 4b-3, coupled with
remat). `residency_all_outer` (input generalized to the unified `sliced_modes`). The remat
pass (4a) is unchanged except that `home_scope` now reads the runtime field.

### 9.3 4b-1: the T1/T2 split and validation

- **T1 -- the correctness change (CCk-gated):** the filter drop + `contracted_modes`
  deletion (steps 1-2). Placement changes. Gate: CCk energies UNCHANGED -- SeQuant's
  `[eval]`/`[eval_ta]` correctness suites green (real contraction results), and the
  cross-repo MPQC CCk validation suite (user-run) is the final energy gate. The dry-run
  witnesses (`[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]`) MAY move -- this is
  the first non-byte-identical phase -- and are RE-BASELINED with justification (unified
  meet corrects the aux under-meet + drops the External/contracted asymmetry), not treated
  as a regression.
- **T2 -- the behavior-neutral de-dup:** retire `seed_residency` (step 3). Because
  `sliced_modes == seed_residency` after T1, `[peak_profile]`/`[placement_remat]`/
  `[lifetime_mask]` figures are UNCHANGED across T2 -- the de-dup is provably inert, which
  is its acceptance gate.

## 10. Phase 4b-2: hash threading + `remat_to_router` (additive), detailed design (2026-08-04)

4b-2 makes the remat placement emittable as `PlacementRouter` overrides, without wiring
anything into runtime (that is 4b-3). Additive; zero production callers of the new emitter.

### 10.1 Value hash vs occurrence key (the load-bearing distinction)

- **Value hash** = `EvalExpr::hash_value()`: the CSE / value identity ("same computed
  result"), batched-slot-BLIND (colors indices by space). This is what folds occurrences
  into ONE remat cell.
- **Occurrence key** = `occurrence_key(node, ctx_modes)` -> `SlotCanonicalizationMetadata`
  (hashed via `SubNetHash`, compared via bliss): the batched-slot-AWARE identity -- the
  Phase-2 `canonicalize_slots` coloring that colors in-scope batched indices DISTINCTLY, so
  `A[i1,i2]` and `A[i1,_]` get DIFFERENT keys. A strictly FINER identity than the value hash;
  effectively "the hash of an occurrence."

The `PlacementRouter` is keyed by the occurrence KEY (per Phase 2). A remat cell is
one-per-VALUE (hash). So a single moved cell emits MULTIPLE overrides -- one per DISTINCT
occurrence key of that value -- all pointing to the SAME `HomeTarget` (the meet home). This
is exactly the partial-overlap case of the Phase-3 negative result (spec `2026-08-03` §3):
`A[i1,i2]` / `A[i1,_]` are two keys sharing the meet home `{i1}`.

### 10.2 The change

1. **`ValueCell` gains `std::size_t hash`** -- `compute_dag_boulevard` already computes it in
   `hash_to_cell` and discards it; just store it. One field, additive; `peak_profile` /
   remat ignore it (byte-identical to them). This is the link from a `RematResult` cell back
   to its forest nodes. We thread the value HASH (cheap, batched-slot-blind), NOT the
   occurrence keys (expensive, computed lazily below).

2. **`remat_to_router(seed_cells, remat_cells, forest) -> PlacementRouter`** (new, in
   `placement_remat.hpp`):
   - Diff `seed_cells` (home = seed `home_scope`) vs `remat_cells` (home = post-spill) by
     `value_id` -> the MOVED set: `{hash -> final home_modes : home_modes != seed home_modes}`.
     Moved-cells-only (unmoved cells' seed home == the runtime's default derivation, so an
     override would be redundant; matches Phase-2's empty-router-is-inert seam).
   - Re-walk the forest (post-order + the same ectx accumulation as `compute_dag_path`); for
     each node whose `hash_value()` is in the MOVED set, compute `occurrence_key(node,
     ectx_modes)` and `router.set_override(key, HomeTarget{final_home_modes, 0})`. The bliss
     canonicalization runs ONLY for moved-cell occurrences -- `compute_dag_boulevard` /
     `peak_profile` are untouched (no per-node bliss cost added to the shared static analysis).
   - Return the populated router.

The runtime depth resolution is free: `HomeTarget` stores the `residency` MODE-SET
(`placement_router.hpp`), and `home_depth(HomeTarget, ectx)` does the rl-walk per live batch
context -- so nothing occurrence-relative is baked into the override (the 3a/3b carry-forward
is satisfied by the router's own design).

### 10.3 Validation

- `remat_to_router` on a constructed forest with a KNOWN moved cell (a shrunk demoted giant)
  emits exactly the expected `set_override`s: one per distinct occurrence key of the moved
  value, each -> `HomeTarget` with the shrunk `home_modes`; unmoved cells emit NONE.
- A value with two DISTINCT-occurrence-key occurrences sharing one (moved) home emits TWO
  overrides -> the SAME `HomeTarget` (the partial-overlap / meet-home case).
- `threshold = inf` (no moves) -> the returned router is EMPTY (nothing moved), i.e. the
  Phase-2 inert seed seam.
- Additive: `[peak_profile]` / existing `[placement_remat]` figures UNCHANGED (the `hash`
  field and the new emitter do not perturb the sizing / sweep).
