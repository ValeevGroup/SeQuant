# Phase 4: threshold-driven O2 spill + the CAT-1 cutover

Status: design, agreed via brainstorming 2026-08-04. Consumes Phase 1 (peak oracle),
Phase 2 (router read+store seam), Phase 3a (`home_scope` seed predictor), Phase 3b
(static `PeakProfile`), all complete on `evaleev/feature/multimode-batched-eval`
(`...f31bda794`). Realizes O2 (spec `2026-08-02-home-scope-placement-design.md` section
7b) and the CUTOVER. **Supersedes** the Phase-3 design doc's roadmap item "Phase 4 --
O2 greedy + CUTOVER ... lands COUPLED behind a flag": there is NO flag; `peak_threshold`
is the sole control (see Section 1).

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
