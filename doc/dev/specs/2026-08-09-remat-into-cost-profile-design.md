# Retiring the replay hold-gate: accurate peak by measurement, budget by placement

_2026-08-09. Branch context: `evaleev/feature/multimode-batched-eval` (SeQuant), paired MPQC batched-CCk work. Depends on the DAG-scope home + replication-priced split just landed (Tasks 3-6 of `2026-08-07-remat-cse-aware-split.md`)._

## Two budgets, not one (the whole source of confusion)

Two distinct, similarly-named knobs govern space, at two different stages. The
C60 witnesses set BOTH to 100 GB, so they read as one thing. They are not.

- **`policy.peak_threshold` -- the DP's batching budget (KEEP).** The DP
  (`optimize`, objective `DenseTimeSpaceBatched`) chooses batching-loop LOCATIONS
  with it. `cost_model.hpp:1727-1730`: walk the chosen schedule tree; at each node
  whose `subtree_peak` exceeds `peak_threshold`, greedily slice a carried external
  mode until the node fits. `subtree_peak` is a **per-node (subtree) proxy** for
  peak -- cheap (one number per node), evaluated across MANY DP candidates, and
  deliberately NOT the global co-resident peak. This is exactly the "crude but
  cheap" per-node budget: it is how the batched schedule is generated, and it must
  stay cheap (see Cost). Fed by `opts.batch_policy.peak_threshold`
  (`optimize.cpp:106`); default `+inf` = no batching.
- **`cfg.max_footprint` -- the replay hold-gate (RETIRE).** Lives in
  `CacheConfig`, used ONLY by `build_dryrun_cache` in the dry-run `cost_profile`.
  It decides which values the REPLAY holds vs recomputes. `grep max_footprint
  core/optimize/` returns nothing -- the DP never sees it. It is a **per-value**
  threshold (not peak-aware) applied AFTER the schedule exists, and it is the
  wrong instrument: per-value not peak, and a SECOND budget that re-judges a
  schedule the DP already shaped for `peak_threshold`.

**Evidence they are independent axes** (from this session's sweeps on
`[.][dryrun-occ-veto]`): sweeping `peak_threshold` 100 -> 40000 GB left unbatched
avoidable flat at ~60%; sweeping `max_footprint` 100 GB -> 10 TB moved it 60% ->
0%. The "60% avoidable at a 40 TB `peak_threshold`" was entirely `max_footprint`
(100 GB) evicting the ~39 TB unbatched giants -- the batching budget never
touched it.

## Summary

Retiring `cfg.max_footprint` splits cleanly into a cheap default and an optional
enforcement step:

- **Prediction (cheap, no remat).** To fix the confound -- report the TRUE peak
  and avoidable of the DP's schedule -- just disable the replay gate
  (`max_footprint = 0`, already the library default) and let the replay MEASURE
  the realized peak. The gate was distorting the *measurement*, not the schedule;
  removing it costs one replay (already paid). If the DP's per-node proxy
  under-counted and the real co-resident peak exceeds B, the measurement now
  honestly says so.
- **Enforcement (optional, remat once).** To make a schedule that actually FITS
  B, run the remat placement pass ONCE on the final fused DAG (home-scope / split
  to fit B), then measure. This is opt-in and is already paid in production: the
  MPQC pre-pass already runs remat once to populate the router.

`policy.peak_threshold` stays throughout -- the DP still generates batching the
same way. The only thing retired is the replay-side per-value gate.

## Background: the current pipeline

`cost_profile(forest, policy, cfg, regime, trace, router, schedule_sink)`
(`cost_profile.hpp`):

1. static per-node cost walk (`model_flops` etc.), placement-blind;
2. `build_dryrun_cache(forest, cfg, regime)` -- replay cache; its lifetime mask +
   the `cfg.max_footprint` gate decide which values persist;
3. `cache.set_placement_router(router)` -- an OPTIONAL router (caller-supplied)
   decides where (scope) a value is homed; null by default;
4. replays each root through `evaluate<Trace::On>` with `make_evaluator(policy)`
   (batched-scratch slicing), tallying recompute-aware `dryrun_flops` /
   `avoidable_flops` and the batched-scratch high-watermark `peak`.

Batching is chosen upstream by `optimize(peak_threshold)` before binarization. The
router is accepted but no dry-run caller populates it from remat; the MPQC
pre-pass does, and the split unit test hand-builds one.

## Design

### A. Confound fix (the cheap default): measure with the gate off

Set `cfg.max_footprint = 0` (its default) and MEASURE. No remat, no router. The
replay executes the DP's schedule faithfully and reports its true realized peak +
avoidable. This alone removes the blunt instrument from the *prediction* path.
Perfect-CSE at large effective memory holds everything and prices to exactly 0
avoidable (already verified: the 1.58e16 floor at `max_footprint` >= ~10 TB).

Caveat this makes explicit: with the gate off and NO placement, the DP's schedule
is measured as-is -- its realized peak may exceed the `peak_threshold` the DP
targeted (per-node `subtree_peak` != global peak). That is the accurate answer
("this schedule needs X GB"), not a bound.

### B. Enforcement (opt-in): remat once on the final DAG

To produce a schedule that FITS B, add a stage between the forest and the replay,
run ONCE:

1. `RematInput in = remat_cells(forest, *cm, block_of)` -- the `ValueCell`s
   (`compute_dag_boulevard`). `block_of` = per-mode block SIZE from
   `policy.batch_target_size` (the same partition the evaluator slices on).
2. `RematResult res = rematerialize_to_budget(in.cells, *cm, block_of,
   in.num_points, B)` -- shrink/split the seed placement until the modelled peak
   fits B (or it reports `RebatchNeeded` / `FactorizationInherent`).
3. `auto router = remat_to_router(in.cells /*seed*/, res.cells /*final*/,
   forest)` -- one DAG-global overlay per moved value.
4. `cost_profile(..., cfg{max_footprint = 0}, &router)` -- replay the placed
   schedule and measure. `res.modeled_recompute` is the placement's forecast
   recompute (report-only; objective stays peak-only, per
   `2026-08-07-remat-cse-aware-split-design.md` section 3).

The recompute the replay then tallies is the recompute the *schedule* incurred
(giants split/re-homed to fit B), not a per-value gate artifact.

## Cost (the reason for the prediction/enforcement split)

The DP explores MANY candidate factorizations/batchings PER TREE, each needing a
space estimate. **Remat must never run inside that loop** -- a greedy re-sweep per
candidate is `O(candidates x remat)`, combinatorial death. So:

- The **DP keeps its cheap per-node `subtree_peak` vs `peak_threshold`** check
  (section "Two budgets"). It is evaluated per-candidate and MUST stay `O(1)` per
  node. Making the DP itself peak-accurate (the accurate model in its inner loop)
  is the thing this ban rules out and is explicitly out of scope.
- **Remat runs ONCE, on the final fused DAG** (path B), not per-candidate. Its
  cost is one greedy re-sweep (`~O(cells^2 x moves)` worst case). This is
  acceptable because (a) it is ALREADY the MPQC production pre-pass cost -- not new
  for a real run; (b) `cost_profile`'s dominant cost is the DP itself
  (`optimize` per term -- the occ-veto scans were minutes, almost all DP), so
  remat-once is very likely small beside it; (c) it is a prediction tool, not a
  hot loop.
- **Make path B opt-in.** Exploratory scans use path A (cheap measure-only); only
  the final "does C60 fit B?" answer pays for remat. Profile remat-once before
  optimizing it.

Residual limitation (honest): remat refines PLACEMENT (home-scope, split) WITHIN
the batching the DP chose -- a shrink block-sizes an existing batch mode, it does
not add or relocate batch loops. So batching-loop LOCATION stays the DP's crude
`peak_threshold` proxy. Fusing the accurate peak model into the DP is a separate,
larger question, out of scope.

## The crux risk (path B only): modelled peak must equal replayed peak

`rematerialize_to_budget` fits B against the **remat peak model**
(`peak_profile_sweep` over cells' liveness intervals + `cell_footprint`). B is
only truly enforced if that model **agrees with the replay's realized `peak`**
(the batched-scratch high-watermark). Two different computations -> they can
diverge, exactly the gap the `2026-08-05-dryrun-wetrun-schedule-equivalence` work
closed between dry and wet. This is path B's main validation:

- **Equivalence check.** For each C60 mode and swept B, assert `replayed_peak <=
  B` (or a bounded, understood slack). Where it exceeds, diff the remat cell
  schedule against the replay's `SCHEDULE_RUN_EVENT` dump (build/fetch per node +
  ctx) -- the same instrument the equivalence work used.
- **Likely divergence sources:** (a) `cell_footprint` sizes at block extent per
  home mode, but the replay's realized slice may differ by tile granularity (the
  16-node aux-loop gap that work found); (b) the sweep prices co-residency by
  static liveness intervals, the replay's `peak` folds batched-scratch
  high-watermarks -- the two residency models must match; (c) the replication
  factor (`ceil(extent/block)`) must match how many times the replay actually
  rebuilds an invariant cell.
- **If they cannot be reconciled cheaply:** keep `max_footprint` as an OPTIONAL
  backstop (not the primary mechanism) and document the residual slack, rather
  than claim a hard bound the replay does not honor.

## Open questions (resolve during implementation)

- **`block_of` source.** Confirm `policy.batch_target_size` is the same partition
  `cell_footprint` / the sweep / the evaluator all use, so remat peak and replay
  peak price the same slices.
- **`build_dryrun_cache` lifetime mask vs the router.** With `max_footprint = 0`
  the footprint gate is off, but the cache still has a lifetime mask (min_repeats,
  cross-occurrence lifetimes). Verify the mask + the router compose (router = scope,
  mask = liveness) with no second implicit budget.
- **Nominal (DP) vs realized peak.** The DP's `subtree_peak` is nominal; the
  realized replay peak was 50 TB at a 40 TB `peak_threshold`. Under path B, remat
  works the realized-side model, so it may need to shrink below the DP's nominal
  target to hit the realized B -- quantify.
- **API shape.** A budget-aware entry (e.g. `cost_profile_to_budget(forest,
  policy, regime, B, ...)`) that runs the path-B remat stage then calls the
  existing `cost_profile` with the derived router + `max_footprint = 0`, leaving
  the low-level `cost_profile` (explicit router / gate) intact for the equivalence
  tests and the MPQC pre-pass.

## Scope boundary

**In scope:** disabling the replay gate on the prediction path (A); the opt-in
remat-once enforcement path (B) and its modelled-vs-replayed peak equivalence
validation; retiring the hardcoded `max_footprint = 1e11` from the C60 witnesses
in favour of a swept single budget B; the corrected C60 peak/avoidable scan as
the committed witness.

**Out of scope:** an accurate peak model inside the DP's inner loop (the cost ban
above); the scope-oriented executor (cross-tree sub-scope CSE), still future work;
fusing `optimize` and remat into one pass; MPQC-side plumbing beyond confirming
the dry-run prediction matches the pre-pass path.

## Testing

1. **Confound gone (path A).** With `max_footprint = 0`, the reported peak +
   avoidable are the schedule's true realized values; results are INDEPENDENT of
   any `max_footprint >= realized peak` -- the direct proof the blunt instrument no
   longer governs the prediction.
2. **Perfect-CSE floor.** At effectively unlimited memory, unbatched avoidable ==
   0 and `dryrun_flops` == the 1.58e16 floor.
3. **Budget enforcement (path B).** For the C60 residual, sweep B; assert the
   REPLAYED `peak <= B` (or documented slack) per batching mode, with the remat
   router driving placement.
4. **Recompute is the schedule's.** `res.modeled_recompute` and the replay's
   `avoidable_flops` agree to within the peak-model slack.
5. **Empty-budget regression.** Existing `cost_profile` callers (no budget, no
   router) stay byte-identical; both paths are opt-in.

## Reproduce / evidence

- Two-knob independence: `peak_threshold` sweep 100 -> 40000 GB leaves unbatched
  avoidable ~60%; `max_footprint` sweep 100 GB -> 10 TB moves it 60% -> 0%.
  `[.][dryrun-occ-veto]` with `SEQUANT_UT_DRYRUN_PEAK_THR_GB` and a `max_footprint`
  env override.
- Perfect-CSE floor: unbatched at `max_footprint` >= ~10 TB -> 0 avoidable,
  `total_flops` = 1.58e16.

## Path B, Task 1 crux probe -- FINDING (2026-08-09)

`[.][dryrun-remat-equiv]` on the C60 aux+occ forest (55 terms) measured the remat
MODELLED peak (`RematResult::profile.peak_bytes`) against the replay REALIZED peak
(`cost_profile(..., &router).peak_bytes`), sweeping B:

| B (GB) | modelled (GB) | realized (GB) | ratio real/mod | status |
|---|---|---|---|---|
| 2000 | 1976.2 | 563.40 | 0.285 | Feasible |
| 1000 |  894.7 | 563.40 | 0.630 | Feasible |
|  500 |  612.4 | 565.03 | 0.923 | RebatchNeeded |

**The models diverge, and the realized peak is essentially INVARIANT to the remat
placement.** The remat model shrinks with B (1976 -> 895 -> 612), but the replay
lands on the byte-identical 563.40 GB at B=2000 and B=1000 -- exactly the occ-veto
gate-off floor -- despite remat producing DIFFERENT placements at those two
budgets. At B=500 remat reports `RebatchNeeded` (its model cannot fit 500 GB)
while the replay actually realizes 565 GB: a FALSE infeasibility.

**Interpretation (two readings, not yet disambiguated).** (a) The remat peak model
(`cell_footprint` co-residency over static liveness intervals) systematically
OVERCOUNTS the replay's batched-scratch high-watermark by up to ~3.5x -- a model
bug to reconcile (Task 1b). (b) On C60 the realized peak is set by the DP's
batched-scratch schedule (~563 GB), and home-scope/split PLACEMENT -- the only
lever remat turns -- cannot move it; remat is then neither necessary nor
sufficient to enforce B, and the DP's `peak_threshold` (batching) is the real
budget lever. The invariance of the realized peak across B=2000/1000 favors (b),
but (a) and (b) can both hold; a `SCHEDULE_RUN_EVENT` diff is needed to apportion.

**Non-uniform tiling caveat (per E.V.).** The probe's `block_of` is a UNIFORM block
size per space (Κ->256, occ->8). Real model tiling is NON-uniform (TA tile
extents vary across a space), so a uniform-block peak model cannot match a
non-uniformly-tiled realized peak exactly even after reconciling (a)/(b). This
bounds how tight any remat-model budget can ever be, and is an independent reason
NOT to claim a hard `peak <= B` from the uniform-block model.

**Consequence for path B.** As written (remat placement enforces a hard `peak <=
B`), path B does NOT hold on C60: placement does not move the realized peak, and
the model that drives it is both pessimistic (reading a) and uniform-block
approximate (tiling caveat). The honest deliverable is path A (shipped: measure
the DP-governed true peak, gate off) plus documenting that the C60 realized peak
is DP-batching-governed. Whether to pursue Task 1b (reconcile the model, keep
remat as a conservative-but-non-binding proxy) or treat placement enforcement as
out of scope for C60 is the open decision.
