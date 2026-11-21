# Perf-first / peak-second single-term optimization objective

**Goal:** Add a perf-first, peak-second selection policy to the single-term
optimizer's peak cost model so that batching-optional factorization *rejects*
FLOPS-catastrophic factorizations (e.g. the 4-PAO AO-basis integral) instead of
*preferring* them for their fully-sliceable peak. Keep the existing
peak-first objectives unchanged (non-breaking).

**Architecture:** Reuse the existing Pareto-frontier + roofline-cache-hybrid
machinery verbatim. Change *only* the two frontier *selectors* to sort by
`(roofline-perf, then peak)` instead of `(peak, then flops)`. `peak_threshold`
shifts from a primary hard feasibility filter to a secondary "slice-down-to"
target. Exposed as a new `ObjectiveFunction` value pair so current behavior is
preserved.

**Tech stack:** SeQuant C++20, `SeQuant/core/optimize/`.

---

## Problem and motivation

The current peak-first objective (`DensePeakSizeBatched`, renamed
`DenseSpaceTimeBatched` below) selects a schedule by `select_root`
(`cost_model.hpp:692`): among root-frontier points with `peak <= peak_threshold`
pick minimum flops (ties by lower peak); if none fit, minimum peak. The
`peak_threshold` is thus a **hard filter over factorizations**, and FLOPS only
breaks ties *among survivors*.

For the CSV-CC particle-particle-ladder (PPL) term this is pathological. With
PAO (mu-tilde) made sliceable, two factorizations compete at
`peak_threshold = 40 GB` (C60, occ=120, PNO=42, mu-tilde=1800, aux=4320):

- **Correct `(gC)(gC)`:** the giant `I(i,i,mu~,K;a<i,i>)` has an *irreducible*
  floor -- its PNO-pair leg `a<i,i>` is not sliceable, so even slicing both
  mu-tilde and aux leaves ~89 GB > 40 GB. **Filtered out.**
- **`g.g -> 4-PAO` `(mu~ mu~|mu~ mu~)`:** every axis is a sliceable PAO;
  slicing all four drives peak to ~34 GB < 40 GB. **Survives.** Its FLOPS
  (~mu-tilde^4 . aux) are astronomically higher, but they never get a vote
  because it is the sole feasible point.

So peak-first *structurally prefers* forming the AO 4-index integral -- the one
thing DF+PNO exists to avoid.

### Evidence (DryRun harness, `tests/unit/test_eval_dryrun.cpp`)

- **aux+PAO sliceable:** the DP forms the 4-PAO monster (`mu~^4`, nominal
  83980 GB) and slices all four axes to 34 GB. Faithful to the objective, not a
  harness artifact -- the cost model is the same one MPQC uses.
- **aux-only (`DRYRUN_AUX_ONLY=1`, aux_ts=35):** no 4-PAO appears, but the
  correct factorization's giants stay mu-tilde-full: largest realized
  `{mu~ mu~ i i} = 373 GB` (carries no aux -> *immune* to aux batching at any
  tile size, explaining why the cluster's `aux_target_size=35` run still OOM'd),
  plus the proper giant `{mu~ a<i,i> K} = 305 GB`. Aux-only avoids the
  pathology but cannot reduce the mu-tilde-full giants -> OOM.

Peak-first is therefore broken both ways: aux-only OOMs on mu-tilde-full giants;
aux+PAO "fixes" peak only by choosing the FLOPS-catastrophic 4-PAO.

### Why the swap fixes it, and why it is safe

Two facts make a pure order-swap correct:

1. **Factorization sets FLOPS; slicing sets PEAK.** The correct `(gC)^2` is
   min-flops; the 4-PAO is max-flops. A perf-first primary rejects the 4-PAO at
   the factorization level, before slicing is even considered.
2. **Slicing is perf-neutral in this model.** In `relax()`
   (`cost_model.hpp:415, 633`) the per-node roofline cost `cflops` is computed
   **once, before the batch/slice loops**, on the **unsliced** footprints
   (`ctx.sz(*, 0)`) and full `flops_of`, then added identically to every
   `(B, Ap)` frontier point. So all slicings of a factorization share one
   `flops`(=roofline) coordinate and differ only in `peak`. A perf-first primary
   therefore leaves slicing entirely to the peak-second selector -- it never
   under-slices. (This retires the earlier worry that roofline would penalize
   slicing's re-reads; the model deliberately does not charge per-batch
   re-read traffic -- see the `relax()` comment "slicing reduces peak, not total
   work".)

`roofline_op_cost` (`cost_model.hpp:281`) is
`max(flops, beta * max(traffic, kappa * flops / sqrt(M/c0)))`: compute-bound ops
reduce to `flops`, bandwidth-bound ops to `beta*Q`. This is retained as the
perf primary precisely because raw FLOPS cannot distinguish the bandwidth-bound
tensor products that occur in PNO-CC; it is exactly the metric that should rank
factorizations.

---

## Design

### Objective naming (`Dense{Primary}{Secondary}`)

The objective names encode the lexicographic optimization order directly:
`Space` = peak/size, `Time` = perf (roofline), and the name lists primary then
secondary. This gives a clean symmetric pair in `ObjectiveFunction`
(`options.hpp:46`):

- `DenseSpaceTime` / `DenseSpaceTimeBatched` -- peak-first, perf-second
  (the current behavior). **Rename** of `DensePeakSize` / `DensePeakSizeBatched`,
  which become deprecated aliases (kept so existing inputs and pinned MPQC do
  not break).
- `DenseTimeSpace` / `DenseTimeSpaceBatched` -- perf-first, peak-second
  (**new**). This is the objective that fixes the C60 pathology.

(`DensePeakSize` alone was ambiguous -- "Peak" can read as either peak *memory*
or peak *performance*; `Space`/`Time` are unambiguous.)

The aliases are implemented as enumerators that share the underlying value of
the renamed constant (`DensePeakSize = DenseSpaceTime`,
`DensePeakSizeBatched = DenseSpaceTimeBatched`), placed AFTER the four primary
enumerators so `DenseSpaceTime` keeps the old numeric value of `DensePeakSize`
(ABI-stable). Because the alias and its target compare equal, every existing
`Metric == DensePeakSize` guard (in `single_term.hpp` and `optimize.cpp`) keeps
working unchanged and catches `DenseSpaceTime`; only the two NEW values
(`DenseTimeSpace{,Batched}`) need fresh arms.

`DenseSpaceTime{,Batched}` keep their current semantics byte-for-byte. The
auto-route (a `DenseTimeSpace` request with an active `is_batchable_index`
becomes `DenseTimeSpaceBatched`) mirrors the existing route in MPQC's
`SeQuantEngine::make_optimize_options` (`sequant_engine.cpp:182-195`).

### Why the selector is a plain lexicographic swap

The Pareto frontier already carries everything the selection needs, and its
structure makes the perf-first rule trivial. `pareto_insert`
(`cost_model.hpp:250`) drops any point dominated in BOTH `(peak, flops)`, so for
two points with the SAME flops the higher-peak one is pruned: the frontier
retains exactly one min-peak point per distinct flops value. Sorted by peak
ascending, flops is therefore strictly descending -- a clean trade-off curve.

Two consequences:

1. **The min-flops frontier point is already the fully-sliced (min-peak)
   realization of the cheapest factorization.** There is no separate "slice
   more" choice to make at selection time -- the DP's `relax()` explored every
   slicing (`Ap` subset) and the frontier collapsed each flops class to its
   lowest peak. So an earlier draft's "among perf-ties, pick least-sliced under
   threshold" steps are unreachable: the frontier never holds two perf-tied
   points at different peaks.
2. **`argmin(flops, then peak)` picks the cheapest factorization AND gets it
   maximally sliced.** For C60 that is `(gC)(gC)` (moderate flops), never the
   4-PAO (flops ~ mu-tilde^4 . aux); and the chosen `(gC)(gC)` point is its
   own min-peak (fully sliced on mu-tilde and aux). Both goals in one key.

So the whole behavioral change is: both `PeakModel` (`cost_model.hpp:303`) and
`PeakBatchedModel` (`cost_model.hpp:486`) gain a `bool perf_first = false` field
(set true for the new objectives), and the selectors sort by `(flops, then
peak)` instead of `(peak, then flops)` when it is set. The Pareto frontier
construction (`pareto_insert`, `relax`, the DP recurrence) is **unchanged**.

**`PeakBatchedModel::select_root` (`cost_model.hpp:692`)** -- perf-first branch:
plain `argmin(flops, ties by peak)` over the root `B=0` frontier `st[root][0]`
(the flops-primary analog of `pareto_best`). `peak_threshold` is NOT consulted:
it can no longer act as a hard filter that forces a FLOPS-catastrophic
factorization for its sliceability (the exact peak-first pathology). The DP
still explores slicing unconditionally, so the chosen point is already
peak-minimized. `peak_threshold` remains a documented knob for the peak-first
objective and is emitted as a diagnostic under perf-first when the selected
peak exceeds it (informational only, no behavior change).

**`PeakModel::reconstruct` (`cost_model.hpp:437-453`)** -- perf-first branch:
plain `argmin(flops, ties by peak)` over the root frontier, bypassing the
`peak_flops_tolerance` epsilon relaxation (a peak-first knob). The non-batched
model has no slicing and no `peak_threshold`.

(`pareto_best` at `cost_model.hpp:263` is used only by the test helper
`reconstructed_batched_peak`, which reconstructs the peak *value* and must stay
peak-min; it is NOT a production factorization selector and is left unchanged.)

### Data flow

`single_term_opt` (`single_term.hpp:65`) constructs the model from the
compile-time `Metric` template parameter via an `if constexpr` chain. It already
selects `PeakModel` vs `PeakBatchedModel` and threads `peak_threshold`,
`volatile_weight`, `roofline`; the two peak arms are broadened to also match the
new `DenseTimeSpace` / `DenseTimeSpaceBatched` values, and each sets
`model.perf_first = (Metric == DenseTimeSpace[Batched])` on the constructed
model. `optimize.cpp:103-125` gains dispatch arms for the two new enum values
(the batched one also routes `out_axes`, mirroring `DenseSpaceTimeBatched`). No
other call site changes: `reconstruct` / `reconstruct_axes` already route
through `select_root`, and the batch-axis annotation path
(`term_batch_axes -> node_batch_axes -> EvalExpr::batch_axes_`) is untouched.

### What does NOT change

- Frontier construction and dominance (`pareto_insert`, `relax`, `sliced_footprints`).
- The roofline cost function, `volatile_weight`, `accumulation_factor`.
- Runtime batched evaluator (`eval.hpp`) -- it consumes the same annotations.
- `DenseSpaceTime`/`DenseSpaceTimeBatched` (formerly `DensePeakSize*`) behavior
  and all their tests -- only the enum name changes, semantics are identical.

---

## MPQC wiring

`SeQuantEngine` (`sequant_engine.{h,cpp}`): extend the
`optimize:objective_function` keyword parser to accept the new canonical strings
`dense_space_time` (-> `DenseSpaceTime`) and `dense_time_space`
(-> `DenseTimeSpace`) alongside `dense_flops` / `dense_size`, and keep
`dense_peak_size` as a deprecated alias for `dense_space_time` so pinned/existing
inputs still parse. Add the auto-route to `DenseTimeSpaceBatched` when
`is_batchable_index` is set (mirroring the existing peak-first route at
`sequant_engine.cpp:182-195`). No change to the `batch:*` keywords:
`peak_threshold`, `aux_target_size`, `pao_target_size` retain their meaning and
values.

CCk input migration is a one-line objective change per input
(`dense_peak_size -> dense_time_space`); no numeric retuning of `peak_threshold`
is required.

---

## Validation (DryRun harness, before any cluster run)

The harness already replays a schedule's predicted sizes with the real cost
model. Add a perf-first path (drive `run_single_term_opt` with the new
objective) and assert, on the real post-transform C60 residual at faithful
config (occ=120, PNO=42, OSV=310, mu-tilde=1800, aux=4320, pao_ts=256,
aux_ts=72, peak_threshold=40 GB):

1. **No 4-PAO:** no `mu~^4` (four-free-PAO) node appears in any summand's
   optimized tree.
2. **Correct factorization:** the PPL giant is the 1-free-PAO
   `I(i,i,mu~,K;a<i,i>)` form (matches the real trace).
3. **Giants slice under threshold:** the largest realized free-mu-tilde
   intermediate is `<= peak_threshold` (the 373 GB `{mu~ mu~ i i}` giant slices
   both mu-tilde to ~7.6 GB; the 305 GB giant slices to <= 40 GB).
4. **A/B vs peak-first:** the same residual under `DenseSpaceTimeBatched`
   reproduces the 4-PAO (regression guard that the swap is what changed the
   outcome).

Only after the harness confirms 1-4 do we consider a C60 cluster run.

---

## Testing

- **Unit (`tests/unit/test_optimize.cpp`):** a minimal PPL-shaped term
  (`g.g` reducible to a 4-PAO vs `(gC)(gC)`) with a small extent regime and a
  finite `peak_threshold`: assert `DenseSpaceTimeBatched` selects the 4-PAO and
  `DenseTimeSpaceBatched` selects `(gC)(gC)`. A second case with
  `peak_threshold=+inf` asserts `DenseTimeSpace` reduces to min-perf-then-peak.
- **Unit (non-batched):** assert `pareto_best`'s perf-first branch returns the
  min-flops (not min-peak) frontier point on a hand-built frontier.
- **Unit (alias):** assert the `dense_peak_size` keyword still parses to
  `DenseSpaceTime` (backward-compat guard).
- **Harness (`test_eval_dryrun.cpp`):** the four C60 assertions above.
- **Regression:** existing `DenseSpaceTime*` (formerly `DensePeakSize*`) tests
  must be byte-unchanged (the new field defaults `perf_first=false`).

---

## Global constraints

- No en-dashes (U+2013) or non-breaking spaces (U+00A0); ASCII hyphens only
  (pre-commit hook rejects otherwise).
- No `Co-Authored-By` trailers in commit messages.
- Cross-repo: this is a SeQuant change; MPQC consumes it via the pin in
  `external/versions.cmake` (`MPQC_TRACKED_SEQUANT_TAG`). The SeQuant change and
  the MPQC wiring land as separate commits/PRs; MPQC repins only after the
  SeQuant objective merges.
- `DenseSpaceTime` / `DenseSpaceTimeBatched` (formerly `DensePeakSize*`)
  semantics and outputs must not change; the old names remain as deprecated
  aliases.
