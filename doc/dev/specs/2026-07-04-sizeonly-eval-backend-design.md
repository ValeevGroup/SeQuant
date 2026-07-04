# Size-only ("Debug") eval backend and PNO-CCSD size harness

_Design spec. Repo: SeQuant. Related: the aux+PAO nested batching feature on
`feature/eval-predicted-peak-trace` (SeQuant) / mpqc4 CSV-CCk._

## Goal

A size-only SeQuant eval backend plus a Catch2 harness that factorizes and
"runs" the PNO-CCSD (CSV-CCk) residual schedule against extent-tracking
tensors (no numeric data, no TiledArray, no MPQC), so we can reproduce
arbitrary size regimes -- including the C60 out-of-memory -- in seconds as a
unit test instead of a multi-hour cluster build-run-debug loop.

## Motivation

The C60 PNO-CCSD run (job 614336) was OUT_OF_MEMORY: aggregate MaxRSS 770 GB
vs a 743 GB node limit. It died building the intermediate

    I(i_3, i_2, mu-tilde_1241, Kappa_2; a_3) = C(...;mu-tilde_1242) * g(mu-tilde_1241, mu-tilde_1242, Kappa_2)

which the trace predicted at ~186 GB. This intermediate carries BOTH a free
PAO index (mu-tilde_1241) and a free aux index (Kappa_2). The batched
evaluator sliced only aux ("over 60 aux batches") and never sliced mu-tilde,
so the intermediate materialized full. This is the exact free-mu-tilde leak
the aux+PAO nested batching feature was meant to close, yet at
`peak_threshold = 40 GB` it did single-mode aux batching only.

Diagnosing this on the cluster is intractable: one build + queue + run to the
CCSD residual is hours. We need a local reproducer that exercises the SAME
factorizer and batched-evaluator code with configurable index sizes.

## Non-goals (YAGNI)

- No numeric correctness (the backend holds no data; exactness is already
  covered by the TiledArray-backed `[eval]` tests).
- No per-pair proto-index instantiation and no real size tables exported from
  MPQC (defeats the arbitrary-regime goal; a moment model suffices).
- No standalone CLI driver in the first cut (a Catch2 case with an editable
  regime is a seconds-long loop). The harness core is a reusable function, so
  a CLI wrapper is a trivial later add if wanted.

## Architecture

Reuse the entire existing pipeline unchanged:

    deserialize -> single_term_opt (peak model + threshold) -> binarize
                -> make_batched_custom_evaluator -> evaluate

Swap in a size-only backend as a fourth backend alongside
`tiledarray`/`btas`/`tapp`. Two facts make this cheap and faithful:

1. The eval trace and the peak high-water mark are already driven by
   `Result::size_in_bytes()`: `CacheManager::note_working_set()` maintains
   `working_set_hwmark_`, printed as `hw=` and documented as the running
   maximum peak. So running the size backend with tracing on reproduces the
   C60 log's `Eval | ... | result=...B | hw=...B` and `BatchGroup` lines
   exactly, instantly, and the realized peak is the final hwmark -- for free.
2. The optimizer's cost model already consumes an `inner_pow(composite, k)`
   moment callback (k-th power mean of a CSV composite sharing its proto set
   with k composites), and a plain index-extent callback. The harness feeds
   both from one config, so the DP-modeled peak and the runtime-realized peak
   are computed from the same size model.

Because the leak is a mismatch in EXTENTS (did slicing happen), not in the
size FORMULA (both moment-aware), the harness faithfully reports whatever the
batched evaluator actually produced -- leak included.

## Components

### Component 1: size-only `Result` backend

Location: `SeQuant/core/eval/backends/sizeonly/{result.hpp, eval_expr.hpp}`
(mirroring the existing `backends/tiledarray` layout).

Two Result types implementing the abstract `Result` interface:

- `ResultSize` (flat / Tensor-of-scalars): carries an ordered list of
  `Index` objects (with proto-index info) and their per-mode extents. No
  buffer.
- `ResultSizeOfSize` (nested / Tensor-of-Tensor, for CSV amplitudes and
  coefficients): carries the outer index list (occ proto indices) and the
  inner index list (PNO/OSV composites).

Virtuals implemented (pure virtuals from `result.hpp`):

- `size_in_bytes()`: moment-aware. Partitions the Result's proto-indexed
  modes by shared proto-index set and sizes each k-group with
  `csv_moment(rank, k)`; non-proto modes contribute their scalar extent;
  multiply by `numeric_size`. Reuses the cost model's `inner_aware_volume` /
  `memsize_counter` helper so the byte formula is IDENTICAL to the DP's, not
  a re-derivation.
- `prod(other, target, DeNest)`: computes the result index set (union minus
  contracted) and its extents; honors `DeNest` for ToT -> T. Returns a new
  `ResultSize`/`ResultSizeOfSize`.
- `sum(other, ...)`: same index set, size-preserving.
- `slice_mode(mode, lo, hi)`: returns a copy with that mode's extent reduced
  to `hi - lo` (this is how a batch slice shrinks a Result).
- `mode_batches(mode, target)`: returns tile-aligned `[lo,hi)` batch ranges
  covering the mode's extent with each batch <= target elements (mirrors
  `mode_batches_of_trange1`, including the `max(target,1)` floor).
- `add_inplace`, `permute`, `symmetrize`, `antisymmetrize`, `mult_by_phase`,
  `adjoint`: size-preserving (no-op on extents, return self-equivalent).
- `type_id()`: a new backend id.

`eval_expr.hpp`: a thin adapter (like the other backends) so the leaf
evaluator can build these Results from an `EvalExpr`, carrying its indices
(including proto structure) and `batch_axes_`.

### Component 2: `SizeRegime` config and leaf yielder

Location: alongside the harness (test-side or a small `sizeonly/regime.hpp`).

`SizeRegime` holds:

- Scalar extents for non-proto spaces: `occ`, `aux` (Kappa), `pao`
  (mu-tilde), and any others that appear (e.g. obs virtual if present).
- A moment table `csv_moment(rank, k)` = mean over pairs of
  (domain(rank))^k, indexed by rank (OSV = 1 proto index, PNO-pair = 2,
  triple = 3) and power k in 1..4. These are INDEPENDENT numbers, not derived
  from a single mean, so both requirements hold:
  - rank-dependent: `<#OSV>` (rank 1, k 1) != `<#PNO>` (rank 2, k 1);
  - moment-dependent: `<#PNO^2>` (rank 2, k 2) != `<#PNO>^2`.
  This is the `average_csv_extent_pow(rank, k)` shape MPQC feeds the
  optimizer. k up to 4 covers the 4-occ PHL terms (CC with a 2-body
  Hamiltonian tops out at 4 co-proto PNO modes).

The same `SizeRegime` is used to build:
- the cost-model callbacks (index extent + `inner_pow`) for `single_term_opt`;
- the leaf yielder `(EvalNode) -> ResultPtr` that materializes each leaf as a
  `ResultSize`/`ResultSizeOfSize`.

so the DP and the runtime see one identical size model.

### Component 3: Catch2 harness/driver

Location: `SeQuant/tests/unit/test_eval_sizeonly.cpp`, with a committed
readable residual fixture under `SeQuant/tests/unit/data/`.

Fixture: the CSV-CCk residual in the human-readable `deserialize` format
(same text form as mpqc4 `csv_eqn_Rs.serialized`), committed as a test data
file. The harness also accepts any inline deserialized expression, so a
single offending term pasted from a trace works.

Driver flow, given a `SizeRegime` + batch config (`peak_threshold`,
`aux_target_size`, `pao_target_size`):

1. deserialize the residual (or a single term);
2. `single_term_opt` with the peak objective + threshold and the regime's
   cost callbacks -> factorized tree with per-node `batch_axes`;
3. `binarize` -> `EvalNode`;
4. `make_batched_custom_evaluator` + `evaluate` against the size backend,
   tracing on;
5. report: DP-modeled peak vs realized hwmark; per-op result sizes; and flag
   any op whose realized result exceeds `peak_threshold` (a leak);
6. assert an invariant: realized peak <= tolerance * peak_threshold, where
   tolerance is a small factor (the plan pins it, e.g. 2x, to absorb batch
   co-residency and tile granularity without admitting a real leak). On the
   C60 regime this assertion FAILS today (reproducing the 186 GB leak) and
   PASSES once the fix lands -- so the same test is both the reproducer and
   the regression lock.

## Data flow

    SizeRegime --+--> cost callbacks (extent, inner_pow) --> single_term_opt --> batch_axes
                 |                                                                    |
                 +--> leaf yielder --> ResultSize/ResultSizeOfSize                    v
                                                     \--> make_batched_custom_evaluator + evaluate
                                                                    |
                                              trace (Eval/BatchGroup/hw)  +  realized peak (hwmark)
                                                                    |
                                              report modeled vs realized + leak flags + assert

## Phasing (drives the implementation plan)

Task 1 is the go/no-go that localizes the C60 bug to factorizer vs runtime;
it is also the smallest useful slice and de-risks the sizing model.

- **Phase 1 (probe / audit): DP-annotation dump on the offending term.**
  Deserialize just the free-mu-tilde PHL 4-occ term (from the C60 trace),
  build the C60 `SizeRegime` (mu-tilde ~1800, Kappa ~4300, PNO moments, occ),
  run `single_term_opt` at 40 GB threshold, and print the per-node
  `batch_axes`. This requires only the cost-model callbacks + a minimal
  size_in_bytes -- not the full runtime backend. It answers directly:
  - if the offending intermediate's node has NO mu-tilde annotation ->
    layer-1 (cost-model/threshold) bug: the DP never asked to slice mu-tilde;
  - if it DOES carry a mu-tilde annotation -> layer-2 (runtime) bug: the
    batched evaluator dropped it.
  Validation gate: the moment sizing must reproduce the real trace numbers
  for this term (the ~1.87 GB aux-batch giant `g(mu-tilde,mu-tilde,Kappa)`
  and the ~186 GB free-mu-tilde intermediate) within a small factor, else the
  regime is mis-modeled and later verdicts are untrustworthy.

- **Phase 2: full size-only Result backend** (Component 1) with moment-aware
  `size_in_bytes`, `slice_mode`, `mode_batches`, `prod`/`sum`/`DeNest`.

- **Phase 3: full harness** (Component 3): deserialize the whole residual,
  run the batched evaluate against the backend with tracing, report modeled
  vs realized peak, flag leaks, assert the peak invariant. This is the
  end-to-end leak reproducer.

- **Phase 4: lock the C60 regime as a regression test** -- failing now,
  passing after the batching fix (whichever layer Phase 1 implicates).

## Validation strategy

Trust the harness only after it reproduces known real numbers. Anchor points
from the C60 log (approximate, per the 8-rank run):

- `g(mu-tilde, mu-tilde, Kappa)` full ~1.87 GB per aux-batch (sliced to 60
  aux batches);
- `I(i_3, i_2, mu-tilde_1241, Kappa_2; a_3)` ~186 GB (free mu-tilde + free
  aux, unsliced).

The size backend, fed the C60 regime, must land within a small factor of
these before its factorization/leak verdicts are believed. (Exact match is
not expected: the regime uses representative moments, not per-pair sizes.)

## File layout

- `SeQuant/core/eval/backends/sizeonly/result.hpp` -- `ResultSize`,
  `ResultSizeOfSize`.
- `SeQuant/core/eval/backends/sizeonly/eval_expr.hpp` -- leaf adapter.
- `SeQuant/core/eval/backends/sizeonly/regime.hpp` -- `SizeRegime` + callback
  builders + leaf yielder (or test-side if it stays test-only).
- `SeQuant/tests/unit/test_eval_sizeonly.cpp` -- probe + harness cases.
- `SeQuant/tests/unit/data/csv_ccsd_residual.txt` -- committed residual
  fixture (readable deserialize form).

## Risks

- **ToT / proto-index nesting fidelity (main risk).** The size model for
  nested (CSV) Results must match how the real backend nests -- the
  outer-occ vs inner-PNO split and the moment-aware sizing -- so
  `prod`+`DeNest`, `slice_mode`, and `mode_batches` behave consistently with
  reality. Mitigated by reusing the cost model's `inner_aware_volume` and by
  the Phase 1 validation gate against real trace numbers.
- **DP/runtime sizing drift.** If the runtime backend and the cost model ever
  computed sizes differently, modeled-vs-realized gaps would be artifacts.
  Mitigated by sharing one `SizeRegime` and one `inner_aware_volume` helper.
- **Fixture staleness.** The committed residual could drift from what MPQC
  generates. Mitigated by keeping it in the readable deserialize form and
  refreshing it from an MPQC dump when the CC/CSV equations change.
