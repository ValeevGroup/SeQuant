# Cost-model replay eval backend and PNO-CCSD cost harness

_Design spec. Repo: SeQuant. Related: the aux+PAO nested batching feature on
`feature/eval-predicted-peak-trace` (SeQuant) / mpqc4 CSV-CCk._

## Goal

A cost-model replay eval backend plus a Catch2 harness that factorizes and
replays the PNO-CCSD (CSV-CCk) residual schedule over `CacheManager` against
thin structural tokens (no numeric data, no TiledArray, no MPQC), reporting
for each evaluation exactly what the optimizer's OWN cost model reports:
projected memory size, FLOPs, and projected execution cost. This lets us
reproduce arbitrary size regimes -- including the C60 out-of-memory -- as a
unit test in seconds instead of a multi-hour cluster build-run-debug loop, and
see the optimizer-projected cost (the `dense_peak_size` frontier) side by side
with the cost actually realized when the schedule is replayed.

## Motivation

The C60 PNO-CCSD run (job 614336) was OUT_OF_MEMORY: aggregate MaxRSS 770 GB
vs a 743 GB node limit. It died building the intermediate

    I(i_3, i_2, mu-tilde_1241, Kappa_2; a_3) = C(...;mu-tilde_1242) * g(mu-tilde_1241, mu-tilde_1242, Kappa_2)

which the trace projected at ~186 GB. This intermediate carries BOTH a free
PAO index (mu-tilde_1241) and a free aux index (Kappa_2). The batched
evaluator sliced only aux ("over 60 aux batches") and never sliced mu-tilde,
so the intermediate materialized full. The optimizer ran
`objective_function: dense_peak_size` at `peak_threshold = 40 GB`, yet the
schedule did single-mode aux batching only.

`dense_peak_size` optimizes on a frontier of two quantities: peak memory and a
projected execution cost (roofline of FLOPs against `machine_balance` /
`fast_mem_elems`). To understand why the frontier search chose an aux-only
schedule -- and whether the runtime then diverged from it -- we want to report
those SAME quantities, per evaluation, on a local replay. Diagnosing this on
the cluster is intractable (hours per iteration).

## Non-goals (YAGNI)

- No numeric correctness (the backend holds no data; exactness is covered by
  the TiledArray-backed `[eval]` tests).
- No per-pair proto-index instantiation and no real size tables from MPQC (a
  moment model suffices; per-pair defeats the arbitrary-regime goal).
- No SparseShape propagation in the first cut. The cost model is a pluggable
  object (see Component 1); a sparse or adversarial model is a later plug-in
  with no change to the replay mechanics.
- No standalone CLI driver in the first cut (a Catch2 case with an editable
  regime is a seconds-long loop). The harness core is reusable, so a CLI is a
  trivial later add.

## Key idea

The quantities we want -- memory, FLOPs, projected execution cost -- are ALL
already produced by one object: the optimizer's cost model. `cost_model.hpp`
defines a CostModel concept with concrete models (`DensePeakSize`,
`DensePeakSizeBatched`, `DenseFLOPs`), each built from the reusable per-op
functors `memsize_counter(idxsz, inner_pow)` and `flops_counter(idxsz,
inner_pow)` plus `roofline_op_cost(flops, block_tiles, ...,  machine_balance,
fast_mem_elems)` (the projected execution cost; pure FLOPs when
`machine_balance == 0`).

So the harness does not re-derive any sizing: it PASSES the cost model object
and queries it per evaluation. The eval machinery (`evaluate` +
`make_batched_custom_evaluator` over `CacheManager`) supplies the SCHEDULE and
the alive-set dynamics (including batching -- a slice alive instead of the
whole); the cost model supplies the numbers. Peak memory is the max over
schedule-time of the summed cost-model memory of the alive set. Because the DP
and the replay use the SAME cost model, a gap between DP-projected and
replay-realized cost is a real DP-vs-runtime bug (the leak), not a modeling
artifact.

## Architecture

Reuse the entire pipeline unchanged:

    deserialize -> single_term_opt (chosen objective + threshold) -> binarize
                -> make_batched_custom_evaluator -> evaluate  (over CacheManager)

Add a cost-model replay backend as a fourth backend alongside
`tiledarray`/`btas`/`tapp`. Two hooks carry the cost model's numbers:

1. Peak memory reuses CacheManager verbatim. `CacheManager::note_working_set`
   maintains `working_set_hwmark_` (printed as `hw=`, documented as the running
   peak), driven by `Result::size_in_bytes()`. The cost-token Result's
   `size_in_bytes()` DELEGATES to the cost model's `memsize`. This is a
   modeling backend: it holds no data, so its only "size" IS the modeled one
   -- not an actual-size abuse but the sole size the token has. So the hwmark
   becomes the cost-model-projected peak, for free.
2. FLOPs and projected execution cost do NOT fit `size_in_bytes`. At each
   `prod`/`sum` the harness queries the SAME cost model object for the op's
   FLOPs and roofline cost, emits them in the trace, and accumulates totals.

## Components

### Component 1: the cost model object

Location: reuse `SeQuant/core/optimize/` cost functors; expose a small bundle
(e.g. `SeQuant/core/eval/backends/costtoken/cost_model_object.hpp`).

A queryable object bundling the optimizer's OWN functors, built from a
`SizeRegime` (Component 3):

- `memsize(index_set) -> bytes` (from `memsize_counter(idxsz, inner_pow)`);
- `flops(op) -> count` (from `flops_counter(idxsz, inner_pow)`);
- `exec_cost(op) -> double` (from `roofline_op_cost(..., machine_balance,
  fast_mem_elems, block_tiles, block_prefactor)`).

Default construction mirrors `DensePeakSize` (dense, moment-aware) so the
replay reports exactly what the C60 `dense_peak_size` optimization used. The
object is the single seam for the size model: a SparseShape-aware or
adversarial cost model plugs in here with NO change to the token or the replay.
Nothing is re-derived -- `memsize_counter`/`flops_counter`/`roofline_op_cost`
are the optimizer's own code.

### Component 2: cost-token Result backend

Location: `SeQuant/core/eval/backends/costtoken/{result.hpp, eval_expr.hpp}`
(mirroring `backends/tiledarray`).

Two Result types implementing the abstract `Result` interface:

- `ResultCostToken` (flat / Tensor-of-scalars): carries an ordered `Index`
  list (with proto info) and per-mode extents. No buffer.
- `ResultCostTokenNested` (Tensor-of-Tensor, for CSV amplitudes/coefficients):
  carries the outer (occ proto) and inner (PNO/OSV composite) index lists.

Virtuals:

- `size_in_bytes()`: returns `cost_model.memsize(this token's current index
  set)`. Delegated to Component 1; no bespoke math in the token. (See
  Architecture hook 1 for the semantics.)
- `prod(other, target, DeNest)`: computes the result index set (union minus
  contracted); honors `DeNest` for ToT -> T. Structural only.
- `sum(other, ...)`: same index set.
- `slice_mode(mode, lo, hi)`: copy with that mode's extent reduced to `hi -
  lo` (how a batch slice shrinks a token; makes its `memsize` drop).
- `mode_batches(mode, target)`: tile-aligned `[lo,hi)` batch ranges,
  `max(target,1)` floor (mirrors `mode_batches_of_trange1`).
- `add_inplace`, `permute`, `symmetrize`, `antisymmetrize`, `mult_by_phase`,
  `adjoint`: size-preserving (structural no-ops).
- `type_id()`: a new backend id.

FLOPs and exec cost are NOT exposed through the token; the harness reports them
per op from the cost model object (Architecture hook 2).

`eval_expr.hpp`: thin adapter so the leaf evaluator builds tokens from an
`EvalExpr`, carrying indices (proto structure) and `batch_axes_`.

### Component 3: `SizeRegime` + harness

Location: `SeQuant/tests/unit/test_eval_costtoken.cpp` plus a committed
readable residual fixture under `SeQuant/tests/unit/data/`. `SizeRegime` may
live test-side or in `backends/costtoken/regime.hpp`.

`SizeRegime` supplies the cost model's parameters:

- scalar extents for non-proto spaces: `occ`, `aux` (Kappa), `pao`
  (mu-tilde);
- a moment table `csv_moment(rank, k)` = mean over pairs of (domain(rank))^k,
  rank in {1 (OSV), 2 (PNO pair), 3 (triple)}, k in 1..4. These are
  INDEPENDENT numbers, so both hold: rank-dependent (`<#OSV>` != `<#PNO>`) and
  moment-dependent (`<#PNO^2>` != `<#PNO>^2`). This is the
  `average_csv_extent_pow(rank, k)` shape MPQC feeds the optimizer; k up to 4
  covers the 4-occ PHL terms (CC with a 2-body Hamiltonian).
- roofline parameters: `machine_balance`, `fast_mem_elems`, block params
  (matching the C60 `optimize` block: `machine_balance 200`, `fast_mem_elems
  1000000`), so the projected exec cost matches the run.

One `SizeRegime` builds the single cost model object used by BOTH
`single_term_opt` (the DP frontier) and the replay (token `size_in_bytes` +
per-op flops/exec). No drift by construction.

Fixture: the CSV-CCk residual in the human-readable `deserialize` format (same
text form as mpqc4 `csv_eqn_Rs.serialized`), committed as test data. The
harness also accepts an inline deserialized expression, so a single offending
term pasted from a trace works.

Driver flow, given `SizeRegime` + batch config (`peak_threshold`,
`aux_target_size`, `pao_target_size`):

1. deserialize the residual (or a single term);
2. `single_term_opt` with the chosen objective + threshold and the regime's
   cost model -> factorized tree with per-node `batch_axes`, and the
   DP-projected frontier point (peak, exec cost);
3. `binarize` -> `EvalNode`;
4. replay: `make_batched_custom_evaluator` + `evaluate` over `CacheManager`
   with the cost-token backend, tracing on;
5. report, per op: cost-model `{memsize, flops, exec_cost}`; and overall:
   realized peak memory (hwmark), total FLOPs, total exec cost; alongside the
   DP-projected frontier point from step 2;
6. flag any op whose realized memsize exceeds `peak_threshold` (a leak);
7. assert: realized peak <= tolerance * peak_threshold (tolerance a small
   factor, e.g. 2x, pinned in the plan). On the C60 regime this FAILS today
   (reproducing the 186 GB leak) and PASSES once the batching fix lands -- so
   the same test is reproducer and regression lock.

## Data flow

    SizeRegime --> CostModel object (memsize, flops, exec_cost)
                     |                                   |
                     v                                   v
              single_term_opt                 replay over CacheManager:
              (DP frontier: peak,             token.size_in_bytes=memsize -> hwmark (peak)
               exec cost) + batch_axes        per-op query -> flops + exec cost
                     |                                   |
                     +------------ report: DP-projected vs realized ----------+
                                   {peak memory, FLOPs, exec cost} + leak flags + assert

## Phasing (drives the implementation plan)

Task 1 is the go/no-go that localizes the C60 bug to factorizer vs runtime,
and de-risks the cost model object.

- **Phase 1 (probe / audit): DP-annotation dump on the offending term.**
  Deserialize just the free-mu-tilde PHL 4-occ term (from the C60 trace),
  build the C60 `SizeRegime` and cost model object, run `single_term_opt` at
  40 GB threshold, and print the per-node `batch_axes`. Needs only the cost
  model object + `single_term_opt` -- not the replay backend. It answers:
  - offending node has NO mu-tilde annotation -> layer-1 (cost-model /
    threshold) bug: the DP never asked to slice mu-tilde;
  - node DOES carry a mu-tilde annotation -> layer-2 (runtime) bug: the batched
    evaluator dropped it.
  Validation gate: the cost model's `memsize` must reproduce the real trace
  numbers for this term (the ~1.87 GB aux-batch `g(mu-tilde,mu-tilde,Kappa)`
  and the ~186 GB free-mu-tilde intermediate) within a small factor, and the
  FLOPs/exec cost should be sane, else the regime is mis-modeled.

- **Phase 2: cost model object + cost-token backend** (Components 1-2), with
  `size_in_bytes` delegating to `memsize`, plus `slice_mode`, `mode_batches`,
  `prod`/`sum`/`DeNest`.

- **Phase 3: full harness** (Component 3): deserialize the whole residual,
  replay over CacheManager, report per-op `{memsize, flops, exec_cost}` +
  realized peak/totals vs the DP-projected frontier, flag leaks, assert the
  peak invariant. The end-to-end leak reproducer.

- **Phase 4: lock the C60 regime as a regression test** -- failing now,
  passing after the batching fix (whichever layer Phase 1 implicates).

## Validation strategy

Trust the harness only after the cost model reproduces known real numbers.
Anchors from the C60 log (approximate, 8-rank run):

- `g(mu-tilde, mu-tilde, Kappa)` ~1.87 GB per aux-batch (sliced to 60 aux
  batches);
- `I(i_3, i_2, mu-tilde_1241, Kappa_2; a_3)` ~186 GB (free mu-tilde + free
  aux, unsliced).

The cost model, fed the C60 regime, must land within a small factor of these
before its factorization/leak verdicts are believed. Exact match is not
expected (representative moments, not per-pair sizes). Sanity-check FLOPs and
exec cost on one known contraction too.

## File layout

- `SeQuant/core/eval/backends/costtoken/cost_model_object.hpp` -- the queryable
  bundle over the optimizer's functors.
- `SeQuant/core/eval/backends/costtoken/result.hpp` -- `ResultCostToken`,
  `ResultCostTokenNested`.
- `SeQuant/core/eval/backends/costtoken/eval_expr.hpp` -- leaf adapter.
- `SeQuant/core/eval/backends/costtoken/regime.hpp` -- `SizeRegime` + builders
  (or test-side if it stays test-only).
- `SeQuant/tests/unit/test_eval_costtoken.cpp` -- probe + harness cases.
- `SeQuant/tests/unit/data/csv_ccsd_residual.txt` -- committed residual
  fixture (readable deserialize form).

## Risks

- **ToT / proto-index nesting fidelity (main risk).** The nested token must
  match how the real backend nests -- outer-occ vs inner-PNO split and the
  moment-aware sizing -- so `prod`+`DeNest`, `slice_mode`, and `mode_batches`
  behave consistently with reality. Mitigated by reusing the optimizer's
  `memsize_counter`/`inner_pow` and by the Phase 1 validation gate against real
  trace numbers.
- **`size_in_bytes` delegation.** Semantically `size_in_bytes` is the actual
  size on real backends; here it returns the modeled `memsize`. Contained: only
  the cost-token backend does this, documented as a modeling backend with no
  data. Real backends are untouched.
- **DP/replay cost drift.** If the DP and the replay ever computed cost
  differently, gaps would be artifacts. Mitigated by building both from one
  cost model object.
- **Fixture staleness.** The committed residual could drift from MPQC's.
  Mitigated by keeping it in the readable deserialize form and refreshing from
  an MPQC dump when the CC/CSV equations change.
