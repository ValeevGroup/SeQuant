# Dry-run cost-profile prediction (SQ + MPQC)

**Goal:** Make the SeQuant DryRun eval backend a faithful, reusable predictor of
a factorized schedule's full cost profile -- peak memory, FLOPs, and roofline
wall-time -- and wire it into MPQC as a predict-then-run pre-flight fed by the
real CSV size moments, so a CCk residual's memory/cost can be predicted from the
same factorized IR that the real run evaluates.

**Architecture:** One dry-run pipeline. The DryRun backend
(`SeQuant/core/eval/backends/dryrun/`) gains a `cost_profile(...)` entry point
that replays the zero-data schedule through the real eval loop + a cache built
from the real cache config, and returns a `CostProfile`. MPQC computes the CSV
domain power means `M_1..M_4`, exposes them in-memory as the same `inner_pow`
provider SeQuant already uses, and (when enabled) calls `cost_profile` on the
real IR before the real evaluation. The single code path differs between SQ tests
and MPQC only in the moment-data source.

**Tech stack:** SeQuant C++20 (`core/eval/`, `core/optimize/`), Catch2; MPQC4
(`chemistry/qc/lcao/cc/`, `.../expression/`). SeQuant is primary; MPQC consumes
it via `FETCHCONTENT_SOURCE_DIR_SEQUANT` (local) and the
`MPQC_TRACKED_SEQUANT_TAG` pin (CI). Cross-repo: SeQuant change lands first, MPQC
repins after.

## Global constraints

- No en-dashes (U+2013) or non-breaking spaces (U+00A0); ASCII hyphens only
  (pre-commit hook rejects both). No `Co-Authored-By` trailers.
- `inner_pow(composite, k)` must return the k-th POWER MEAN
  `M_k = (mean_over_pairs d^k)^(1/k)`, NOT the raw moment `mean(d^k)`, so that a
  group of `k` composites sharing a proto-set contributes
  `prod M_k = M_k^k = mean(d^k)` = the true block-sparse volume. `M_1` = today's
  printed average (clean extension).
- Predict-then-run is non-fatal: a dry-run failure logs and continues to the
  real evaluation.
- The dry-run and the real run share the SAME factorized IR, the SAME cache
  config (max_footprint, batchable-axis veto, persistence), and the SAME
  `inner_pow`; only the leaf evaluator differs (zero-data vs real tiles).
- clang-format every changed C/C++ file (`/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i`).

## The moment interface (b)

`inner_pow` is the single moment provider for BOTH the DP cost model and the
runtime Result sizing. It is a `std::function<double(Index const&, size_t k)>`
returning `M_k` for the index's composite arity (1 proto -> OSV domain, 2 protos
-> PNO domain). SeQuant already threads `inner_pow` into the optimizer
(`OptimizeOptions::inner_pow`) and the DryRun `SizeRegime`
(`csv_pno_moment[k]` / `csv_osv_moment[k]`); today the harness populates it with
a CONSTANT (M_k = the mean, independent of k). This design replaces the constant
with real `M_1..M_4`.

- **MPQC (producer):** the CSV/PNO construction already computes each pair's PNO
  domain size and each orbital's OSV domain size to print the average. Extend
  that reduction to accumulate `sum(d)`, `sum(d^2)`, `sum(d^3)`, `sum(d^4)` and
  emit `M_k = (sum(d^k)/N)^(1/k)` for k=1..4 (PNO and OSV), replacing/extending
  the single "Average ... per pair" line. Expose the eight numbers in-memory as
  an `inner_pow` closure for the dry-run.
- **SeQuant (consumer):** `inner_pow(ix, k)` looks up `M_k` (k clamped to
  [1,4]; for k>4 fall back to `M_4`, which over-estimates slightly -- documented
  and adequate, twin composites use k<=2). SQ standalone tests hardcode the C60
  `M_1..M_4` (read from an MPQC run's printout) into `df_regime`.

## Components

### a1. Verify (do not "fix") moment-aware runtime sizing (SeQuant)

The DryRun `Result::size_in_bytes()` already routes through the SAME `inner_pow`
the DP uses: `size_in_bytes() = cm_->memsize(indices_, overrides_)`
(`result.hpp:254,384`), and `memsize` runs `inner_aware_volume`
(`single_term_detail.hpp:82`), which for a group of `k` composites sharing a
proto-set multiplies `inner_pow(c, k)` over members. So a result with `k`
proto-sharing PNO legs over one occ-pair sizes as `occ^N * M_k^k = occ^N *
<d^k>` -- the true block-sparse volume -- PROVIDED `inner_pow` returns the power
mean `M_k`. This is correct today: `df_regime` sets `csv_pno_moment[k] = pno`
(the power mean of a CONSTANT domain), so a 4-PNO intermediate
`{a_1<i,i>..a_4<i,i>}` sizes as `occ^2 * pno^4` (= 358 GB at C60 constant
extents) -- the genuine size of the CC doubles PPL `W`, NOT an artifact.

Therefore a1 is a VERIFICATION task, not a fix: add a unit test pinning that a
k-composite result sizes as `occ^N * M_k^k` for a NON-constant `M_k`, so the
power-mean contract is locked before (b) supplies real, heavy-tailed moments.
Real moments have `M_k > M_1` for `k>1` (Jensen), so composite intermediates
size LARGER than the constant-domain model, not smaller; the dry-run's job is to
quantify that. (The earlier "twin-PNO occ^2*PNO^4 sizing artifact" framing was a
misdiagnosis: the sizing was always moment-aware; only the moment VALUES were a
constant stand-in.)

### a2. Faithful dry-run cache (SeQuant)

Build the dry-run cache with the gated factory
`cache_manager(nodes, is_volatile, min_repeats, footprint_of, max_footprint, is_batchable_index)`
(`cache_manager.hpp:477`), using the real run's values: `max_footprint` from the
config (e.g. 1e11), `footprint_of` = the moment-aware node size (a1),
`is_batchable_index` = the batch policy's predicate (so free-mu-tilde/K giants
are NOT cached whole, matching the real run), `is_volatile` from the volatile
leaf policy. `reset()` after each summand (persistent entries survive), so the
cache models cross-term persistence without accumulating per-term NP scratch.

### a3. Scratch-fold global peak (SeQuant)

Thread a shared peak accumulator into `make_batched_scratch`
(`eval.hpp:1142`) so each batch scratch's `working_set_hwmark` is folded into one
global peak (the outer cache misses batched transients today -- proven in
`[dryrun-perf]`: 0.2 GB outer vs 38.9 GB scratch). The accumulator is the true
peak of the whole replay. Because the batched evaluator is shared, this also
makes the REAL run's reported peak faithful (a free bonus enabling
predicted-vs-actual comparison). The peak sink is an opt-in parameter (a pointer,
null = current behavior) so non-dry-run callers are byte-unchanged unless they
pass one.

**Scope of the predicted peak.** It is the single-logical working-set high
watermark (whole-DistArray element sizes, replicated accounting -- the same
quantity the real run's `working_set_hwmark` reports on each rank), NOT per-rank
RSS nor the SLURM step-cgroup aggregate that actually triggers OOM. Relating the
logical peak to node memory requires the run's distribution/replication model,
which is out of scope here; predicted-vs-actual compares logical peak to logical
peak.

### a4. CostProfile API (SeQuant)

```
struct CostProfile {
  double peak_bytes = 0;   // global scratch-folded working-set high watermark
  double flops = 0;        // summed unweighted contraction FLOPs
  double exec_cost = 0;    // summed roofline op cost (machine-balance aware)
  std::size_t n_ops = 0;   // contraction nodes evaluated
};
CostProfile cost_profile(forest, cache_config, inner_pow, leaf_evaluator,
                         optional trace_stream);
```

`cost_profile` builds the gated cache (a2) with the peak sink (a3), replays each
tree zero-data (moment-aware sizing a1), aggregates peak/flops/exec_cost from the
per-op accounting, and (if a trace stream is given) writes the per-op trace.
`[dryrun-perf]` and `[dryrun-trace]` switch to this API; the ad-hoc cache setup,
manual hwmark reads, and gap-diagnosis prints in those tests are replaced by the
`CostProfile` return (their assertions become: twin-PNO sizes as occ^2*M_2^2,
peak in a realistic band, perf-first flops <= peak-first flops).

### b. MPQC CSV moment print + provider

In the CSV/PNO construction (`pao_to_pno_mp2.ipp` and/or the CSV builder that
prints "Average number of PNOs per pair"), extend the per-pair reduction to
`M_1..M_4` for PNO and OSV, print them, and build an `inner_pow` closure over
the eight numbers for the dry-run to consume.

### c. MPQC dry-run hook

Keyword `sequant:eval:dry_run` (bool, default false) parsed by `SeQuantEngine`.
When true, `process_for_evaluation` (`sequant.cpp`), immediately before the real
residual eval, calls SeQuant's `cost_profile` on the ALREADY-built factorized IR
(the per-summand optimized eval forest) with: the in-memory `M_1..M_4`
`inner_pow` (b), the real cache config (max_footprint, batch policy, volatile
predicate), and the zero-data DryRun leaf evaluator. It logs the `CostProfile`
(peak/flops/exec_cost/n_ops) and optionally writes the per-op trace to a file
(`sequant:eval:dry_run_trace` path). Then the real evaluation proceeds unchanged
(predict-then-run). The dry-run builds its OWN cache (from the same config), so
it does not perturb the real run's cache.

## Data flow

CSV domains -> `M_1..M_4` (print + in-memory `inner_pow`) -> CCk builds the
factorized IR (unchanged) -> if `dry_run`: `cost_profile(IR, real cache config,
M_k inner_pow, zero-data leaf)` -> log CostProfile (+ optional trace) -> real
eval runs. Predicted vs actual: both sides read the now-faithful scratch-folded
peak (a3).

## Error handling

- Dry-run replay throws -> catch, log "dry-run failed: <what>", continue to the
  real eval (predict-then-run is never blocked by prediction).
- No moment provider (dry-run driven without moments) -> constant-moment fallback
  with a one-line warning (current behavior; sizes are then indicative only).
- `k > 4` requested from `inner_pow` -> return `M_4` (documented slight
  over-estimate; twin composites need only k<=2).

## Testing

- **SeQuant unit:** (1) `inner_pow` power-mean helper: `M_k` from a small size
  list equals `(mean d^k)^(1/k)`. (2) moment-aware Result sizing: a twin-PNO
  result sizes as `occ^2*M_2^2`, not `occ^2*M_1^4`. (3) gated-cache dry-run peak:
  a free-batchable-axis giant is NOT cached (footprint gate + veto), so the
  peak matches the sliced size. (4) scratch-fold: a batched contraction's global
  peak >= its per-batch scratch hwmark.
- **SeQuant harness:** `[dryrun-perf]` / `[dryrun-trace]` use `cost_profile`
  with real C60 `M_1..M_4`; assert twin-PNO sizing fixed and the whole-residual
  peak lands in a realistic band (order ~1 TB or below, not 11.9 TB).
- **MPQC unit:** `M_1..M_4` computed correctly from a synthetic domain-size list;
  `sequant:eval:dry_run` parses.
- **MPQC validation:** a small CCk input with `dry_run=true` logs a CostProfile
  and still completes the real run; a coarse predicted-vs-actual peak check
  (same order of magnitude).

## Sequencing and cross-repo

Implement b -> a -> c: (b) MPQC moment print is standalone and gives real
numbers to hardcode; (a1-a4) SeQuant realism consumes them; (c) MPQC hook drives
the SeQuant API. The SeQuant changes (a, and the `inner_pow` contract) land and
are pushed first; MPQC repins `MPQC_TRACKED_SEQUANT_TAG` (CURRENT + PREVIOUS) in
its own commit before/with the MPQC wiring (b, c). Local MPQC builds see SeQuant
via `FETCHCONTENT_SOURCE_DIR_SEQUANT`; the repin is for CI/clean builds.
