# Dry-run cost-profile prediction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the SeQuant DryRun eval backend a faithful predictor of a factorized schedule's full cost profile (peak memory, FLOPs, roofline cost), and wire it into MPQC as a predict-then-run pre-flight fed by the real CSV domain moments.

**Architecture:** SeQuant gains a `cost_profile(...)` entry point that replays the zero-data eval forest through the real eval loop plus a cache built from the real cache config, folding batched-scratch high-watermarks into one global peak, and returns a `CostProfile`. MPQC computes the CSV domain power means `M_1..M_4`, exposes them in-memory as the `inner_pow` provider SeQuant already threads through the optimizer and DryRun `SizeRegime`, and (when enabled) calls `cost_profile` on the real IR before the real evaluation. SeQuant lands first; MPQC repins and wires.

**Tech Stack:** SeQuant C++20 (`core/eval/backends/dryrun/`, `core/optimize/`), Catch2; MPQC4 (`chemistry/qc/lcao/mbpt/`, `.../expression/`). MPQC consumes SeQuant via `FETCHCONTENT_SOURCE_DIR_SEQUANT` (local) and `MPQC_TRACKED_SEQUANT_TAG` (CI).

## Global Constraints

- No en-dashes (U+2013) or non-breaking spaces (U+00A0); ASCII hyphens only (a pre-commit hook rejects both). No `Co-Authored-By` trailers in commits.
- `inner_pow(composite, k)` MUST return the k-th POWER MEAN `M_k = (mean_over_pairs d^k)^(1/k)`, NOT the raw moment `mean(d^k)`. The optimizer's `inner_aware_volume` (`SeQuant/core/optimize/single_term_detail.hpp:87-92`) multiplies `inner_pow(c,k)` once per member of a k-composite group, so the group product is `M_k^k = mean(d^k)` and `outer_nocc^N * M_k^k = Sum_pairs d^k` = the true block-sparse volume. `M_1` = today's printed average. At k=1, `M_1 = mean(d)`, matching the current single-average print exactly.
- The dry-run and the real run share the SAME factorized IR, the SAME cache config (`max_footprint`, batchable-axis veto, volatile predicate, `min_repeats`), and the SAME `inner_pow`; only the leaf evaluator differs (zero-data DryRun vs real tiles).
- Predict-then-run is non-fatal: a dry-run failure logs and continues to the real evaluation.
- The predicted peak is the single-logical working-set high watermark (whole-DistArray element sizes, replicated accounting), NOT per-rank RSS nor the SLURM step-cgroup aggregate. Predicted-vs-actual compares logical peak to logical peak.
- clang-format every changed C/C++ file: `/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i <file>`. Do NOT format `.tpp` files.
- SeQuant threads its unit-test binaries internally; NEVER run two SeQuant unit-test binaries concurrently. No parallel MPQC runs. Use `cmake-build-release` for any MPQC validation run.
- Commit after each task; push only the SeQuant branch when all SeQuant tasks are green (Task 6 gates the MPQC repin on the pushed SeQuant SHA).

Reference spec: `SeQuant/doc/dev/specs/2026-07-05-dryrun-cost-profile-prediction-design.md`.

---

## Repo A: SeQuant (Tasks 1-5)

Branch: continue on the existing `feature/eval-predicted-peak-trace` (SeQuant). Source tree: `/Users/efv/code/SeQuant`. Build dir: `cmake-build-release` (create with the project's standard cmake invocation if absent). Unit-test binary target: `unit_tests-sequant`. Run a single test by tag with `unit_tests-sequant "<tag>"`.

### Task 1: Lock the power-mean sizing contract (a1 verification)

The DryRun sizing is ALREADY moment-aware (`Result::size_in_bytes()` -> `cm_->memsize` -> `inner_aware_volume`). This task does NOT change sizing; it pins the power-mean contract with a test so later real-moment data cannot silently regress to raw moments, and documents the contract at the `SizeRegime` definition.

**Files:**
- Test: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add one `TEST_CASE`)
- Modify (doc only): `SeQuant/core/eval/backends/dryrun/size_regime.hpp:18-52` (comment above `csv_pno_moment`)

**Interfaces:**
- Consumes: `sequant::eval::dryrun::SizeRegime` with fields `space_extent`, `csv_pno_moment[5]`, `csv_osv_moment[5]`, method `inner_pow(Index const&, size_t) const`, and `inner_pow_fn()`. The optimizer helper `sequant::opt::detail::inner_aware_volume(tot_idxs, ixex, inner_pow)` (already exists).
- Produces: nothing new; a locked invariant.

- [ ] **Step 1: Write the failing test**

Add to `test_eval_dryrun.cpp`. This builds a `SizeRegime` whose PNO moment array holds a NON-constant power-mean sequence (a heavy-tailed domain: `M_1<M_2<M_3<M_4`), then sizes a synthetic 2-composite (twin-PNO) result over an occ-pair and asserts the size equals `occ^2 * M_2^2`, NOT `occ^2 * M_1^2` (the flat-average under-count) and NOT `occ^2 * M_1^4` (the old misdiagnosed "artifact").

```cpp
TEST_CASE("dryrun power-mean sizing contract", "[dryrun][sizing]") {
  using namespace sequant;
  using sequant::eval::dryrun::SizeRegime;

  // Heavy-tailed PNO domain: power means strictly increasing in k.
  // (For any non-degenerate distribution, M_1 < M_2 < M_3 < M_4 by Jensen.)
  SizeRegime r;
  const std::size_t occ = 10;
  r.space_extent[L"i"] = occ;  // occupied index extent
  const std::array<double, 5> Mk{1.0, 30.0, 40.0, 55.0, 75.0};  // M_0..M_4
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = Mk[k];
    r.csv_osv_moment[k] = Mk[k];
  }
  auto ip = r.inner_pow_fn();

  // A composite PNO index a<i_1,i_2> has a 2-proto (PNO) domain.
  Index i1(L"i_1"), i2(L"i_2");
  Index a1 = Index(L"a_1", {i1, i2});
  Index a2 = Index(L"a_2", {i1, i2});

  // Two composites sharing the SAME proto-set => a k=2 group.
  // inner_aware_volume multiplies ip(c,2) once per member => M_2^2.
  opt::detail::tot_indices tot;
  tot.outer = {i1, i2};              // occ^2 outer
  tot.inner = {a1, a2};              // twin PNO composites
  auto ixex = [&](Index const& ix) {
    return static_cast<double>(r.extent(ix));
  };
  const double vol = opt::detail::inner_aware_volume(tot, ixex, ip);

  const double expected = double(occ) * double(occ) * Mk[2] * Mk[2];  // occ^2 * M_2^2
  CHECK(vol == Catch::Approx(expected));
  // Guard against the two wrong models:
  CHECK(vol != Catch::Approx(double(occ) * occ * Mk[1] * Mk[1]));      // flat-average under-count
  CHECK(vol != Catch::Approx(double(occ) * occ * std::pow(Mk[1], 4))); // old "occ^2*PNO^4 artifact"
}
```

Before writing, open `size_regime.hpp` and `single_term_detail.hpp` to confirm the exact spellings of `tot_indices`, its `outer`/`inner` members, the `Index` proto-index constructor, and the `inner_aware_volume` namespace path; adjust the `#include`s and qualified names in the test to match. If `tot_indices` is not default-constructible with public `outer`/`inner`, build it via the existing `tot_indices(...)` splitter used elsewhere in this file (grep `tot_indices` in the test for a live example) rather than aggregate-initializing.

- [ ] **Step 2: Run test to verify it fails**

Run: `unit_tests-sequant "[sizing]"`
Expected: FAIL to COMPILE if a name is wrong (fix names), or PASS immediately if the contract already holds. If it PASSES on first run, that is the desired end state for this task, but you MUST still confirm it fails when the contract is broken: temporarily change the test's `expected` to `occ*occ*Mk[1]*Mk[1]` and confirm FAIL, then revert. Record that you did this in the commit body.

- [ ] **Step 3: Document the contract**

Edit the comment block directly above `csv_pno_moment` in `size_regime.hpp` to state verbatim:

```cpp
  // csv_pno_moment[k] / csv_osv_moment[k] hold the k-th POWER MEAN
  // M_k = (mean_over_pairs d^k)^(1/k) of the per-pair PNO / per-orbital OSV
  // domain size d, for k in [1,4] (index 0 is unused, set to 1). inner_pow()
  // returns M_k so that inner_aware_volume's per-member product over a
  // k-composite group is M_k^k = mean(d^k), and outer_nocc^N * M_k^k equals
  // the true block-sparse volume Sum_pairs d^k. Do NOT store raw moments
  // mean(d^k) here: that would over-count k-composite groups by a further
  // power of k. For a constant domain d, M_k = d for all k.
```

- [ ] **Step 4: Run test to verify it passes**

Run: `unit_tests-sequant "[sizing]"`
Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/tests/unit/test_eval_dryrun.cpp SeQuant/core/eval/backends/dryrun/size_regime.hpp
git add SeQuant/tests/unit/test_eval_dryrun.cpp SeQuant/core/eval/backends/dryrun/size_regime.hpp
git commit -m "dryrun: lock power-mean sizing contract with a regression test"
```

---

### Task 2: Faithful gated dry-run cache with per-term reset (a2)

Build the dry-run cache with the gated `cache_manager` factory so free-batchable-axis giants are not cached whole (matching the real run), and `reset()` between summands so persistent entries survive but per-term NP scratch does not accumulate across terms.

**Files:**
- Modify: `SeQuant/core/eval/backends/dryrun/` (the dry-run driver that today builds a cache; grep for `cache_manager(` under `backends/dryrun/`). If the driver is inside the test harness rather than a library header, the cache-build helper you add lives in a new small header `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` created in Task 4; for this task add the helper to that header's eventual home but expose it now as a free function.
- Test: `SeQuant/tests/unit/test_eval_dryrun.cpp`

**Interfaces:**
- Consumes: the gated factory `sequant::eval::cache_manager(nodes, is_volatile, min_repeats, footprint_of, max_footprint, is_batchable_index)` (`SeQuant/core/eval/cache_manager.hpp:477`); the DryRun node-size functor (moment-aware, from Task 1's `inner_pow`).
- Produces: `struct CacheConfig { double max_footprint; std::size_t min_repeats; std::function<bool(EvalNode const&)> is_volatile; std::function<bool(Index const&)> is_batchable_index; };` (define in `cost_profile.hpp`, consumed by Task 4). `is_volatile`/`is_batchable_index` element types must match the gated factory's actual parameter types; confirm them by reading `cache_manager.hpp:461-490` and copy the exact signatures.

- [ ] **Step 1: Write the failing test**

Add a `TEST_CASE("dryrun gated cache vetoes free-batchable giant", "[dryrun][cache]")` that:
1. Builds a `SizeRegime` (reuse `df_regime` from this file) and an eval forest containing one node whose result carries a free batchable axis (a `mu~`/`K`-carrying giant) plus a small cacheable intermediate.
2. Builds the cache via the gated factory with `max_footprint = 1e11`, an `is_batchable_index` that returns true for the giant's free axis, `footprint_of` = the moment-aware node size, `min_repeats = 1`.
3. Asserts the giant is NOT resident in the cache after evaluation (query the cache's residency the same way `[dryrun-perf]` inspects it today; grep `cache` in the test for the accessor), while the small intermediate IS cached when repeated.

Model the forest construction on the existing `[dryrun-perf]` test in the same file (it already optimizes and replays the giant term). Reuse its term selection rather than authoring a new expression.

- [ ] **Step 2: Run test to verify it fails**

Run: `unit_tests-sequant "[cache]"`
Expected: FAIL (giant currently cached, or the gated factory not yet wired into the dry-run path).

- [ ] **Step 3: Wire the gated factory + per-term reset**

In the dry-run cache build path, replace the simple `cache_manager(nodes, min_repeats)` call with the gated overload, passing `CacheConfig`'s fields. After each summand's replay, call `cache.reset()` (persistent entries survive by the factory's persistence rule; NP scratch is dropped). Read `cache_manager.hpp:461-490` to confirm `reset()` semantics (it must keep persistent/`min_repeats`-qualified entries and drop transients); if `reset()` clears everything, use the narrower per-term-scratch clear the class exposes, or scope the transient cache per summand and keep a separate persistent cache across summands (mirror how the real eval loop in `eval.hpp` scopes persistence across summands - grep `reset(` in `eval.hpp`).

- [ ] **Step 4: Run test to verify it passes**

Run: `unit_tests-sequant "[cache]"`
Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i <changed .hpp/.cpp files>
git add -A && git commit -m "dryrun: gated cache with batchable-axis veto and per-term reset"
```

---

### Task 3: Scratch-fold global peak accumulator (a3)

The batched evaluator runs each batch in a separate scratch cache (`make_batched_scratch`, `SeQuant/core/eval/eval.hpp:1142`), so batched-inner transients are invisible to the outer cache's `working_set_hwmark` (proven: 0.2 GB outer vs 38.9 GB scratch). Thread an opt-in shared peak sink so each scratch's high-watermark folds into one global peak. Null sink = byte-unchanged behavior for the real run.

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (around `make_batched_scratch`, line ~1142, and the batched contraction loop that calls it)
- Modify: `SeQuant/core/eval/cache_manager.hpp` (expose `working_set_hwmark()` if not already public; it is per the summary)
- Test: `SeQuant/tests/unit/test_eval_dryrun.cpp`

**Interfaces:**
- Consumes: `CacheManager::working_set_hwmark()` (existing).
- Produces: an optional peak sink parameter. Define `using PeakSink = std::atomic<double>*;` (or a plain `double*` guarded by the eval's existing serialization - confirm the batched loop's threading model in `eval.hpp` before choosing; if batches run concurrently use `std::atomic<double>` with a `fetch_max` CAS loop, else a plain `double*` max-assign). Add the sink as a trailing defaulted parameter (`PeakSink peak = nullptr`) to `make_batched_scratch` and to the batched-eval entry it feeds, so all existing callers compile unchanged.

- [ ] **Step 1: Write the failing test**

Add `TEST_CASE("dryrun scratch-fold captures batched peak", "[dryrun][peak]")`:
1. Reuse the `[dryrun-perf]` giant-term forest and objective that produces a batched (K-sliced) contraction.
2. Evaluate once with a `PeakSink` passed in; capture `global_peak`.
3. Also read the outer cache's `working_set_hwmark()` (as today).
4. Assert `global_peak >= outer_hwmark` and, for this specifically-batched term, `global_peak` is strictly greater (the batched transient exceeds the outer residency): `CHECK(global_peak > outer_hwmark * 2.0)` (the observed gap is ~195x; a 2x floor is a safe, non-flaky assertion).

- [ ] **Step 2: Run test to verify it fails**

Run: `unit_tests-sequant "[peak]"`
Expected: FAIL to compile (no `PeakSink` parameter yet) or FAIL the assertion (peak not folded).

- [ ] **Step 3: Implement the fold**

In `make_batched_scratch` and its caller: after each batch's scratch cache finishes, read `scratch.working_set_hwmark()` and fold it into `*peak` (max). If batches run concurrently, use a `fetch_max` loop:
```cpp
if (peak) {
  double cur = peak->load(std::memory_order_relaxed);
  const double cand = scratch.working_set_hwmark();
  while (cand > cur &&
         !peak->compare_exchange_weak(cur, cand, std::memory_order_relaxed)) {}
}
```
Also fold the OUTER cache's `working_set_hwmark()` into the same sink at end of replay, so `*peak = max(outer, all scratch)` = the true whole-replay peak. Confirm the exact place the outer replay completes (grep the function that owns the outer `CacheManager` and drives all summands).

- [ ] **Step 4: Run test to verify it passes**

Run: `unit_tests-sequant "[peak]"`
Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/core/eval/eval.hpp SeQuant/core/eval/cache_manager.hpp \
  SeQuant/tests/unit/test_eval_dryrun.cpp
git add -A && git commit -m "eval: opt-in scratch-fold peak sink for faithful batched-replay peak"
```

---

### Task 4: CostProfile struct and cost_profile() API (a4)

Wrap Tasks 2-3 in a single reusable entry point that both SeQuant tests and MPQC call.

**Files:**
- Create: `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` (declarations + inline template body, header-only per SeQuant convention; if it must cross a TU boundary add `cost_profile.cpp` and register in `SeQuant/core/CMakeLists.txt` under the eval sources)
- Modify: `SeQuant/core/CMakeLists.txt` only if a `.cpp` is added
- Test: `SeQuant/tests/unit/test_eval_dryrun.cpp`

**Interfaces:**
- Consumes: `CacheConfig` (Task 2), `PeakSink` (Task 3), the DryRun leaf evaluator and zero-data `Result`, `inner_pow` (`std::function<double(Index const&, std::size_t)>`).
- Produces:
```cpp
namespace sequant::eval::dryrun {

struct CostProfile {
  double peak_bytes = 0;   // global scratch-folded working-set high watermark
  double flops = 0;        // summed unweighted contraction FLOPs
  double exec_cost = 0;    // summed roofline op cost (machine-balance aware)
  std::size_t n_ops = 0;   // contraction nodes evaluated
};

/// Replays a factorized eval forest zero-data through the real eval loop with a
/// cache built from `cfg`, folding batched-scratch high-watermarks into
/// CostProfile::peak_bytes. Optionally writes a per-op trace to `trace`.
/// \param forest        per-summand optimized eval forest (the real IR)
/// \param cfg           real cache config (max_footprint, veto, volatile, repeats)
/// \param inner_pow     M_k power-mean provider (see size_regime.hpp contract)
/// \param leaf_eval     zero-data DryRun leaf evaluator
/// \param trace         optional per-op trace sink (nullptr = no trace)
CostProfile cost_profile(
    EvalForest const& forest, CacheConfig const& cfg,
    std::function<double(Index const&, std::size_t)> inner_pow,
    LeafEvaluator const& leaf_eval, std::wostream* trace = nullptr);

}  // namespace sequant::eval::dryrun
```
Replace `EvalForest`, `LeafEvaluator` with the concrete types the existing `[dryrun-perf]` replay uses (grep the test for the forest container type and the leaf functor type; use those exact spellings so MPQC and the tests share them).

- [ ] **Step 1: Write the failing test**

Add `TEST_CASE("cost_profile returns peak/flops/exec/n_ops", "[dryrun][cost_profile]")`:
1. Build the `df_regime` with real C60 `M_1..M_4` (hardcode the numbers you will later read from an MPQC printout; until Task 7 exists, use the constant `M_k = 42` PNO / `310` OSV that `df_regime(1800,4320,120,42.0,310.0)` already encodes, extended so `csv_pno_moment[k]` differs per k once real numbers land - for now constant is fine, the test asserts structure not magnitude).
2. Build `CacheConfig` (max_footprint 1e11, veto on the giant's free axis, min_repeats 1).
3. Call `cost_profile(forest, cfg, regime.inner_pow_fn(), leaf_eval)`.
4. Assert: `cp.n_ops > 0`, `cp.flops > 0`, `cp.peak_bytes > 0`, and `cp.peak_bytes` equals the Task-3 folded peak for the same forest (call both, compare).

- [ ] **Step 2: Run test to verify it fails**

Run: `unit_tests-sequant "[cost_profile]"`
Expected: FAIL (no such symbol).

- [ ] **Step 3: Implement cost_profile**

In `cost_profile.hpp`: build the gated cache (Task 2 helper), allocate a `PeakSink`, loop over summands replaying zero-data through the real eval loop (reuse the exact replay call the `[dryrun-perf]` test performs), accumulate `flops`/`exec_cost`/`n_ops` from the per-op accounting the DryRun backend already produces (grep the test for where it reads per-op flop/cost today), set `peak_bytes` from the sink, and if `trace` is non-null pass it into the replay's trace hook. Return the `CostProfile`.

- [ ] **Step 4: Run test to verify it passes**

Run: `unit_tests-sequant "[cost_profile]"`
Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/cost_profile.hpp \
  SeQuant/tests/unit/test_eval_dryrun.cpp
git add -A && git commit -m "dryrun: CostProfile struct and cost_profile() reusable entry point"
```

---

### Task 5: Retrofit [dryrun-perf] and [dryrun-trace] onto cost_profile (a4)

Replace the ad-hoc cache setup, manual hwmark reads, and gap-diagnosis prints in the two harness tests with `cost_profile`, so there is one code path.

**Files:**
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (`[dryrun-perf]` and `[.][dryrun-trace]` cases)

**Interfaces:**
- Consumes: `cost_profile` (Task 4).
- Produces: nothing new.

- [ ] **Step 1: Update [dryrun-perf]**

Replace its manual cache build + replay + hwmark read with a `cost_profile` call under each objective (perf-first, peak-first). Keep the objective comparison. New assertions:
- `perf.flops <= peak.flops` (perf-first minimizes flops).
- Both `peak_bytes` land in a realistic band for the constant-moment C60 giant: `CHECK(perf.peak_bytes < 1e12); CHECK(perf.peak_bytes > 1e11);` (perf-first forms the ~358 GB 4-PNO `W`; band 100 GB..1 TB is safe and non-flaky). Document in a comment that with real heavy-tailed moments this rises.

- [ ] **Step 2: Update [.][dryrun-trace]**

Replace its manual whole-residual replay with a `cost_profile(forest, cfg, inner_pow, leaf, &trace_file)` call, keeping the env knobs (`DRYRUN_OBJ`, `DRYRUN_TRACE_FILE`, `DRYRUN_MAX_TERMS`). Assert the returned `CostProfile.peak_bytes` is finite and the trace file is non-empty.

- [ ] **Step 3: Run both tests**

Run: `unit_tests-sequant "[dryrun-perf]"` then `unit_tests-sequant "[dryrun-trace]"` (serially; never concurrent).
Expected: PASS (`[dryrun-trace]` is hidden `[.]`, invoke by explicit tag).

- [ ] **Step 4: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i SeQuant/tests/unit/test_eval_dryrun.cpp
git add SeQuant/tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: route [dryrun-perf]/[dryrun-trace] through cost_profile"
```

- [ ] **Step 5: Push SeQuant branch**

```bash
git -C /Users/efv/code/SeQuant push
git -C /Users/efv/code/SeQuant rev-parse HEAD   # record SHA for Task 6
```

---

## Repo B: MPQC (Tasks 6-8)

Source tree: `/Users/efv/code/mpqc4`. Branch: continue on `feature/eval-predicted-peak-trace`. Build dir: `cmake-build-release`. Executable target: `mpqc`.

### Task 6: Repin SeQuant

**Files:**
- Modify: `external/versions.cmake:13-14`

**Interfaces:**
- Consumes: the pushed SeQuant SHA from Task 5 Step 5.
- Produces: an updated pin so CI/clean builds see the new API.

- [ ] **Step 1: Update the pin**

Set `MPQC_TRACKED_SEQUANT_PREVIOUS_TAG` to the current `MPQC_TRACKED_SEQUANT_TAG` value (`91dd48a2dbb9617a3812141aa01c0e3de43b7d2e`), and set `MPQC_TRACKED_SEQUANT_TAG` to the SHA recorded in Task 5.

- [ ] **Step 2: Commit (pin change alone, per CLAUDE.md)**

```bash
git add external/versions.cmake
git commit -m "external: repin SeQuant to <new-sha> (dry-run cost_profile API)"
```

---

### Task 7: MPQC CSV moment print + inner_pow provider (b)

Extend the existing CC-pair PNO average reduction to `M_1..M_4` power means, print them, and build an `inner_pow` closure the dry-run hook (Task 8) consumes. OSV moments are computed where per-pair data exists; where only a total is available, emit `M_1` and set `M_2..M_4 = M_1` with a one-line caveat (documented limitation, not a placeholder - PNO composites dominate the C60 giant).

**Files:**
- Modify: `src/mpqc/chemistry/qc/lcao/mbpt/pao_to_pno_mp2.ipp:570-597` (PNO moment accumulation + print) and `:1474-1494` (OSV)
- Modify: the CSV-holder struct that already carries per-pair PNO/OSV counts to store the eight moments and expose an `inner_pow` closure (grep the return type of the function containing line 570; the moments live next to the CSV maps handed to CCk). If there is no natural holder, add a small `struct CsvMoments { std::array<double,5> pno{1,1,1,1,1}; std::array<double,5> osv{1,1,1,1,1}; std::function<double(...)> inner_pow() const; };` next to the CSV data type and thread it to where CCk builds its `SeQuantEngine`.
- Test: `tests/unit/` new file `test_csv_moments.cpp` (add to `tests/unit/CMakeLists.txt`)

**Interfaces:**
- Consumes: the per-pair `cc_pnos` accumulation loop at `pao_to_pno_mp2.ipp:574-585`.
- Produces: `CsvMoments{ std::array<double,5> pno, osv; }` with `M_k = (sum(d^k)/N)^(1/k)`, and `inner_pow(Index const&, size_t k) -> double` returning `pno[k]` for 2-proto composites, `osv[k]` for 1-proto, `pow(extent,k)` otherwise. The `inner_pow` closure signature and composite-arity dispatch MUST match `SizeRegime::inner_pow` (`size_regime.hpp:36-45`) exactly so the two are interchangeable.

- [ ] **Step 1: Write the failing test**

`tests/unit/test_csv_moments.cpp`: a pure function `csv_power_means(std::vector<std::size_t> const& d) -> std::array<double,5>` computing `M_k = (mean d^k)^(1/k)` for k=0..4 (M_0 := 1). Test with a small heavy-tailed list, e.g. `{1, 2, 2, 40}`:
```cpp
TEST_CASE("csv_power_means", "[csv][moments]") {
  const std::vector<std::size_t> d{1, 2, 2, 40};
  const auto M = mpqc::lcao::csv_power_means(d);
  const double N = 4;
  auto mean_pow = [&](int k) {
    double s = 0; for (auto x : d) s += std::pow(double(x), k); return s / N;
  };
  CHECK(M[1] == Catch::Approx(mean_pow(1)));                       // M_1 = mean(d)
  CHECK(M[2] == Catch::Approx(std::pow(mean_pow(2), 1.0 / 2)));
  CHECK(M[4] == Catch::Approx(std::pow(mean_pow(4), 1.0 / 4)));
  CHECK(M[4] > M[1]);                                              // heavy tail: M_4 > M_1
}
```
Add `test_csv_moments.cpp` to `tests/unit/CMakeLists.txt`.

- [ ] **Step 2: Run test to verify it fails**

Run: `ctest --test-dir cmake-build-release -R csv_moments` (or run the unit binary directly). Expected: FAIL (undefined `csv_power_means`).

- [ ] **Step 3: Implement `csv_power_means` and wire the reduction**

- Add `csv_power_means` (free function, in an `mbpt` header included by `pao_to_pno_mp2.ipp`, e.g. a small `csv_moments.h`; declare in `mpqc::lcao`).
- In the CC-pair loop (`:574-585`), replace the single `cc_pnos += pno` with four accumulators `cc_pno_pk[k] += std::pow(double(pno), k)` for k=1..4 (keep `cc_pnos` for the existing average line, or derive M_1 from it). `world.gop.sum` each. After reduction compute `M_k = pow(cc_pno_pk[k]/cc_npairs, 1.0/k)`.
- Replace/extend the print at `:594-597`:
```cpp
ExEnv::out0() << "Average number of PAOs per pair: " << avg_paos << "\n";
ExEnv::out0() << "PNO domain power means M_1..M_4 per pair: "
              << M_pno[1] << " " << M_pno[2] << " " << M_pno[3] << " "
              << M_pno[4] << "\n";
```
(Keep the "Average number of PNOs per pair" line = `M_pno[1]` for backward-compatible log grepping.)
- Store `M_pno[1..4]` into the `CsvMoments` on the CSV holder; build the `inner_pow` closure. For OSV (`:1474-1494`): if per-pair OSV domain sizes are reachable from `pno_builder`, accumulate the same four power sums; otherwise set `M_osv[1] = nosvs_total/numIJ_total` and `M_osv[2..4] = M_osv[1]`, and print `"OSV domain power means (M_2..M_4 approximated by M_1: per-pair OSV sizes not exposed)"`.

- [ ] **Step 4: Run test to verify it passes**

Run: `ctest --test-dir cmake-build-release -R csv_moments`. Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  src/mpqc/chemistry/qc/lcao/mbpt/pao_to_pno_mp2.ipp \
  src/mpqc/chemistry/qc/lcao/mbpt/csv_moments.h \
  tests/unit/test_csv_moments.cpp
git add -A && git commit -m "mbpt: print CSV PNO/OSV domain power means M_1..M_4 and expose inner_pow"
```

---

### Task 8: MPQC dry-run hook (predict-then-run) (c)

Add `sequant:eval:dry_run` (bool, default false). When true, immediately before the real residual evaluation, call SeQuant's `cost_profile` on the already-built factorized IR with the in-memory `M_1..M_4` provider (Task 7), the real cache config, and the zero-data DryRun leaf evaluator; log the `CostProfile`; then run the real evaluation unchanged.

**Files:**
- Modify: `src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp` (parse the keyword, ~line 78-90 block) and `sequant_engine.h` (add a `bool dry_run_` field + accessor; default in `SeQuantEngineDefaults`)
- Modify: `src/mpqc/chemistry/qc/lcao/expression/sequant.cpp` (`process_for_evaluation`, the function that holds the optimized forest and drives the real eval; grep `optimize`/`eval` near line 274-380) to call `cost_profile` when `engine.dry_run()` is set
- Test: MPQC validation input under `tests/validation/reference/inputs/` (small CCk with `dry_run: true`)

**Interfaces:**
- Consumes: `sequant::eval::dryrun::cost_profile` (Task 4), `CsvMoments::inner_pow()` (Task 7), the engine's existing cache config (`eval:cache:max_footprint`, `eval:cache:min_repeats`), the optimized eval forest already built in `process_for_evaluation`.
- Produces: a logged `CostProfile` line; no change to the real result.

- [ ] **Step 1: Parse the keyword**

In `sequant_engine.cpp` alongside the `eval:cache:*` parses (line ~78-90):
```cpp
dry_run_ = kv.value<bool>("eval:dry_run", d.dry_run);
if (kv.exists("eval:dry_run_trace"))
  dry_run_trace_path_ = kv.value<std::string>("eval:dry_run_trace", std::string{});
```
Add `bool dry_run_ = false; std::string dry_run_trace_path_;` and accessors `bool dry_run() const; std::string const& dry_run_trace_path() const;` to `sequant_engine.h`; add `bool dry_run = false;` to `SeQuantEngineDefaults`. Document the two keys in the engine's Doxygen key table.

- [ ] **Step 2: Write the failing validation test**

Copy the smallest existing CCk PNO input from `tests/validation/reference/inputs/` to a new `<name>-dryrun.json`, add `"dry_run": true` to its `sequant:eval` block, and its reference output (tighter precision than input) identical to the base input's energy (the dry-run must NOT change the result). Register the pair per `tests/validation/`'s CMake/JSON harness. Run it; expect FAIL (dry_run not yet honored -> either unknown-key error or no CostProfile in the log, depending on how you assert). Assert via the harness that the run completes with the reference energy AND the log contains `"dry-run CostProfile"`.

- [ ] **Step 3: Wire the predict-then-run call**

In `process_for_evaluation`, immediately before the real per-summand eval, guard `if (engine.dry_run())` and:
```cpp
try {
  auto cfg = sequant::eval::dryrun::CacheConfig{
      engine.cache_max_footprint(), engine.cache_min_repeats(),
      /*is_volatile=*/engine.volatile_predicate(),
      /*is_batchable_index=*/engine.batchable_predicate()};
  std::wofstream trace;
  std::wostream* tp = nullptr;
  if (!engine.dry_run_trace_path().empty()) {
    trace.open(engine.dry_run_trace_path());
    tp = &trace;
  }
  const auto cp = sequant::eval::dryrun::cost_profile(
      forest, cfg, csv_moments.inner_pow(), make_dryrun_leaf_evaluator(), tp);
  ExEnv::out0() << "dry-run CostProfile: peak=" << cp.peak_bytes / 1e9
                << " GB flops=" << cp.flops << " exec=" << cp.exec_cost
                << " n_ops=" << cp.n_ops << "\n";
} catch (std::exception const& e) {
  ExEnv::out0() << "dry-run failed (continuing to real eval): " << e.what()
                << "\n";
}
```
Replace `engine.volatile_predicate()`, `engine.batchable_predicate()`, `make_dryrun_leaf_evaluator()`, `csv_moments`, and `forest` with the real accessors/locals in scope at that point (grep the function for the forest variable name and how the real cache config is currently assembled; reuse those exact expressions so dry-run and real share config). The DryRun leaf evaluator is the zero-data one from `backends/dryrun/`; if MPQC has no helper, add a one-line factory next to the eval call.

- [ ] **Step 4: Run the validation test**

Run: `ctest --test-dir cmake-build-release -R <name>-dryrun`. Expected: PASS (real energy matches reference; log carries the CostProfile line). Serialize; no other MPQC run concurrently.

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i \
  src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp \
  src/mpqc/chemistry/qc/lcao/expression/sequant_engine.h \
  src/mpqc/chemistry/qc/lcao/expression/sequant.cpp
git add -A && git commit -m "sequant: sequant:eval:dry_run predict-then-run cost-profile pre-flight"
```

---

## Self-review notes (author)

- **Spec coverage:** a1 -> Task 1; a2 -> Task 2; a3 -> Task 3; a4 -> Tasks 4-5; b -> Task 7; c -> Task 8; cross-repo repin -> Task 6. Data flow (moments -> IR -> cost_profile -> log -> real eval) is realized by Tasks 7+8. Error handling (dry-run non-fatal, k>4 clamp, no-moment fallback) is in the `inner_pow` closure (Task 7) and the try/catch (Task 8 Step 3).
- **Sequencing:** b (Task 7) is standalone but the plan orders SeQuant first (Tasks 1-5) because Task 4's API is what Task 8 calls and Task 6 gates the repin on the pushed SeQuant SHA; Task 7 (MPQC moment print) has no SeQuant dependency and MAY be done any time after Task 6, but is listed after it to keep one linear branch history. This matches the spec's "b -> a -> c" intent while respecting the cross-repo push/repin gate.
- **Known soft spots the implementer must resolve by reading code (flagged inline, not placeholders):** exact `tot_indices` construction (Task 1), `reset()` persistence semantics (Task 2), batched-loop threading model for the peak sink (Task 3), the concrete `EvalForest`/`LeafEvaluator` type spellings (Task 4), the CSV-holder that carries moments to CCk and the OSV per-pair reachability (Task 7), and the in-scope forest/config/leaf expressions in `process_for_evaluation` (Task 8). Each is named with the file/line to read.
