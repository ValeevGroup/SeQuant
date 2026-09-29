# Perf-first / peak-second objective (DenseTimeSpace) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a perf-first / peak-second single-term optimization objective
(`DenseTimeSpace` / `DenseTimeSpaceBatched`) that selects a factorization by
roofline-perf first and peak second, so the optimizer never prefers a
FLOPS-catastrophic factorization (the C60 4-PAO integral) merely because it is
sliceable below `peak_threshold`.

**Architecture:** Reuse the existing Pareto-frontier + roofline cost machinery
verbatim. Add a `bool perf_first` flag to `PeakModel` / `PeakBatchedModel`; when
set, the root-frontier selector sorts by `(flops, then peak)` instead of
`(peak, then flops)`. Rename `DensePeakSize{,Batched}` to `DenseSpaceTime{,Batched}`
(keeping the old names as deprecated aliases sharing the same enum values) and
add the two new perf-first values. Wire the new objective through
`single_term_opt`, `optimize.cpp`, and MPQC's `SeQuantEngine`.

**Tech Stack:** SeQuant C++20 (`SeQuant/core/optimize/`), Catch2 unit tests,
MPQC4 (`src/mpqc/.../sequant_engine.cpp`). SeQuant is the primary repo; MPQC
consumes it via `FETCHCONTENT_SOURCE_DIR_SEQUANT` (local edits are live) and via
the `MPQC_TRACKED_SEQUANT_TAG` pin (CI/clean builds; repin last).

## Global Constraints

- No en-dashes (U+2013) or non-breaking spaces (U+00A0) anywhere; ASCII hyphens
  only. A pre-commit hook rejects both. (This applies to source, comments, and
  commit messages.)
- No `Co-Authored-By:` trailers in commit messages.
- `DenseSpaceTime` / `DenseSpaceTimeBatched` (formerly `DensePeakSize*`) must
  keep byte-for-byte identical behavior; the old names remain as deprecated
  aliases with the same underlying enum values.
- `perf_first` defaults to `false`; every existing model construction and test
  must be unaffected.
- The frontier construction (`pareto_insert`, `relax`, `sliced_footprints`, the
  DP recurrence) does NOT change. `pareto_best` (test helper) does NOT change.
- Naming: `Dense{Primary}{Secondary}`, `Space` = peak/size, `Time` = perf.
- clang-format every changed C/C++ file before committing:
  `/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i <file>`.
- SeQuant unit test build/run (from a configured build dir, e.g.
  `cmake-build-release`): `cmake --build cmake-build-release --target unit_tests-sequant`
  then `./cmake-build-release/tests/unit/unit_tests-sequant "<test name>"`.
- SeQuant threads internally: never run more than one unit-test binary at once.

---

### Task 1: Rename the enum and add the two perf-first values

**Files:**
- Modify: `SeQuant/core/optimize/options.hpp:46-51` (the `ObjectiveFunction` enum)
- Modify: `SeQuant/core/optimize/options.hpp:19-45` (the doc comment above it)
- Test: `tests/unit/test_optimize.cpp` (a new static-assert-style `SECTION`)

**Interfaces:**
- Produces: `ObjectiveFunction::DenseSpaceTime` (== old `DensePeakSize`, value 2),
  `DenseSpaceTimeBatched` (== old `DensePeakSizeBatched`, value 3),
  `DenseTimeSpace` (new, value 4), `DenseTimeSpaceBatched` (new, value 5), plus
  deprecated aliases `DensePeakSize = DenseSpaceTime`,
  `DensePeakSizeBatched = DenseSpaceTimeBatched`.

- [ ] **Step 1: Write the failing test**

Add this `SECTION` inside the existing `TEST_CASE("optimize", "[optimize]")`
block in `tests/unit/test_optimize.cpp` (place it right after the opening of the
test case's body, near the other structural sections):

```cpp
    SECTION("ObjectiveFunction aliases share enum values") {
      using sequant::ObjectiveFunction;
      // The deprecated names are aliases (same underlying value), so all
      // existing `== DensePeakSize` guards keep catching the renamed constant.
      STATIC_REQUIRE(ObjectiveFunction::DensePeakSize ==
                     ObjectiveFunction::DenseSpaceTime);
      STATIC_REQUIRE(ObjectiveFunction::DensePeakSizeBatched ==
                     ObjectiveFunction::DenseSpaceTimeBatched);
      // The new perf-first values are distinct from the peak-first ones.
      STATIC_REQUIRE(ObjectiveFunction::DenseTimeSpace !=
                     ObjectiveFunction::DenseSpaceTime);
      STATIC_REQUIRE(ObjectiveFunction::DenseTimeSpaceBatched !=
                     ObjectiveFunction::DenseSpaceTimeBatched);
      // ABI: DenseSpaceTime keeps DensePeakSize's old numeric value (2).
      STATIC_REQUIRE(static_cast<int>(ObjectiveFunction::DenseSpaceTime) == 2);
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build cmake-build-release --target unit_tests-sequant`
Expected: COMPILE FAILURE (`DenseSpaceTime` / `DenseTimeSpace` are not members
of `ObjectiveFunction`).

- [ ] **Step 3: Implement the enum rename + new values**

Replace the enum at `SeQuant/core/optimize/options.hpp:46-51` with:

```cpp
enum class ObjectiveFunction {
  DenseFLOPs,
  DenseSize,
  /// Peak-first, perf-second: minimize peak memory, break ties by roofline
  /// perf. (Formerly `DensePeakSize`.)
  DenseSpaceTime,
  /// Batched variant of `DenseSpaceTime`. (Formerly `DensePeakSizeBatched`.)
  DenseSpaceTimeBatched,
  /// Perf-first, peak-second: minimize roofline perf, break ties by peak. This
  /// never prefers a FLOPS-catastrophic factorization for its sliceability.
  DenseTimeSpace,
  /// Batched variant of `DenseTimeSpace`.
  DenseTimeSpaceBatched,
  /// Deprecated aliases (same underlying values as the renamed constants above,
  /// placed AFTER them so `DenseSpaceTime` keeps `DensePeakSize`'s old value).
  /// Kept so existing code and JSON inputs ("dense_peak_size") keep compiling.
  DensePeakSize = DenseSpaceTime,
  DensePeakSizeBatched = DenseSpaceTimeBatched
};
```

Then, in the doc comment above (`options.hpp:19-45`), update the two bullet
lines that describe `DensePeakSize` / `DensePeakSizeBatched` to name
`DenseSpaceTime` / `DenseSpaceTimeBatched` (noting the old names are deprecated
aliases), and add two bullets describing `DenseTimeSpace` /
`DenseTimeSpaceBatched` as the perf-first / peak-second duals. Keep it factual;
do not restate the whole design.

- [ ] **Step 4: Run test to verify it passes**

Run: `cmake --build cmake-build-release --target unit_tests-sequant && ./cmake-build-release/tests/unit/unit_tests-sequant "optimize"`
Expected: PASS (the new `SECTION` compiles and its `STATIC_REQUIRE`s hold).

- [ ] **Step 5: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i SeQuant/core/optimize/options.hpp tests/unit/test_optimize.cpp
git add SeQuant/core/optimize/options.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: rename DensePeakSize -> DenseSpaceTime, add DenseTimeSpace values"
```

---

### Task 2: Add `perf_first` and the perf-first selectors

**Files:**
- Modify: `SeQuant/core/optimize/cost_model.hpp:303-326` (add `perf_first` to `PeakModel`)
- Modify: `SeQuant/core/optimize/cost_model.hpp:437-453` (perf-first branch in `PeakModel::reconstruct`)
- Modify: `SeQuant/core/optimize/cost_model.hpp:486-524` (add `perf_first` to `PeakBatchedModel`)
- Modify: `SeQuant/core/optimize/cost_model.hpp:692-719` (perf-first branch in `PeakBatchedModel::select_root`)

**Interfaces:**
- Consumes: nothing new.
- Produces: `PeakModel::perf_first` and `PeakBatchedModel::perf_first` (both
  `bool`, default `false`), added as the LAST data member of each struct so
  existing positional aggregate initializers are unaffected. When `perf_first`
  is set, `PeakModel::reconstruct` and `PeakBatchedModel::select_root` select by
  `(flops, then peak)`.

- [ ] **Step 1: Add the `perf_first` member to both models**

In `PeakModel` (after the `peak_flops_tolerance` member at
`cost_model.hpp:326`, before the `struct FrontPoint` declaration):

```cpp
  /// Perf-first / peak-second selection: when true, `reconstruct` selects the
  /// root-frontier point by (flops, then peak) instead of (peak, then flops),
  /// bypassing `peak_flops_tolerance`. Default false = peak-first (unchanged).
  bool perf_first = false;
```

In `PeakBatchedModel` (after the `numeric_size` member at
`cost_model.hpp:524`, before the `struct BFrontPoint` declaration):

```cpp
  /// Perf-first / peak-second selection: when true, `select_root` selects the
  /// root-frontier point by (flops, then peak) and does NOT consult
  /// `peak_threshold` as a feasibility gate (it can no longer force a
  /// FLOPS-catastrophic factorization for its sliceability). Default false =
  /// peak-first threshold-gated selection (unchanged).
  bool perf_first = false;
```

- [ ] **Step 2: Write the failing test**

Extend the existing `TEST_CASE("threshold gates batching: ...",
"[optimize][threshold]")` in `tests/unit/test_optimize.cpp` (ends at line ~2264,
right after `CHECK_FALSE(integral_at(2, 1.0));`). Add a perf-first variant of
`integral_at` and the threshold-insensitivity checks. Insert BEFORE the closing
`}` of that `TEST_CASE`:

```cpp
  // Perf-first (DenseTimeSpaceBatched): min-flops is primary and peak_threshold
  // is NOT a feasibility gate, so the persistent Kappa-free integral (the
  // min-(replay-weighted)-flops choice) is formed REGARDLESS of the threshold.
  // Contrast the peak-first checks above: at a near-zero threshold peak-first
  // fell back to min-peak (CHECK_FALSE), whereas perf-first still forms it.
  auto integral_at_perf = [&](std::size_t Kb, double peak_threshold) -> bool {
    std::function<std::size_t(Index const&)> bts = [Kb](Index const&) {
      return Kb;
    };
    CostParams cost{is_t, 100.0, 0.0};
    cost.peak_threshold = peak_threshold;
    auto res = opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
        prod, idxsz, false, cost, is_batch, bts, ip);
    return forms_integral(res);
  };
  // Threshold-insensitive: forms the integral at +inf AND at near-zero.
  CHECK(integral_at_perf(/*Kb=*/2, inf));
  CHECK(integral_at_perf(/*Kb=*/2, /*peak_threshold(bytes)=*/1.0));
```

- [ ] **Step 3: Run test to verify it fails**

Run: `cmake --build cmake-build-release --target unit_tests-sequant && ./cmake-build-release/tests/unit/unit_tests-sequant "threshold gates batching*"`
Expected: at Step 2 the build currently FAILS to compile (no
`DenseTimeSpaceBatched` arm in `single_term_opt` yet -> `static_assert` fires or
the arm is unreachable). If it compiles because Task 3 is not yet done, the
second `CHECK(integral_at_perf(2, 1.0))` FAILS (perf-first not yet implemented;
the value routes to the peak-first path and falls back to min-peak). Either way
this task is not green until Steps 4 AND Task 3 are done. (If ordering makes the
build fail here, that is an expected RED; proceed to Step 4, then complete Task
3, then re-run.)

- [ ] **Step 4: Implement the perf-first branch in `PeakBatchedModel::select_root`**

At the TOP of `select_root` (`cost_model.hpp:692`), immediately after the
`auto peak_bytes = ...;` lambda (line 698) and before `int best = -1;`, insert:

```cpp
    if (perf_first) {
      // Perf-first / peak-second: min flops, ties by lower peak. The frontier
      // keeps one min-peak point per distinct flops value (pareto_insert prunes
      // equal-flops higher-peak points), so this both picks the cheapest
      // factorization and takes its fully-sliced (min-peak) realization.
      // peak_threshold is deliberately NOT consulted here.
      int pbest = 0;
      for (int i = 1; i < static_cast<int>(rootf.size()); ++i)
        if (rootf[i].flops < rootf[pbest].flops ||
            (rootf[i].flops == rootf[pbest].flops &&
             rootf[i].peak < rootf[pbest].peak))
          pbest = i;
      return pbest;
    }
```

- [ ] **Step 5: Implement the perf-first branch in `PeakModel::reconstruct`**

At the TOP of `reconstruct` (`cost_model.hpp:437`), immediately after
`size_t const full = st.size() - 1;` (line 439) and before the
`// epsilon-tolerant selection` comment, insert:

```cpp
    if (perf_first) {
      // Perf-first / peak-second (non-batched): min flops, ties by lower peak,
      // bypassing the peak_flops_tolerance epsilon band (a peak-first knob).
      auto const& rootf = st[full];
      int pbest = 0;
      for (int i = 1; i < static_cast<int>(rootf.size()); ++i)
        if (rootf[i].flops < rootf[pbest].flops ||
            (rootf[i].flops == rootf[pbest].flops &&
             rootf[i].peak < rootf[pbest].peak))
          pbest = i;
      // Reuse the existing back-pointer walk with the chosen root index.
      std::function<EvalSequence(size_t, int)> pbuild =
          [&](size_t n, int idx) -> EvalSequence {
        if (std::popcount(n) == 1)
          return EvalSequence{static_cast<int>(std::countr_zero(n))};
        FrontPoint const& fp = st[n][idx];
        size_t const fs = fp.lp_first ? fp.lp : fp.rp;
        int const fi = fp.lp_first ? fp.lp_idx : fp.rp_idx;
        size_t const ss = fp.lp_first ? fp.rp : fp.lp;
        int const si = fp.lp_first ? fp.rp_idx : fp.lp_idx;
        EvalSequence s = pbuild(fs, fi);
        EvalSequence b = pbuild(ss, si);
        s.insert(s.end(), b.begin(), b.end());
        s.push_back(-1);
        return s;
      };
      return pbuild(full, pbest);
    }
```

- [ ] **Step 6: Run the tests to verify they pass (after Task 3 is also done)**

This task's test cannot go green until Task 3 wires `DenseTimeSpaceBatched`
through `single_term_opt`. Complete Task 3, then run:
`./cmake-build-release/tests/unit/unit_tests-sequant "threshold gates batching*"`
Expected: PASS (both `integral_at_perf` checks true).

- [ ] **Step 7: clang-format and commit** (after Task 3 compiles)

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i SeQuant/core/optimize/cost_model.hpp
git add SeQuant/core/optimize/cost_model.hpp
git commit -m "optimize: add perf_first (flops, then peak) selectors to peak models"
```

---

### Task 3: Wire the new objectives through `single_term_opt` and `optimize.cpp`

**Files:**
- Modify: `SeQuant/core/optimize/single_term.hpp:98-183` (the `if constexpr` chain)
- Modify: `SeQuant/core/optimize/optimize.cpp:102-126` (the dispatch lambda)

**Interfaces:**
- Consumes: `PeakModel::perf_first`, `PeakBatchedModel::perf_first` (Task 2);
  the four enum values (Task 1).
- Produces: `single_term_opt<DenseTimeSpace>` and
  `single_term_opt<DenseTimeSpaceBatched>` route to the peak models with
  `perf_first = true`; `optimize()` dispatches both new objectives.

- [ ] **Step 1: Broaden the non-batched peak arm in `single_term.hpp`**

Change the arm header at `single_term.hpp:98` from:

```cpp
  if constexpr (Metric == ObjectiveFunction::DensePeakSize) {
```
to:
```cpp
  if constexpr (Metric == ObjectiveFunction::DenseSpaceTime ||
                Metric == ObjectiveFunction::DenseTimeSpace) {
```

Then replace the `return run_single_term_opt(PeakModel{...}, network, tidxs);`
statement (lines 109-114) with a named local that sets `perf_first`:

```cpp
    PeakModel model{idxsz,
                    inner_pow,
                    is_volatile_leaf,
                    volatile_weight,
                    roofline.machine_balance,
                    roofline.fast_mem_elems,
                    roofline.block_tiles,
                    roofline.block_prefactor,
                    peak_flops_tolerance};
    model.perf_first = (Metric == ObjectiveFunction::DenseTimeSpace);
    return run_single_term_opt(model, network, tidxs);
```

- [ ] **Step 2: Broaden the batched peak arm in `single_term.hpp`**

Change the arm header at `single_term.hpp:115` from:

```cpp
  } else if constexpr (Metric == ObjectiveFunction::DensePeakSizeBatched) {
```
to:
```cpp
  } else if constexpr (Metric == ObjectiveFunction::DenseSpaceTimeBatched ||
                       Metric == ObjectiveFunction::DenseTimeSpaceBatched) {
```

Then, immediately AFTER the `PeakBatchedModel model{...};` aggregate
initialization (i.e. after line 135's `cost.peak_threshold};`), insert:

```cpp
    model.perf_first = (Metric == ObjectiveFunction::DenseTimeSpaceBatched);
```

- [ ] **Step 3: Update the trailing `static_assert` message**

In the final `else` branch (`single_term.hpp:166-168`), update the message text
to list the current objective names (informational only):

```cpp
    static_assert(Metric == ObjectiveFunction::DenseSize,
                  "Only DenseFLOPs, DenseSize, DenseSpaceTime, "
                  "DenseSpaceTimeBatched, DenseTimeSpace, and "
                  "DenseTimeSpaceBatched ObjectiveFunction supported.");
```

- [ ] **Step 4: Add dispatch arms in `optimize.cpp`**

Replace the body of the `run` lambda at `optimize.cpp:102-126` (from the first
`if (opts.objective_function == ObjectiveFunction::DenseFLOPs)` through the final
batched `return`) with:

```cpp
  auto run = [&]() -> ExprPtr {
    if (opts.objective_function == ObjectiveFunction::DenseFLOPs)
      return opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.batch_policy.is_batchable_index,
          opts.batch_policy.batch_target_size, opts.inner_pow);
    if (opts.objective_function == ObjectiveFunction::DenseSize)
      return opt::single_term_opt<ObjectiveFunction::DenseSize>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.batch_policy.is_batchable_index,
          opts.batch_policy.batch_target_size, opts.inner_pow);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTime)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTime>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.batch_policy.is_batchable_index,
          opts.batch_policy.batch_target_size, opts.inner_pow);
    if (opts.objective_function == ObjectiveFunction::DenseTimeSpace)
      return opt::single_term_opt<ObjectiveFunction::DenseTimeSpace>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.batch_policy.is_batchable_index,
          opts.batch_policy.batch_target_size, opts.inner_pow);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTimeBatched)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTimeBatched>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.batch_policy.is_batchable_index,
          opts.batch_policy.batch_target_size, opts.inner_pow,
          opts.batch_policy.persistent_only,
          opts.term_batch_axes ? &node_axes : nullptr);
    SEQUANT_ASSERT(opts.objective_function ==
                   ObjectiveFunction::DenseTimeSpaceBatched);
    return opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
        prod, opts.idx_to_extent, subnet_cse, cost,
        opts.batch_policy.is_batchable_index,
        opts.batch_policy.batch_target_size, opts.inner_pow,
        opts.batch_policy.persistent_only,
        opts.term_batch_axes ? &node_axes : nullptr);
  };
```

(Note: `DenseSpaceTime` == `DensePeakSize` and `DenseSpaceTimeBatched` ==
`DensePeakSizeBatched` as enum values, so these arms also cover the deprecated
names.)

- [ ] **Step 5: Build and run the full optimize + threshold tests**

Run: `cmake --build cmake-build-release --target unit_tests-sequant && ./cmake-build-release/tests/unit/unit_tests-sequant "[optimize]"`
Expected: PASS, including `optimize` (Task 1 section), `threshold gates
batching*` (Task 2 perf-first checks), and all pre-existing `[optimize]` cases
(regression: peak-first paths unchanged).

- [ ] **Step 6: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i SeQuant/core/optimize/single_term.hpp SeQuant/core/optimize/optimize.cpp
git add SeQuant/core/optimize/single_term.hpp SeQuant/core/optimize/optimize.cpp
git commit -m "optimize: route DenseTimeSpace{,Batched} to perf-first peak models"
```

---

### Task 4: Perf-first characterization test (perf-first matches flops-optimal factorization)

**Files:**
- Modify: `tests/unit/test_optimize.cpp` (extend the `[optimize][osv]` "PPL term:
  form-W vs fold-t" TEST_CASE at line ~1943, which already has `idxsz`,
  `is_batch`, `bts`, `ip`, `is_t`, `wflops`, `rep`, `prod`, `flops_w`)

**Interfaces:**
- Consumes: `DenseTimeSpaceBatched` (Task 3), the in-test `wflops` / `rep` /
  `flops_w` helpers already defined in that TEST_CASE.

- [ ] **Step 1: Write the failing test**

Immediately after the existing `CHECK(pbat_w_vw1 > flops_w);` line (the last
assertion of the "form-W vs fold-t" TEST_CASE, ~line 1941), insert:

```cpp
  // Perf-first (DenseTimeSpaceBatched) is flops-primary, so it forms the
  // flop-optimal W REGARDLESS of peak_flops_tolerance and regardless of a tight
  // peak_threshold: its replay-weighted flop count equals the FLOPs optimum.
  // (Peak-first at strict tolerance did NOT: peak_w_strict > flops_w above.)
  double tbat_w = rep(
      L"TimeSpaceB/vw100:",
      opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
          prod, idxsz, false, CostParams{is_t, 100.0, 0.0}, is_batch, bts, ip));
  CHECK(tbat_w == flops_w);
  // Even with a near-zero peak budget, perf-first stays flop-optimal (threshold
  // is not a feasibility gate), unlike peak-first's min-peak fallback.
  CostParams tight{is_t, 100.0, 0.0};
  tight.peak_threshold = 1.0;
  double tbat_w_tight =
      rep(L"TimeSpaceB/tight:",
          opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
              prod, idxsz, false, tight, is_batch, bts, ip));
  CHECK(tbat_w_tight == flops_w);
```

- [ ] **Step 2: Run test to verify it fails, then passes**

Run: `./cmake-build-release/tests/unit/unit_tests-sequant "*form-W vs fold-t*"`
(the TEST_CASE name may differ; use `--list-tests "[osv]"` to find the exact
name, then quote it). Expected BEFORE Task 3: compile/route failure. AFTER Tasks
2-3: PASS (both `CHECK(... == flops_w)` hold).

If `tbat_w != flops_w` unexpectedly, do NOT weaken the assertion. Investigate:
dump `render_tree` of the perf-first result vs the `FLOPs` result; a mismatch
means the perf-first selector is not reaching the flop-optimal factorization
(re-check Task 2 Step 4/5 and that `machine_balance == 0` here so
`roofline_op_cost` reduces to raw flops).

- [ ] **Step 3: clang-format and commit**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i tests/unit/test_optimize.cpp
git add tests/unit/test_optimize.cpp
git commit -m "optimize: test DenseTimeSpaceBatched matches the flops-optimal factorization"
```

---

### Task 5: MPQC SeQuantEngine keyword + auto-route

**Files:**
- Modify: `src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp:51-69` (parser)
- Modify: `src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp:179-196` (auto-route)
- Modify: `src/mpqc/chemistry/qc/lcao/expression/sequant_engine.h:87` (keyword doc)

(All paths in the mpqc4 repo: `/Users/efv/code/mpqc4`.)

**Interfaces:**
- Consumes: `sequant::ObjectiveFunction::DenseTimeSpace{,Batched}` (Task 1),
  live in local MPQC builds via `FETCHCONTENT_SOURCE_DIR_SEQUANT`.
- Produces: `optimize:objective_function = "dense_time_space"` selects the
  perf-first objective and auto-routes to `DenseTimeSpaceBatched` when a batch
  policy is active. `"dense_peak_size"` and a new `"dense_space_time"` both map
  to the peak-first `DenseSpaceTime`.

- [ ] **Step 1: Extend the parser**

In `sequant_engine.cpp`, replace the parse branches (lines 58-68) with:

```cpp
    if (of == "dense_flops")
      objective_function_ = sequant::ObjectiveFunction::DenseFLOPs;
    else if (of == "dense_size")
      objective_function_ = sequant::ObjectiveFunction::DenseSize;
    else if (of == "dense_space_time" || of == "dense_peak_size")
      objective_function_ = sequant::ObjectiveFunction::DenseSpaceTime;
    else if (of == "dense_time_space")
      objective_function_ = sequant::ObjectiveFunction::DenseTimeSpace;
    else
      throw InputError(
          "invalid \"sequant:optimize:objective_function\" (must be "
          "\"dense_flops\", \"dense_size\", \"dense_space_time\" (alias "
          "\"dense_peak_size\"), or \"dense_time_space\")",
          __FILE__, __LINE__, "optimize:objective_function", of.c_str());
```

(The default-echo string at lines 51-57 stays as-is: the default is
`DenseFLOPs`, and the else-arm's `"dense_peak_size"` fallback remains a valid
alias for `DenseSpaceTime`.)

- [ ] **Step 2: Extend the auto-route**

Replace the auto-route block (lines 182-196) with:

```cpp
  if (objective_function_ == sequant::ObjectiveFunction::DenseSpaceTime ||
      objective_function_ == sequant::ObjectiveFunction::DenseTimeSpace) {
    if (!opts.idx_to_extent)
      throw ProgrammingError(
          "peak objectives require index extents (idx_to_extent on the "
          "EvalContext)",
          __FILE__, __LINE__);
    if (opts.CSE.subnet)
      throw InputError(
          "sequant:optimize:objective_function peak objectives are "
          "incompatible with sequant:optimize:cse:subnet=true",
          __FILE__, __LINE__, "optimize:cse:subnet", "true");
    if (opts.batch_policy.is_batchable_index)
      opts.objective_function =
          (objective_function_ == sequant::ObjectiveFunction::DenseTimeSpace)
              ? sequant::ObjectiveFunction::DenseTimeSpaceBatched
              : sequant::ObjectiveFunction::DenseSpaceTimeBatched;
  }
```

- [ ] **Step 3: Update the keyword doc comment**

In `sequant_engine.h:87`, change the `dense_flops, dense_size, or
dense_peak_size` list to `dense_flops, dense_size, dense_space_time (alias
dense_peak_size; peak-first), or dense_time_space (perf-first)`, keeping the
rest of that Doxygen row intact and within the existing style.

- [ ] **Step 4: Build MPQC and sanity-check the parser end-to-end**

Build the mpqc target against the locally-edited SeQuant (release build dir):

Run: `cmake --build /Users/efv/code/mpqc4/cmake-build-release --target mpqc 2>&1 | tail -30`
Expected: links cleanly (the enum change + wiring compile in MPQC too).

There is no dedicated parser unit test; verify via an existing CCk validation
input. Pick the smallest CCk JSON that sets `optimize:objective_function`, copy
it, change that value to `dense_time_space`, and run it once (serialized, single
rank) to confirm it parses and completes without an InputError:

Run: `MAD_NUM_THREADS=6 mpiexec -n 1 /Users/efv/code/mpqc4/cmake-build-release/src/bin/mpqc/mpqc <copied-input>.json 2>&1 | tail -40`
Expected: no `invalid "sequant:optimize:objective_function"` error; the job runs
the CCk solver. (A full numeric validation is not required here; Task 6 does the
C60 dry-run. This step only confirms the keyword plumbs through.)

- [ ] **Step 5: clang-format and commit (mpqc4 repo)**

```bash
cd /Users/efv/code/mpqc4
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp src/mpqc/chemistry/qc/lcao/expression/sequant_engine.h
git add src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp src/mpqc/chemistry/qc/lcao/expression/sequant_engine.h
git commit -m "SeQuantEngine: accept dense_time_space objective, auto-route to batched"
```

---

### Task 6: DryRun harness validation on C60 (perf-first vs peak-first)

**Files:**
- Modify: `tests/unit/test_eval_dryrun.cpp` (the `[dryrun-eval]` / `[dryrun-df]`
  reproducers; currently uncommitted WIP on branch
  `feature/eval-predicted-peak-trace`, driving `DensePeakSizeBatched` on the
  faithful C60 config)

**Interfaces:**
- Consumes: `DenseTimeSpaceBatched` (Task 3). The harness already builds the real
  C60 residual, the real config (occ=120, pno=42, osv=310, mu-tilde=1800,
  aux=4320, pao_ts=256, aux_ts=72, peak_threshold=40 GB), and the top-down
  realized-peak walk that reports free-mu-tilde escapees and node dumps.

This is a characterization/observation task (not a red-green unit test): its job
is to EMPIRICALLY confirm the fix before any cluster run. Do not weaken it into
a tautology; record the actual numbers.

- [ ] **Step 1: Add a perf-first pass to the harness**

In the C60 dryrun reproducer that currently runs the batched optimize with
`DensePeakSizeBatched` (via the `OptimizeOptions` it builds), add a SECOND run
of the identical config with `opts.objective_function =
ObjectiveFunction::DenseTimeSpace` (the engine/opt path auto-routes to
`DenseTimeSpaceBatched` because `is_batchable_index` is set; if the harness calls
`opt::single_term_opt<...>` directly rather than through `optimize()`, use
`ObjectiveFunction::DenseTimeSpaceBatched` explicitly). Keep the existing
peak-first run so the two can be compared side by side.

For each of the two runs, emit (the harness already has the machinery):
- whether any node is a four-free-PAO (`mu~^4`) contraction (the 4-PAO integral);
- the LARGEST realized free-mu-tilde intermediate (GB), after ancestor slicing;
- the count of intermediates whose free mu-tilde escaped slicing.

- [ ] **Step 2: Build and run the harness (serialized; SeQuant threads internally)**

Run: `cmake --build cmake-build-release --target unit_tests-sequant && ./cmake-build-release/tests/unit/unit_tests-sequant "[dryrun-eval]" "[dryrun-df]" 2>&1 | tail -60`
Expected: the perf-first pass reports NO `mu~^4` node and a largest realized
free-mu-tilde intermediate substantially below the peak-first pass's (and below
the ~40 GB budget if the fix is complete). Record both passes' numbers in the
commit message and in the ledger.

- [ ] **Step 3: Interpret the result and gate the cluster decision**

- If perf-first shows NO 4-PAO and the largest realized intermediate is well
  under the peak-first result (and under available node memory): the design is
  validated in-harness. STOP here; the C60 cluster run is a separate,
  user-initiated step (do not launch it).
- If perf-first STILL shows a large free-mu-tilde giant: this is a genuine
  finding, NOT a test to loosen. It means the min-flops factorization itself
  carries an irreducible large peak (the dual concern noted in the spec). Record
  the exact giant (indices + realized GB) and surface it as a design question
  for the user (a `peak_threshold`-aware secondary among near-min-flops points
  may then be warranted). Use `superpowers:systematic-debugging` before
  proposing any follow-up change.

- [ ] **Step 4: clang-format and commit the harness apparatus**

```bash
/opt/homebrew/opt/llvm@17/bin/clang-format --style=file -i tests/unit/test_eval_dryrun.cpp
git add tests/unit/test_eval_dryrun.cpp
git commit -m "eval: dry-run C60 perf-first vs peak-first factorization comparison"
```

---

### Task 7: Repin MPQC to the updated SeQuant (CI/clean builds)

**Files:**
- Modify: `/Users/efv/code/mpqc4/external/versions.cmake` (`MPQC_TRACKED_SEQUANT_TAG`
  and its `PREVIOUS_*` companion)

**Interfaces:**
- Consumes: the merged/pushed SeQuant commits from Tasks 1-4 and 6.

This task is LAST and is gated on the SeQuant changes being pushed to a commit
CI can fetch. Local MPQC builds already see the edits via
`FETCHCONTENT_SOURCE_DIR_SEQUANT`; the pin only matters for CI/clean builds.

- [ ] **Step 1: Push the SeQuant branch and get its head SHA**

```bash
cd /Users/efv/code/SeQuant
git push   # push feature/eval-predicted-peak-trace (or the agreed branch)
git rev-parse HEAD
```

- [ ] **Step 2: Update the pin**

In `/Users/efv/code/mpqc4/external/versions.cmake`, set
`MPQC_TRACKED_SEQUANT_TAG` to the new SHA and move the old value into the
corresponding `MPQC_TRACKED_PREVIOUS_SEQUANT_TAG` (match the file's existing
CURRENT/PREVIOUS convention exactly; read the surrounding lines first).

- [ ] **Step 3: Commit the pin on its own**

```bash
cd /Users/efv/code/mpqc4
git add external/versions.cmake
git commit -m "external: repin SeQuant to <new-sha7> (DenseTimeSpace objective)"
```

(Keep this pin bump in its own commit, separate from the Task 5 wiring, per the
repo's commit rules.)

---

## Self-Review

**1. Spec coverage:**
- Enum rename + aliases + two new values -> Task 1. (spec "Objective naming")
- `perf_first` field + `select_root` perf-first + `PeakModel::reconstruct`
  perf-first -> Task 2. (spec "Why the selector is a plain lexicographic swap")
- `single_term_opt` arms + `optimize.cpp` dispatch + `perf_first` set from
  `Metric` -> Task 3. (spec "Data flow")
- MPQC keyword `dense_time_space` + auto-route + doc -> Task 5. (spec "MPQC wiring")
- C60 DryRun validation, no-4-PAO, realized-peak comparison -> Task 6. (spec
  "Validation")
- `DenseSpaceTime*` unchanged / alias parse guard -> Task 1 (static-assert) +
  Task 3 Step 5 (regression run of full `[optimize]`).
- Repin -> Task 7. (spec "Global constraints" cross-repo note)
- `pareto_best` unchanged, frontier unchanged -> stated in Global Constraints;
  no task touches them (correct).

Gap check: the spec's unit-test bullet "assert `pareto_best`'s perf-first branch
returns the min-flops point" is intentionally DROPPED -- analysis showed
`pareto_best` is a test helper that must stay peak-min and is not on the
perf-first path, so there is no perf-first `pareto_best` branch to test. The
non-batched perf-first path is instead covered structurally by Task 2 Step 5 +
the full `[optimize]` regression run. This is a correction to the spec's
testing list, not a coverage gap.

**2. Placeholder scan:** No "TBD"/"handle edge cases"/prose-only steps; every
code step shows the exact code. Task 4 Step 2 flags that the exact TEST_CASE
name must be looked up with `--list-tests` (a concrete command), not a
placeholder.

**3. Type consistency:** `perf_first` is `bool` in both models (Task 2) and set
identically in both `single_term.hpp` arms (Task 3). `DenseTimeSpaceBatched` /
`DenseTimeSpace` spelled consistently across Tasks 1/3/4/5/6. `CostParams`
aggregate order (`is_volatile_leaf, volatile_weight, footprint_weight,
peak_flops_tolerance, roofline, accumulation_factor, peak_threshold`) matches
`options.hpp:93-117` and the `.peak_threshold =` field-assignment idiom used in
the existing `[optimize][threshold]` test. `PeakModel` / `PeakBatchedModel`
aggregate-init field orders in Task 3 match `cost_model.hpp` (perf_first added
LAST in each, so positional inits are unaffected).
