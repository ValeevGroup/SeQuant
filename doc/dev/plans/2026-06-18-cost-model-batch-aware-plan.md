# Batch-aware cost model - Phase 1: peak-memory objective

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a peak-memory single-term-optimization objective (`ObjectiveFunction::DensePeakSize`) and prove its dynamic program finds the true minimum-peak contraction tree by validating against an independent brute-force oracle.

**Architecture:** A new DP variant (`single_term_opt_peak_impl`) computes, per tensor subset, the minimum peak memory to evaluate that subtree via the pebbling recurrence (which chooses the cheaper of the two child evaluation orders). Result footprints `S[subset]` are precomputed once and shared with an independent test oracle that enumerates all contraction sequences and simulates memory. The existing additive DP (`single_term_opt_impl`) is untouched, so `DenseFLOPs`/`DenseSize` behavior cannot regress.

**Tech Stack:** C++20, SeQuant `core/optimize`, Catch2 unit tests (`tests/unit/test_optimize.cpp`).

## Global Constraints

- Style: Google-based `.clang-format`, 80-column, 2-space indent, no tabs. Run `clang-format --style=file -i <file>` before each commit.
- No en-dashes (U+2013) or non-breaking spaces anywhere (pre-commit hook `forbid-en-dashes` rejects them) - use ASCII hyphen.
- No `Co-Authored-By` trailers in commit messages.
- Branch: `feature/cost-model-batch-aware` (already created, based on `feature/optimize-options-terminology`).
- Phase 1 is SeQuant-only; no mpqc changes. `DensePeakSizeBatched`, the model-is-evaluator unification, `OptimizeOptions` -> `CostModel` refactor, and mpqc wiring are later phases (see spec `doc/dev/specs/2026-06-18-cost-model-batch-aware-design.md`).
- Phase 1 does NOT support `subnet_cse` for the peak objective (spec risk item); the peak path asserts `subnet_cse == false`.

## File structure

- `core/optimize/options.hpp` - add `DensePeakSize` to the `ObjectiveFunction` enum.
- `core/optimize/single_term.hpp` - add `subset_footprints()` (shared primitive), `single_term_opt_peak_impl()` (the peak DP), and dispatch in the detail `single_term_opt<Metric>`.
- `core/optimize/optimize.cpp` - dispatch `DensePeakSize` in `opt_pure_product`.
- `tests/unit/test_optimize.cpp` - add the brute-force oracle helper and new SECTIONs.

Build/run for every test step:
```bash
cmake --build build-test --target unit_tests-sequant -j
./build-test/tests/unit/unit_tests-sequant "[optimize]"
```

---

### Task 1: `subset_footprints` shared primitive

**Files:**
- Modify: `core/optimize/single_term.hpp` (add helper near `footprint_counter`, ~line 116)
- Test: `tests/unit/test_optimize.cpp` (new SECTION inside `TEST_CASE("optimize")`)

**Interfaces:**
- Produces: `std::vector<double> sequant::opt::detail::subset_footprints(TensorNetwork const& network, meta::range_of<Index> auto const& tidxs, IdxToSz&& idxsz)` - returns `S`, sized `2^nt`, where `S[n]` is the product of extents of the open indices of subset `n` (0 for the empty subset and any scalar result). Reuses `init_results` to obtain per-subset open indices.

- [ ] **Step 1: Write the failing test**

Add inside `TEST_CASE("optimize", "[optimize]")` after the existing sections:

```cpp
SECTION("subset_footprints") {
  using namespace sequant;
  // i,j occ (size 2); a,b virt (size 4). Tensors: T0=g{a;i}, T1=g{b;j}.
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto g1 = parse_expr(L"g{a1;i1}", Symmetry::nonsymm);
  auto g2 = parse_expr(L"g{a2;i2}", Symmetry::nonsymm);
  TensorNetwork net{std::vector<ExprPtr>{g1, g2}};
  container::svector<Index> targets;  // fully contracted-to-scalar? keep a1,a2,i1,i2 open
  auto S = opt::detail::subset_footprints(net, targets, idxsz);
  REQUIRE(S.size() == 4u);
  REQUIRE(S[0] == 0.0);                 // empty subset
  // singleton {T0}: open indices a1,i1 -> 4*2 = 8
  REQUIRE(S[0b01] == Catch::Approx(8.0));
  REQUIRE(S[0b10] == Catch::Approx(8.0));
  // full {T0,T1}: open a1,i1,a2,i2 -> 4*2*4*2 = 64
  REQUIRE(S[0b11] == Catch::Approx(64.0));
}
```

- [ ] **Step 2: Run to verify it fails**

Run the build/run commands above.
Expected: compile error - `subset_footprints` not declared.

- [ ] **Step 3: Implement `subset_footprints`**

Add to `core/optimize/single_term.hpp` in `namespace sequant::opt::detail`, after `footprint_counter` (~line 116). It mirrors `init_results`'s open-index computation:

```cpp
/// \brief Footprint (dense element count) of every subset's result tensor.
///
/// `S[n]` is the product of extents of the open indices of subset `n`
/// (those remaining after contracting the tensors in `n`, given `tidxs` as the
/// final target indices). `S[0]` (empty subset) and any scalar result are 0.
/// Shared by the peak DP and its tests so both agree on per-subset sizes.
template <typename TIdxs, typename IdxToSz>
container::vector<double> subset_footprints(TensorNetwork const& network,
                                            TIdxs const& tidxs,
                                            IdxToSz&& idxsz) {
  container::vector<OptRes> results((size_t{1} << network.tensors().size()));
  init_results(network, tidxs, results);
  auto fp = footprint_counter(std::forward<IdxToSz>(idxsz));
  container::vector<double> S(results.size(), 0.0);
  for (size_t n = 0; n < results.size(); ++n)
    S[n] = (n == 0) ? 0.0 : fp(results[n].indices);
  return S;
}
```

- [ ] **Step 4: Run to verify it passes**

Run the build/run commands. Expected: PASS.

- [ ] **Step 5: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/single_term.hpp
git add core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: add subset_footprints primitive for peak costing"
```

---

### Task 2: brute-force peak oracle (test helper)

**Files:**
- Test: `tests/unit/test_optimize.cpp` (free helper above `TEST_CASE`, plus a SECTION)

**Interfaces:**
- Produces (test-local): `double brute_force_min_peak(std::vector<double> const& S, size_t nt)` - the minimum over all valid pairwise contraction sequences of the sequence's peak memory, where a step that forms subset `m` from live subsets `a,b` has instantaneous memory `S[m] + sum of S over all currently-live subsets`. Independent of the DP recurrence (simulates memory directly), so it validates the DP end-to-end.

- [ ] **Step 1: Write the failing test**

Add a free function above `TEST_CASE("optimize")` in `test_optimize.cpp`:

```cpp
// Minimum peak memory over ALL pairwise contraction sequences of `nt` leaves,
// using S[subset] (subset_footprints) for sizes. Independent oracle for the
// peak DP: enumerates schedules and simulates memory; no recurrence assumed.
static double brute_force_min_peak(std::vector<double> const& S, size_t nt) {
  double const full = static_cast<double>((size_t{1} << nt) - 1);
  // `live` is the set of subset-masks currently materialized (a partition of
  // the full set). Recurse over every pair to merge.
  std::function<double(std::vector<size_t>)> rec =
      [&](std::vector<size_t> live) -> double {
    if (live.size() == 1) return 0.0;  // nothing more to compute
    double live_sum = 0.0;
    for (auto m : live) live_sum += S[m];
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i < live.size(); ++i)
      for (size_t j = i + 1; j < live.size(); ++j) {
        size_t merged = live[i] | live[j];
        // instantaneous peak when forming `merged`: all live results plus the
        // new result momentarily co-resident.
        double step_peak = live_sum + S[merged];
        std::vector<size_t> next;
        next.reserve(live.size() - 1);
        for (size_t k = 0; k < live.size(); ++k)
          if (k != i && k != j) next.push_back(live[k]);
        next.push_back(merged);
        best = std::min(best, std::max(step_peak, rec(next)));
      }
    return best;
  };
  std::vector<size_t> leaves;
  for (size_t b = 0; b < nt; ++b) leaves.push_back(size_t{1} << b);
  (void)full;
  return rec(std::move(leaves));
}

TEST_CASE("brute_force_min_peak oracle", "[optimize]") {
  // 3 leaves, hand-checkable sizes by subset mask:
  //   S[001]=2 S[010]=2 S[100]=2 (leaves)
  //   S[011]=1 S[101]=1 S[110]=8 S[111]=1 (pair/full results)
  std::vector<double> S(8, 0.0);
  S[0b001] = S[0b010] = S[0b100] = 2.0;
  S[0b011] = 1.0; S[0b101] = 1.0; S[0b110] = 8.0; S[0b111] = 1.0;
  // Best schedule avoids ever materializing the size-8 pair {010,100}.
  // Contract {001,010}->{011}(1): live {011(1),100(2)} step_peak=2+2+1=5;
  //   then {011,100}->{111}(1): live_sum=1+2=3, +1 = 4. Peak = 5.
  // The schedule that forms {110}=8 has a step >= 8. So min peak = 5.
  REQUIRE(brute_force_min_peak(S, 3) == Catch::Approx(5.0));
}
```

Add includes at the top of the file if missing: `#include <functional>`, `#include <limits>`, `#include <vector>`, `#include <algorithm>`.

- [ ] **Step 2: Run to verify it fails, then passes**

Run build/run. Expected: this is a pure test helper, so it should compile and PASS immediately once the helper is added. If it FAILS, fix the helper until the hand-computed `5.0` matches (debugging the oracle now is the point - it is the source of truth).

- [ ] **Step 3: Commit**

```bash
git add tests/unit/test_optimize.cpp
git commit -m "optimize: add brute-force min-peak oracle for tests"
```

---

### Task 3: peak DP (`single_term_opt_peak_impl`) and `DensePeakSize`

**Files:**
- Modify: `core/optimize/options.hpp:16` (enum)
- Modify: `core/optimize/single_term.hpp` (new impl + dispatch in detail `single_term_opt`)
- Test: `tests/unit/test_optimize.cpp` (SECTION matching DP to oracle)

**Interfaces:**
- Consumes: `subset_footprints` (Task 1).
- Produces: `EvalSequence sequant::opt::detail::single_term_opt_peak_impl(TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz)` - optimal min-peak contraction order; and `double sequant::opt::detail::peak_cost(TensorNetwork const&, TIdxs const&, IdxToSz&&)` returning the achieved `peak[full]` (for tests). `ObjectiveFunction::DensePeakSize` selects this path in `single_term_opt<Metric>`.

- [ ] **Step 1: Add the enum value**

In `core/optimize/options.hpp:16`:

```cpp
enum class ObjectiveFunction { DenseFLOPs, DenseSize, DensePeakSize };
```

- [ ] **Step 2: Write the failing test (DP peak == oracle)**

Add a SECTION to `TEST_CASE("optimize")`:

```cpp
SECTION("DensePeakSize DP matches brute-force oracle") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // A 4-tensor chain whose intermediates differ in size by contraction order.
  auto t0 = parse_expr(L"g{a1;i1}", Symmetry::nonsymm);
  auto t1 = parse_expr(L"g{a1;a2}", Symmetry::nonsymm);
  auto t2 = parse_expr(L"g{a2;a3}", Symmetry::nonsymm);
  auto t3 = parse_expr(L"g{a3;i2}", Symmetry::nonsymm);
  TensorNetwork net{std::vector<ExprPtr>{t0, t1, t2, t3}};
  container::svector<Index> targets;  // i1,i2 left open
  auto S = opt::detail::subset_footprints(net, targets, idxsz);
  double oracle = brute_force_min_peak(S, 4);
  double dp = opt::detail::peak_cost(net, targets, idxsz);
  REQUIRE(dp == Catch::Approx(oracle));
}
```

- [ ] **Step 3: Run to verify it fails**

Expected: compile error - `peak_cost` / `single_term_opt_peak_impl` not declared.

- [ ] **Step 4: Implement the peak DP**

Add to `core/optimize/single_term.hpp` in `namespace sequant::opt::detail`, after `single_term_opt_impl`. The pebbling recurrence with the order choice and back-pointers:

```cpp
/// Per-subset state for the peak DP.
struct PeakRes {
  double peak = std::numeric_limits<double>::max();  // min peak to build subtree
  size_t lp = 0, rp = 0;     // winning bipartition (0 for singletons)
  bool lp_first = true;      // winning evaluation order
};

/// Minimum peak memory to evaluate the whole network, plus the order.
/// Fills `pr[n].peak` for every subset via:
///   peak[n] = min over (bipartition lp|rp, order) of
///             max( peak[first], S[first] + peak[second],
///                  S[lp] + S[rp] + S[n] )
/// S[n] is build-independent, and the recurrence is monotone in child peaks,
/// so the per-subset minimum is globally optimal (optimal substructure).
template <typename TIdxs, typename IdxToSz>
container::vector<PeakRes> peak_dp(TensorNetwork const& network,
                                   TIdxs const& tidxs, IdxToSz&& idxsz,
                                   container::vector<double> const& S) {
  auto const nt = network.tensors().size();
  container::vector<PeakRes> pr(size_t{1} << nt);
  for (size_t n = 0; n < pr.size(); ++n) {
    if (std::popcount(n) == 0) {
      pr[n].peak = 0.0;
      continue;
    }
    if (std::popcount(n) == 1) {
      pr[n].peak = S[n];  // a leaf, resident at its own size
      continue;
    }
    for (auto&& [lp, rp] : bits::bipartitions(n)) {
      if (lp == 0 || rp == 0) continue;
      double const both = S[lp] + S[rp] + S[n];
      double const lp_first =
          std::max({pr[lp].peak, S[lp] + pr[rp].peak, both});
      double const rp_first =
          std::max({pr[rp].peak, S[rp] + pr[lp].peak, both});
      double const cand = std::min(lp_first, rp_first);
      if (cand < pr[n].peak) {
        pr[n].peak = cand;
        pr[n].lp = lp;
        pr[n].rp = rp;
        pr[n].lp_first = (lp_first <= rp_first);
      }
    }
  }
  return pr;
}

template <typename TIdxs, typename IdxToSz>
double peak_cost(TensorNetwork const& network, TIdxs const& tidxs,
                 IdxToSz&& idxsz) {
  auto S = subset_footprints(network, tidxs, idxsz);
  auto pr = peak_dp(network, tidxs, idxsz, S);
  return pr.back().peak;
}

/// Reconstruct the EvalSequence from the peak DP back-pointers, honoring the
/// chosen evaluation order at each node (the lower-peak child is emitted last).
template <typename TIdxs, typename IdxToSz>
EvalSequence single_term_opt_peak_impl(TensorNetwork const& network,
                                       TIdxs const& tidxs, IdxToSz&& idxsz) {
  using ranges::views::concat;
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  auto S = subset_footprints(network, tidxs, idxsz);
  auto pr = peak_dp(network, tidxs, idxsz, S);
  // bottom-up emit: a subset's sequence is (first-child seq)(second-child seq)-1
  container::vector<EvalSequence> seq(pr.size());
  for (size_t n = 0; n < pr.size(); ++n) {
    if (std::popcount(n) == 1) {
      seq[n] = EvalSequence{static_cast<int>(std::countr_zero(n))};
    } else if (std::popcount(n) >= 2) {
      auto const& a = pr[n].lp_first ? seq[pr[n].lp] : seq[pr[n].rp];
      auto const& b = pr[n].lp_first ? seq[pr[n].rp] : seq[pr[n].lp];
      seq[n] = concat(a, b) | ranges::to<EvalSequence>;
      seq[n].push_back(-1);
    }
  }
  return seq.back();
}
```

- [ ] **Step 5: Dispatch `DensePeakSize` in `single_term_opt<Metric>`**

In the detail `single_term_opt<ObjectiveFunction Metric>` (single_term.hpp ~410), add before the existing `if constexpr`:

```cpp
  if constexpr (Metric == ObjectiveFunction::DensePeakSize) {
    SEQUANT_ASSERT(!subnet_cse &&
                   "subnet_cse not supported with DensePeakSize (Phase 1)");
    (void)is_volatile_leaf;
    (void)volatile_weight;
    (void)footprint_weight;
    return single_term_opt_peak_impl(network, tidxs, idxsz);
  } else if constexpr (Metric == ObjectiveFunction::DenseFLOPs) {
```

(Convert the existing `if constexpr (Metric == ObjectiveFunction::DenseFLOPs)` into the `else if constexpr` above, and update the trailing `static_assert` to also accept `DensePeakSize`.)

- [ ] **Step 6: Run to verify the test passes**

Run build/run. Expected: `DensePeakSize DP matches brute-force oracle` PASSES. Also run the existing `[optimize]` sections to confirm no regression.

- [ ] **Step 7: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/single_term.hpp core/optimize/options.hpp
git add core/optimize/single_term.hpp core/optimize/options.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: add DensePeakSize peak-memory DP (validated vs oracle)"
```

---

### Task 4: random-network regression + public `optimize()` exposure

**Files:**
- Modify: `core/optimize/optimize.cpp:41-44` (dispatch in `opt_pure_product`)
- Test: `tests/unit/test_optimize.cpp` (randomized regression + a DenseFLOPs-vs-DensePeakSize divergence case)

**Interfaces:**
- Consumes: `single_term_opt_peak_impl`, `peak_cost`, `brute_force_min_peak`.
- Produces: `optimize(expr, {.objective_function = ObjectiveFunction::DensePeakSize, .idx_to_extent = ...})` returns a min-peak-ordered product.

- [ ] **Step 1: Write the failing regression test**

Add a SECTION sweeping several small networks (n = 3..5) with varied space sizes, asserting DP == oracle each time:

```cpp
SECTION("DensePeakSize matches oracle over a battery of small networks") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  std::vector<std::vector<std::wstring>> nets = {
    {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"},
    {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;a3}", L"g{a3;i2}"},
    {L"g{i1;a1}", L"t{a1,a2;i1,i2}", L"g{a2;i2}"},
    {L"g{a1,a2;i1,i2}", L"t{a1;i1}", L"t{a2;i2}"},
  };
  for (auto const& spec : nets) {
    std::vector<ExprPtr> ts;
    for (auto const& s : spec) ts.push_back(parse_expr(s, Symmetry::nonsymm));
    TensorNetwork net{ts};
    container::svector<Index> targets;
    auto S = opt::detail::subset_footprints(net, targets, idxsz);
    double oracle = brute_force_min_peak(S, ts.size());
    double dp = opt::detail::peak_cost(net, targets, idxsz);
    REQUIRE(dp == Catch::Approx(oracle));
  }
}

SECTION("DensePeakSize can pick a lower-peak tree than DenseFLOPs") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // Construct a product where the flop-optimal order materializes a large
  // intermediate the peak-optimal order avoids; assert the peak objective's
  // achieved peak is <= the flop objective's tree peak.
  auto prod = parse_expr(
      L"g{a1;i1} * g{a1;a2} * g{a2;a3} * g{a3;i2}",
      Symmetry::nonsymm)->as<Product>();
  auto flops_tree = opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
      prod, idxsz, /*subnet_cse=*/false);
  auto peak_tree = opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
      prod, idxsz, /*subnet_cse=*/false);
  // Both must evaluate to equivalent expressions (same leaf set).
  REQUIRE(count_tensor_leaves(flops_tree) == count_tensor_leaves(peak_tree));
  // (Stronger numeric peak comparison added once peak_cost is exposed on trees.)
}
```

- [ ] **Step 2: Run to verify it fails**

Expected: compile/link error - `single_term_opt<DensePeakSize>(Product, ...)` reachable, but `optimize()` public dispatch for `DensePeakSize` missing; the regression SECTION should compile but the public path is exercised next.

- [ ] **Step 3: Dispatch in `opt_pure_product`**

In `core/optimize/optimize.cpp`, replace the `SEQUANT_ASSERT(... == DenseSize)` tail (~41-44) with a switch over all three:

```cpp
  if (opts.objective_function == ObjectiveFunction::DenseSize)
    return opt::single_term_opt<ObjectiveFunction::DenseSize>(
        prod, opts.idx_to_extent, subnet_cse, opts.is_volatile_leaf,
        opts.volatile_weight, opts.footprint_weight);
  SEQUANT_ASSERT(opts.objective_function == ObjectiveFunction::DensePeakSize);
  return opt::single_term_opt<ObjectiveFunction::DensePeakSize>(
      prod, opts.idx_to_extent, subnet_cse, opts.is_volatile_leaf,
      opts.volatile_weight, opts.footprint_weight);
```

(Keep the existing `DenseFLOPs` branch above it.)

- [ ] **Step 4: Run to verify it passes**

Run build/run. Expected: both new SECTIONs PASS; all existing `[optimize]` sections still PASS.

- [ ] **Step 5: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/optimize.cpp tests/unit/test_optimize.cpp
git add core/optimize/optimize.cpp tests/unit/test_optimize.cpp
git commit -m "optimize: expose DensePeakSize via optimize(); add peak regression battery"
```

---

## Self-review notes

- **Spec coverage (Phase 1 subset):** peak objective (spec 5) - Tasks 3-4; optimal-substructure validation (spec 5.2, 10) - Tasks 2-3 (independent oracle); order-dependence (spec 5.1) - encoded in `peak_dp`'s `min(lp_first, rp_first)`. Out of Phase 1 by design: footprint-as-monomial / batching (spec 6), CostModel abstraction (spec 3-4), model-is-evaluator (spec 3B), mpqc wiring (spec 8.2) - later plans.
- **`subnet_cse`:** asserted off for the peak path per the spec risk item; the additive path keeps CSE.
- **No regressions:** `single_term_opt_impl` (additive) is untouched; `DenseFLOPs`/`DenseSize` paths unchanged.
- **Naming consistency:** `subset_footprints`, `peak_dp`, `peak_cost`, `single_term_opt_peak_impl`, `PeakRes`, `ObjectiveFunction::DensePeakSize` used identically across tasks.
- **Known follow-ups to fold into Phase 2 planning:** expose `peak_cost` on an arbitrary (already-built) tree so the divergence test can assert a strict numeric peak reduction; confirm `parse_expr` index-space sizes via the test `Context` clone pattern already used at the top of `TEST_CASE("optimize")` (copy that setup into the new SECTIONs if `approximate_size()` is not configured by default).
