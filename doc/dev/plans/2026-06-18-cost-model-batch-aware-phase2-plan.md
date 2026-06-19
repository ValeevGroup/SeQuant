# Batch-aware cost model - Phase 2: DensePeakSizeBatched

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a batchability-aware peak-memory objective (`ObjectiveFunction::DensePeakSizeBatched`) that prices an intermediate carrying a free batchable (aux) index at its *sliced* footprint when it sits at/under a persistent batch-index contraction - validated against the Phase-1 peak DP (slice mode) and an independent batch-aware brute-force oracle (full mode).

**Architecture:** A two-mode dynamic program. `peak_slice[n]` = the all-co-resident peak of subtree `n` evaluated entirely sliced (= the Phase-1 DP run with batch-index extents replaced by `min(extent, batch_target_size)`). `peak_full[n]` = the peak at full extents, except a child that is a **batchable frontier** (a batch index is internal to it AND its subtree is persistent) may instead contribute its `peak_slice`. The objective allows the root itself to be batched. Phase 1's `single_term_opt_impl` / `peak_dp` are untouched.

**Tech Stack:** C++20, SeQuant `core/optimize`, Catch2 (`tests/unit/test_optimize.cpp`).

## Global Constraints

- Style: Google `.clang-format`, 80-col, 2-space, no tabs. `clang-format --style=file -i <file>` before each commit.
- No en-dashes (U+2013) or non-breaking spaces anywhere (pre-commit hook rejects them); ASCII hyphen only.
- No `Co-Authored-By` trailers.
- Branch: `feature/cost-model-batch-aware` (continue on it; Phase 1 already committed here).
- Memory model is **all-co-resident (model A)**, identical to Phase 1: every live tensor counts (intermediates + resident input leaves). Reuse Phase 1's `subset_footprints` and `peak_cost`; do NOT reimplement them.
- Phase 2 is SeQuant-only; no mpqc changes. Batchability inputs ride on `OptimizeOptions` (`is_batchable_index`, `batch_target_size`); persistence reuses the existing `is_volatile_leaf`. The `CostModel` abstraction and the model-is-evaluator unification remain later phases.
- `DensePeakSizeBatched` does not support `subnet_cse` (same Phase-1 limitation); the path asserts `subnet_cse == false`.
- Reference: spec `doc/dev/specs/2026-06-18-cost-model-batch-aware-design.md` section 6.

## File structure

- `core/optimize/options.hpp` - add `DensePeakSizeBatched` to `ObjectiveFunction`; add `is_batchable_index` (predicate) and `batch_target_size` (size_t) to `OptimizeOptions`.
- `core/optimize/single_term.hpp` - add `batched_extent` (extent wrapper), `subset_batchable_internal`, `leaf_volatile_mask`, `PeakBatchedRes`, `peak_dp_batched`, `peak_cost_batched`, `single_term_opt_peak_batched_impl`, and a 4th `if constexpr` arm in the detail `single_term_opt`.
- `core/optimize/optimize.cpp` - `DensePeakSizeBatched` branch in `opt_pure_product`.
- `tests/unit/test_optimize.cpp` - the batch-aware oracle and SECTIONs.

Build/run for every test step:
```bash
cmake --build build-test --target unit_tests-sequant -j
./build-test/tests/unit/unit_tests-sequant "[optimize]"
```

---

### Task 1: Batchability input tables

**Files:**
- Modify: `core/optimize/single_term.hpp` (helpers near `subset_footprints`)
- Modify: `core/optimize/options.hpp` (enum + fields)
- Test: `tests/unit/test_optimize.cpp` (SECTION)

**Interfaces:**
- Produces:
  - `auto detail::batched_extent(IdxToSz idxsz, std::function<bool(Index const&)> is_batchable, std::size_t batch)` -> a callable `(Index)->size_t` returning `is_batchable(ix) ? min(idxsz(ix), batch) : idxsz(ix)` (captures by value).
  - `container::vector<char> detail::subset_batchable_internal(TensorNetwork const& net, TIdxs const& tidxs, std::function<bool(Index const&)> const& is_batchable)` -> `internal[n] != 0` iff some batchable index appears in a tensor of subset `n` but NOT in `n`'s open indices (`init_results`).
  - `std::size_t detail::leaf_volatile_mask(TensorNetwork const& net, std::function<bool(Tensor const&)> const& is_volatile_leaf)` -> bit `i` set iff tensor `i` is a volatile leaf; 0 if predicate empty (mirrors the inline logic in `single_term_opt_impl`).
  - `OptimizeOptions::is_batchable_index` (`std::function<bool(Index const&)>`, default empty), `OptimizeOptions::batch_target_size` (`std::size_t`, default 0); `ObjectiveFunction::DensePeakSizeBatched`.

- [ ] **Step 1: Add the enum value and OptimizeOptions fields**

`options.hpp`: `enum class ObjectiveFunction { DenseFLOPs, DenseSize, DensePeakSize, DensePeakSizeBatched };` and document the new enumerator (peak with batchable-index slicing under a persistent batch-contraction). Add fields with doc comments:
```cpp
  /// Predicate marking an Index as living in a batchable space the runtime
  /// evaluator slices over (e.g. the DF/RI auxiliary space; = the eval cache's
  /// accept_aux). Only consulted by ObjectiveFunction::DensePeakSizeBatched.
  std::function<bool(Index const&)> is_batchable_index = {};

  /// Slice size for batchable indices: such an index contributes
  /// min(extent, batch_target_size) to the sliced footprint. 0 disables the
  /// batched discount. Only consulted by DensePeakSizeBatched.
  std::size_t batch_target_size = 0;
```

- [ ] **Step 2: Write the failing test**

Add inside the cloned-context block of `TEST_CASE("optimize")` (where index-space sizes are set, as Tasks 1-4 of Phase 1 did):
```cpp
SECTION("batchability input tables") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  // Treat the aux/density-fitting space as batchable. Adjust the space tag to
  // whatever the test context uses for a DF/aux index; the point is a predicate
  // that is true on exactly the aux index of `g{...;...;K}` integrals.
  auto is_batchable = [](Index const& ix) -> bool {
    return ix.space().base_key() == L"F";  // aux/fitting space tag in this ctx
  };
  // g{a1;i1;F1} carries an aux index F1; contracting two such over F1 makes F1
  // internal to that pair.
  auto t0 = deserialize(L"g{a1;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
  auto t1 = deserialize(L"g{a1;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
  TensorNetwork net{std::vector<ExprPtr>{t0, t1}};
  container::svector<Index> targets;
  auto internal = opt::detail::subset_batchable_internal(net, targets,
                                                         is_batchable);
  REQUIRE(internal.size() == 4u);
  REQUIRE(internal[0b01] == 0);  // singleton: F1 still open, not internal
  REQUIRE(internal[0b10] == 0);
  REQUIRE(internal[0b11] != 0);  // pair contracts F1 -> internal
  // batched_extent shrinks only the aux index.
  auto be = opt::detail::batched_extent(idxsz, is_batchable, 1);
  Index F1 = t0->as<Tensor>().aux().at(0);
  REQUIRE(be(F1) == 1u);                      // batchable -> min(size,1)
  Index a1 = t0->as<Tensor>().bra().at(0);
  REQUIRE(be(a1) == idxsz(a1));               // non-batchable unchanged
}
```
NOTE: confirm the aux-space tag/predicate against the test context; if `g{..;..;F}` aux indices use a different space key, use that. The asserted *structure* (singleton not internal, pair internal; aux shrinks, non-aux does not) is what matters.

- [ ] **Step 3: Run to verify it fails**

Expected: compile error - `subset_batchable_internal` / `batched_extent` not declared.

- [ ] **Step 4: Implement the helpers**

Add to `core/optimize/single_term.hpp` in `namespace sequant::opt::detail`, after `subset_footprints`:
```cpp
/// Wrap an extent provider so a batchable index reports min(extent, batch).
template <typename IdxToSz>
auto batched_extent(IdxToSz idxsz, std::function<bool(Index const&)> is_batchable,
                    std::size_t batch) {
  return [idxsz = std::move(idxsz), is_batchable = std::move(is_batchable),
          batch](Index const& ix) -> std::size_t {
    std::size_t const e = idxsz(ix);
    return (is_batchable && is_batchable(ix)) ? std::min(e, batch) : e;
  };
}

/// internal[n] != 0 iff a batchable index is contracted *within* subset n: it
/// appears in some tensor of n but not in n's open indices. Such an n is a
/// candidate batch frontier (persistence is checked separately).
template <typename TIdxs>
container::vector<char> subset_batchable_internal(
    TensorNetwork const& network, TIdxs const& tidxs,
    std::function<bool(Index const&)> const& is_batchable) {
  auto const nt = network.tensors().size();
  container::vector<OptRes> results(size_t{1} << nt);
  init_results(network, tidxs, results);
  // Per-tensor index sets (bra+ket+aux).
  container::vector<container::vector<Index>> tix(nt);
  {
    size_t i = 0;
    for (auto&& t : network.tensors()) {
      auto tp = std::dynamic_pointer_cast<Tensor>(t);
      for (auto&& ix : ranges::views::concat(tp->bra(), tp->ket(), tp->aux()))
        tix[i].push_back(ix);
      ++i;
    }
  }
  container::vector<char> internal(results.size(), 0);
  for (size_t n = 1; n < results.size(); ++n) {
    if (std::popcount(n) < 2) continue;  // singletons contract nothing
    auto const& open = results[n].indices;
    bool found = false;
    for (size_t b = 0; b < nt && !found; ++b) {
      if (!(n & (size_t{1} << b))) continue;
      for (auto&& ix : tix[b]) {
        if (is_batchable && is_batchable(ix) &&
            ranges::find(open, ix) == ranges::end(open)) {
          found = true;
          break;
        }
      }
    }
    internal[n] = found ? 1 : 0;
  }
  return internal;
}

/// Bit i set iff network tensor i is a volatile leaf. Mirrors the inline mask
/// built in single_term_opt_impl; 0 when the predicate is empty.
inline std::size_t leaf_volatile_mask(
    TensorNetwork const& network,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  std::size_t mask = 0;
  if (!is_volatile_leaf) return 0;
  std::size_t i = 0;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    if (tp && is_volatile_leaf(*tp)) mask |= (std::size_t{1} << i);
    ++i;
  }
  return mask;
}
```

- [ ] **Step 5: Run to verify it passes**

Run build/run. Expected: PASS (new SECTION) and all existing `[optimize]` SECTIONs still pass.

- [ ] **Step 6: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/single_term.hpp core/optimize/options.hpp tests/unit/test_optimize.cpp
git add core/optimize/single_term.hpp core/optimize/options.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: batchability input tables (extent wrapper, internal, volatile mask)"
```

---

### Task 2: Batch-aware brute-force oracle

**Files:**
- Test: `tests/unit/test_optimize.cpp` (free helper + hand-checked TEST_CASE)

**Interfaces:**
- Produces (test-local): `double batched_min_peak(std::vector<double> const& S_full, std::vector<double> const& S_slc, std::vector<char> const& is_frontier, size_t nt)` - the minimum, over all full binary trees x evaluation orders x subsets-of-eligible-frontiers-to-batch, of the model-A peak where a live subset `m` is sized `S_slc[m]` iff `m` is a STRICT subset of some batched frontier, else `S_full[m]`. `is_frontier[n]` is the precomputed (batch-internal AND persistent) flag.

This is the Phase-2 trust anchor; its correctness is established by the hand-checked TEST_CASE, so derive that value independently before trusting the helper.

- [ ] **Step 1: Write the helper + failing test**

Add above `TEST_CASE("optimize")`:
```cpp
// Model-A peak of a single contraction tree under a fixed batch-set, minimized
// over child-evaluation orders, computed by direct memory simulation (NOT the DP
// recurrence). `batched` is the set of frontier subsets chosen to be sliced; a
// live subset m is sliced iff m is a strict subset of some batched frontier.
static double batched_tree_peak(size_t n, std::vector<double> const& S_full,
                                std::vector<double> const& S_slc,
                                std::set<size_t> const& batched);

// recursive: best peak to evaluate subtree `n` given the batch-set, trying both
// child orders. Returns the subtree's internal peak (model A). Sizing of a
// subset m uses sliced(m): S_slc[m] if m strictly inside a batched frontier.
static double batched_subtree(size_t n, size_t lp_choice_unused,
                              std::vector<double> const& S_full,
                              std::vector<double> const& S_slc,
                              std::set<size_t> const& batched);

// size of a live subset under the batch-set.
static double sz(size_t m, std::vector<double> const& S_full,
                 std::vector<double> const& S_slc,
                 std::set<size_t> const& batched) {
  for (size_t F : batched)
    if (m != F && (m & F) == m) return S_slc[m];  // m strict subset of F
  return S_full[m];
}

// L (leaf-sum) of subset n under the batch-set.
static double Lsz(size_t n, size_t nt, std::vector<double> const& S_full,
                  std::vector<double> const& S_slc,
                  std::set<size_t> const& batched) {
  double s = 0.0;
  for (size_t b = 0; b < nt; ++b)
    if (n & (size_t{1} << b)) s += sz(size_t{1} << b, S_full, S_slc, batched);
  return s;
}

static double batched_subtree(size_t n, std::vector<double> const& S_full,
                              std::vector<double> const& S_slc,
                              std::set<size_t> const& batched, size_t nt) {
  if (std::popcount(n) == 1) return sz(n, S_full, S_slc, batched);
  double best = std::numeric_limits<double>::max();
  // enumerate bipartitions of n into (lp, rp), lp the lower-bit-containing half
  for (size_t lp = (n - 1) & n; lp; lp = (lp - 1) & n) {
    size_t rp = n ^ lp;
    if (lp > rp) continue;  // each unordered pair once
    double pl = batched_subtree(lp, S_full, S_slc, batched, nt);
    double pr = batched_subtree(rp, S_full, S_slc, batched, nt);
    double both = sz(lp, S_full, S_slc, batched) + sz(rp, S_full, S_slc, batched) +
                  sz(n, S_full, S_slc, batched);
    double lp_first = std::max({Lsz(rp, nt, S_full, S_slc, batched) + pl,
                                sz(lp, S_full, S_slc, batched) + pr, both});
    double rp_first = std::max({Lsz(lp, nt, S_full, S_slc, batched) + pr,
                                sz(rp, S_full, S_slc, batched) + pl, both});
    best = std::min(best, std::min(lp_first, rp_first));
  }
  return best;
}

static double batched_min_peak(std::vector<double> const& S_full,
                               std::vector<double> const& S_slc,
                               std::vector<char> const& is_frontier, size_t nt) {
  size_t const full = (size_t{1} << nt) - 1;
  // eligible frontiers among all subsets
  std::vector<size_t> elig;
  for (size_t n = 1; n <= full; ++n)
    if (is_frontier[n]) elig.push_back(n);
  double best = std::numeric_limits<double>::max();
  // enumerate every subset of eligible frontiers to batch
  for (size_t mask = 0; mask < (size_t{1} << elig.size()); ++mask) {
    std::set<size_t> batched;
    for (size_t k = 0; k < elig.size(); ++k)
      if (mask & (size_t{1} << k)) batched.insert(elig[k]);
    best = std::min(best, batched_subtree(full, S_full, S_slc, batched, nt));
  }
  return best;
}

TEST_CASE("batched_min_peak oracle", "[optimize]") {
  // 2 leaves T0,T1 each size 4, sharing an aux index contracted in the pair.
  //   S_full[01]=S_full[10]=4 (leaves); S_full[11]=2 (aux-free result).
  //   sliced sizes: S_slc[01]=S_slc[10]=2 (aux halved); S_slc[11]=2 (no aux).
  //   is_frontier[11]=1 (aux internal + persistent), else 0.
  std::vector<double> Sf(4), Ss(4);
  Sf[0b01] = Sf[0b10] = 4; Sf[0b11] = 2;
  Ss[0b01] = Ss[0b10] = 2; Ss[0b11] = 2;
  std::vector<char> fr(4, 0); fr[0b11] = 1;
  // Not batched: peak = 4+4+2 = 10. Batched {11}: leaves sliced to 2 each ->
  // 2+2+2 = 6. min = 6.
  REQUIRE(batched_min_peak(Sf, Ss, fr, 2) == Catch::Approx(6.0));
}
```
Add includes if missing: `<set>`, `<functional>`, `<limits>`, `<vector>`, `<algorithm>`.

- [ ] **Step 2: Run, debug the oracle until the hand value holds**

Run build/run. The oracle is the trust anchor: independently confirm the 6.0 by hand (un-batched 10 vs batched 6) before accepting. If your correct implementation yields a different number, the brief's value is wrong - report it with your derivation rather than editing the assertion blindly.

- [ ] **Step 3: Commit**

```bash
clang-format --style=file -i tests/unit/test_optimize.cpp
git add tests/unit/test_optimize.cpp
git commit -m "optimize: batch-aware brute-force min-peak oracle for tests"
```

---

### Task 3: Two-mode batched DP

**Files:**
- Modify: `core/optimize/single_term.hpp` (`PeakBatchedRes`, `peak_dp_batched`, `peak_cost_batched`, dispatch)
- Test: `tests/unit/test_optimize.cpp` (slice-mode==Phase1 + full-mode==oracle SECTIONs)

**Interfaces:**
- Consumes: Task 1 helpers, Phase-1 `subset_footprints`/`peak_cost`/`peak_dp`.
- Produces:
  - `struct detail::PeakBatchedRes { double peak_full, peak_slice; size_t lp_full, rp_full; bool lp_first_full; size_t lp_slice, rp_slice; bool lp_first_slice; };`
  - `container::vector<PeakBatchedRes> detail::peak_dp_batched(net, tidxs, idxsz, is_batchable, batch, volatile_mask)`.
  - `double detail::peak_cost_batched(net, tidxs, idxsz, is_batchable, batch, is_volatile_leaf)` -> the objective `cc(root)` (root may itself be batched).

- [ ] **Step 1: Write the failing tests**

Two SECTIONs (inside the sized-context block):
```cpp
SECTION("DensePeakSizeBatched slice mode equals Phase-1 peak at batched extents") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) -> bool {
    return ix.space().base_key() == L"F";  // aux tag (match Task 1)
  };
  std::size_t const batch = 1;
  auto t0 = deserialize(L"g{a1;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
  auto t1 = deserialize(L"g{a2;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
  auto t2 = deserialize(L"g{a2;i2;F2}", {.def_perm_symm = Symmetry::Nonsymm});
  TensorNetwork net{std::vector<ExprPtr>{t0, t1, t2}};
  container::svector<Index> targets;
  auto vmask = opt::detail::leaf_volatile_mask(net, {});
  auto pr = opt::detail::peak_dp_batched(net, targets, idxsz, is_batchable,
                                         batch, vmask);
  // slice mode of the full set == Phase-1 peak_cost run with batched extents.
  auto be = opt::detail::batched_extent(idxsz, is_batchable, batch);
  double phase1_batched = opt::detail::peak_cost(net, targets, be);
  REQUIRE(pr.back().peak_slice == Catch::Approx(phase1_batched));
}

SECTION("DensePeakSizeBatched full mode matches batch-aware oracle") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) -> bool {
    return ix.space().base_key() == L"F";
  };
  std::size_t const batch = 1;
  std::vector<std::vector<std::wstring>> nets = {
    {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},
    {L"g{a1;i1;F1}", L"g{a1;a2;F1}", L"g{a2;i2;F2}", L"g{a3;i2;F2}"},
  };
  for (auto const& spec : nets) {
    std::vector<ExprPtr> ts;
    for (auto const& s : spec)
      ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
    TensorNetwork net{ts};
    container::svector<Index> targets;
    auto vmask = opt::detail::leaf_volatile_mask(net, {});
    auto be = opt::detail::batched_extent(idxsz, is_batchable, batch);
    auto S_full_c = opt::detail::subset_footprints(net, targets, idxsz);
    auto S_slc_c = opt::detail::subset_footprints(net, targets, be);
    auto internal = opt::detail::subset_batchable_internal(net, targets,
                                                           is_batchable);
    std::vector<double> S_full(S_full_c.begin(), S_full_c.end());
    std::vector<double> S_slc(S_slc_c.begin(), S_slc_c.end());
    std::vector<char> fr(internal.size());
    for (size_t n = 0; n < internal.size(); ++n)
      fr[n] = (internal[n] && ((vmask & n) == 0)) ? 1 : 0;
    double oracle = batched_min_peak(S_full, S_slc, fr, ts.size());
    double dp = opt::detail::peak_cost_batched(net, targets, idxsz,
                                               is_batchable, batch, {});
    REQUIRE(dp == Catch::Approx(oracle));
  }
}
```

- [ ] **Step 2: Run to verify it fails**

Expected: compile error - `peak_dp_batched`/`peak_cost_batched` not declared.

- [ ] **Step 3: Implement the two-mode DP**

Add to `single_term.hpp` (detail), after `peak_cost`:
```cpp
struct PeakBatchedRes {
  double peak_full = std::numeric_limits<double>::max();
  double peak_slice = std::numeric_limits<double>::max();
  size_t lp_full = 0, rp_full = 0;
  bool lp_first_full = true;
  size_t lp_slice = 0, rp_slice = 0;
  bool lp_first_slice = true;
};

/// Two-mode batched peak DP (spec 6.2). peak_slice[n] = Phase-1 peak with
/// batch-index extents sliced; peak_full[n] = peak at full extents with a
/// batchable-frontier child contributing min(full, slice). A subset c is a
/// batchable frontier iff a batch index is internal to it (internal[c]) and it
/// is persistent ((volatile_mask & c) == 0).
template <typename TIdxs, typename IdxToSz>
container::vector<PeakBatchedRes> peak_dp_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::size_t volatile_mask) {
  auto const nt = network.tensors().size();
  auto S_full = subset_footprints(network, tidxs, idxsz);
  auto S_slc = subset_footprints(
      network, tidxs, batched_extent(idxsz, is_batchable, batch));
  auto internal = subset_batchable_internal(network, tidxs, is_batchable);
  container::vector<double> L_full(S_full.size(), 0.0), L_slc(S_slc.size(), 0.0);
  for (size_t n = 0; n < S_full.size(); ++n)
    for (size_t b = 0; b < nt; ++b)
      if (n & (size_t{1} << b)) {
        L_full[n] += S_full[size_t{1} << b];
        L_slc[n] += S_slc[size_t{1} << b];
      }
  auto frontier = [&](size_t c) {
    return internal[c] && ((volatile_mask & c) == 0);
  };
  container::vector<PeakBatchedRes> pr(size_t{1} << nt);
  for (size_t n = 0; n < pr.size(); ++n) {
    if (std::popcount(n) == 0) {
      pr[n].peak_full = pr[n].peak_slice = 0.0;
      continue;
    }
    if (std::popcount(n) == 1) {
      pr[n].peak_full = S_full[n];
      pr[n].peak_slice = S_slc[n];
      continue;
    }
    for (auto&& [lp, rp] : bits::bipartitions(n)) {
      if (lp == 0 || rp == 0) continue;
      // --- slice mode: Phase-1 recurrence with sliced sizes ---
      double both_s = S_slc[lp] + S_slc[rp] + S_slc[n];
      double lpf_s =
          std::max({L_slc[rp] + pr[lp].peak_slice, S_slc[lp] + pr[rp].peak_slice,
                    both_s});
      double rpf_s =
          std::max({L_slc[lp] + pr[rp].peak_slice, S_slc[rp] + pr[lp].peak_slice,
                    both_s});
      double cand_s = std::min(lpf_s, rpf_s);
      if (cand_s < pr[n].peak_slice) {
        pr[n].peak_slice = cand_s;
        pr[n].lp_slice = lp;
        pr[n].rp_slice = rp;
        pr[n].lp_first_slice = (lpf_s <= rpf_s);
      }
      // --- full mode: frontier substitution on each child ---
      double cc_lp = frontier(lp)
                         ? std::min(pr[lp].peak_full, pr[lp].peak_slice)
                         : pr[lp].peak_full;
      double cc_rp = frontier(rp)
                         ? std::min(pr[rp].peak_full, pr[rp].peak_slice)
                         : pr[rp].peak_full;
      double both_f = S_full[lp] + S_full[rp] + S_full[n];
      double lpf_f = std::max({L_full[rp] + cc_lp, S_full[lp] + cc_rp, both_f});
      double rpf_f = std::max({L_full[lp] + cc_rp, S_full[rp] + cc_lp, both_f});
      double cand_f = std::min(lpf_f, rpf_f);
      if (cand_f < pr[n].peak_full) {
        pr[n].peak_full = cand_f;
        pr[n].lp_full = lp;
        pr[n].rp_full = rp;
        pr[n].lp_first_full = (lpf_f <= rpf_f);
      }
    }
  }
  return pr;
}

/// Objective value: peak_full[root], or its sliced form if the root itself is a
/// batchable frontier (the whole term sliced over the batch index).
template <typename TIdxs, typename IdxToSz>
double peak_cost_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  auto const vmask = leaf_volatile_mask(network, is_volatile_leaf);
  auto pr = peak_dp_batched(network, tidxs, idxsz, is_batchable, batch, vmask);
  auto internal = subset_batchable_internal(network, tidxs, is_batchable);
  size_t const root = pr.size() - 1;
  bool root_frontier = internal[root] && ((vmask & root) == 0);
  return root_frontier ? std::min(pr[root].peak_full, pr[root].peak_slice)
                       : pr[root].peak_full;
}
```

- [ ] **Step 4: Dispatch the new objective in `single_term_opt`**

In the detail `single_term_opt<Metric>` (single_term.hpp), add a 4th arm BEFORE the `DensePeakSize` arm (or alongside): when `Metric == DensePeakSizeBatched`, assert `!subnet_cse`, compute `vmask` from `is_volatile_leaf`, and `return single_term_opt_peak_batched_impl(...)` (implemented in Task 4). For this task, to make the test SECTIONs reachable, it is enough that `peak_dp_batched`/`peak_cost_batched` compile and are called directly by the tests; the `single_term_opt` arm can be added here or in Task 4 - put a placeholder arm that static_asserts unreachable only if you do not wire it now. (Recommended: wire the arm fully in Task 4 with reconstruction; Task 3 tests call `peak_cost_batched` directly.)

- [ ] **Step 5: Run to verify it passes**

Run build/run. Expected: both new SECTIONs PASS; all existing `[optimize]` SECTIONs still pass. If the full-mode-vs-oracle SECTION fails, the divergence localizes the bug (DP vs the independent oracle) - debug the DP recurrence, not the oracle (Task 2's oracle is hand-anchored).

- [ ] **Step 6: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git add core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: two-mode DensePeakSizeBatched DP (slice==Phase1, full==oracle)"
```

---

### Task 4: Reconstruction, dispatch, public API

**Files:**
- Modify: `core/optimize/single_term.hpp` (`single_term_opt_peak_batched_impl` + `single_term_opt` arm)
- Modify: `core/optimize/optimize.cpp` (`opt_pure_product` branch)
- Test: `tests/unit/test_optimize.cpp` (reconstruction-achieves-optimum + public API)

**Interfaces:**
- Produces: `EvalSequence detail::single_term_opt_peak_batched_impl(net, tidxs, idxsz, is_batchable, batch, is_volatile_leaf)` - the contraction order realizing `peak_cost_batched`; `optimize(expr, {objective_function=DensePeakSizeBatched, is_batchable_index, batch_target_size, ...})` works.

- [ ] **Step 1: Write the failing tests**

```cpp
SECTION("DensePeakSizeBatched reconstructed sequence achieves the optimum") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  auto is_batchable = [](Index const& ix) -> bool {
    return ix.space().base_key() == L"F";
  };
  std::size_t const batch = 1;
  // peak_of_sequence_batched: model-A peak of an EvalSequence, with the chosen
  // batch-set inferred as "every persistent batch-internal subtree is sliced"
  // -- i.e. simulate at full extents but size a subset sliced iff it is a strict
  // subset of a frontier subset that the reconstruction batched. Simplest exact
  // check: compare peak_cost_batched against the batch-aware oracle (already in
  // Task 3) AND assert the reconstructed sequence's leaf set is intact.
  std::vector<std::wstring> spec = {L"g{a1;i1;F1}", L"g{a2;i1;F1}",
                                    L"g{a2;i2;F2}"};
  std::vector<ExprPtr> ts;
  for (auto const& s : spec)
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;
  auto seq = opt::detail::single_term_opt_peak_batched_impl(
      net, targets, idxsz, is_batchable, batch, {});
  // valid postorder over all leaves: count of non-negative tokens == nt, count
  // of -1 == nt-1.
  size_t leaves = 0, merges = 0;
  for (int tok : seq) (tok >= 0 ? leaves : merges)++;
  REQUIRE(leaves == ts.size());
  REQUIRE(merges == ts.size() - 1);
}

SECTION("optimize() public API dispatches DensePeakSizeBatched") {
  using namespace sequant;
  auto expr = deserialize(L"g{a1;i1;F1} * g{a2;i1;F1} * g{a2;i2;F2}",
                          {.def_perm_symm = Symmetry::Nonsymm});
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = [](Index const& ix) -> std::size_t {
    return ix.space().approximate_size();
  };
  opts.is_batchable_index = [](Index const& ix) -> bool {
    return ix.space().base_key() == L"F";
  };
  opts.batch_target_size = 1;
  auto optimized = optimize(expr, opts);
  REQUIRE(optimized);
  REQUIRE(count_tensor_leaves(optimized) == 3u);
}
```

- [ ] **Step 2: Run to verify it fails**

Expected: the reconstruction SECTION fails to compile (`single_term_opt_peak_batched_impl` missing); the public-API SECTION fails at runtime (`opt_pure_product` `SEQUANT_ASSERT(objective_function == DensePeakSize)` fires for the unhandled `DensePeakSizeBatched`).

- [ ] **Step 3: Implement reconstruction**

Add to `single_term.hpp` (detail). It mirrors `single_term_opt_peak_impl` but walks two-mode back-pointers: each subset is reconstructed in the mode chosen for it (a frontier child taken in slice mode reconstructs from its slice back-pointers, else full).
```cpp
template <typename TIdxs, typename IdxToSz>
EvalSequence single_term_opt_peak_batched_impl(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  using ranges::views::concat;
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  auto const vmask = leaf_volatile_mask(network, is_volatile_leaf);
  auto pr = peak_dp_batched(network, tidxs, idxsz, is_batchable, batch, vmask);
  auto internal = subset_batchable_internal(network, tidxs, is_batchable);
  auto frontier = [&](size_t c) {
    return internal[c] && ((vmask & c) == 0);
  };
  // emit subset n in a given mode (false=full, true=slice); a frontier child is
  // emitted in min(full,slice) mode.
  container::vector<EvalSequence> seq_full(pr.size()), seq_slice(pr.size());
  std::function<EvalSequence const&(size_t, bool)> emit =
      [&](size_t n, bool slice) -> EvalSequence const& {
    auto& slot = slice ? seq_slice[n] : seq_full[n];
    if (!slot.empty() || std::popcount(n) == 0) return slot;
    if (std::popcount(n) == 1) {
      slot = EvalSequence{static_cast<int>(std::countr_zero(n))};
      return slot;
    }
    size_t lp = slice ? pr[n].lp_slice : pr[n].lp_full;
    size_t rp = slice ? pr[n].rp_slice : pr[n].rp_full;
    bool lp_first = slice ? pr[n].lp_first_slice : pr[n].lp_first_full;
    // child mode: in slice mode everything is sliced; in full mode a frontier
    // child uses whichever of full/slice was cheaper.
    auto child_mode = [&](size_t c) {
      if (slice) return true;
      return frontier(c) && (pr[c].peak_slice < pr[c].peak_full);
    };
    auto const& a = lp_first ? emit(lp, child_mode(lp)) : emit(rp, child_mode(rp));
    auto const& b = lp_first ? emit(rp, child_mode(rp)) : emit(lp, child_mode(lp));
    slot = concat(a, b) | ranges::to<EvalSequence>;
    slot.push_back(-1);
    return slot;
  };
  size_t const root = pr.size() - 1;
  bool root_slice = frontier(root) && (pr[root].peak_slice < pr[root].peak_full);
  return emit(root, root_slice);
}
```
NOTE: `emit` references `seq_full`/`seq_slice` slots by reference while recursing and assigning; ensure the recursive `emit(...)` calls complete (populating child slots) BEFORE `concat` reads them - the code above sequences the `emit` calls into `a`/`b` locals first, then concats, which is correct. If returning references through `std::function` proves fragile, return `EvalSequence` by value and memoize in the vectors explicitly.

Wire the `single_term_opt` arm (if not already): in the detail `single_term_opt<Metric>`,
```cpp
  if constexpr (Metric == ObjectiveFunction::DensePeakSizeBatched) {
    SEQUANT_ASSERT(!subnet_cse &&
                   "subnet_cse not supported with DensePeakSizeBatched");
    (void)volatile_weight;
    (void)footprint_weight;
    return single_term_opt_peak_batched_impl(network, tidxs, idxsz,
                                             is_batchable_index,
                                             batch_target_size, is_volatile_leaf);
  } else if constexpr ( ... existing arms ... )
```
This requires the detail `single_term_opt` to receive `is_batchable_index` and `batch_target_size`. Thread two new parameters through the detail `single_term_opt` (both overloads) and the public `single_term_opt(Product, ...)`; default them to `{}` / `0` so existing callers are unaffected. Update the `static_assert` message to list all four objectives.

- [ ] **Step 4: Dispatch in `opt_pure_product`**

`core/optimize/optimize.cpp`: add, before the final `SEQUANT_ASSERT`,
```cpp
  if (opts.objective_function == ObjectiveFunction::DensePeakSizeBatched)
    return opt::single_term_opt<ObjectiveFunction::DensePeakSizeBatched>(
        prod, opts.idx_to_extent, subnet_cse, opts.is_volatile_leaf,
        opts.volatile_weight, opts.footprint_weight, opts.is_batchable_index,
        opts.batch_target_size);
```
and pass the two new trailing args to the other three `single_term_opt` calls (they ignore them). Keep the final `SEQUANT_ASSERT` exhaustive.

- [ ] **Step 5: Run to verify it passes**

Run build/run. Expected: both new SECTIONs pass; all prior `[optimize]` SECTIONs pass.

- [ ] **Step 6: clang-format and commit**

```bash
clang-format --style=file -i core/optimize/single_term.hpp core/optimize/optimize.cpp tests/unit/test_optimize.cpp
git add core/optimize/single_term.hpp core/optimize/optimize.cpp tests/unit/test_optimize.cpp
git commit -m "optimize: reconstruct + dispatch DensePeakSizeBatched via optimize()"
```

---

## Self-review notes

- **Spec coverage:** monomial footprint (spec 6.1) -> Task 1 `batched_extent` + slice footprints; two-mode recurrence (6.2) -> Task 3 `peak_dp_batched`; frontier gate (6.2) -> Task 1 `subset_batchable_internal` + persistence; root-batchable objective -> `peak_cost_batched`; validation (6.4) -> Task 3 slice==Phase1 and full==oracle SECTIONs; plumbing (6.5) -> Task 1 OptimizeOptions fields, Task 4 dispatch.
- **Trust anchors:** slice mode is tied to the *already-validated* Phase-1 `peak_cost`; full mode to the Task-2 oracle (hand-anchored at 6.0). Neither is circular.
- **No regression:** Phase-1 `peak_dp`/`single_term_opt_impl`/`peak_cost` and `DenseFLOPs`/`DenseSize`/`DensePeakSize` paths are untouched; only a new arm + two defaulted params are added.
- **Open risk to flag to reviewer:** (a) the `is_batchable` predicate in tests keys on an aux space tag (`base_key() == L"F"`) that must match the test context's DF/aux index - if `g{..;..;F}` aux indices are not in such a space, pick the correct predicate; the asserted structure is what matters. (b) The oracle enumerates `2^(#eligible frontiers)` batch-sets; keep test networks small (nt <= 4, few aux contractions). (c) Confirm `single_term_opt_peak_batched_impl`'s reference-returning memoized `emit` is sound under the compiler; fall back to by-value memoization if not.
- **Deferred:** the `peak_of_sequence_batched` numeric reconstruction check is intentionally replaced by a structural check (valid postorder, intact leaf set) plus the `peak_cost_batched == oracle` equality from Task 3; a full numeric "reconstructed sequence achieves the batched optimum" simulator is a nice-to-have for a later hardening pass (it needs the batch-set the reconstruction chose, which is derivable from the per-node mode decisions).
