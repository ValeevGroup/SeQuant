# Batch-aware cost model - Phase 2: DensePeakSizeBatched (per-index multi-mode)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `ObjectiveFunction::DensePeakSizeBatched` - a peak-memory objective that prices each batchable (aux) index as sliceable **independently**: an intermediate is charged its sliced footprint over exactly the batchable indices whose persistent contraction frontiers sit above it.

**Architecture:** A multi-mode DP `peak[n][B]` over subsets `n` and sliced-sets `B` (subsets of the term's distinct batchable indices `K`, `m = |K|` small). Sizes come from `2^m` Phase-1 footprint tables (one per `B`). At each contraction, the batchable indices contracted there (`A`) may, if the subtree is persistent, be added to the slice context for that subtree. Objective = `peak[root][empty]`. Phase 1's `peak_dp`/`peak_cost`/`single_term_opt_impl` are untouched. See spec `doc/dev/specs/2026-06-18-cost-model-batch-aware-design.md` section 6.

**Tech Stack:** C++20, SeQuant `core/optimize`, Catch2 (`tests/unit/test_optimize.cpp`).

## Global Constraints

- Style: Google `.clang-format`, 80-col, 2-space, no tabs. `clang-format --style=file -i <file>` before each commit.
- No en-dashes (U+2013) or non-breaking spaces (pre-commit hook rejects them); ASCII hyphen only.
- No `Co-Authored-By` trailers.
- Branch: `feature/cost-model-batch-aware` (continue on it).
- Memory model is **all-co-resident (model A)**, as in Phase 1. Reuse `subset_footprints`, `peak_cost`, `init_results`, `bits::bipartitions`; do not reimplement.
- SeQuant-only; no mpqc. Batchability rides on `OptimizeOptions` (`is_batchable_index` predicate, single shared `batch_target_size`); persistence reuses `is_volatile_leaf`. `CostModel` abstraction and eval-side unification are later phases.
- `DensePeakSizeBatched` does not support `subnet_cse`; the path asserts `subnet_cse == false`.
- Keep `m` (distinct batchable indices) and `nt` small in tests; DP is `2^nt * 2^m` states and the oracle is `tree x order x 2^(contracted aux per node)`.

## File structure

- `core/optimize/options.hpp` - `DensePeakSizeBatched` enumerator; `OptimizeOptions::{is_batchable_index, batch_target_size}`.
- `core/optimize/single_term.hpp` - `batchable_index_list`, `sliced_footprints`, `subset_open_aux` (or reuse open indices), `leaf_volatile_mask`, `BatchedRes`, `peak_dp_batched`, `peak_cost_batched`, `single_term_opt_peak_batched_impl`, dispatch arm.
- `core/optimize/optimize.cpp` - `opt_pure_product` branch.
- `tests/unit/test_optimize.cpp` - per-index oracle + SECTIONs.

Build/run each test step:
```bash
cmake --build build-test --target unit_tests-sequant -j
./build-test/tests/unit/unit_tests-sequant "[optimize]"
```

---

### Task 1: Per-index batchability input tables

**Files:** `core/optimize/options.hpp`, `core/optimize/single_term.hpp`, test.

**Interfaces (produces, all in `sequant::opt::detail`):**
- `container::vector<Index> batchable_index_list(TensorNetwork const&, std::function<bool(Index const&)> const& is_batchable)` - the distinct batchable indices in canonical (sorted-by-appearance, deduped) order; index `i` is bit `i` of a sliced-set `B`.
- `container::vector<container::vector<double>> sliced_footprints(net, tidxs, idxsz, is_batchable, batch, aux_list)` - `tables[B]` is `subset_footprints` with extents wrapped so that an index equal to `aux_list[k]` reports `min(extent, batch)` iff bit `k` of `B` is set, for every `B` in `[0, 2^m)`.
- `std::size_t leaf_volatile_mask(net, is_volatile_leaf)` - bit `i` set iff tensor `i` is a volatile leaf (mirrors `single_term_opt_impl`'s inline mask; 0 if predicate empty).
- per-subset open indices via `init_results` (already available) to compute, per bipartition, the contracted batchable set.
- `OptimizeOptions::is_batchable_index` (`std::function<bool(Index const&)>` = {}), `OptimizeOptions::batch_target_size` (`std::size_t` = 0), `ObjectiveFunction::DensePeakSizeBatched`.

- [ ] **Step 1: Enum + OptimizeOptions fields** (options.hpp). Add `DensePeakSizeBatched`; document it (per-index sliceable peak). Add:
```cpp
  /// Predicate marking an Index as living in a batchable space the runtime
  /// slices over (e.g. DF/RI aux; = the eval cache's accept_aux). Each distinct
  /// batchable index is sliced independently. Only consulted by
  /// ObjectiveFunction::DensePeakSizeBatched.
  std::function<bool(Index const&)> is_batchable_index = {};
  /// Shared slice size: a sliced batchable index contributes
  /// min(extent, batch_target_size). 0 disables the batched discount. Only
  /// consulted by DensePeakSizeBatched.
  std::size_t batch_target_size = 0;
```

- [ ] **Step 2: Failing test** (inside the sized-context block of `TEST_CASE("optimize")`). Use the aux/DF space tag the context provides for `g{..;..;F}` aux indices (confirm against the context; the asserted structure is the contract):
```cpp
SECTION("per-index batchability tables") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) {
    return ix.space().base_key() == L"F";  // aux/fitting space tag in this ctx
  };
  auto t0 = deserialize(L"g{a1;i1;F1}", {.def_perm_symm = Symmetry::Nonsymm});
  auto t1 = deserialize(L"g{a2;i1;F2}", {.def_perm_symm = Symmetry::Nonsymm});
  TensorNetwork net{std::vector<ExprPtr>{t0, t1}};
  container::svector<Index> targets;
  auto aux = opt::detail::batchable_index_list(net, is_batchable);
  REQUIRE(aux.size() == 2u);  // F1, F2 distinct
  auto tables = opt::detail::sliced_footprints(net, targets, idxsz,
                                               is_batchable, 1, aux);
  REQUIRE(tables.size() == 4u);                 // 2^2 sliced-sets
  // B=00 (none sliced) is the full footprint; B=11 (both) the all-sliced one.
  REQUIRE(tables[0b00][0b11] > tables[0b11][0b11]);  // full > all-sliced result
  // slicing only F1 (bit 0) shrinks the F1-leaf but not the F2-leaf.
  size_t f1bit = 0;  // aux[0]==F1 by appearance order
  REQUIRE(tables[size_t{1} << f1bit][0b01] < tables[0b00][0b01]);
}
```

- [ ] **Step 3: Run to verify it fails** (undeclared `batchable_index_list`/`sliced_footprints`).

- [ ] **Step 4: Implement** (single_term.hpp, `detail`):
```cpp
template <typename TIdxs>
container::vector<Index> batchable_index_list(
    TensorNetwork const& network, TIdxs const&,
    std::function<bool(Index const&)> const& is_batchable) {
  container::vector<Index> aux;
  if (!is_batchable) return aux;
  for (auto&& t : network.tensors()) {
    auto tp = std::dynamic_pointer_cast<Tensor>(t);
    for (auto&& ix :
         ranges::views::concat(tp->bra(), tp->ket(), tp->aux()))
      if (is_batchable(ix) && ranges::find(aux, ix) == ranges::end(aux))
        aux.push_back(ix);
  }
  return aux;
}

// tables[B][n] = subset_footprints with aux_list[k] sliced iff bit k of B set.
template <typename TIdxs, typename IdxToSz>
container::vector<container::vector<double>> sliced_footprints(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable, std::size_t batch,
    container::vector<Index> const& aux_list) {
  std::size_t const m = aux_list.size();
  container::vector<container::vector<double>> tables(std::size_t{1} << m);
  for (std::size_t B = 0; B < tables.size(); ++B) {
    auto extent = [&, B](Index const& ix) -> std::size_t {
      std::size_t e = idxsz(ix);
      if (is_batchable && is_batchable(ix)) {
        auto it = ranges::find(aux_list, ix);
        if (it != ranges::end(aux_list)) {
          std::size_t k = static_cast<std::size_t>(it - ranges::begin(aux_list));
          if (B & (std::size_t{1} << k)) return std::min(e, batch);
        }
      }
      return e;
    };
    tables[B] = subset_footprints(network, tidxs, extent);
  }
  return tables;
}
```
plus `leaf_volatile_mask` (identical body to the Phase-1 version mentioned in the prior plan; mirrors `single_term_opt_impl`'s mask).

- [ ] **Step 5: Run to verify it passes.** All existing `[optimize]` SECTIONs still pass.

- [ ] **Step 6: clang-format + commit**
```bash
clang-format --style=file -i core/optimize/single_term.hpp core/optimize/options.hpp tests/unit/test_optimize.cpp
git add core/optimize/single_term.hpp core/optimize/options.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: per-index batchability tables (aux list, sliced footprints)"
```

---

### Task 2: Per-index batch-aware oracle

**Files:** test only.

**Interfaces (test-local):** `double batched_min_peak(tables, contracted_at, persistent, nt, m)` where `tables[B][n]` are the sliced footprints (Task 1), `contracted_at[n]` packs (per subset, over all bipartitions is wrong - instead) the batchable-index bitmask **internal to `n`** = present in `n`'s tensors but not in `n`'s open indices (a superset used only to prune; the actual per-bipartition contracted set is derived inside), and `persistent[n] = (volatile_mask & n)==0`. The oracle enumerates trees x orders x, at each persistent contraction, which subset of the indices contracted *there* to batch; simulates the model-A peak where a live subset `mm` is sized `tables[ctx][mm]` with `ctx` = union of batched indices over frontiers strictly containing `mm`. (See spec 6.4.)

Implementation note: the cleanest oracle is the recursive form that threads a context `B` (bitmask of currently-sliced aux), mirroring spec 6.2 but enumerating (no memoization) and computing peak by the model-A pebble at the current context:
```cpp
// open-aux helper: bitmask of aux_list indices that are open indices of subset s
// (precomputed by the test from init_results + aux_list), passed in as open_aux[s].
static double oracle_rec(size_t n, size_t B,
    std::vector<std::vector<double>> const& T,        // T[ctx][subset]
    std::vector<size_t> const& open_aux, std::vector<char> const& persistent,
    size_t nt) {
  if (std::popcount(n) == 1) return T[B & open_aux[n]][n];
  auto sz = [&](size_t s, size_t ctx){ return T[ctx & open_aux[s]][s]; };
  auto Lof = [&](size_t s, size_t ctx){
    double r=0; for (size_t b=0;b<nt;b++) if (s&(size_t{1}<<b)) r+=sz(size_t{1}<<b,ctx);
    return r; };
  double best = std::numeric_limits<double>::max();
  for (size_t lp=(n-1)&n; lp; lp=(lp-1)&n) {
    size_t rp = n^lp; if (lp>rp) continue;
    // aux contracted at THIS node = open_aux[lp]|open_aux[rp] minus open_aux[n]
    size_t Acand = persistent[n] ? ((open_aux[lp]|open_aux[rp]) & ~open_aux[n]) : 0;
    for (size_t Ap = Acand; ; Ap = (Ap-1)&Acand) {   // every subset of Acand
      size_t C = B | Ap;
      double pl = oracle_rec(lp, C, T, open_aux, persistent, nt);
      double pr = oracle_rec(rp, C, T, open_aux, persistent, nt);
      double both = sz(lp,C)+sz(rp,C)+sz(n,B);
      double lpf = std::max({Lof(rp,C)+pl, sz(lp,C)+pr, both});
      double rpf = std::max({Lof(lp,C)+pr, sz(rp,C)+pl, both});
      best = std::min(best, std::min(lpf,rpf));
      if (Ap==0) break;
    }
  }
  return best;
}
```
This recursion threads the slice context (the independent, structural element) and computes peak by the pebble. It is "independent of the production DP" in that it does NOT memoize and re-derives every tree+context, but it shares the pebble formula. To ALSO guard the pebble itself, Task 4 adds a memory-simulation reconstruction check; and the `m=1` all-sliced corner is tied to the Phase-1 oracle in Task 3. NOTE the contract: `best = oracle_rec(full, 0, ...)`.

- [ ] **Step 1: Write the helper + hand-checked test.** Provide `batched_min_peak` wrapping `oracle_rec(full, 0, ...)`. Hand-checked case (2 leaves, one shared aux F contracted in the pair):
```cpp
TEST_CASE("batched_min_peak oracle (per-index)", "[optimize]") {
  // 2 leaves sharing aux F (m=1). open_aux: leaves have F open (bit0=1);
  // the pair has F internal (bit0=0). persistent everywhere.
  // T[ctx][subset]: full (ctx=0) leaves 4, pair result 2; sliced (ctx=1) leaves 2.
  std::vector<std::vector<double>> T(2, std::vector<double>(4, 0.0));
  T[0][0b01]=T[0][0b10]=4; T[0][0b11]=2;       // full
  T[1][0b01]=T[1][0b10]=2; T[1][0b11]=2;       // F sliced -> leaves halve
  std::vector<size_t> open_aux(4,0);
  open_aux[0b01]=1; open_aux[0b10]=1; open_aux[0b11]=0;  // F open in leaves, internal in pair
  std::vector<char> persistent(4,1);
  REQUIRE(batched_min_peak(T, open_aux, persistent, /*nt=*/2) == Catch::Approx(6.0));
  // not batched: 4+4+2=10; batch F at the pair: leaves sliced 2 -> 2+2+2=6. min 6.
}
```
Independently re-derive 6.0 before trusting. Add `<set>`/`<vector>`/`<algorithm>`/`<limits>` includes if needed.

- [ ] **Step 2: Run; debug oracle to the hand value.** If a correct impl yields a different number, report it with the hand derivation.

- [ ] **Step 3: clang-format + commit**
```bash
git add tests/unit/test_optimize.cpp && git commit -m "optimize: per-index batch-aware min-peak oracle for tests"
```

---

### Task 3: Multi-mode batched DP

**Files:** `core/optimize/single_term.hpp`, test.

**Interfaces:**
- `struct BatchedRes { double peak = inf; size_t lp=0, rp=0; bool lp_first=true; size_t aprime=0; };`
- `container::vector<BatchedRes> peak_dp_batched(net, tidxs, idxsz, is_batchable, batch, volatile_mask)` - flat table indexed `n*(2^m)+B`; `peak[n][B]` per spec 6.2.
- `double peak_cost_batched(net, tidxs, idxsz, is_batchable, batch, is_volatile_leaf)` - returns `peak[root][0]`.

- [ ] **Step 1: Failing tests** (sized-context block). Two SECTIONs:
```cpp
SECTION("DensePeakSizeBatched all-sliced corner equals Phase-1 batched peak") {
  using namespace sequant;
  auto idxsz = [](Index const& ix){ return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix){ return ix.space().base_key()==L"F"; };
  std::size_t const batch = 1;
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"})
    ts.push_back(deserialize(s, {.def_perm_symm = Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;
  auto aux = opt::detail::batchable_index_list(net, is_batchable);
  std::size_t const m = aux.size();
  auto pr = opt::detail::peak_dp_batched(net, targets, idxsz, is_batchable, batch,
                                         opt::detail::leaf_volatile_mask(net, {}));
  size_t root = (size_t{1} << ts.size()) - 1;
  size_t allK = (size_t{1} << m) - 1;
  double dp_allsliced = pr[root * (size_t{1} << m) + allK].peak;
  auto be = opt::detail::batched_extent(idxsz, is_batchable, batch);  // from prior plan, or inline a wrapper slicing all aux
  double phase1 = opt::detail::peak_cost(net, targets, be);
  REQUIRE(dp_allsliced == Catch::Approx(phase1));
}

SECTION("DensePeakSizeBatched objective matches per-index oracle") {
  using namespace sequant;
  auto idxsz = [](Index const& ix){ return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix){ return ix.space().base_key()==L"F"; };
  std::size_t const batch = 1;
  std::vector<std::vector<std::wstring>> nets = {
    {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},          // shared F1
    {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"},          // two distinct aux
  };
  for (auto const& spec : nets) {
    std::vector<ExprPtr> ts;
    for (auto s : spec) ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
    TensorNetwork net{ts}; container::svector<Index> targets;
    auto aux = opt::detail::batchable_index_list(net, is_batchable);
    std::size_t const m = aux.size();
    auto vmask = opt::detail::leaf_volatile_mask(net, {});
    auto tables = opt::detail::sliced_footprints(net, targets, idxsz,
                                                 is_batchable, batch, aux);
    // build open_aux[s] from init_results open indices vs aux list (test-side helper)
    std::vector<size_t> open_aux = build_open_aux(net, targets, aux);  // see note
    std::vector<std::vector<double>> T(tables.begin(), tables.end());
    std::vector<char> persistent(size_t{1} << ts.size());
    for (size_t n=0;n<persistent.size();++n) persistent[n]=((vmask&n)==0)?1:0;
    double oracle = batched_min_peak(T, open_aux, persistent, ts.size());
    double dp = opt::detail::peak_cost_batched(net, targets, idxsz,
                                               is_batchable, batch, {});
    REQUIRE(dp == Catch::Approx(oracle));
  }
}
```
`build_open_aux(net, targets, aux)`: a small test helper that runs the same open-index computation (`init_results`) and, per subset, sets bit `k` iff `aux[k]` is in that subset's open indices. (The implementer may expose a `detail::subset_open_aux` and reuse it in both the DP and the test.)

- [ ] **Step 2: Run to verify it fails** (undeclared `peak_dp_batched`).

- [ ] **Step 3: Implement the multi-mode DP** per spec 6.2. Key structure:
```cpp
// precompute: aux=batchable_index_list; m; tables=sliced_footprints; open_aux[n]
//   (bitmask of aux open in n) and Ltab[B][n] (leaf-sum of tables[B]).
// pr sized (2^nt)*(2^m), indexed idx(n,B)=n*(2^m)+B.
// for n in increasing order:
//   if popcount(n)==0: all B -> peak 0
//   if popcount(n)==1: pr[idx(n,B)].peak = tables[B & open_aux[n]][n]
//   else for each B in [0,2^m):
//     for each bipartition (lp,rp):
//       Acand = persistent(n) ? ((open_aux[lp]|open_aux[rp]) & ~open_aux[n]) : 0
//       for Ap = every subset of Acand:
//         C = B | Ap
//         sz(s,ctx)=tables[ctx & open_aux[s]][s];  L(s,ctx)=Ltab over ctx
//         both = sz(lp,C)+sz(rp,C)+sz(n,B)
//         lpf = max(L(rp,C)+pr[idx(lp,C)].peak, sz(lp,C)+pr[idx(rp,C)].peak, both)
//         rpf = max(L(lp,C)+pr[idx(rp,C)].peak, sz(rp,C)+pr[idx(lp,C)].peak, both)
//         cand = min(lpf,rpf); if < pr[idx(n,B)].peak: record peak,lp,rp,lp_first,aprime=Ap
// peak_cost_batched: vmask -> peak_dp_batched -> pr[idx(root,0)].peak
```
Implement `subset_open_aux(net, tidxs, aux_list)` returning `open_aux[n]`. Use the `&open_aux[s]` mask in `sz`/`L` so `tables` is always indexed by a context restricted to what is actually open in the sized subset (mirrors the oracle).

- [ ] **Step 4: Run to verify it passes.** Both SECTIONs green; all prior `[optimize]` SECTIONs pass. A failure of the oracle SECTION localizes a DP bug (oracle is hand-anchored).

- [ ] **Step 5: clang-format + commit**
```bash
git add core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: multi-mode (per-index) DensePeakSizeBatched DP"
```

---

### Task 4: Reconstruction (with full numeric check), dispatch, public API

**Files:** `core/optimize/single_term.hpp`, `core/optimize/optimize.cpp`, test.

**Interfaces:** `EvalSequence single_term_opt_peak_batched_impl(net, tidxs, idxsz, is_batchable, batch, is_volatile_leaf)`; `optimize(expr, {DensePeakSizeBatched, is_batchable_index, batch_target_size, ...})`.

- [ ] **Step 1: Failing tests.** Include the **full numeric reconstruction check**: a tree-walk simulator that recomputes the chosen reconstruction's model-A peak by direct memory simulation (independent of the DP's max/+) and asserts it equals `peak_cost_batched`.
```cpp
SECTION("DensePeakSizeBatched reconstruction achieves the optimum (numeric)") {
  using namespace sequant;
  auto idxsz = [](Index const& ix){ return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix){ return ix.space().base_key()==L"F"; };
  std::size_t const batch = 1;
  for (auto const& spec : std::vector<std::vector<std::wstring>>{
         {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},
         {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"}}) {
    std::vector<ExprPtr> ts;
    for (auto s : spec) ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
    TensorNetwork net{ts}; container::svector<Index> targets;
    // recompute the chosen tree's peak by simulation over the pr back-pointers:
    double recon = opt::detail::reconstructed_batched_peak(
        net, targets, idxsz, is_batchable, batch, {});
    double dp = opt::detail::peak_cost_batched(net, targets, idxsz,
                                               is_batchable, batch, {});
    REQUIRE(recon == Catch::Approx(dp));
  }
}
SECTION("optimize() public API dispatches DensePeakSizeBatched") {
  using namespace sequant;
  auto expr = deserialize(L"g{a1;i1;F1} * g{a2;i1;F1} * g{a2;i2;F2}",
                          {.def_perm_symm = Symmetry::Nonsymm});
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = [](Index const& ix){ return ix.space().approximate_size(); };
  opts.is_batchable_index = [](Index const& ix){ return ix.space().base_key()==L"F"; };
  opts.batch_target_size = 1;
  auto optimized = optimize(expr, opts);
  REQUIRE(optimized);
  REQUIRE(count_tensor_leaves(optimized) == 3u);
}
```

- [ ] **Step 2: Run to verify it fails** (`reconstructed_batched_peak`/`single_term_opt_peak_batched_impl` missing; public path hits the exhaustiveness `SEQUANT_ASSERT`).

- [ ] **Step 3: Implement reconstruction + the numeric checker.**
  - `reconstructed_batched_peak`: run `peak_dp_batched`, then recursively walk the chosen back-pointers from `(root, B=0)` - at `(n,B)` read `aprime`, set `C = B | aprime`, recurse children at context `C` in `lp_first` order - and compute the subtree peak by **simulation** (model-A: evaluate first child to completion holding nothing extra, then the second holding the first's result, summing all-co-resident sizes at the active context, sizes from `tables[ctx & open_aux[.]][.]`). Return the root peak. This re-derives peak independent of the DP recurrence, so `recon == peak_cost_batched` validates both the DP and the chosen tree.
  - `single_term_opt_peak_batched_impl`: same back-pointer walk, emitting the `EvalSequence` (`concat(first,second), -1`) in `lp_first` order; child context `C = B | aprime`; start at `(root, 0)`.
- [ ] **Step 4: Dispatch.** Add the `DensePeakSizeBatched` arm to the detail `single_term_opt` (thread two new defaulted params `is_batchable_index`, `batch_target_size` through both `single_term_opt` overloads; assert `!subnet_cse`; update the `static_assert` to list all four). Add the `opt_pure_product` branch passing `opts.is_batchable_index`, `opts.batch_target_size`; pass the two new trailing args (ignored) to the other three `single_term_opt` calls; keep the final `SEQUANT_ASSERT` exhaustive.
- [ ] **Step 5: Run to verify it passes.** All SECTIONs green.
- [ ] **Step 6: clang-format + commit**
```bash
git add core/optimize/single_term.hpp core/optimize/optimize.cpp tests/unit/test_optimize.cpp
git commit -m "optimize: reconstruct (numeric-checked) + dispatch DensePeakSizeBatched"
```

---

## Self-review notes

- **Spec coverage:** monomial/per-B tables (6.1) -> Task 1 `sliced_footprints`; multi-mode `peak[n][B]` + per-index frontier `A'` (6.2) -> Task 3; objective `peak[root][empty]` -> `peak_cost_batched`; validation (6.4) -> all-sliced==Phase-1 (Task 3), objective==per-index oracle (Task 3), full numeric reconstruction (Task 4); plumbing (6.5) -> Task 1 fields + Task 4 dispatch.
- **Trust anchors:** the all-sliced corner ties to the proven Phase-1 `peak_cost`; the objective ties to the hand-anchored per-index oracle (6.0); the reconstruction ties to an independent memory simulation. The oracle threads the slice context structurally (the new element) though it shares the pebble; the reconstruction's memory simulation is the independent guard on the pebble itself.
- **No regression:** Phase-1 `peak_dp`/`peak_cost`/`single_term_opt_impl` and `DenseFLOPs`/`DenseSize`/`DensePeakSize` untouched; only a new arm + two defaulted params.
- **Risks to flag to reviewer:** (a) the aux predicate `base_key()==L"F"` must match the test context's DF/aux space; the implementer pins it, the asserted structure is fixed. (b) DP/oracle are `2^nt * 2^m` / `tree x 2^(aux-per-node)`; keep `nt<=4`, `m<=2` in tests. (c) `subset_open_aux` must be shared between DP, oracle (`build_open_aux`), and reconstruction so all three agree on what is open. (d) index identity for the aux list relies on `Index` equality/`find`; confirm two `F1` occurrences compare equal and `F1 != F2`.
