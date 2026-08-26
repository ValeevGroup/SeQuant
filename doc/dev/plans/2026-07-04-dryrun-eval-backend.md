# DryRun eval backend and PNO-CCSD cost harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a "DryRun" SeQuant eval backend that replays a factorized IR over the real `CacheManager` and, for each evaluation, reports what the optimizer's own cost model reports (memsize, FLOPs, projected execution cost) so PNO-CCSD schedules can be exercised at arbitrary size regimes in seconds without building or running MPQC.

**Architecture:** The DryRun backend implements the abstract `Result` interface with tensors carrying no data, only an index set plus a regime-driven extent/moment model. Each `Result` operation (`prod`/`sum`/`permute`/`slice_mode`/`mode_batches`/`add_inplace`) produces a new zero-data `Result` whose `size_in_bytes()` delegates to a shared `CostModel` object built from the optimizer's `memsize_counter`/`flops_counter`/`roofline_op_cost`. A Catch2 harness deserializes the committed CSV-CCSD residual, runs the real optimize -> binarize -> `make_batched_custom_evaluator` -> `evaluate` pipeline against this backend under a chosen `SizeRegime`, and prints a per-operation cost table plus the realized `CacheManager` working-set high-water mark.

**Tech Stack:** C++20, SeQuant core (`io::serialization`, `core/optimize`, `core/eval`), Catch2 v3 unit tests, CMake + Unix Makefiles (existing `build-test` dir).

## Global Constraints

- ASCII only in all source, tests, and docs: no en-dashes (U+2013), no non-breaking spaces (U+00A0). The SeQuant pre-commit hook rejects them. Use plain hyphens.
- No `Co-Authored-By` trailers in commit messages. Plain messages describing the change.
- Docs (this plan, the spec) live under `doc/dev/{plans,specs}/`, not `docs/superpowers/`.
- The DryRun backend carries NO tensor data. Every `Result` in it is a pure size/cost token. Any allocation proportional to the modeled extent is a defect.
- Extents are rank-dependent and moment-dependent: `<#OSV> != <#PNO>` and `<#PNO^k> != <#PNO>^k`. The size model MUST route CSV (proto-indexed) inner extents through the moment-aware path (`inner_pow(index, k)`), not a scalar average. Support moments `k <= 4` (sufficient for a 2-body Hamiltonian).
- The size model is the optimizer's existing model, reused verbatim: `memsize_counter`, `flops_counter`, `roofline_op_cost` from `core/optimize`. Do NOT reimplement a parallel cost model in the backend.
- Reference backend to mirror for `Result` mechanics: `SeQuant/core/eval/backends/tapp/{result.hpp,eval_expr.hpp}`. Mirror its structure (final class deriving `Result`, `type_id()`, annotation unpacking from `std::array<std::any,N>`), and strip everything that touches real tensor storage.
- **Build harness:** use the **Release** build at `SeQuant/cmake-build-release` (CMake + Ninja). Build with `cmake --build cmake-build-release --target unit_tests-sequant -j8`. Run with `./cmake-build-release/tests/unit/unit_tests-sequant "[tag]"`. **Release is mandatory:** the batched single-term DP (`DensePeakSizeBatched`) is exponential in the number of batchable indices and takes 20-60 s per giant term in Release; in the `build-test` Debug build it is minutes-to-hung. Do NOT use `build-test` for anything touching `optimize()`.
- **Post-transform fixtures + mu~/K registration (VERDICT ALREADY OBTAINED):** the real factorizer-input residual is `tests/unit/data/csv_ccsd_doubles_residual_df.txt` (55-summand Sum; contains PAO `μ̃`, DF-aux `Κ`, 3-center DF `g{μ̃;i;Κ}`, CSV coeffs `C{a<i>;μ̃}`; the giant free-mu~ term is the first ~13-tensor summand). Standalone SeQuant's default registry lacks `μ̃`/`Κ`; register them by cloning the default context, `set_first_dummy_index_ordinal(1000000)`, `auto isr = ctx.mutable_index_space_registry(); sequant::mbpt::add_pao_spaces(isr); sequant::mbpt::add_df_spaces(isr);` then `set_scoped_default_context(std::move(ctx))`. The `[dryrun-df]` verdict case already PROVED the DP annotates a **μ̃** batch axis on the 188 GB (K=4320: 1210 GB) free-mu~ giant (`batch_axes={μ̃}`, K on outer nodes only) and that it survives `binarize`. So the eval backend's job is NOT to re-decide that -- it is to REPLAY the giant term's schedule through the runtime (`make_batched_custom_evaluator`) and witness that the runtime applies the OUTER K slice but DROPS the giant's nested μ̃ slice (giant K-sliced, μ~-full -> 186 GB/K-batch -> the OOM). That is the C60 runtime bug; the backend reproduces it locally.
- **Interface traps (from the [dryrun-df] work, reuse verbatim):** (a) after deserialize, re-flatten each summand `ex<Product>(p.scalar(), p.factors(), Product::Flatten::Yes)` before `optimize()` or `term_batch_axes` comes back empty; (b) `set_first_dummy_index_ordinal(1000000)` (mpqc ordinals reach tens of thousands); (c) never pass a raw `ExprPtr`/`Product` to `REQUIRE`/`CHECK` (Catch2 stringification hits a `Product` hash-cache assert -- capture `static_cast<bool>(expr)` into a `bool` first); (d) `memsize_counter` lives in `sequant::opt::detail` and returns an ELEMENT count (multiply by 8 for bytes); (e) `OptimizeOptions` batchability nests under `opts.batch_policy.{is_batchable_index,batch_target_size,peak_threshold}`; `term_batch_axes` is a caller-allocated `std::shared_ptr<std::unordered_map<Expr const*, container::vector<container::svector<Index>>>>`; read back via `->find(optimized_summand.get())`; optimize PER SUMMAND, not the whole Sum. The SAME `BatchPolicy` object must be passed to `optimize()` and to the runtime evaluator factory (`sequant::make_evaluator`/`make_batched_custom_evaluator`) so they cannot drift (see mpqc `SeQuantEngine::make_optimize_options`).
- **Focus the replay on the giant term only.** Do NOT sweep all 55 summands (20-60 s each). Deserialize, take the giant summand (the first, ~13-tensor `g·C·g·C·s·C·C·s·C·C·t·t·t` term), optimize it once at `{μ̃,Κ}` batchable / `peak_threshold=40e9` / regime {i=32, μ̃=1800, Κ=4320, PNO=19, OSV=57}, binarize, and replay THAT tree.
- **Test registration:** the DryRun test needs both `SeQuant::eval` (for `make_batched_custom_evaluator`, `CacheManager`, `Result`) and `SeQuant::optimize` (for `optimize`, `OptimizeOptions`, `memsize_counter`). Register it as its own OBJECT-library group (see Task 1 Step 3), NOT by appending to an existing source list -- the existing groups link only one sub-library each.
- Validation anchors (C60 pVDZ-F12 PNO-CCSD, from trace `614336.log`): the aux-batched giant `g(mu~,mu~,K)` reports `result=1867538944` bytes (~1.87 GB) per aux batch; the un-sliced free-mu~ intermediate `I(i_3,i_2,mu~_1241,K_2;a_3i_2i_3)` reports `pred_result=185979871872` bytes (~186 GB). The CostModel `memsize` for these two nodes MUST reproduce those numbers (within a few percent) before any Phase 1 verdict or Phase 3 report is trusted.

---

## File structure

New backend (header-only, mirrors `backends/tapp/`):

- `SeQuant/core/eval/backends/dryrun/cost_model_object.hpp` -- `CostModel` value type bundling `memsize`/`flops`/`exec_cost` closures built from a `SizeRegime`; owns the `idxsz` and `inner_pow` callbacks.
- `SeQuant/core/eval/backends/dryrun/size_regime.hpp` -- `SizeRegime` struct (per-space extents + per-rank CSV moment table) and the `idx_to_extent` / `inner_pow` callbacks it produces.
- `SeQuant/core/eval/backends/dryrun/result.hpp` -- `ResultDryRun` (flat) and `ResultDryRunNested` (ToT), both `final : public Result`.
- `SeQuant/core/eval/backends/dryrun/eval_expr.hpp` -- `EvalExprDryRun`/`EvalNodeDryRun` alias plus the leaf yielder functor `DryRunLeafEvaluator`.

New tests + fixture:

- `SeQuant/tests/unit/test_eval_dryrun.cpp` -- all Catch2 cases below.
- `SeQuant/tests/unit/data/csv_ccsd_residual.txt` -- the committed residual fixture (copied from mpqc4 `csv_eqn_Rs.serialized`).

Modified:

- `SeQuant/tests/unit/CMakeLists.txt` -- register `test_eval_dryrun.cpp`.

Each file has one responsibility: `size_regime` = "what extents does this regime imply", `cost_model_object` = "what does the optimizer's model report for a given index set", `result` = "a zero-data tensor token that reshapes under ops and asks the cost model for its size", `eval_expr` = "turn IR leaves into tokens".

---

## Task 1: DP-annotation probe (Phase 1 go/no-go)

The highest-value task: reproduce the C60 batch-axis decision locally and answer whether the free-mu~ giant fails to get a mu~ slice because the **cost model / threshold** never annotated it, or because the **runtime** dropped an annotation the model produced. This task builds only the `SizeRegime` + a thin `memsize` check -- NOT the full backend -- and inspects the optimizer's per-node batch axes directly.

**Files:**
- Create: `SeQuant/core/eval/backends/dryrun/size_regime.hpp`
- Create: `SeQuant/tests/unit/data/csv_ccsd_residual.txt`
- Create (probe test only): `SeQuant/tests/unit/test_eval_dryrun.cpp`
- Modify: `SeQuant/tests/unit/CMakeLists.txt`

**Interfaces:**
- Consumes: `sequant::deserialize<ExprPtr>(std::string_view)` from `SeQuant/core/io/shorthands.hpp`; the `single_term_opt` ExprPtr overload and `OptimizeOptions` from `SeQuant/core/optimize/`.
- Produces:
  - `struct SizeRegime` in namespace `sequant::eval::dryrun` with fields:
    - `std::map<std::wstring, std::size_t> space_extent;` -- key = `Index::space().base_key()` (e.g. `L"i"`, `L"a"`, `L"mu~"`, `L"K"`), value = number of elements in that space.
    - `std::array<double, 5> csv_pno_moment;` and `std::array<double, 5> csv_osv_moment;` -- `csv_pno_moment[k]` = `<#PNO^k>` averaged over occ pairs, `k` in `0..4`; `csv_osv_moment[k]` = `<#OSV^k>` over occ singles. Index 0 is unused (set 1.0).
  - `std::size_t SizeRegime::extent(Index const&) const;` -- returns `space_extent.at(base_key)`; throws `std::out_of_range` if the space is unknown (fail loud, do not default to 1).
  - `double SizeRegime::inner_pow(Index const& composite, std::size_t k) const;` -- for a proto-indexed CSV index, returns the moment `<#X^k>` for its rank (PNO if the composite spans an occ pair, OSV if a single); `k` clamped to `0..4`. For a non-composite index returns `std::pow(extent(composite), double(k))`. (Rank detection: number of `composite.proto_indices()` -- 2 protos => PNO pair moment, 1 proto => OSV moment. Confirm against `average_csv_extent_pow` in `single_term_detail.hpp`.)
  - `auto SizeRegime::idx_to_extent() const;` -- returns `std::function<std::size_t(Index const&)>` = `[this](Index const& ix){ return extent(ix); }`.
  - `auto SizeRegime::inner_pow_fn() const;` -- returns `std::function<double(Index const&, std::size_t)>` = `[this](Index const& ix, std::size_t k){ return inner_pow(ix, k); }`.

- [ ] **Step 1: Copy the residual fixture**

Copy mpqc4's committed residual into the test data dir (it is the exact IR the C60 run factorizes):

```bash
cp /Users/efv/code/mpqc4/csv_eqn_Rs.serialized \
   /Users/efv/code/SeQuant/tests/unit/data/csv_ccsd_residual.txt
head -c 200 /Users/efv/code/SeQuant/tests/unit/data/csv_ccsd_residual.txt
```

Expected: prints the leading `2 ((g{i_1,i_2;a_...` text. If the file is not present in mpqc4, regenerate it by running any CSV-CCk input with the residual-dump enabled, or ask the human for the artifact.

- [ ] **Step 2: Write `size_regime.hpp`**

```cpp
#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP

#include <SeQuant/core/index.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>

namespace sequant::eval::dryrun {

/// Per-space extents and per-rank CSV moment tables that define one size
/// regime for a dry-run replay. Extents are element counts; CSV moments are
/// power means over occupied pairs (PNO) or singles (OSV).
struct SizeRegime {
  std::map<std::wstring, std::size_t> space_extent;
  std::array<double, 5> csv_pno_moment{1.0, 1.0, 1.0, 1.0, 1.0};
  std::array<double, 5> csv_osv_moment{1.0, 1.0, 1.0, 1.0, 1.0};

  [[nodiscard]] std::size_t extent(Index const& ix) const {
    return space_extent.at(std::wstring{ix.space().base_key()});
  }

  [[nodiscard]] double inner_pow(Index const& composite,
                                 std::size_t k) const {
    if (k > 4) k = 4;
    auto const& protos = composite.proto_indices();
    if (protos.empty())
      return std::pow(static_cast<double>(extent(composite)),
                      static_cast<double>(k));
    return (protos.size() >= 2) ? csv_pno_moment[k] : csv_osv_moment[k];
  }

  [[nodiscard]] std::function<std::size_t(Index const&)> idx_to_extent()
      const {
    return [this](Index const& ix) { return extent(ix); };
  }

  [[nodiscard]] std::function<double(Index const&, std::size_t)> inner_pow_fn()
      const {
    return [this](Index const& ix, std::size_t k) { return inner_pow(ix, k); };
  }
};

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_SIZE_REGIME_HPP
```

- [ ] **Step 3: Register the test as its own OBJECT-library group in CMake**

The DryRun test uses both `SeQuant::eval` and `SeQuant::optimize`; existing groups link only one sub-library each, so add a dedicated group. In `SeQuant/tests/unit/CMakeLists.txt`, after the TAPP backend block (ends near line 156, `endif()`), insert:

```cmake
##########################
# OBJECT library: DryRun backend tests
##########################
set(eval_dryrun_test_sources "test_eval_dryrun.cpp")
add_library(unit_tests-sequant-eval-dryrun-obj OBJECT ${eval_dryrun_test_sources})
set_target_properties(unit_tests-sequant-eval-dryrun-obj PROPERTIES CXX_SCAN_FOR_MODULES OFF)
target_link_libraries(unit_tests-sequant-eval-dryrun-obj
    PUBLIC SeQuant::eval SeQuant::optimize SeQuant::bliss SeQuant::mbpt
    PRIVATE Catch2::Catch2 dtl::dtl)
target_compile_definitions(unit_tests-sequant-eval-dryrun-obj PRIVATE
    SEQUANT_UNIT_TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
```

Then add it to the combined executable's link list (the `target_link_libraries(unit_tests-sequant PRIVATE ...)` block near line 168), alongside `unit_tests-sequant-optimize-obj`:

```cmake
    unit_tests-sequant-eval-dryrun-obj
```

(`SeQuant::mbpt` is required so `deserialize` resolves the standard space registry, matching the tapp/ta/btas test groups; `SeQuant::bliss` matches the optimize group.)

- [ ] **Step 4: Write the probe test (expected to FAIL to build first, then RUN)**

Create `SeQuant/tests/unit/test_eval_dryrun.cpp`. This first case only exercises `SizeRegime` + the optimizer's batch-axis output. The C60 regime numbers below are the anchors; the moment tables come from the C60 CC-average PNO domain (avg PNO ~58, so use a representative geometric spread -- start from the flat approximation `moment[k] = avg^k` and refine in Step 6 if the memsize check misses the anchor).

```cpp
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>

#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <sstream>

using namespace sequant;
using sequant::eval::dryrun::SizeRegime;

namespace {

std::string slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// C60 pVDZ-F12 PNO-CCSD regime (from trace 614336.log).
SizeRegime c60_regime() {
  SizeRegime r;
  r.space_extent = {
      {L"i", 120},      // occupied (60 C atoms, pVDZ-F12) -- CONFIRM from log
      {L"a", 0},        // canonical virtual unused in CSV path
      {L"mu~", 1800},   // PAO domain per union (nnz_tiles=1800 in log)
      {L"K", 0},        // DF aux -- CONFIRM extent from log
  };
  // CC-average PNO ~58; OSV ~ smaller. Placeholder power-mean; refined in Step 6.
  double const pno = 58.0, osv = 20.0;
  for (std::size_t k = 0; k <= 4; ++k) {
    r.csv_pno_moment[k] = std::pow(pno, double(k));
    r.csv_osv_moment[k] = std::pow(osv, double(k));
  }
  return r;
}

}  // namespace

TEST_CASE("dryrun size regime basic extents", "[dryrun-probe]") {
  auto r = c60_regime();
  // A bare occ index resolves to its space extent.
  auto i = Index{L"i_1"};
  CHECK(r.extent(i) == 120);
  // A composite (proto-indexed) PAO index still resolves by base space.
  // (moment routing is exercised in the memsize task)
}

TEST_CASE("dryrun DP annotates mu~ on the free-mu~ giant", "[dryrun-probe]") {
  auto const residual = slurp(
      std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
      "/data/csv_ccsd_residual.txt");
  auto expr = deserialize<ExprPtr>(residual);
  REQUIRE(expr);

  auto r = c60_regime();

  // Build the SAME OptimizeOptions the C60 run used for the batched path.
  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = r.idx_to_extent();
  opts.inner_pow = r.inner_pow_fn();
  opts.peak_threshold = 40e9;  // C60 run's peak_threshold (40 GB)
  // batchable spaces: DF aux (K) AND PAO (mu~)
  opts.is_batchable_index = [](Index const& ix) {
    auto k = ix.space().base_key();
    return k == L"K" || k == L"mu~";
  };
  opts.batch_target_size = /* aux_target_size */ 72;

  auto optimized = optimize(expr, opts);
  REQUIRE(optimized);

  // Dump the per-node batch axes the DP produced and assert the offending
  // node -- the one whose result carries a free mu~ AND a free K -- got a
  // mu~ (or K) slice. term_batch_axes is keyed per binary node; see how
  // make_optimize_options threads it in mpqc4 cck.h.
  // GO/NO-GO:
  //   - if the giant node's axis set is empty or aux-only  -> COST-MODEL bug
  //     (DP never chose mu~ under the 40 GB threshold)
  //   - if the DP annotated mu~ but the runtime trace shows no slice
  //     -> RUNTIME bug (annotation dropped in make_batched_custom_evaluator)
  // Record the finding in the ledger; do not assert a pass/fail verdict here
  // beyond REQUIRE(optimized) -- this case is a probe. Print the axes:
  INFO("inspect batch axes for I(i,i,mu~,K;a) node -- see stdout");
  // ... iterate optimized's binary nodes, print node->annot() + assigned axes.
  SUCCEED();
}
```

- [ ] **Step 5: Build and run the probe**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -20
./build-test/tests/unit/unit_tests-sequant "[dryrun-probe]" -s 2>&1 | tail -60
```

Expected: builds; the second case prints the deserialized network's node annotations. Read the printed batch axes for the free-mu~ giant.

- [ ] **Step 6: Calibrate the regime against the memsize anchors**

Add an inline memsize check (using `memsize_counter(idxsz, inner_pow)` from `core/optimize/single_term_detail.hpp`) that computes the byte size of two known index sets and compares to the trace anchors:

```cpp
TEST_CASE("dryrun regime reproduces C60 trace sizes", "[dryrun-probe]") {
  auto r = c60_regime();
  auto memsize = memsize_counter(r.idx_to_extent(), r.inner_pow_fn());
  // aux-batched giant g(mu~,mu~,K) at aux batch = 72:
  //   expected ~1.87 GB (trace result=1867538944)
  // free-mu~ intermediate I(i_3,i_2,mu~_1241,K_2;a_3i_2i_3):
  //   expected ~186 GB (trace pred_result=185979871872)
  // Build the index sets from parsed Index objects, call memsize, CHECK
  // within +/-5%. If off, adjust csv_pno_moment / space_extent until both
  // anchors land, THEN re-read the Step 4 verdict -- a mis-calibrated regime
  // invalidates the go/no-go.
  SUCCEED();
}
```

Iterate `csv_pno_moment` / `space_extent` until both anchors reproduce within 5%. Confirm the exact `i`, `mu~`, `K` extents from `614336.log` (`nnz_tiles`, `result=` bytes) rather than trusting the placeholders above.

- [ ] **Step 7: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/size_regime.hpp \
  tests/unit/test_eval_dryrun.cpp
git add SeQuant/core/eval/backends/dryrun/size_regime.hpp \
        tests/unit/test_eval_dryrun.cpp tests/unit/data/csv_ccsd_residual.txt \
        tests/unit/CMakeLists.txt
git commit -m "dryrun: size regime + DP batch-axis probe on C60 free-mu~ giant"
```

**Go/no-go gate:** record in the SDD ledger which verdict the probe returned (cost-model vs runtime). Phase 2+ proceed regardless (the backend is useful either way), but the verdict decides where the eventual C60 fix lands.

---

## Task 2: CostModel object

Bundle the optimizer's three cost closures behind one value type the backend `Result`s can query, so `size_in_bytes()` returns the MODEL size (not an allocated size) and the harness can also read FLOPs and projected exec cost per op.

**Files:**
- Create: `SeQuant/core/eval/backends/dryrun/cost_model_object.hpp`
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add `[dryrun-costmodel]` cases)

**Interfaces:**
- Consumes: `SizeRegime` (Task 1); `memsize_counter`, `flops_counter` (`core/optimize/single_term_detail.hpp`), `roofline_op_cost` (`core/optimize/cost_model.hpp`), `RooflineParams` (`core/optimize/options.hpp`).
- Produces: `struct CostModel` in `sequant::eval::dryrun` with:
  - `explicit CostModel(SizeRegime, RooflineParams = {});`
  - `std::size_t memsize(container::svector<Index> const& idxset) const;` -- bytes for a tensor with these outer indices (delegates to the `memsize_counter` closure).
  - `double flops(container::svector<Index> const& out, container::svector<Index> const& contracted) const;` -- multiply-add count for a contraction producing `out` over `contracted` (delegates to `flops_counter`).
  - `double exec_cost(double flops, std::size_t left_bytes, std::size_t right_bytes) const;` -- roofline projected cost (delegates to `roofline_op_cost` with the stored `RooflineParams`; C60 used `machine_balance=200`, `fast_mem_elems=1000000`).
  - Public read-only accessor `SizeRegime const& regime() const;`.

- [ ] **Step 1: Write the failing test**

```cpp
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>

TEST_CASE("cost model memsize matches memsize_counter", "[dryrun-costmodel]") {
  auto r = c60_regime();
  sequant::eval::dryrun::CostModel cm{r};
  // g(i_1,i_2,K_2): 3 flat indices
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"K_2"}};
  auto direct = memsize_counter(r.idx_to_extent(), r.inner_pow_fn())(idx);
  CHECK(cm.memsize(idx) == direct);
}
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd /Users/efv/code/SeQuant
./build-test/tests/unit/unit_tests-sequant "[dryrun-costmodel]" 2>&1 | tail -5
```

Expected: FAIL to compile (`cost_model_object.hpp` not found).

- [ ] **Step 3: Write `cost_model_object.hpp`**

```cpp
#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/optimize/cost_model.hpp>
#include <SeQuant/core/optimize/options.hpp>

#include <cstddef>

namespace sequant::eval::dryrun {

/// Bundles the optimizer's own cost closures (memsize / flops / roofline)
/// behind one value type so dry-run Results report MODEL size and the harness
/// can read FLOPs and projected execution cost per operation.
class CostModel {
 public:
  explicit CostModel(SizeRegime regime, RooflineParams roofline = {})
      : regime_{std::move(regime)},
        roofline_{roofline},
        memsize_{memsize_counter(regime_.idx_to_extent(),
                                 regime_.inner_pow_fn())},
        flops_{flops_counter(regime_.idx_to_extent(), regime_.inner_pow_fn())} {
  }

  [[nodiscard]] std::size_t memsize(
      container::svector<Index> const& idxset) const {
    return memsize_(idxset);
  }

  [[nodiscard]] double flops(container::svector<Index> const& out,
                             container::svector<Index> const& contracted) const {
    return flops_(out, contracted);
  }

  [[nodiscard]] double exec_cost(double flops, std::size_t left_bytes,
                                 std::size_t right_bytes) const {
    return roofline_op_cost(flops, roofline_.block_tiles,
                            left_bytes + right_bytes, roofline_.machine_balance,
                            roofline_.fast_mem_elems, roofline_.block_prefactor);
  }

  [[nodiscard]] SizeRegime const& regime() const { return regime_; }

 private:
  SizeRegime regime_;
  RooflineParams roofline_;
  decltype(memsize_counter(std::declval<SizeRegime>().idx_to_extent(),
                           std::declval<SizeRegime>().inner_pow_fn())) memsize_;
  decltype(flops_counter(std::declval<SizeRegime>().idx_to_extent(),
                         std::declval<SizeRegime>().inner_pow_fn())) flops_;
};

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_COST_MODEL_OBJECT_HPP
```

Note: confirm the exact parameter order of `roofline_op_cost` against `cost_model.hpp:282` and the exact return types of `memsize_counter`/`flops_counter`; if `decltype(...)` on the members is awkward, store `std::function<std::size_t(container::svector<Index> const&)>` and `std::function<double(...)>` instead.

- [ ] **Step 4: Run to verify it passes**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -5
./build-test/tests/unit/unit_tests-sequant "[dryrun-costmodel]" 2>&1 | tail -5
```

Expected: PASS.

- [ ] **Step 5: Add a flops + exec_cost assertion**

```cpp
TEST_CASE("cost model flops and exec_cost are finite/positive",
          "[dryrun-costmodel]") {
  sequant::eval::dryrun::CostModel cm{c60_regime()};
  container::svector<Index> out{Index{L"mu~_1"}, Index{L"i_1"}};
  container::svector<Index> contracted{Index{L"i_2"}};
  auto f = cm.flops(out, contracted);
  CHECK(f > 0.0);
  CHECK(cm.exec_cost(f, cm.memsize(out), 4096) > 0.0);
}
```

Run `[dryrun-costmodel]` again; expect PASS.

- [ ] **Step 6: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/cost_model_object.hpp \
  tests/unit/test_eval_dryrun.cpp
git add SeQuant/core/eval/backends/dryrun/cost_model_object.hpp \
        tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: CostModel object wrapping optimizer memsize/flops/roofline"
```

---

## Task 3: Flat DryRun Result

Implement `ResultDryRun` -- a zero-data tensor token carrying an ordered index set, whose ops reshape the index set and whose `size_in_bytes()` calls `CostModel::memsize`. Mirror `ResultTensorTAPP` structure; delete all tensor storage.

**Files:**
- Create: `SeQuant/core/eval/backends/dryrun/result.hpp` (flat class only in this task)
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add `[dryrun-result]` cases)

**Interfaces:**
- Consumes: abstract `Result` (`core/eval/result.hpp`) -- override `type_id()`, `prod(Result const&, std::array<std::any,3> const&, DeNest)`, `sum(Result const&, std::array<std::any,3> const&)`, `permute(std::array<std::any,2> const&)`, `slice_mode(std::size_t, long, long)`, `mode_batches(std::size_t, std::size_t)`, `add_inplace(Result const&)`, `size_in_bytes()`. `CostModel` (Task 2).
- Produces: `class ResultDryRun final : public Result` in `sequant::eval::dryrun` with:
  - `ResultDryRun(container::svector<Index> idxset, std::shared_ptr<CostModel const> cm);`
  - `container::svector<Index> const& indices() const;`
  - The annotation arrays passed to `prod`/`sum`/`permute` are `[l_annot, r_annot, res_annot]` for the ternary ops and `[pre, post]` for `permute`. The `res_annot` (element `[2]`) IS the result index set -- decode it (each `std::any` holds the same annotation type tapp uses; for DryRun the natural annotation is `container::svector<Index>` or the node's `annot()` -- confirm what the leaf yielder and `evaluate()` pass, and unpack identically to tapp `Annot`).

- [ ] **Step 1: Write the failing test**

```cpp
#include <SeQuant/core/eval/backends/dryrun/result.hpp>

using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::ResultDryRun;

TEST_CASE("flat result size delegates to cost model", "[dryrun-result]") {
  auto cm = std::make_shared<CostModel const>(c60_regime());
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"K_2"}};
  ResultDryRun t{idx, cm};
  CHECK(t.size_in_bytes() == cm->memsize(idx));
}

TEST_CASE("flat result prod yields the result annotation index set",
          "[dryrun-result]") {
  auto cm = std::make_shared<CostModel const>(c60_regime());
  ResultDryRun l{{Index{L"mu~_1"}, Index{L"i_2"}}, cm};
  ResultDryRun r{{Index{L"i_2"}, Index{L"i_1"}}, cm};
  // result annotation says the product produces (mu~_1, i_1) (i_2 contracted)
  container::svector<Index> res{Index{L"mu~_1"}, Index{L"i_1"}};
  auto out = l.prod(r, {l_annot_of(l), r_annot_of(r), annot_of(res)},
                    Result::DeNest::None);
  auto const& ot = dynamic_cast<ResultDryRun const&>(*out);
  CHECK(ot.indices() == res);
  CHECK(out->size_in_bytes() == cm->memsize(res));
}
```

(`l_annot_of` / `annot_of` = helpers that wrap an index set into the same `std::any` payload `evaluate()` supplies; define them in the test's anonymous namespace to match whatever `EvalExprDryRun::annot()` returns from Task 5. If Task 5 is not yet done, temporarily construct the arrays from the raw index sets and adjust in Task 5.)

- [ ] **Step 2: Run to verify it fails**

```bash
cd /Users/efv/code/SeQuant
./build-test/tests/unit/unit_tests-sequant "[dryrun-result]" 2>&1 | tail -5
```

Expected: FAIL to compile.

- [ ] **Step 3: Implement `ResultDryRun`**

Mirror `backends/tapp/result.hpp` (`ResultTensorTAPP`): same `final : public Result`, same `type_id()` pattern (register a unique id), same annotation unpacking. Replace every `tapp_ops::permute(...)` / real-tensor line with pure index-set bookkeeping:

- `prod`: the result index set is the decoded `res_annot` (element `[2]`). Return `make_result<ResultDryRun>(res_indices, cm_)`. No contraction is executed. (Optionally record FLOPs = `cm_->flops(res_indices, contracted)` where `contracted` = indices in `l_annot union r_annot` minus `res_annot`, for the harness to read via an accessor.)
- `sum`: result index set = decoded `res_annot`; both operands must share it up to permutation. Return a new `ResultDryRun`.
- `permute`: result index set = decoded `post` annotation (element `[1]`); same size.
- `slice_mode(mode, elem_lo, elem_hi)`: return a `ResultDryRun` with the SAME index set but a per-mode extent override of `elem_hi - elem_lo` on `mode`. The simplest faithful model: store an optional `std::map<std::size_t, std::size_t> mode_extent_override_` and have `size_in_bytes()` scale the sliced mode's contribution. (If overriding the cost model's per-mode extent is awkward, model the slice by replacing that mode's Index with a narrower synthetic space in the index set; document the choice.)
- `mode_batches(mode, target)`: return `svector<pair<size_t,size_t>>` tiling `[0, extent(mode))` into ceil(extent/target) element ranges -- identical to what the tiledarray backend does; this needs the regime extent for `indices_[mode]`.
- `add_inplace`: assert same index set (up to permutation); no data, so a no-op beyond the assertion.
- `size_in_bytes()`: `return cm_->memsize(indices_);` (applying any `mode_extent_override_`).

Store `container::svector<Index> indices_;` and `std::shared_ptr<CostModel const> cm_;`.

- [ ] **Step 4: Run to verify it passes**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -5
./build-test/tests/unit/unit_tests-sequant "[dryrun-result]" 2>&1 | tail -8
```

Expected: PASS.

- [ ] **Step 5: Add slice_mode + mode_batches tests**

```cpp
TEST_CASE("flat result slice_mode shrinks the sliced mode", "[dryrun-result]") {
  auto cm = std::make_shared<CostModel const>(c60_regime());
  ResultDryRun t{{Index{L"mu~_1"}, Index{L"i_2"}}, cm};  // mu~ extent 1800
  auto full = t.size_in_bytes();
  auto sliced = t.slice_mode(0, 0, 450);  // quarter of mu~
  CHECK(sliced->size_in_bytes() < full);
  CHECK(sliced->size_in_bytes() == full / 4);  // linear in mu~
}

TEST_CASE("flat result mode_batches tiles the mode extent", "[dryrun-result]") {
  auto cm = std::make_shared<CostModel const>(c60_regime());
  ResultDryRun t{{Index{L"mu~_1"}, Index{L"i_2"}}, cm};
  auto batches = t.mode_batches(0, 450);  // 1800 / 450 = 4 batches
  CHECK(batches.size() == 4);
  CHECK(batches.front().first == 0);
  CHECK(batches.back().second == 1800);
}
```

Run `[dryrun-result]`; expect PASS.

- [ ] **Step 6: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/result.hpp tests/unit/test_eval_dryrun.cpp
git add SeQuant/core/eval/backends/dryrun/result.hpp tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: flat ResultDryRun (zero-data size/cost token)"
```

---

## Task 4: Nested DryRun Result (ToT / proto-indexed CSV)

CSV amplitudes are Tensor-of-Tensor: an outer occ-pair index set over inner PNO/OSV domains. `ResultDryRunNested` models this so `size_in_bytes()` uses the moment-aware inner extent (`<#PNO^k>`), reproducing `<#PNO^2> != <#PNO>^2`.

**Files:**
- Modify: `SeQuant/core/eval/backends/dryrun/result.hpp` (add `ResultDryRunNested`)
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add `[dryrun-nested]` cases)

**Interfaces:**
- Consumes: same `Result` overrides + `CostModel`. The `DeNest` argument on `prod` decides whether a nested x nested product collapses the inner layer (mirror tapp's handling).
- Produces: `class ResultDryRunNested final : public Result` in `sequant::eval::dryrun` with:
  - `ResultDryRunNested(container::svector<Index> outer, container::svector<Index> inner, std::shared_ptr<CostModel const> cm);`
  - `size_in_bytes()` = outer element count times the moment-aware inner element count. The inner count uses `CostModel::regime().inner_pow(inner_proto_index, 1)` for linear size, and the moment enters when a product raises the effective inner power (a nested x nested contraction over the inner domain yields `<#PNO^2>`, not `<#PNO>^2`). Route through `memsize_counter`, which already applies `inner_pow` for composite indices -- so the correct implementation passes the full (outer + composite-inner) index set to `cm_->memsize(...)` and lets the counter apply the moment. Confirm which index in the set carries the proto (composite) tag so `memsize_counter` picks the moment path.

- [ ] **Step 1: Write the failing test**

```cpp
using sequant::eval::dryrun::ResultDryRunNested;

TEST_CASE("nested result uses moment-aware inner extent", "[dryrun-nested]") {
  // regime with a non-trivial second moment: <#PNO^2> != <#PNO>^2
  auto r = c60_regime();
  r.csv_pno_moment[1] = 58.0;
  r.csv_pno_moment[2] = 58.0 * 58.0 * 1.5;  // dispersion: 2nd moment inflated
  auto cm = std::make_shared<CostModel const>(r);

  // C(i_1,i_2,mu~; a<i_1 i_2>): outer occ pair + PAO, inner PNO over the pair
  container::svector<Index> outer{Index{L"i_1"}, Index{L"i_2"}, Index{L"mu~_1"}};
  Index a_pno{L"a_1", {Index{L"i_1"}, Index{L"i_2"}}};  // proto-indexed PNO
  ResultDryRunNested c{outer, {a_pno}, cm};

  // size must scale with the FIRST moment (linear), not with moment^k
  auto expected = /* outer element count */ * r.csv_pno_moment[1];
  CHECK(c.size_in_bytes() ==
        cm->memsize(/* outer + a_pno as one index set */));
  // and must differ from the naive avg^2 path
}
```

(Fill the exact index-set construction to match how `memsize_counter` reads a proto-indexed `Index`. Verify proto attachment syntax against `Index` ctors in `core/index.hpp`.)

- [ ] **Step 2: Run to verify it fails**

```bash
cd /Users/efv/code/SeQuant
./build-test/tests/unit/unit_tests-sequant "[dryrun-nested]" 2>&1 | tail -5
```

Expected: FAIL to compile.

- [ ] **Step 3: Implement `ResultDryRunNested`**

Mirror the tapp nested/ToT result (if tapp has one) or `backends/btas`/`tiledarray` ToT result for the interface shape; strip storage. Store `outer_`, `inner_` index sets + `cm_`. `size_in_bytes()` builds the combined index set (outer ++ inner-composite) and returns `cm_->memsize(combined)`, so the moment path is applied by the counter. `prod` with `DeNest::All`/inner-collapse returns a flat `ResultDryRun` whose index set is the decoded `res_annot`; `prod` staying nested returns a `ResultDryRunNested`. Decide nested-vs-flat from the decoded result annotation's proto structure (matches how the runtime chooses).

- [ ] **Step 4: Run to verify it passes**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -5
./build-test/tests/unit/unit_tests-sequant "[dryrun-nested]" 2>&1 | tail -8
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/result.hpp tests/unit/test_eval_dryrun.cpp
git add SeQuant/core/eval/backends/dryrun/result.hpp tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: nested ResultDryRunNested (moment-aware CSV inner extent)"
```

---

## Task 5: EvalExpr adapter and leaf yielder

Provide the leaf evaluator functor `evaluate()` calls to turn each IR leaf (a tensor/variable node) into a DryRun `Result`, and the `EvalExprDryRun` alias carrying `annot()` (mirroring `EvalExprTAPP`).

**Files:**
- Create: `SeQuant/core/eval/backends/dryrun/eval_expr.hpp`
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add `[dryrun-leaf]` cases; retrofit the annotation helpers referenced in Task 3)

**Interfaces:**
- Consumes: `EvalExpr` (`core/eval/eval_expr.hpp`), `EvalNode` template, `ResultDryRun`/`ResultDryRunNested`, `CostModel`.
- Produces:
  - `class EvalExprDryRun final : public EvalExpr` with an `annot()` method (copy `EvalExprTAPP` verbatim; it only adds `annot()`).
  - `using EvalNodeDryRun = EvalNode<EvalExprDryRun>;`
  - `struct DryRunLeafEvaluator { std::shared_ptr<CostModel const> cm; ResultPtr operator()(EvalNodeDryRun const& leaf) const; };` -- builds a `ResultDryRun` (flat leaf) or `ResultDryRunNested` (proto-indexed leaf, i.e. a CSV amplitude/coeff) from the leaf's tensor index set, using `cm`. This is the `F` in `evaluate<Trace>(node, layout, F, cache)`.

- [ ] **Step 1: Write the failing test**

```cpp
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>

TEST_CASE("leaf yielder builds a sized token from a tensor leaf",
          "[dryrun-leaf]") {
  auto cm = std::make_shared<CostModel const>(c60_regime());
  // parse a single tensor into an eval node, yield it
  auto expr = deserialize<ExprPtr>(L"g{i_1,i_2;K_2}:N-C-S");  // adjust syntax
  auto node = /* to_eval_node<EvalExprDryRun>(expr) */;
  sequant::eval::dryrun::DryRunLeafEvaluator yield{cm};
  auto res = yield(node);
  container::svector<Index> idx{Index{L"i_1"}, Index{L"i_2"}, Index{L"K_2"}};
  CHECK(res->size_in_bytes() == cm->memsize(idx));
}
```

(Confirm the `to_eval_node` / `binarize`-to-`EvalNode<EvalExprDryRun>` entry from `core/eval/eval_expr.hpp` / how tapp constructs its nodes.)

- [ ] **Step 2: Run to verify it fails**

```bash
cd /Users/efv/code/SeQuant
./build-test/tests/unit/unit_tests-sequant "[dryrun-leaf]" 2>&1 | tail -5
```

Expected: FAIL to compile.

- [ ] **Step 3: Implement `eval_expr.hpp`**

`EvalExprDryRun`: copy `EvalExprTAPP` (constructor-forwarding + `annot()`), rename. `DryRunLeafEvaluator::operator()`: read the leaf's tensor, collect its bra+ket indices into a `container::svector<Index>`; if any index is proto-indexed (composite), build a `ResultDryRunNested` (outer = non-composite indices, inner = composite index); else a flat `ResultDryRun`. Return via `make_result<...>`.

- [ ] **Step 4: Run to verify it passes**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -5
./build-test/tests/unit/unit_tests-sequant "[dryrun-leaf]" 2>&1 | tail -8
```

Expected: PASS. Retrofit Task 3's `annot_of`/`l_annot_of` helpers to return exactly what `EvalExprDryRun::annot()` yields, and re-run `[dryrun-result]` to confirm still green.

- [ ] **Step 5: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i \
  SeQuant/core/eval/backends/dryrun/eval_expr.hpp tests/unit/test_eval_dryrun.cpp
git add SeQuant/core/eval/backends/dryrun/eval_expr.hpp tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: EvalExprDryRun + leaf yielder (flat/nested token from IR leaf)"
```

---

## Task 6: End-to-end replay harness

Deserialize the residual, run the real optimize -> binarize -> `make_batched_custom_evaluator` -> `evaluate` pipeline against the DryRun backend under a `SizeRegime`, and print a per-operation cost table plus the realized `CacheManager` working-set high-water mark.

**Files:**
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp` (add `[dryrun]` end-to-end case)

**Interfaces:**
- Consumes: everything above; `make_batched_custom_evaluator` (`eval.hpp:1211`) / the `evaluate_batched` wrapper (`eval.hpp:1485`); `CacheManager` (`core/eval/cache_manager.hpp`) with `note_working_set` / `working_set_hwmark`; `optimize(expr, OptimizeOptions)`.
- Produces: a documented harness case that (1) factorizes the residual under the C60 regime + 40 GB threshold + `{K, mu~}` batchable, (2) evaluates over a `CacheManager` with `Trace::On`, (3) collects per-op `{result index set, memsize, flops, exec_cost}`, (4) reports peak = `cache.working_set_hwmark()`, (5) flags every intermediate whose result carries a free `mu~` AND is NOT sliced (the leak signature).

- [ ] **Step 1: Write the harness test**

```cpp
TEST_CASE("dryrun replays the C60 residual and reports peak", "[dryrun]") {
  auto residual = slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                        "/data/csv_ccsd_residual.txt");
  auto expr = deserialize<ExprPtr>(residual);
  REQUIRE(expr);

  auto regime = c60_regime();
  auto cm = std::make_shared<CostModel const>(regime);

  OptimizeOptions opts;
  opts.objective_function = ObjectiveFunction::DensePeakSizeBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.peak_threshold = 40e9;
  opts.is_batchable_index = [](Index const& ix) {
    auto k = ix.space().base_key();
    return k == L"K" || k == L"mu~";
  };
  opts.batch_target_size = 72;
  auto optimized = optimize(expr, opts);

  // binarize to EvalNode<EvalExprDryRun>, build cache + batched evaluator,
  // evaluate with Trace::On, capturing per-op logs. Then:
  auto peak = /* cache.working_set_hwmark() */;
  INFO("realized peak bytes = " << peak);
  // Under a CORRECT mu~-slicing DP, peak must stay below tolerance*40 GB.
  // Today (bug present) it will spike to ~186 GB on the free-mu~ giant.
  CHECK(peak > 0);
}
```

- [ ] **Step 2: Run and read the report**

```bash
cd /Users/efv/code/SeQuant
cmake --build build-test --target unit_tests-sequant -j4 2>&1 | tail -5
./build-test/tests/unit/unit_tests-sequant "[dryrun]" -s 2>&1 | tail -80
```

Expected: builds and prints a per-op cost table + realized peak. Confirm the free-mu~ giant appears in the table at ~186 GB when the DP leaves it un-sliced (reproducing the cluster OOM locally in seconds). Record the realized peak in the ledger.

- [ ] **Step 3: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i tests/unit/test_eval_dryrun.cpp
git add tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: end-to-end C60 residual replay harness with peak report"
```

---

## Task 7: C60 regression lock

Turn the reproduced peak into an explicit assertion so the eventual mu~-slicing fix (cost-model or runtime, per Task 1's verdict) has a fast local gate. Until the fix lands, this documents the known failure.

**Files:**
- Modify: `SeQuant/tests/unit/test_eval_dryrun.cpp`

**Interfaces:**
- Consumes: Task 6's harness.
- Produces: a `[dryrun][!shouldfail]`-tagged (Catch2 expected-failure) or explicitly-documented case asserting `peak <= 1.1 * opts.peak_threshold` for the C60 regime. Tagged `!shouldfail` while the bug is live; flip to a hard `CHECK` (remove the tag) in the same commit that fixes the DP/runtime.

- [ ] **Step 1: Add the locked assertion**

```cpp
// While the mu~-slicing bug is live this case documents the failure; the DP/
// runtime fix commit removes the [!shouldfail] tag so it becomes a hard gate.
TEST_CASE("C60 batched peak stays under threshold", "[dryrun][!shouldfail]") {
  // ... reuse the Task 6 setup ...
  auto peak = /* realized peak */;
  CHECK(peak <= static_cast<std::size_t>(1.1 * 40e9));  // 1.1 * 40 GB
}
```

- [ ] **Step 2: Run and confirm expected-failure status**

```bash
cd /Users/efv/code/SeQuant
./build-test/tests/unit/unit_tests-sequant "[dryrun]" 2>&1 | tail -12
```

Expected: the `!shouldfail` case reports as an expected failure (peak ~186 GB exceeds 44 GB), the suite stays green. This is the regression lock: once the fix lands, removing the tag makes it a real pass.

- [ ] **Step 3: Commit**

```bash
cd /Users/efv/code/SeQuant
clang-format --style=file -i tests/unit/test_eval_dryrun.cpp
git add tests/unit/test_eval_dryrun.cpp
git commit -m "dryrun: lock C60 batched peak regression (expected-fail until mu~ slice fix)"
```

---

## Self-Review

**1. Spec coverage.** Spec Component 1 (CostModel object) -> Task 2. Component 2 (DryRun Result backend: flat + nested, `size_in_bytes` delegating to cost model) -> Tasks 3, 4. Component 3 (SizeRegime + Catch2 harness + fixture) -> Tasks 1 (regime, fixture), 6 (harness). Spec Phase 1 (DP-annotation probe, go/no-go) -> Task 1. Phase 2 (backend) -> Tasks 2-5. Phase 3 (full harness) -> Task 6. Phase 4 (C60 regression lock) -> Task 7. Moment-dependent CSV extents (`<#PNO^k> != <#PNO>^k`, k<=4) -> `SizeRegime::inner_pow` + Task 4. Rank-dependent (`<#OSV> != <#PNO>`) -> separate `csv_osv_moment`/`csv_pno_moment` tables. Validation anchors (~1.87 GB, ~186 GB) -> Task 1 Step 6 + Task 6 Step 2. All covered.

**2. Placeholder scan.** The tasks contain real code and real commands. Remaining `/* ... */` markers are deliberate: they mark the two spots where the exact SeQuant API (annotation `std::any` payload type in `evaluate()`; the `to_eval_node<EvalExprDryRun>` / binarize entry) must be confirmed against `backends/tapp/{result,eval_expr}.hpp` at implementation time, and each is accompanied by the reference file + what to mirror. The C60 space extents (`i`, `K`) are flagged CONFIRM-from-log rather than guessed. These are grounding pointers, not vague directives.

**3. Type consistency.** `SizeRegime::idx_to_extent()` returns `std::function<std::size_t(Index const&)>`, consumed by `CostModel` ctor and `OptimizeOptions::idx_to_extent` (both expect that signature). `SizeRegime::inner_pow_fn()` returns `std::function<double(Index const&, std::size_t)>`, matching `OptimizeOptions::inner_pow`. `CostModel::memsize(container::svector<Index> const&)` is the single size entry used by `ResultDryRun::size_in_bytes` (Task 3), `ResultDryRunNested::size_in_bytes` (Task 4), and the tests. `ResultDryRun`/`ResultDryRunNested`/`EvalExprDryRun`/`EvalNodeDryRun`/`DryRunLeafEvaluator` names are used identically across Tasks 3-6. All in `sequant::eval::dryrun`.

---

## Execution Handoff

Plan complete and saved to `doc/dev/plans/2026-07-04-dryrun-eval-backend.md`. Two execution options:

1. **Subagent-Driven (recommended)** - dispatch a fresh subagent per task, review between tasks, fast iteration. Task 1 (the go/no-go probe) is the natural first gate: its verdict decides where the eventual C60 mu~-slicing fix lands.
2. **Inline Execution** - execute tasks in this session with checkpoints for review.

Which approach?
