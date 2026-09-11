# Cost-model Phase 3: generic single-term-opt driver + CostModel types

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Factor the duplicated subset-lattice + bipartition orchestration of SeQuant's four single-term-optimization objectives into one compile-time generic driver `run_single_term_opt<Model>`, with each objective expressed as a `CostModel` type, **without changing any optimization result**.

**Architecture:** A generic driver owns the subset lattice + bipartition loop; each objective is a model type with associated `State`/`Context` and `leaf/init/relax/finalize/reconstruct`. Production dispatch (`opt_pure_product`, detail `single_term_opt<Metric>`) routes through the models; the existing standalone DP/cost functions stay as independent reference oracles, and per-objective equivalence tests prove the model path reproduces them.

**Tech Stack:** C++20, SeQuant `core/optimize`, Catch2 (`tests/unit/test_optimize.cpp`).

## Global Constraints

- Style: Google `.clang-format`, 80-col, 2-space, no tabs. `clang-format --style=file -i <file>` before each commit.
- No en-dashes (U+2013) or non-breaking spaces (pre-commit hook rejects them); ASCII hyphen only.
- No `Co-Authored-By` trailers.
- Branch: `feature/cost-model-batch-aware` (continue on it).
- **Behavior-preserving (hard):** no existing `[optimize]` assertion may change. The full `[optimize]` suite passing after each task is the primary gate.
- **Scope:** cost side only. NO evaluator face (`evaluate<Backend>`), NO mpqc, NO removal of the existing standalone DP/cost functions (`single_term_opt_impl`, `peak_dp`/`peak_cost`, `single_term_opt_peak_impl`, `peak_dp_batched`/`peak_cost_batched`, `single_term_opt_peak_batched_impl`) - those remain as reference oracles; their removal is a deliberate fast-follow. The brute-force oracles (`brute_force_min_peak`, `batched_min_peak`) also stay.
- Reference: spec `doc/dev/specs/2026-06-18-cost-model-phase3-generic-driver-design.md`.
- Models are compile-time TYPES; the `ObjectiveFunction` enum stays as the runtime selector.

## File structure

- `core/optimize/cost_model.hpp` (new) - the `CostModel` concept, `solve_single_term`/`run_single_term_opt` driver, and the four built-in model types (`AdditiveModel`, `PeakModel`, `PeakBatchedModel`). Includes `single_term.hpp` for the reused helpers.
- `core/optimize/single_term.hpp` - unchanged helpers (`subset_footprints`, `sliced_footprints`, `subset_open_aux`, `build_subnet_metadata`, `init_results`, the cost/footprint counters, the standalone DP/cost functions kept as oracles). The detail `single_term_opt<Metric>` arms are rerouted to build a model + call the driver.
- `core/optimize/optimize.cpp` - `opt_pure_product` unchanged in behavior (still dispatches the enum), now via the model path through `single_term_opt<Metric>`.
- `tests/unit/test_optimize.cpp` - existing tests unchanged; add per-objective equivalence SECTIONs + the concept/custom-model checks.

Build/run each test step:
```bash
cmake --build build-test --target unit_tests-sequant -j
./build-test/tests/unit/unit_tests-sequant "[optimize]"
```

---

### Task 1: Generic driver + `CostModel` concept + `AdditiveModel`

**Files:**
- Create: `core/optimize/cost_model.hpp`
- Modify: `core/optimize/single_term.hpp` (reroute the `DenseFLOPs`/`DenseSize` arms of detail `single_term_opt<Metric>`)
- Test: `tests/unit/test_optimize.cpp`

**Interfaces:**
- Produces: `sequant::opt::detail::CostModel` concept; `template<class M> container::vector<typename M::State> solve_single_term(M const&, TensorNetwork const&, TIdxs const&)`; `template<class M> EvalSequence run_single_term_opt(M const&, TensorNetwork const&, TIdxs const&)`; `struct AdditiveModel<CostFn>` with nested `State`/`Context`. The detail `single_term_opt<DenseFLOPs|DenseSize>` returns `run_single_term_opt(AdditiveModel{...}, ...)`.

- [ ] **Step 1: Write the failing equivalence test**

Add a SECTION to `TEST_CASE("optimize")` (inside the sized-context block). It asserts the model-via-driver yields the SAME `EvalSequence` as the existing `single_term_opt_impl` for FLOPs and Size, including the CSE, volatile-weight, and footprint-weight paths:
```cpp
SECTION("AdditiveModel via driver == single_term_opt_impl") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  std::vector<std::wstring> spec = {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;a3}",
                                    L"g{a3;i2}"};
  std::vector<ExprPtr> ts;
  for (auto s : spec) ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
  TensorNetwork net{ts};
  container::svector<Index> targets;
  for (bool cse : {false, true})
    for (double fw : {0.0, 1000.0}) {
      // OLD path:
      auto old_flops = opt::detail::single_term_opt_impl(
          net, targets, opt::detail::flops_counter(idxsz), cse,
          /*volatile_mask=*/0u, /*volatile_weight=*/1.0,
          opt::detail::footprint_counter(idxsz), fw);
      // NEW path via the model+driver:
      opt::detail::AdditiveModel model{
          opt::detail::flops_counter(idxsz), opt::detail::footprint_counter(idxsz),
          /*volatile_mask=*/0u, /*volatile_weight=*/1.0, fw, cse};
      auto neww = opt::detail::run_single_term_opt(model, net, targets);
      REQUIRE(neww == old_flops);
    }
}
```
(`EvalSequence` is a `container::vector<int>`; `==` compares element-wise.)

- [ ] **Step 2: Run to verify it fails** (compile error: `AdditiveModel`/`run_single_term_opt` undeclared).

- [ ] **Step 3: Write `cost_model.hpp` - the concept + driver + `AdditiveModel`**

Create `core/optimize/cost_model.hpp`. The driver (the spec sketch, split into `solve` + `run`):
```cpp
#ifndef SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP
#define SEQUANT_CORE_OPTIMIZE_COST_MODEL_HPP

#include <SeQuant/core/optimize/single_term.hpp>  // helpers + EvalSequence + OptRes

#include <bit>

namespace sequant::opt::detail {

/// Fills the per-subset State table bottom-up via the model's hooks.
template <class Model, typename TIdxs>
container::vector<typename Model::State> solve_single_term(
    Model const& m, TensorNetwork const& network, TIdxs const& tidxs,
    typename Model::Context& ctx) {
  auto const nt = network.tensors().size();
  container::vector<typename Model::State> st(size_t{1} << nt);
  for (size_t n = 1; n < st.size(); ++n) {
    if (std::popcount(n) == 1) {
      st[n] = m.leaf(ctx, n);
      continue;
    }
    typename Model::State acc = m.init(ctx, n);
    for (auto&& [lp, rp] : bits::bipartitions(n))
      if (lp != 0 && rp != 0) m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);
    st[n] = std::move(acc);
    m.finalize(ctx, n, st);
  }
  return st;
}

/// Generic single-term optimization: build context, solve, reconstruct.
template <class Model, typename TIdxs>
EvalSequence run_single_term_opt(Model const& m, TensorNetwork const& network,
                                 TIdxs const& tidxs) {
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  typename Model::Context ctx = m.build_context(network, tidxs);
  auto st = solve_single_term(m, network, tidxs, ctx);
  return m.reconstruct(ctx, st);
}

}  // namespace sequant::opt::detail

#endif
```
Then `AdditiveModel` in the same file. Its `State`/`Context`/methods relocate the body of `single_term_opt_impl` (single_term.hpp ~389-460): `Context` holds `container::vector<OptRes> results` (from `init_results`, providing per-subset open indices + the `ops`/`lp`/`rp`/`subnets` scratch), `cost_fn`, `footprint_fn`, `footprint_weight`, `volatile_mask`, `volatile_weight`, the bool `subnet_cse`, and (when CSE) `meta_ids` + `unique_meta_costs` from `build_subnet_metadata`. `State` mirrors the relevant `OptRes` fields the DP mutates: `{double ops; size_t lp, rp; container::vector<size_t> subnets;}`.
- `build_context`: `init_results(network, tidxs, results)`; if `subnet_cse` run `build_subnet_metadata`.
- `leaf(ctx,n)`: `{0.0, 0, 0, {}}` (ops 0; the singleton sequence is implicit in reconstruct).
- `init(ctx,n)`: `{max, 0, 0, {}}`.
- `relax(ctx,n,lp,rp,lp_st,rp_st,acc)`: the body of the bipartition loop (single_term.hpp ~304-353) - compute `w`, `fp`, `new_cost` (CSE vs non-CSE branch using `unique_meta_costs`/`combined_subnets`), and if `< acc.ops` update `acc.ops/lp/rp/subnets`.
- `finalize(ctx,n,st)`: the post-loop CSE bookkeeping (~356-371) updating `ctx.unique_meta_costs[ctx.meta_ids[n]]` and `st[n].subnets`; no-op when `!subnet_cse`.
- `reconstruct(ctx,st)`: relocate the sequence assembly (~374-383) - walk `st[*].lp/rp` from the full set, emitting `(lseq, rseq)` ordered by `lseq[0] < rseq[0]`, then `-1`; build a `seq` table over subsets bottom-up exactly as the old code does its final pass.

IMPORTANT: copy the existing logic verbatim (same `w`, `fp`, CSE accounting, tie-break) so results match bit-for-bit; the equivalence test guards this.

- [ ] **Step 4: Reroute the additive arms of `single_term_opt<Metric>`**

In `single_term.hpp` detail `single_term_opt<Metric>` (~931-996), replace the `DenseFLOPs` and `DenseSize` arm bodies so each builds an `AdditiveModel` (with `flops_counter`/`memsize_counter`, the `volatile_mask` computed as today via `leaf_volatile_mask` or the existing inline mask, `volatile_weight`, `footprint_counter`, `footprint_weight`, `subnet_cse`) and returns `run_single_term_opt(model, network, tidxs)`. Include `cost_model.hpp`. Leave the `DensePeakSize`/`DensePeakSizeBatched` arms untouched (Tasks 2-3).

- [ ] **Step 5: Run to verify the equivalence test AND the full suite pass**

Run build/run. Expected: the new SECTION passes AND all existing `[optimize]` assertions still pass (the public-API DenseFLOPs/DenseSize tests now exercise `AdditiveModel`).

- [ ] **Step 6: clang-format + commit**
```bash
clang-format --style=file -i core/optimize/cost_model.hpp core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git add core/optimize/cost_model.hpp core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git commit -m "optimize: generic run_single_term_opt driver + AdditiveModel (FLOPs/Size)"
```

---

### Task 2: `PeakModel`

**Files:** `core/optimize/cost_model.hpp`, `core/optimize/single_term.hpp` (reroute `DensePeakSize` arm), test.

**Interfaces:**
- Consumes: the driver + concept (Task 1).
- Produces: `struct PeakModel` (nested `State`/`Context`); detail `single_term_opt<DensePeakSize>` returns `run_single_term_opt(PeakModel{...}, ...)`.

- [ ] **Step 1: Failing equivalence test**
```cpp
SECTION("PeakModel via driver == single_term_opt_peak_impl") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  for (auto const& spec : std::vector<std::vector<std::wstring>>{
         {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"},
         {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;a3}", L"g{a3;i2}"}}) {
    std::vector<ExprPtr> ts;
    for (auto s : spec) ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
    TensorNetwork net{ts}; container::svector<Index> targets;
    auto old_seq = opt::detail::single_term_opt_peak_impl(net, targets, idxsz);
    auto new_seq = opt::detail::run_single_term_opt(
        opt::detail::PeakModel{idxsz}, net, targets);
    REQUIRE(new_seq == old_seq);
  }
}
```

- [ ] **Step 2: Run to verify it fails** (`PeakModel` undeclared).

- [ ] **Step 3: Implement `PeakModel`** in `cost_model.hpp`. `Context` holds `S` (`subset_footprints`) and `L` (per-subset leaf-sum, as `peak_dp` precomputes). `State` = `{double peak; size_t lp, rp; bool lp_first;}`. `build_context` asserts `!subnet_cse` (PeakModel takes no CSE; if a `subnet_cse` flag is threaded, assert it false). `leaf` = `{S[n],0,0,true}`. `init` = `{max,0,0,true}`. `relax` = the `peak_dp` bipartition body (single_term.hpp ~516-560): both/lp_first/rp_first/cand, update `acc` + `lp_first`. `reconstruct` = `single_term_opt_peak_impl`'s emission (~737-770). Relocate verbatim.

- [ ] **Step 4: Reroute the `DensePeakSize` arm** of detail `single_term_opt<Metric>` to `return run_single_term_opt(PeakModel{idxsz}, network, tidxs);` (keep the `SEQUANT_ASSERT(!subnet_cse)`).

- [ ] **Step 5: Run - equivalence SECTION + full suite green** (the public-API DensePeakSize tests + the reconstruction/oracle SECTIONs that go through the public path now exercise `PeakModel`; the direct `peak_cost`/oracle SECTIONs still test the kept standalone functions).

- [ ] **Step 6: clang-format + commit**
```bash
clang-format --style=file -i core/optimize/cost_model.hpp core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git add -A && git commit -m "optimize: PeakModel via generic driver (DensePeakSize)"
```

---

### Task 3: `PeakBatchedModel`

**Files:** `core/optimize/cost_model.hpp`, `core/optimize/single_term.hpp` (reroute `DensePeakSizeBatched` arm), test.

**Interfaces:**
- Consumes: driver + concept.
- Produces: `struct PeakBatchedModel` (nested `State = container::vector<BatchedRes>`, `Context`); detail `single_term_opt<DensePeakSizeBatched>` returns `run_single_term_opt(PeakBatchedModel{...}, ...)`.

- [ ] **Step 1: Failing equivalence test** (shared-aux AND two-distinct-aux nets; reuse the `F`-space setup from the existing batched SECTIONs):
```cpp
SECTION("PeakBatchedModel via driver == single_term_opt_peak_batched_impl") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix) { return ix.space().base_key()==L"F"; };
  std::size_t const batch = 1;
  for (auto const& spec : std::vector<std::vector<std::wstring>>{
         {L"g{a1;i1;F1}", L"g{a2;i1;F1}", L"g{a2;i2;F2}"},
         {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"}}) {
    std::vector<ExprPtr> ts;
    for (auto s : spec) ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
    TensorNetwork net{ts}; container::svector<Index> targets;
    auto old_seq = opt::detail::single_term_opt_peak_batched_impl(
        net, targets, idxsz, is_batchable, batch, {});
    auto new_seq = opt::detail::run_single_term_opt(
        opt::detail::PeakBatchedModel{idxsz, is_batchable, batch,
                                      /*is_volatile_leaf=*/{}},
        net, targets);
    REQUIRE(new_seq == old_seq);
  }
}
```
(`F` must be registered in the cloned context as the existing batched SECTIONs do.)

- [ ] **Step 2: Run to verify it fails** (`PeakBatchedModel` undeclared).

- [ ] **Step 3: Implement `PeakBatchedModel`** in `cost_model.hpp`. Members: `idxsz`, `is_batchable`, `batch`, `is_volatile_leaf`. `Context` holds `aux=batchable_index_list`, `m`, `tables=sliced_footprints`, per-`B` `L`, `open_aux=subset_open_aux`, `volatile_mask=leaf_volatile_mask`, `nB=1<<m`. `State` = `container::vector<BatchedRes>` (size `nB`). `build_context` asserts `!subnet_cse`. `leaf(ctx,n)`: per `B`, `tables[B & open_aux[n]][n]`. `init(ctx,n)`: per `B`, `max`. `relax`: the `peak_dp_batched` body (single_term.hpp ~632-720) over `B` and `A'`, updating `acc[B]` + `aprime`. `reconstruct`: the `single_term_opt_peak_batched_impl` back-pointer walk from `(root,0)` (~779-900). Relocate verbatim; the only change is reading children from the driver's `st` table instead of an internal one.

- [ ] **Step 4: Reroute the `DensePeakSizeBatched` arm** to `return run_single_term_opt(PeakBatchedModel{idxsz, is_batchable_index, batch_target_size, is_volatile_leaf}, network, tidxs);` (keep `SEQUANT_ASSERT(!subnet_cse)`).

- [ ] **Step 5: Run - equivalence SECTION + full suite green** (public-API + reconstruction batched SECTIONs now exercise `PeakBatchedModel`; the direct `peak_cost_batched`/oracle SECTIONs still test the kept standalone functions).

- [ ] **Step 6: clang-format + commit**
```bash
clang-format --style=file -i core/optimize/cost_model.hpp core/optimize/single_term.hpp tests/unit/test_optimize.cpp
git add -A && git commit -m "optimize: PeakBatchedModel via generic driver (DensePeakSizeBatched)"
```

---

### Task 4: Concept conformance + custom-model extension point

**Files:** `core/optimize/cost_model.hpp` (the `CostModel` concept definition, if expressed as a C++20 `concept`), test.

**Interfaces:**
- Consumes: the driver + the four models.
- Produces: a `template<class M> concept CostModel = ...;` (a compile-checkable concept) and a `static_assert` that the four built-ins satisfy it; a tiny custom model demonstrating `run_single_term_opt<Custom>` directly.

- [ ] **Step 1: Define the `CostModel` concept** in `cost_model.hpp`:
```cpp
template <class M>
concept CostModel = requires(M const& m, typename M::Context& ctx,
                             TensorNetwork const& net,
                             container::svector<Index> const& tidxs, size_t n,
                             container::vector<typename M::State>& st,
                             typename M::State& acc) {
  { m.build_context(net, tidxs) } -> std::same_as<typename M::Context>;
  { m.leaf(ctx, n) } -> std::same_as<typename M::State>;
  { m.init(ctx, n) } -> std::same_as<typename M::State>;
  m.relax(ctx, n, n, n, acc, acc, acc);
  m.finalize(ctx, n, st);
  { m.reconstruct(ctx, st) } -> std::same_as<EvalSequence>;
};
```
(Adjust the exact `requires` expressions to the real method signatures from Tasks 1-3; the point is a compile-time conformance gate.)

- [ ] **Step 2: Failing test - conformance + a custom model**
```cpp
SECTION("CostModel concept conformance + custom model") {
  using namespace sequant;
  auto idxsz = [](Index const& ix) { return ix.space().approximate_size(); };
  static_assert(opt::detail::CostModel<
                opt::detail::AdditiveModel<decltype(opt::detail::flops_counter(idxsz))>>);
  static_assert(opt::detail::CostModel<opt::detail::PeakModel<decltype(idxsz)>>);
  // a trivial custom model: gross-size additive with a doubled footprint
  // weight, driven directly (proves the public extension point).
  opt::detail::AdditiveModel custom{
      opt::detail::memsize_counter(idxsz), opt::detail::footprint_counter(idxsz),
      0u, 1.0, /*footprint_weight=*/2.0, /*subnet_cse=*/false};
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;i1}", L"g{a1;a2}", L"g{a2;i2}"})
    ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
  TensorNetwork net{ts}; container::svector<Index> targets;
  auto seq = opt::detail::run_single_term_opt(custom, net, targets);
  size_t leaves = 0, merges = 0;
  for (int t : seq) (t >= 0 ? leaves : merges)++;
  REQUIRE(leaves == 3u);
  REQUIRE(merges == 2u);
}
```
(If `AdditiveModel`/`PeakModel` are class templates, spell the template args to match Tasks 1-2; if non-template, drop the `<...>`.)

- [ ] **Step 3: Make the conformance compile** - adjust the `concept` and/or the models so the `static_assert`s hold; no behavior change.

- [ ] **Step 4: Run - the new SECTION + full suite green.**

- [ ] **Step 5: clang-format + commit**
```bash
clang-format --style=file -i core/optimize/cost_model.hpp tests/unit/test_optimize.cpp
git add -A && git commit -m "optimize: CostModel concept + custom-model extension-point test"
```

---

## Self-review notes

- **Spec coverage:** generic driver + concept (spec 2) -> Task 1; the four models (spec 3) -> Tasks 1-3 (`AdditiveModel` x2, `PeakModel`, `PeakBatchedModel`); CSE via `finalize` (spec 3) -> Task 1 `AdditiveModel::finalize`; dispatch transition (spec 4) -> the per-task arm rerouting; scope/no-evaluator/no-mpqc (spec 5) -> Global Constraints; testing (spec 6: regression net + concept + custom) -> existing suite green each task + Task 4. Spec 7 file split -> `cost_model.hpp` (Task 1). Deliberately deferred per Global Constraints: removal of the old standalone DP/cost functions (spec 4's "replaced by" is realized for *dispatch*; the functions remain as oracles).
- **Behavior preservation:** every model method is a verbatim relocation of an existing loop body (exact source line ranges named); the equivalence SECTIONs assert `model-via-driver == old-impl` per objective, and the full existing suite must stay green - no assertion changes.
- **Type consistency:** `State`/`Context` names and the `leaf/init/relax/finalize/reconstruct` signatures are used identically across Tasks 1-4; `run_single_term_opt`/`solve_single_term` names fixed in Task 1.
- **Risk to flag to reviewer:** (a) whether `AdditiveModel`/`PeakModel` end up class templates (over the cost-fn / idxsz type) - the concept `static_assert`s and the dispatch sites must spell matching template args; (b) `volatile_mask` is computed in the additive arm exactly as today (reuse the existing inline computation or `leaf_volatile_mask`) so the equivalence holds; (c) `reconstruct` for additive must reproduce the `lseq[0] < rseq[0]` tie-break to match sequences bit-for-bit.
