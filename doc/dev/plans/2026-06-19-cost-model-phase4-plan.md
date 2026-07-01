# Cost-model Phase 4: shared BatchPolicy (optimizer + eval)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Define the batchability policy once - a `BatchPolicy{is_batchable_index, batch_target_size (per-index), is_volatile_leaf}` - consumed by BOTH the optimizer (`OptimizeOptions` embeds it) and the runtime batched evaluator (a thin `make_evaluator` adapter), so the two can no longer drift; and generalize `batch_target_size` from a scalar to a per-index function.

**Architecture:** Stage A (SeQuant): generalize `batch_target_size` scalar->`std::function<size_t(Index)>` through the batched path + `make_batched_custom_evaluator`; introduce `BatchPolicy` (neutral `core/batch_policy.hpp`); embed it in `OptimizeOptions`; add the eval-layer `make_evaluator` adapter (lifts the Tensor volatile predicate to EvalNode). Stage B (mpqc): construct one `BatchPolicy`, feed both the optimizer and the eval cache, delete the duplicated definitions. Behavior-preserving throughout.

**Tech Stack:** C++20, SeQuant `core/optimize` + `core/eval`, Catch2 (`tests/unit/test_optimize.cpp`); mpqc `src/mpqc/chemistry/qc/lcao/{expression/sequant_engine.cpp,cc/cck.ipp}`.

## Global Constraints

- Style: Google `.clang-format`, 80-col, 2-space, no tabs. `clang-format --style=file -i <file>` before each commit (SeQuant). mpqc uses the same.
- No en-dashes (U+2013) / non-breaking spaces (SeQuant pre-commit hook rejects); ASCII hyphen only. No `Co-Authored-By`.
- Branch: `feature/cost-model-batch-aware` (continue; Phases 1-3 already here). mpqc work is on the corresponding mpqc4 working tree (SeQuant is consumed live via FetchContent override, so Stage-A changes are visible to mpqc builds).
- **Behavior-preserving (hard):** SeQuant - no existing `[optimize]` assertion changes (Stage A). mpqc - a CSV-CCk reference run yields the SAME energy (Stage B).
- **Known internal window:** Stage A changes `make_batched_custom_evaluator`'s `target_batch_size` parameter (scalar->function), which `cck.ipp` calls directly; mpqc will not compile between Stage A and Stage B. This is closed by Stage B on the same branch; do NOT build mpqc between A and B (Stage A gates are SeQuant unit tests only).
- Reference: spec `doc/dev/specs/2026-06-19-cost-model-phase4-batch-policy-design.md`.
- SeQuant source is under `SeQuant/core/...`; tests under `tests/unit/`. Build/run SeQuant tests:
  ```bash
  cd /Users/efv/code/SeQuant && cmake --build build-test --target unit_tests-sequant -j
  ./build-test/tests/unit/unit_tests-sequant "[optimize]"
  ```

## File structure

- `SeQuant/core/batch_policy.hpp` (new) - the `BatchPolicy` struct.
- `SeQuant/core/optimize/options.hpp` - `batch_target_size` field type change (A1), then embed `BatchPolicy` (A2).
- `SeQuant/core/optimize/single_term.hpp` + `cost_model.hpp` - thread the per-index `batch_target_size` (A1); read `opts.batch_policy` (A2).
- `SeQuant/core/eval/eval.hpp` - `make_batched_custom_evaluator` param scalar->function (A1); add `make_evaluator` adapter (A3).
- `tests/unit/test_optimize.cpp` - update batched tests to the function form (A1); add adapter + per-index tests (A1/A3).
- mpqc `sequant_engine.cpp` + `cck.ipp` - construct + feed one policy, delete dup (B1).

---

### Task A1: Generalize `batch_target_size` scalar -> per-index function

**Files:** `SeQuant/core/optimize/{options.hpp,single_term.hpp,cost_model.hpp}`, `SeQuant/core/eval/eval.hpp`, `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Produces: every batched-path API that took `std::size_t batch` now takes `std::function<std::size_t(Index const&)> batch_target_size`; slicing applies `min(extent(ix), batch_target_size(ix))`. `make_batched_custom_evaluator`'s `target_batch_size` becomes the same function type, used as `target_batch_size(*K)` / `target_batch_size(*Kk)` at its two `mode_batches` sites.

- [ ] **Step 1: Write a failing per-index test**

Add a SECTION to `TEST_CASE("optimize")` (sized-context block, `F` aux registered as the existing batched SECTIONs do) asserting per-index slice sizes are honored: build a two-distinct-aux net, give `F1` batch size 1 and `F2` batch size 2 via a lambda, and assert `peak_cost_batched` differs from the all-size-1 result (proving the function is consulted per index, not a single scalar):
```cpp
SECTION("per-index batch_target_size honored") {
  using namespace sequant;
  auto idxsz = [](Index const& ix){ return ix.space().approximate_size(); };
  auto is_batchable = [](Index const& ix){ return ix.space().base_key()==L"F"; };
  std::vector<ExprPtr> ts;
  for (auto s : {L"g{a1;i1;F1}", L"g{a2;i1;F2}", L"g{a2;i2;F2}"})
    ts.push_back(deserialize(s, {.def_perm_symm=Symmetry::Nonsymm}));
  TensorNetwork net{ts}; container::svector<Index> targets;
  auto all1   = [](Index const&) -> std::size_t { return 1; };
  auto mixed  = [](Index const& ix) -> std::size_t {
    return ix.label().find(L"F1") != std::wstring::npos ? 1 : 2; };
  double c_all1  = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                  is_batchable, all1, {});
  double c_mixed = opt::detail::peak_cost_batched(net, targets, idxsz,
                                                  is_batchable, mixed, {});
  REQUIRE(c_mixed >= c_all1);          // larger F2 slice -> no smaller peak
  REQUIRE(c_mixed != c_all1);          // the per-index size is actually used
}
```
(Adjust the `mixed` discriminator to whatever reliably distinguishes `F1` from `F2` in this Index API; the point is two different sizes by space/label.)

- [ ] **Step 2: Run to verify it fails** (compile error: `peak_cost_batched` takes `std::size_t`, not a lambda).

- [ ] **Step 3: Change the batched-path signatures to the function**

In `single_term.hpp`: change the `batch` parameter from `std::size_t` to `std::function<std::size_t(Index const&)> batch_target_size` in `sliced_footprints`, `peak_dp_batched`, `peak_cost_batched`, `single_term_opt_peak_batched_impl`, and the `batched_extent` wrapper. In `batched_extent`, slice as `is_batchable(ix) ? std::min(idxsz(ix), batch_target_size(ix)) : idxsz(ix)`. Thread the function unchanged through the other functions (they only forward it into `sliced_footprints`/`batched_extent`).

In `cost_model.hpp`: change `PeakBatchedModel`'s `batch` member type to the function; `build_context` forwards it into `sliced_footprints`.

- [ ] **Step 4: Change `make_batched_custom_evaluator`**

In `eval.hpp`, change `make_batched_custom_evaluator`'s `target_batch_size` parameter from `std::size_t` to `std::function<std::size_t(Index const&)>`. At the trigger `mode_batches` (the `auto const batches = le(leaf->first)->mode_batches(leaf->second, target_batch_size)` line), pass `target_batch_size(*K)`. At the replay-group candidate probe (`... mode_batches(lk->second, target_batch_size) != batches`), pass `target_batch_size(*Kk)` (the candidate's own axis), so group-compatibility still compares like with like.

- [ ] **Step 5: Update the OptimizeOptions field + dispatch arm + existing tests**

In `options.hpp`: change `OptimizeOptions::batch_target_size` from `std::size_t` to `std::function<std::size_t(Index const&)>` (default `{}`). In `single_term.hpp` the `DensePeakSizeBatched` arm passes `opts.batch_target_size` (now the function) into the model - no other change. Update every existing batched test that passed `batch=1` (a `std::size_t`) to a constant lambda `[](Index const&){ return std::size_t{1}; }` (search `test_optimize.cpp` for the batched SECTIONs from Phase 2/3). Update any SeQuant eval unit test that calls `make_batched_custom_evaluator` with a scalar likewise.

- [ ] **Step 6: Run - new SECTION + full `[optimize]` suite green**

The constant-lambda substitution reproduces the prior scalar results (behavior-preserving); the new per-index SECTION passes.

- [ ] **Step 7: clang-format + commit**
```bash
clang-format --style=file -i SeQuant/core/optimize/single_term.hpp SeQuant/core/optimize/cost_model.hpp SeQuant/core/optimize/options.hpp SeQuant/core/eval/eval.hpp tests/unit/test_optimize.cpp
git add -A && git commit -m "optimize/eval: batch_target_size scalar -> per-index function"
```

---

### Task A2: `BatchPolicy` struct + embed in `OptimizeOptions`

**Files:** `SeQuant/core/batch_policy.hpp` (new), `SeQuant/core/optimize/options.hpp`, `SeQuant/core/optimize/{single_term.hpp,cost_model.hpp}`, `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Produces: `sequant::BatchPolicy{ std::function<bool(Index const&)> is_batchable_index; std::function<std::size_t(Index const&)> batch_target_size; std::function<bool(Tensor const&)> is_volatile_leaf; }`. `OptimizeOptions` gains `BatchPolicy batch_policy;` and DROPS the loose `is_batchable_index`, `batch_target_size`, `is_volatile_leaf` fields (keeps `volatile_weight`). Models/arms read `opts.batch_policy.*`.

- [ ] **Step 1: Create `batch_policy.hpp`**
```cpp
#ifndef SEQUANT_CORE_BATCH_POLICY_HPP
#define SEQUANT_CORE_BATCH_POLICY_HPP
#include <cstddef>
#include <functional>
namespace sequant {
class Index;
class Tensor;
/// One batchability policy shared by the single-term optimizer and the runtime
/// batched evaluator (see make_evaluator). All three predicates default empty.
struct BatchPolicy {
  std::function<bool(Index const&)> is_batchable_index = {};
  std::function<std::size_t(Index const&)> batch_target_size = {};
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
};
}  // namespace sequant
#endif
```

- [ ] **Step 2: Failing test - models read from `batch_policy`**

Update the existing batched public-API SECTION (the one that builds `OptimizeOptions` with `objective_function = DensePeakSizeBatched`) to set `opts.batch_policy = {is_batchable, batch_lambda, {}}` instead of the three loose fields, and assert it still returns a valid 3-leaf product. This SECTION must fail to compile first (the loose fields removed in Step 3 don't exist), then pass.

- [ ] **Step 3: Embed `BatchPolicy` in `OptimizeOptions`**

`options.hpp`: `#include <SeQuant/core/batch_policy.hpp>`; remove the three loose fields; add `BatchPolicy batch_policy = {};`. Keep `volatile_weight`.

- [ ] **Step 4: Rewire models + dispatch arms**

In `single_term.hpp` `single_term_opt<Metric>`: the `DensePeakSizeBatched` arm builds `PeakBatchedModel{idxsz, opts... }` reading `opts.batch_policy.is_batchable_index` / `opts.batch_policy.batch_target_size` / `opts.batch_policy.is_volatile_leaf`; the `DenseFLOPs`/`DenseSize` arms read `opts.batch_policy.is_volatile_leaf` for the volatile mask (where they previously read `opts.is_volatile_leaf`); `volatile_weight` unchanged. (Whatever signatures `single_term_opt<Metric>` uses to receive these - parameters or `opts` - reroute the reads to `opts.batch_policy`.) `opt_pure_product` (optimize.cpp) passes `opts` through unchanged.

- [ ] **Step 5: Update remaining tests**

Any test or call site that set the loose `OptimizeOptions::{is_batchable_index,batch_target_size,is_volatile_leaf}` now sets `opts.batch_policy.{...}`. The direct `opt::detail::*` calls (which take the predicates as explicit args, not via OptimizeOptions) are unaffected.

- [ ] **Step 6: Run - full `[optimize]` suite green.**

- [ ] **Step 7: clang-format + commit**
```bash
clang-format --style=file -i SeQuant/core/batch_policy.hpp SeQuant/core/optimize/options.hpp SeQuant/core/optimize/single_term.hpp SeQuant/core/optimize/cost_model.hpp tests/unit/test_optimize.cpp
git add -A && git commit -m "optimize: introduce BatchPolicy, embed in OptimizeOptions"
```

---

### Task A3: eval-layer `make_evaluator` adapter

**Files:** `SeQuant/core/eval/eval.hpp`, `tests/unit/` (an eval test file - use the existing eval test TU; if `[optimize]` is the only convenient place, add it there with the eval includes).

**Interfaces:**
- Produces: `template<class F, class ScopeGuardFactory> auto sequant::make_evaluator(BatchPolicy const& policy, F yielder, ScopeGuardFactory make_scope_guard = {})` returning the same custom-evaluator callable `make_batched_custom_evaluator` produces, built from `policy`.

- [ ] **Step 1: Failing test - adapter == hand-built evaluator**

Add a test (in the eval test TU that already exercises `make_batched_custom_evaluator`; mirror its fixture) that builds two evaluators on the same small batched network - one via `make_evaluator(policy, yielder, guard)` and one via `make_batched_custom_evaluator(yielder, const_size_fn, policy.is_batchable_index, guard, is_volatile_node)` (the hand lift) - drives both through the cache on the same node, and asserts identical results (same `ResultPtr` value / tensor). If the existing eval tests don't have a ready batched fixture, assert the weaker but real property: `make_evaluator(policy, yielder)` returns a non-null callable that, set on the cache, evaluates a known batched node to the same value as the non-batched path.

- [ ] **Step 2: Run to verify it fails** (`make_evaluator` undeclared).

- [ ] **Step 3: Implement `make_evaluator`** in `eval.hpp`, after `make_batched_custom_evaluator`:
```cpp
template <class F, class ScopeGuardFactory = /* the same default as
          make_batched_custom_evaluator */>
[[nodiscard]] auto make_evaluator(BatchPolicy const& policy, F yielder,
                                  ScopeGuardFactory make_scope_guard = {}) {
  auto is_volatile_node = [p = policy.is_volatile_leaf](auto const& n) {
    return n.leaf() && n->is_tensor() && p(n->as_tensor());
  };
  return make_batched_custom_evaluator(
      std::move(yielder), policy.batch_target_size, policy.is_batchable_index,
      std::move(make_scope_guard), std::move(is_volatile_node));
}
```
`#include <SeQuant/core/batch_policy.hpp>` in `eval.hpp`. Use the same `ScopeGuardFactory` default type `make_batched_custom_evaluator` uses (so the no-arg form matches).

- [ ] **Step 4: Run - the new test + full suite green.**

- [ ] **Step 5: clang-format + commit**
```bash
clang-format --style=file -i SeQuant/core/eval/eval.hpp tests/unit/<eval-test>.cpp
git add -A && git commit -m "eval: make_evaluator(BatchPolicy) adapter over make_batched_custom_evaluator"
```

---

### Task B1: mpqc - construct one BatchPolicy, feed both, delete dup

**Files:** mpqc `src/mpqc/chemistry/qc/lcao/cc/cck.ipp`, `src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp` (and `sequant_engine.h` if the policy is threaded via a member/context).

**Interfaces:**
- Consumes: `sequant::BatchPolicy`, `sequant::make_evaluator` (Stage A). `OptimizeOptions::batch_policy`.

- [ ] **Step 1: Construct one `BatchPolicy` in the CSV-CCk setup**

In `cck.ipp` where `csv_batch_aux_target_size_ > 0` is handled, build:
```cpp
auto const aux_space = to_sequant_space(SPIndex::Type::dfbs);
sequant::BatchPolicy policy;
policy.is_batchable_index = [aux_space](sequant::Index const& ix) {
  return ix.space() == aux_space;
};
policy.batch_target_size =
    [n = csv_batch_aux_target_size_](sequant::Index const&) { return n; };
policy.is_volatile_leaf = [](sequant::Tensor const& t) {
  return t.label() == L"t";
};
```

- [ ] **Step 2: Feed the eval cache via `make_evaluator`**

Replace the `cache.set_custom_evaluator(make_batched_custom_evaluator(yielder, csv_batch_aux_target_size_, accept_aux, make_scope_guard, is_volatile))` call with `cache.set_custom_evaluator(sequant::make_evaluator(policy, yielder, make_scope_guard))`. Delete the now-dead `accept_aux` and the EvalNode-based `is_volatile` lambdas.

- [ ] **Step 3: Feed the optimizer via `batch_policy`**

Thread `policy` to where `SeQuantEngine::make_optimize_options` builds `OptimizeOptions`, and set `opts.batch_policy = policy`. Use the construction site that reaches both this and Step 2 with one `policy` (the `EvalContext` `SeQuantEngine::optimize` takes, or a `CCk` member the engine reads). Remove the separate `is_batchable_index`/`batch_target_size`/`is_volatile_leaf` assignments in `sequant_engine.cpp` (now from `batch_policy`).

- [ ] **Step 4: Build mpqc + run a CSV-CCk reference**

```bash
cd /Users/efv/code/mpqc4 && cmake --build cmake-build-release --target mpqc -j
```
Run a CSV-CCk reference calculation (e.g. the `he10` CSV-CCk batched validation input used in this work, or `tests/validation/reference/inputs/*csv-cck*batched*`), confirming the correlation energy matches the value from before this change (the pre-Phase-4 baseline). Use a small input so the run is quick; serialize (one mpqc process at a time).

- [ ] **Step 5: Commit (mpqc repo)**
```bash
cd /Users/efv/code/mpqc4
clang-format --style=file -i src/mpqc/chemistry/qc/lcao/cc/cck.ipp src/mpqc/chemistry/qc/lcao/expression/sequant_engine.cpp
git add -A && git commit -m "cck/SeQuantEngine: single BatchPolicy feeds optimizer + eval cache"
```
(Follow mpqc commit conventions: no Co-Authored-By trailers.)

---

## Self-review notes

- **Spec coverage:** BatchPolicy struct (spec 2) -> A2; per-index batch_target_size (spec 1,2) -> A1; eval adapter + volatile lift (spec 3) -> A3; OptimizeOptions embed + model rewire (spec 4) -> A2; mpqc one-policy feed-both + delete-dup (spec 5) -> B1; two-stage scope (spec 6) -> A* then B1; testing (spec 7) -> A1 per-index test + A3 adapter test + suite-green gates + B1 CSV-CCk validation. The spec-8 `mode_batches` two-site caveat is handled explicitly in A1 Step 4.
- **Behavior preservation:** A1 is a type generalization proven by the constant-lambda substitution keeping the full suite green; A2/A3 are additive/rewire proven by the suite; B1 by the CSV-CCk energy match.
- **Type consistency:** `batch_target_size : std::function<std::size_t(Index const&)>` and `is_batchable_index : std::function<bool(Index const&)>` and `is_volatile_leaf : std::function<bool(Tensor const&)>` used identically across A1/A2/A3/B1; `BatchPolicy`/`make_evaluator` names fixed in A2/A3.
- **Risks to flag to reviewer:** (a) the A->B mpqc-compile window (stated in Global Constraints); (b) the `make_batched_custom_evaluator` second `mode_batches` site (replay-group probe) must use the candidate axis `*Kk`, not the trigger axis, or group compatibility breaks; (c) the eval test fixture availability (A3 falls back to a weaker non-null+same-value assertion if no batched eval fixture exists); (d) the mpqc threading point (EvalContext vs member) - pick one construction site reaching both consumers.
