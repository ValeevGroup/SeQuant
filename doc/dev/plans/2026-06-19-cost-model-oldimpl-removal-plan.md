# CostModel old-impl removal (full DRY) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the now-dead standalone single-term-DP driver/table functions left behind by the Phase 3 CostModel refactor, routing the kept cost-only oracle helpers through the models (full DRY), so each DP recurrence exists in exactly one place.

**Architecture:** Phase 3 made every production single-term-optimization arm route through `run_single_term_opt(Model)`, leaving three reconstruction-driver functions (`single_term_opt_impl`, `single_term_opt_peak_impl`, `single_term_opt_peak_batched_impl`) and two DP table-builders (`peak_dp`, `peak_dp_batched`) referenced only by tests and by the cost-only oracle helpers (`peak_cost`, `peak_cost_batched`, `reconstructed_batched_peak`). The chosen scope (option B) is full DRY: rewire the three oracle helpers to delegate to the models (`PeakModel`/`PeakBatchedModel` via `solve_single_term`), which requires moving them from `single_term.hpp` into `cost_model.hpp` (the models live there, and `cost_model.hpp` already `#include`s `single_term.hpp`); then delete all five dead functions plus the `PeakRes` struct. The independent regression net becomes the brute-force oracles (`brute_force_min_peak`, `batched_min_peak`) plus the reconstruction-simulation checks — the model-vs-old-impl equivalence SECTIONs (refactor scaffolding) are deleted.

**Tech Stack:** C++20, SeQuant `sequant::opt::detail` namespace, Catch2 v3 unit tests (`tests/unit/test_optimize.cpp`), CMake+Ninja build.

## Global Constraints

- Behavior-preserving: the full `[optimize]` Catch2 suite stays green with no test-number changes (the kept oracle helpers reproduce their prior values bit-for-bit, because the models are the Phase-3 behavior-preserving extractions of the deleted DPs). — copied from the Phase-3 spec hard constraint.
- No `Co-Authored-By` / tooling-attribution trailers in commits (matches SeQuant + mpqc `git log`).
- Keep the independent oracles: `brute_force_min_peak` and `batched_min_peak` (the brute-force simulations in `test_optimize.cpp`) are NOT touched — they remain the independent cross-checks.
- All work on branch `feature/cost-model-batch-aware` (the Phases 1-4 branch); do not create a new branch.
- `peak_cost` / `peak_cost_batched` / `reconstructed_batched_peak` keep their exact current signatures (only their bodies and their home header change), so test call sites need no signature edits.

---

### Task 1: Rewire the cost-only oracle helpers to delegate to the models

Move `peak_cost`, `peak_cost_batched`, and `reconstructed_batched_peak` out of `single_term.hpp` into `cost_model.hpp` (placed AFTER the model definitions, before the `CostModel` concept / `#endif`), rewiring each body to build the matching model, run `solve_single_term`, and read the root cell — instead of calling `peak_dp` / `peak_dp_batched`. This makes the models the single source of each recurrence while preserving the helpers' signatures and values.

**Files:**
- Modify: `SeQuant/core/optimize/single_term.hpp` — delete the bodies of `peak_cost` (lines ~558-569), `peak_cost_batched` (lines ~712-741), `reconstructed_batched_peak` (lines ~827-end of that function). Leave `peak_dp` / `peak_dp_batched` / the `*_impl` functions in place for now (Task 3 deletes them); they are still referenced by the not-yet-rewired tests.
- Modify: `SeQuant/core/optimize/cost_model.hpp` — add the three rewired free functions after `PeakBatchedModel` (so they can name the models) and before `concept CostModel`.
- Test: `tests/unit/test_optimize.cpp` (no edits this task — call sites keep working because the signatures are unchanged and the test TU already includes `cost_model.hpp`).

**Interfaces:**
- Consumes (from `cost_model.hpp`, already present): `solve_single_term(Model const&, TensorNetwork const&, TIdxs const&, typename Model::Context&) -> container::vector<Model::State>`; `PeakModel{IdxToSz}` with `State{double peak; size_t lp,rp; bool lp_first;}` and `Context build_context(net,tidxs)`; `PeakBatchedModel{IdxToSz, is_batchable, batch, is_volatile_leaf}` with `using State = container::vector<BatchedRes>` and `Context` exposing `double sz(size_t s, size_t ctx) const` and `double Lof(size_t s, size_t ctx) const`, plus members `aux`, `m`, `nB`, `nt`, `open_aux`.
- Consumes (from `single_term.hpp`, stay): `BatchedRes{double peak; size_t lp,rp; bool lp_first; size_t aprime;}`, `bits::bipartitions`.
- Produces (signatures UNCHANGED, new home `cost_model.hpp`): `double peak_cost(TensorNetwork const&, TIdxs const&, IdxToSz&&)`; `double peak_cost_batched(TensorNetwork const&, TIdxs const&, IdxToSz&&, std::function<bool(Index const&)> const&, std::function<std::size_t(Index const&)> const&, std::function<bool(Tensor const&)> const&)`; `double reconstructed_batched_peak(/* same 6 params as peak_cost_batched */)`.

- [ ] **Step 1: Verify the suite is green before touching anything (baseline)**

Run (build dir already configured at `/Users/efv/code/SeQuant/build`; if absent, configure with `cmake -S /Users/efv/code/SeQuant -B /Users/efv/code/SeQuant/build -G Ninja`):
```bash
cmake --build /Users/efv/code/SeQuant/build --target unit_tests-sequant 2>&1 | tail -5
/Users/efv/code/SeQuant/build/tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -8
```
Expected: build succeeds; all `[optimize]` assertions pass (`All tests passed` line lists the SECTION count).

- [ ] **Step 2: Add the rewired `peak_cost` to `cost_model.hpp`**

Insert after `PeakBatchedModel`'s closing `};` (currently line ~450) and before `concept CostModel`:
```cpp
/// \brief Achieved minimum peak memory (the DensePeakSize objective value) for
/// the whole network under its optimal order. Builds \ref PeakModel, runs the
/// generic driver's \ref solve_single_term, and returns the root subset's peak.
/// Used by tests to compare against the brute-force oracle.
template <typename TIdxs, typename IdxToSz>
double peak_cost(TensorNetwork const& network, TIdxs const& tidxs,
                 IdxToSz&& idxsz) {
  PeakModel<std::decay_t<IdxToSz>> model{std::forward<IdxToSz>(idxsz)};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  return st.back().peak;  // root subset == full set == last element
}
```

- [ ] **Step 3: Add the rewired `peak_cost_batched` to `cost_model.hpp`**

Immediately below `peak_cost`:
```cpp
/// \brief Achieved minimum batched peak memory: the peak[root][B=0] objective
/// of the multi-mode batched DP. Builds \ref PeakBatchedModel, runs
/// \ref solve_single_term, and returns the root subset's B=0 peak.
template <typename TIdxs, typename IdxToSz>
double peak_cost_batched(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{
      std::forward<IdxToSz>(idxsz), is_batchable, batch_target_size,
      is_volatile_leaf};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  return st.back()[0].peak;  // root subset's B=0 cell
}
```

- [ ] **Step 4: Add the rewired `reconstructed_batched_peak` to `cost_model.hpp`**

Immediately below `peak_cost_batched`. Build the model, solve, then preserve the existing simulation walk VERBATIM from the current `reconstructed_batched_peak` body in `single_term.hpp` (the `build` lambda that walks back-pointers and the live-set memory simulation, lines ~868 to the function's closing `}`), swapping only the data sources: where the old body reads `at(n, B)` (its local `pr[n*nB+B]`) read `st[n][B]`; where it calls its local `sz`/`Lof` lambdas use `ctx.sz`/`ctx.Lof`; delete the old local setup (`vmask`, `pr`, `aux`, `nB`, `tables`, `open_aux`, the `at`/`sz`/`Lof` lambdas) since `ctx` and `st` now supply them. The function signature is unchanged.
```cpp
/// \brief Independent memory-simulation recomputation of the chosen batched
/// reconstruction's model-A peak. Builds \ref PeakBatchedModel, runs
/// \ref solve_single_term for the back-pointer table, and recomputes the
/// subtree peak by direct memory simulation (NOT the DP's max/+ formula). Must
/// EQUAL \ref peak_cost_batched; a mismatch signals a DP/reconstruction bug.
template <typename TIdxs, typename IdxToSz>
double reconstructed_batched_peak(
    TensorNetwork const& network, TIdxs const& tidxs, IdxToSz&& idxsz,
    std::function<bool(Index const&)> const& is_batchable,
    std::function<std::size_t(Index const&)> const& batch_target_size,
    std::function<bool(Tensor const&)> const& is_volatile_leaf) {
  PeakBatchedModel<std::decay_t<IdxToSz>> model{
      std::forward<IdxToSz>(idxsz), is_batchable, batch_target_size,
      is_volatile_leaf};
  auto ctx = model.build_context(network, tidxs);
  auto st = solve_single_term(model, network, tidxs, ctx);
  auto const nt = network.tensors().size();
  // <-- preserve the existing simulation walk verbatim here, reading st[n][B]
  //     for back-pointers and ctx.sz / ctx.Lof for sizes (see step text). The
  //     walk starts at (root = (1<<nt)-1, B=0) and returns the simulated peak.
}
```
Note: `ctx.sz`/`ctx.Lof` are `const` methods (Task 0 of Phase 4 added them); `ctx` may be `const auto&` here. If the existing simulation used `m` (number of aux) or `aux`, read `ctx.m` / `ctx.aux`.

- [ ] **Step 5: Delete the three helper bodies from `single_term.hpp`**

Remove the entire definitions of `peak_cost` (the `template ... double peak_cost(...) { ... }` block, ~lines 558-569), `peak_cost_batched` (~lines 712-741), and `reconstructed_batched_peak` (~lines 827 through its closing brace). Leave `peak_dp`, `peak_dp_batched`, `single_term_opt_impl`, `single_term_opt_peak_impl`, `single_term_opt_peak_batched_impl` in place (Task 3 removes them). Do NOT remove `subset_open_aux`, `BatchedRes`, `subset_footprints`, `sliced_footprints` — they are reused by the models.

- [ ] **Step 6: Build and run the suite — values must be unchanged**

Run:
```bash
clang-format --style=file -i /Users/efv/code/SeQuant/SeQuant/core/optimize/cost_model.hpp /Users/efv/code/SeQuant/SeQuant/core/optimize/single_term.hpp
cmake --build /Users/efv/code/SeQuant/build --target unit_tests-sequant 2>&1 | tail -5
/Users/efv/code/SeQuant/build/tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -8
```
Expected: clean build; `[optimize]` all pass (same SECTION/assertion counts as Step 1). The peak/batched oracle SECTIONs (`DensePeakSize DP matches brute-force oracle`, `DensePeakSizeBatched objective matches per-index oracle`, `DensePeakSizeBatched reconstruction achieves the optimum (numeric)`, etc.) now run through the models and still match their independent oracles.

- [ ] **Step 7: Commit**

```bash
cd /Users/efv/code/SeQuant
git add SeQuant/core/optimize/cost_model.hpp SeQuant/core/optimize/single_term.hpp
git commit -m "optimize: route peak_cost/peak_cost_batched/reconstructed_batched_peak through the models"
```

---

### Task 2: Rewire the tests off the soon-to-be-removed functions

Delete the three model-vs-old-impl equivalence SECTIONs (pure refactor scaffolding) and rewire the two remaining tests that use a soon-removed function as their production source, so that after this task NO test references `single_term_opt_impl`, `single_term_opt_peak_impl`, `single_term_opt_peak_batched_impl`, `peak_dp`, or `peak_dp_batched`.

**Files:**
- Modify: `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Consumes: `opt::detail::run_single_term_opt(Model const&, net, targets) -> EvalSequence`; `opt::detail::PeakModel{idxsz}`; `opt::detail::PeakBatchedModel{idxsz, is_batchable, batch_fn, is_volatile_leaf}`; `opt::detail::solve_single_term(model, net, targets, ctx)`; `PeakBatchedModel::Context` (from `model.build_context(net, targets)`) and `State = container::vector<BatchedRes>` (so `st[n][B].peak`); `opt::detail::peak_cost` (now from `cost_model.hpp`, signature unchanged).
- Produces: a test file with zero references to the five doomed symbols.

- [ ] **Step 1: Delete the three equivalence SECTIONs**

Remove these three SECTIONs in full (they assert `model-via-driver == old standalone impl`, which becomes meaningless once the old impls are deleted; the models' correctness is covered by the brute-force-oracle and reconstruction SECTIONs):
- `SECTION("AdditiveModel via driver == single_term_opt_impl")` (~lines 541-583).
- `SECTION("PeakModel via driver == single_term_opt_peak_impl")` (~lines 905-924).
- `SECTION("PeakBatchedModel via driver == single_term_opt_peak_batched_impl")` (~lines 926-956).

- [ ] **Step 2: Rewire the `single_term_opt_peak_impl` use in the reconstruction SECTION**

In `SECTION("DensePeakSize reconstructed sequence achieves the DP optimum")`, replace the line:
```cpp
        auto seq = opt::detail::single_term_opt_peak_impl(net, targets, idxsz);
```
with:
```cpp
        auto seq = opt::detail::run_single_term_opt(
            opt::detail::PeakModel{idxsz}, net, targets);
```
(`peak_of_sequence(seq, S, ts.size()) == peak_cost(...)` still holds: the model's reconstruct is the Phase-3 extraction of `single_term_opt_peak_impl`.)

- [ ] **Step 3: Rewire the direct `peak_dp_batched` use in the all-sliced-corner SECTION**

In `SECTION("DensePeakSizeBatched all-sliced corner equals Phase-1 batched peak")`, replace the block:
```cpp
      auto pr = opt::detail::peak_dp_batched(
          net, targets, idxsz, is_batchable, batch_fn,
          opt::detail::leaf_volatile_mask(net, {}));
      size_t root = (size_t{1} << ts.size()) - 1;
      size_t allK = (size_t{1} << m) - 1;
      double dp_allsliced = pr[root * (size_t{1} << m) + allK].peak;
```
with a solve through the model (its `State` is the per-subset `[B]`-vector of `BatchedRes`, so `st[root][allK].peak` is the same cell `pr[root*nB+allK].peak` was):
```cpp
      opt::detail::PeakBatchedModel model{idxsz, is_batchable, batch_fn,
                                          /*is_volatile_leaf=*/{}};
      auto ctx = model.build_context(net, targets);
      auto st = opt::detail::solve_single_term(model, net, targets, ctx);
      size_t root = (size_t{1} << ts.size()) - 1;
      size_t allK = (size_t{1} << m) - 1;
      double dp_allsliced = st[root][allK].peak;
```
(The downstream `phase1 = opt::detail::peak_cost(net, targets, be)` and `REQUIRE(dp_allsliced == phase1)` are unchanged.)

- [ ] **Step 4: Confirm no test still references the doomed symbols**

Run:
```bash
grep -nE '\b(single_term_opt_impl|single_term_opt_peak_impl|single_term_opt_peak_batched_impl|peak_dp|peak_dp_batched)\b' /Users/efv/code/SeQuant/tests/unit/test_optimize.cpp
```
Expected: NO output (empty). If any line prints, rewire it as above before continuing.

- [ ] **Step 5: Build and run the suite**

Run:
```bash
clang-format --style=file -i /Users/efv/code/SeQuant/tests/unit/test_optimize.cpp
cmake --build /Users/efv/code/SeQuant/build --target unit_tests-sequant 2>&1 | tail -5
/Users/efv/code/SeQuant/build/tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -8
```
Expected: clean build; all `[optimize]` pass (three fewer SECTIONs than Task 1 Step 1; all kept SECTIONs still pass).

- [ ] **Step 6: Commit**

```bash
cd /Users/efv/code/SeQuant
git add tests/unit/test_optimize.cpp
git commit -m "optimize tests: drop model-vs-old-impl equivalence SECTIONs, route remaining uses through the models"
```

---

### Task 3: Delete the dead standalone functions and the `PeakRes` struct

With no production or test caller remaining, delete the five dead functions and the now-unused `PeakRes` struct from `single_term.hpp`, and fix any dangling Doxygen `\ref` to them.

**Files:**
- Modify: `SeQuant/core/optimize/single_term.hpp`.

**Interfaces:**
- Consumes: nothing new.
- Produces: a `single_term.hpp` with the five functions and `PeakRes` removed; the surviving helpers (`subset_footprints`, `sliced_footprints`, `subset_open_aux`, `init_results`, `build_subnet_metadata`, the counters, `batchable_index_list`, `leaf_volatile_mask`, `BatchedRes`, `bits::bipartitions`) untouched.

- [ ] **Step 1: Confirm no production caller remains**

Run:
```bash
grep -rnE '\b(single_term_opt_impl|single_term_opt_peak_impl|single_term_opt_peak_batched_impl|peak_dp|peak_dp_batched)\b' /Users/efv/code/SeQuant/SeQuant | grep -vE '^\S+:[0-9]+:\s*(///|//|\*)'
```
Expected: only the definition lines in `single_term.hpp` (and possibly a `\ref` inside a comment, which is fine — Step 3 fixes comments). No call sites in other `.hpp`/`.cpp`. If a real (non-comment) call site appears outside the definitions, STOP and report — the prior tasks missed a caller.

- [ ] **Step 2: Delete the five function definitions and `PeakRes`**

Remove from `single_term.hpp`:
- `struct PeakRes { ... };` (~lines 500-506) and its preceding `/// Per-subset state for the peak DP.` doc comment.
- `peak_dp(...)` (the `template ... container::vector<PeakRes> peak_dp(...) { ... }` block, ~lines 508-556, including its doc comment).
- `peak_dp_batched(...)` (~lines 603-710, the doc comment + body; keep `struct BatchedRes` which precedes it and is reused by `PeakBatchedModel`).
- `single_term_opt_impl(...)` (the additive reconstruction driver, ~line 393 through its closing brace, including its doc comment).
- `single_term_opt_peak_impl(...)` (~lines 743-768).
- `single_term_opt_peak_batched_impl(...)` (~lines 770-825).

Be careful to KEEP `struct BatchedRes` (~lines 603-610) and `subset_open_aux` (~lines 585-601) — both are consumed by `PeakBatchedModel`.

- [ ] **Step 3: Fix dangling `\ref` in comments**

Run to find comment references to the removed names:
```bash
grep -nE '\\ref (single_term_opt_impl|single_term_opt_peak_impl|single_term_opt_peak_batched_impl|peak_dp|peak_dp_batched)|single_term_opt_impl' /Users/efv/code/SeQuant/SeQuant/core/optimize/single_term.hpp /Users/efv/code/SeQuant/SeQuant/core/optimize/cost_model.hpp
```
For each hit (e.g. the line ~235 comment `construction inside \ref single_term_opt_impl` and the model doc comments in `cost_model.hpp` that say `Reproduces \ref single_term_opt_impl` / `\ref peak_dp` / `\ref single_term_opt_peak_impl` / `\ref single_term_opt_peak_batched_impl` / `\ref peak_dp_batched`), reword to refer to the model itself rather than the deleted symbol — e.g. change `Reproduces \ref single_term_opt_peak_impl / \ref peak_dp exactly, factored into...` to `Implements the all-co-resident pebble-game DP, factored into...`. Do not leave a `\ref` pointing at a deleted symbol (it breaks the Doxygen build).

- [ ] **Step 4: Build the full SeQuant library + tests and run the suite**

Run:
```bash
clang-format --style=file -i /Users/efv/code/SeQuant/SeQuant/core/optimize/single_term.hpp /Users/efv/code/SeQuant/SeQuant/core/optimize/cost_model.hpp
cmake --build /Users/efv/code/SeQuant/build --target unit_tests-sequant 2>&1 | tail -8
/Users/efv/code/SeQuant/build/tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -8
```
Expected: clean build (no "undeclared identifier" for any removed name); all `[optimize]` pass.

- [ ] **Step 5: Commit**

```bash
cd /Users/efv/code/SeQuant
git add SeQuant/core/optimize/single_term.hpp SeQuant/core/optimize/cost_model.hpp
git commit -m "optimize: delete dead standalone single-term DP drivers/tables (peak_dp, *_impl, PeakRes)"
```

---

### Task 4: Negative-direction `CostModel` concept test (deferred follow-up #4)

Strengthen the concept-conformance SECTION with a negative check proving the `CostModel` concept actually rejects a non-conforming type (today only positive `static_assert`s exist, so a vacuously-true concept would pass unnoticed).

**Files:**
- Modify: `tests/unit/test_optimize.cpp` — inside `SECTION("CostModel concept conformance + custom model")` (~line 1009), after the existing positive `static_assert`s.

**Interfaces:**
- Consumes: `opt::detail::CostModel<M>` concept.
- Produces: a compile-time negative assertion.

- [ ] **Step 1: Add the negative `static_assert`**

After the positive `static_assert`s (and before / alongside the custom-model runtime check) in that SECTION, add a local incomplete type that is missing the required members and assert the concept rejects it:
```cpp
      // Negative direction: a type lacking State/Context/the six methods must
      // NOT satisfy CostModel (guards against a vacuously-true concept).
      struct NotAModel {};
      static_assert(!opt::detail::CostModel<NotAModel>);
```
(A local `struct` inside a Catch2 SECTION is a function-local type; that is fine for a `static_assert` referencing it in the same scope.)

- [ ] **Step 2: Build and run the SECTION**

Run:
```bash
cmake --build /Users/efv/code/SeQuant/build --target unit_tests-sequant 2>&1 | tail -5
/Users/efv/code/SeQuant/build/tests/unit/unit_tests-sequant "[optimize]" 2>&1 | tail -6
```
Expected: clean build (the negative `static_assert` compiles, i.e. the concept is genuinely `false` for `NotAModel`); `[optimize]` all pass. If the build fails on that `static_assert`, the concept is broken (accepts a non-model) — report it.

- [ ] **Step 3: Commit**

```bash
cd /Users/efv/code/SeQuant
git add tests/unit/test_optimize.cpp
git commit -m "optimize tests: add negative-direction CostModel concept check"
```

---

## Self-Review

**Spec coverage (vs. the option-B scope decision + Phase-3 spec §7):**
- "Remove the per-objective `*_impl` driver functions" → Task 3 deletes `single_term_opt_impl`, `single_term_opt_peak_impl`, `single_term_opt_peak_batched_impl`. ✓
- "Full DRY: also remove `peak_dp`/`peak_dp_batched`, route `peak_cost`/`peak_cost_batched` through the models" → Task 1 rewires the oracles to the models; Task 3 deletes both table-builders. ✓
- "Keep the brute-force oracles" → `brute_force_min_peak`/`batched_min_peak` are untouched (Global Constraints; not referenced by any task's removals). ✓
- "Delete the equivalence tests; rewire reconstruction tests to the model path" → Task 2. ✓
- Deferred follow-up #4 (negative concept assert) → Task 4. ✓

**Placeholder scan:** The only intentional "preserve verbatim" is Task 1 Step 4's simulation walk — the engineer reads the existing body in `single_term.hpp` and performs the named mechanical data-source swap (`at(n,B)`→`st[n][B]`, local `sz`/`Lof`→`ctx.sz`/`ctx.Lof`); the surrounding signature, model build, and solve are spelled out. No "TBD"/"handle edge cases"/"similar to" placeholders elsewhere.

**Type consistency:** `peak_cost`/`peak_cost_batched`/`reconstructed_batched_peak` signatures are copied unchanged from the current source; `PeakModel`/`PeakBatchedModel` ctor field order matches `cost_model.hpp` (`{idxsz}` and `{idxsz, is_batchable, batch, is_volatile_leaf}`); `State` access `st.back().peak` (PeakModel), `st.back()[0].peak` / `st[root][allK].peak` (PeakBatchedModel, `State = vector<BatchedRes>`) matches the model definitions; `ctx.sz`/`ctx.Lof` match the Phase-4 `PeakBatchedModel::Context` methods. Line numbers are approximate (`~`) and the engineer locates by symbol name.
