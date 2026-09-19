# External batch loops: runtime consumed_upscope fix -- implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the batched evaluator from hoisting loop-local external-carrying
intermediates to full extent above their external loop, so the C60 external-occ
replay peak drops from 5860.877 GB toward the factorizer-modeled ~35 GB.

**Architecture:** Implements fix (a) of
`doc/dev/specs/2026-07-27-external-loop-scope-semantics-design.md`. The
`scope_level`-over-`nest` computation (spec Sec. 2.1) ALREADY LANDED this session
(`cost_model.hpp:1922-1943`), yet the peak is unchanged -- proof (spec Sec. 1.2)
that correcting `scope_level` is necessary-but-not-sufficient. The remaining,
load-bearing fix is: the DP emits a `consumed_upscope` bit (spec Sec. 2.3), and
`hoist_invariants` excludes from its full-extent real-cache hoist any node that
carries a batched external mode and is NOT `consumed_upscope` (loop-local) -- so
that node is rebuilt sliced per external block instead of stored whole. The
scatter's full destination is already reached only by the scatter trigger (the
outermost carrier == a `consumed_upscope` node); Task 3 confirms this by trace and
adds a scatter gate ONLY if the trace shows a loop-local node scattering full.

**Tech Stack:** C++20, SeQuant, Catch2, the dry-run eval backend
(`core/eval/backends/dryrun/`). **Build/run in `cmake-build-release`** -- the
factorizer is too slow in Debug; the dry-run witness replays the real C60 forest
in seconds in Release. The eval unit tests compile in Release (the older
"debug-only" note is stale).

## Global Constraints

- **`scope_level` = innermost CARRIED loop** (external or contracted), already
  emitted over the combined `nest` (`cost_model.hpp:1922-1943`); this plan does
  NOT change it. `scope_level == -1` = carries no batched index; the sentinel
  `kNoBatchScopeLevel` (`INT_MIN`) = no order-aware emit (OFF path).
- **`consumed_upscope` (new bit).** A node is `consumed_upscope` iff it is the
  outermost carrier of a batched external mode -- its parent does not carry that
  mode (spec Sec. 2.3). For a single term the term's external output (the scatter
  target) is the only `consumed_upscope` node; every external-carrying
  intermediate is loop-local (`consumed_upscope == false`).
- **Full/real-cache hoist only for a node consumed up-scope.** A node carrying a
  batched external mode with `consumed_upscope == false` must NOT receive a
  full-extent real-cache entry anywhere; it is rebuilt sliced per block. A pure
  contracted intermediate (carries no external) is UNAFFECTED -- its hoisting is
  the legitimate contracted-path behavior and must not regress.
- **OFF-path byte-identical.** With `order_aware_recompute == false` OR
  `batch_spectator_indices == false`, no external loop is emitted, every node
  keeps `scope_level == kNoBatchScopeLevel`, `consumed_upscope` stays `false` and
  is never consulted (the hoist's existing `sl != kNoBatchScopeLevel` guard
  already excludes such nodes), and nothing here engages -- byte-identical.
- **Exactness.** Converged energy / replayed result unchanged; the batched-scratch
  peak reflects only loop-local block slices plus the genuine outputs' full
  results.
- **Process:** branch `evaleev/feature/suppress-heuristic-fallback`; commit there;
  ONE task per commit; NO push / PR / remote / repin. No `Co-Authored-By`
  trailers. ASCII hyphens only (pre-commit rejects U+2013 / U+00A0). clang-format
  is a pre-commit hook; re-add + re-commit if it reformats. Cap builds at `-j6`.

---

## File structure

- `core/eval/node_batch_annotation.hpp` -- `NodeBatchAnnotation`: gains
  `bool consumed_upscope = false;` (Task 1).
- `core/eval/eval_expr.{hpp,cpp}` -- `EvalExpr`: gains `batch_consumed_upscope()`
  accessor / `set_batch_consumed_upscope()` / member, stamped in binarize
  alongside `scope_level` (`eval_expr.cpp:601-605`) (Task 1).
- `core/optimize/cost_model.hpp` -- `reconstruct_batched_modes::build`
  (`:1793-1972`): compute `consumed_upscope` from the parent's carried batched-
  external set threaded down the recursion; set `ann.consumed_upscope` under
  `order_aware_recompute` (Task 2).
- `core/eval/eval.hpp` -- `hoist_invariants::collect` (`:1475-1487`): exclude a
  node carrying a batched external mode with `consumed_upscope == false` from the
  hoist target set; confirm/gate the scatter full-dest (`:1564-1650`) (Task 3).
- `tests/unit/test_optimize.cpp` (`[loop-tree]`) -- emit-level `consumed_upscope`
  test (Task 2).
- `tests/unit/test_eval_ta.cpp` (`[eval]`) -- runtime block-extent + slice-exact
  test (Task 3).
- `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`) -- C60 peak gate
  (Task 4).

---

## Task 1: Plumb the `consumed_upscope` bit (annotation -> EvalExpr, inert)

Pure plumbing: add the field, accessor, and binarize stamp. No emit logic yet
(the DP leaves it `false`), so this task is behavior-neutral everywhere.

**Files:** `core/eval/node_batch_annotation.hpp`, `core/eval/eval_expr.hpp`,
`core/eval/eval_expr.cpp`, `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Produces: `NodeBatchAnnotation::consumed_upscope` (default `false`);
  `EvalExpr::batch_consumed_upscope() -> bool` / `set_batch_consumed_upscope(bool)`;
  the binarize stamp copying `ann.consumed_upscope` onto the node.

- [ ] **Step 1: Write the failing test** (`test_optimize.cpp`, `[loop-tree]`).
  A minimal check that the accessor exists and defaults `false` through binarize:

```cpp
TEST_CASE("consumed_upscope default is false", "[loop-tree]") {
  using namespace sequant;
  // Any binarized single contraction; no order-aware emit -> default false.
  auto expr = parse_expr(L"A{i1;a1} * B{a1;i2}", Symmetry::nonsymm);
  auto node = binarize(expr);  // default BinarizationOptions
  // Root is a contraction node; the accessor must exist and default false.
  REQUIRE(node->batch_consumed_upscope() == false);
}
```
  (Reuse whatever `binarize(...)` / `parse_expr(...)` idiom the file already
  uses -- match a neighbouring `[loop-tree]` test's construction exactly.)
- [ ] **Step 2: Run to confirm it fails.**
  `./cmake-build-release/tests/unit/unit_tests-sequant "[loop-tree]"` -- FAILS to
  compile (`batch_consumed_upscope` undeclared).
- [ ] **Step 3: Implement the plumbing.**
  - `node_batch_annotation.hpp`, in `struct NodeBatchAnnotation` after
    `effective_count`:
    ```cpp
    /// True iff this node is the outermost carrier of a batched external mode
    /// (its parent does not carry that mode) -- the scatter target / a genuine
    /// output read from an enclosing (up-scope) cache. A loop-local
    /// external-carrying intermediate is false and must NOT be hoisted full.
    bool consumed_upscope = false;
    ```
  - `eval_expr.hpp`: after `batch_effective_count()` (`:317-319`) add
    ```cpp
    [[nodiscard]] bool batch_consumed_upscope() const noexcept {
      return batch_consumed_upscope_;
    }
    ```
    after `set_batch_effective_count` (`:330-332`) add
    ```cpp
    void set_batch_consumed_upscope(bool v) noexcept {
      batch_consumed_upscope_ = v;
    }
    ```
    and in the members after `batch_effective_count_` (`:356`) add
    ```cpp
    /// See \c batch_consumed_upscope.
    bool batch_consumed_upscope_ = false;
    ```
  - `eval_expr.cpp`, in the binarize node-axes stamp after
    `result.set_batch_effective_count(ann.effective_count);` (`:605`) add
    ```cpp
    result.set_batch_consumed_upscope(ann.consumed_upscope);
    ```
- [ ] **Step 4: Run to confirm it passes.** Then no-regression:
  `./cmake-build-release/tests/unit/unit_tests-sequant "[optimize][eval]"` --
  counts unchanged (this task adds one defaulted field).
- [ ] **Step 5: clang-format touched files + Commit.**

---

## Task 2: Emit `consumed_upscope` in the DP `build` walk

Set the bit per node: `consumed_upscope` iff the node carries a batched external
mode its PARENT does not carry (spec Sec. 2.3). Thread the parent's carried
batched-external set down the `build` recursion.

**Files:** `core/optimize/cost_model.hpp` (`reconstruct_batched_modes`,
`build` `:1793-1972`), `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Consumes: Task 1's `NodeBatchAnnotation::consumed_upscope`.
- Produces: `ann.consumed_upscope` set (under `order_aware_recompute` only);
  `false` on the OFF path.

- [ ] **Step 1: Write the failing test** (`test_optimize.cpp`, `[loop-tree]`).
  Build the network `R{i1,i2} = M{i1,i2,i3} * B{i3}` (an external pair `i1,i2`
  free on the result, contracted `i3`), configure `is_batchable_external_index =
  {i1,i2}`, `is_batchable_contracted_index = {i3}`,
  `batch_spectator_indices = true`, `order_aware_recompute = true`, run
  `reconstruct_batched_modes`. Locate the ROOT `R` (carries `i1,i2`, is the
  output) and the intermediate `M` (carries `i1,i2`, its parent `R` also carries
  them). Assert:
  ```cpp
  REQUIRE(root_ann.consumed_upscope == true);   // outermost carrier
  REQUIRE(M_ann.consumed_upscope == false);      // parent carries i1,i2 -> loop-local
  ```
  Reuse the network-build + `model.order_aware_recompute = true` idiom already in
  `test_optimize.cpp` (search for `order_aware_recompute` and the
  `is_batchable_external_index` setter used by the `[loop-tree]` /
  `[role-filter]` tests) and the RPN/annotation-inspection idiom those tests use
  to read `node_axes`.
- [ ] **Step 2: Run to confirm it fails.**
  `./cmake-build-release/tests/unit/unit_tests-sequant "[loop-tree]"` -- FAILS
  (`consumed_upscope` is `false` on the root; not yet computed).
- [ ] **Step 3: Implement.** In `reconstruct_batched_modes`, before `build`,
  compute the batched-external mask once:
  ```cpp
  // Modes that are batched externals: adopted external placements with >1 batch.
  std::size_t batched_ext_mask = 0;
  for (auto const& [nn, mask] : placed_at_node) batched_ext_mask |= mask;
  batched_ext_mask |= chosen_seed_mask;  // legacy root-seed path
  {
    std::size_t filtered = 0;
    for (std::size_t k = 0; k < ctx.m; ++k)
      if ((batched_ext_mask & (std::size_t{1} << k)) && ctx.nbatches[k] > 1.0 &&
          is_external_mode(ctx, k))
        filtered |= (std::size_t{1} << k);
    batched_ext_mask = filtered;
  }
  ```
  Extend `build`'s signature with the parent's carried batched-external set:
  change the `std::function<...>` type and lambda to
  `build(std::size_t n, std::size_t B, int idx,
         container::svector<std::uint8_t> nest, std::size_t parent_ext)`, root
  call `build(root, 0, best, {}, 0)`. Inside, after `ann.axes` is set and within
  the `if (order_aware_recompute)` block, compute:
  ```cpp
  std::size_t const n_ext = ctx.open_modes[n] & batched_ext_mask;
  ann.consumed_upscope = (n_ext & ~parent_ext) != 0;  // n carries an ext its parent does not
  ```
  and recurse into children with `parent_ext = n_ext`:
  ```cpp
  EvalSequence s = build(fs, C, fi, child_nest, n_ext);
  EvalSequence b = build(ss, C, si, child_nest, n_ext);
  ```
  Guard the whole computation under `order_aware_recompute` (already the
  enclosing block) so the OFF path leaves `consumed_upscope == false`. Do NOT set
  it in the leaf early-return (`std::popcount(n) == 1`) -- leaves are not hoist
  targets.
- [ ] **Step 4: Run to confirm it passes.** Then no-regression:
  `./cmake-build-release/tests/unit/unit_tests-sequant "[optimize][role-filter][loop-tree]"`
  -- OFF-path counts unchanged; record counts.
- [ ] **Step 5: clang-format + Commit.**

---

## Task 3: Runtime -- exclude loop-local external intermediates from the full hoist

The pinned root cause: `hoist_invariants` collects a loop-local external-carrying
descendant (`sl != kNoBatchScopeLevel && sl < depth`) and builds it via
`evaluate(d, le)` at FULL external extent (the `le` at hoist time has not sliced
the external mode -- the scatter loop has not started), storing it in an ancestor
cache read whole by every block. Evidence: the summand-46 operand
`I(i2,i1,μ̃;a1<i1,i2>)` is `Cache|Access`'d at 8.7 GB full BEFORE
`BatchScatter|Begin` (investigation ledger). Fix: never hoist such a node.

**Files:** `core/eval/eval.hpp` (`hoist_invariants::collect` `:1475-1487`; the
external scatter branch `:1564-1650`), `tests/unit/test_eval_ta.cpp`.

**Interfaces:**
- Consumes: Task 2's `batch_consumed_upscope()`.
- Produces: a runtime where an external-carrying loop-local intermediate is
  realized at block-external extent, never stored full.

- [ ] **Step 1: Confirm the site + whether a scatter gate is also needed.**
  Read `hoist_invariants` (`:1471-1540`) and the external scatter branch
  (`:1564-1650`). Confirm (a) `collect` (`:1475-1487`) is the only place a
  descending node is built as a full unit via `evaluate(d, le)` (`:1535`), and
  (b) whether the scatter's `pre_sized_zeros_over_mode` (`:1634`) can be reached
  by a NON-`consumed_upscope` node: a loop-local node sharing the SAME external
  mode as an enclosing scatter has that mode already sliced by `le_g`
  (`:1616-1622`) so it yields a single batch and is skipped (`:1573-1576`
  invariant) -- meaning only the outermost carrier (a `consumed_upscope` node)
  scatters. Record the finding in the commit message. If the read shows a
  loop-local node CAN reach `pre_sized_zeros_over_mode` with an unsliced external,
  add the scatter gate in Step 3 (below); otherwise the hoist-exclusion alone
  suffices.
- [ ] **Step 2: Write the failing test** (`test_eval_ta.cpp`, `[eval]`). Build an
  external-batched network `R{e} = M{e,c} * B{c}` whose intermediate `M{e,c}`
  carries the external `e` and is contracted into `R` within the `e`-loop
  (a giant-shaped node). Configure external `e` batched (block < full extent),
  `order_aware_recompute = true`, `batch_spectator_indices = true`. Replay through
  the batched evaluator and assert:
  ```cpp
  // 1. slice-exact: batched result equals the unbatched reference exactly.
  REQUIRE(norm(batched_R - reference_R) == Approx(0.0).margin(1e-12));
  // 2. M is never realized at full-e extent: no hoisted full entry for M.
  //    (Use the dry-run op trace / a size probe: the max realized extent of e
  //    on M is the block size, not the full range.)
  REQUIRE(max_realized_e_extent_of_M <= block_size);
  ```
  Use the `[eval]` dry-run / TA replay idiom already in `test_eval_ta.cpp`
  (search for an existing external-batched or `batch_spectator_indices` test to
  copy the harness; if none, model the network on the
  `test_optimize.cpp` `[loop-tree]` builder and evaluate it via the same path the
  `[eval]` suite uses).
- [ ] **Step 3: Implement the hoist-exclusion.** The loop-local giants carry the
  external mode via D1's `stamp_carrying_descendants`, so every carrying descendant
  has an `External` entry in `batched_here()` (investigation ledger: the giants are
  stamped `batched_here={i_1:EXT i_2:EXT}`). The predicate is therefore exact
  without a new emitted bit: `has-External-in-batched_here && !consumed_upscope`.
  In `hoist_invariants::collect` (`:1475-1487`), before the hoist `if`, add:
  ```cpp
  auto has_ext = [](node_t const& x) {
    for (auto const& [ix, knd] : x->batched_here())
      if (knd == BatchModeType::External) return true;
    return false;
  };
  ```
  and extend the guard so a loop-local external-carrying node is not collected (it
  is then descended into via the existing `self(self, ...)` fall-through and
  rebuilt per block):
  ```cpp
  int const sl = n->batch_scope_level();
  bool const ext_loop_local = has_ext(n) && !n->batch_consumed_upscope();
  if (sl != kNoBatchScopeLevel && sl < static_cast<int>(depth) &&
      !subtree_any(n, is_volatile) && !ext_loop_local) {
    if (std::none_of(targets.begin(), targets.end(),
                     [&](node_t const* p) { return eq(*p, n); }))
      targets.push_back(&n);
    return;  // built as a unit -- do not descend into it
  }
  self(self, n.left());
  self(self, n.right());
  ```
  A pure contracted intermediate has no `External` entry -> `has_ext == false` ->
  `ext_loop_local == false` -> its hoisting is UNCHANGED. The scatter trigger `R`
  is `consumed_upscope == true` -> not excluded (and is never collected as a
  descendant anyway).
  If Step 1 found a loop-local node can reach the scatter's
  `pre_sized_zeros_over_mode` (`:1634`) with an unsliced external, ALSO gate it so
  a non-`consumed_upscope` trigger returns its block partial instead of a full
  destination; otherwise leave the scatter branch unchanged (record the finding).
- [ ] **Step 4: Run to confirm it passes** + reference-equality. No-regression:
  `./cmake-build-release/tests/unit/unit_tests-sequant "[eval]~[tot]"` and
  `"[cache_manager]"` unchanged on the OFF path.
- [ ] **Step 5: clang-format + Commit** (record the Step 1 finding -- hoist-only
  vs hoist+scatter -- in the message).

---

## Task 4: GATE -- C60 external-occ peak drops on the dry-run witness

**Files:** `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`; external-occ
config: occ in the external role, `batch_spectator_indices = true`,
`order_aware_recompute = true`, `peak_threshold = 100 GB`).

**Interfaces:** Consumes Tasks 1-3.

- [ ] **Step 1: Run the witness** (full term count):
  `SEQUANT_UT_DRYRUN_NTERMS=<full> ./cmake-build-release/tests/unit/unit_tests-sequant "[dryrun-occ-veto]"`.
  Assert the summand-46 operand/giant is realized at BLOCK occ (no `result=`
  near the 2930 GB full giant / 8.7 GB full operand), and the batched-scratch
  peak (`cost_profile(...).peak_bytes`) drops from 5860.877 GB toward the
  factorizer-modeled ~35 GB. Record before/after peak.
- [ ] **Step 2: No-regression.**
  `./cmake-build-release/tests/unit/unit_tests-sequant "[optimize][eval][cache_manager][dryrun-objective]"`
  -- OFF-path byte-identical; the `peak_threshold=+inf` heuristic path unchanged.
- [ ] **Step 3: Record** before/after peak + the giant's realized extent in the
  commit message. **Commit.**

---

## Task 5: MPQC enablement (handoff -- separate session)

**Files:** MPQC `src/mpqc/chemistry/qc/lcao/cc/csv_batch_policy.h`,
`mpqc4/external/versions.cmake` (own commit).

- [ ] **Step 1:** Ensure MPQC sets BOTH `order_aware_recompute` and
  `batch_spectator_indices` on the CCk `BatchPolicy` for the external-occ path,
  with occ in the EXTERNAL role (`is_batchable_external_index`). Local CSV/PNO-CCSD
  sanity: converged energy unchanged vs unbatched.
- [ ] **Step 2:** Repin `MPQC_TRACKED_SEQUANT_TAG` (own commit, update
  `PREVIOUS_*`), paired with the MPQC code that requires it. **Do not push /
  repin unilaterally -- hand to the user.**
- [ ] **Step 3:** Hand off Owl C60 occ+aux end-to-end (separate SLURM session,
  refresh the TA pin): giants bounded under budget, energy unchanged.

---

## Notes for the executor

- **The scope_level fix is already in** (`cost_model.hpp:1922-1943`); do NOT
  re-implement it. This plan adds the `consumed_upscope` bit + the hoist-exclusion
  ONLY.
- **The exclusion is narrow.** Only a node that CARRIES a batched external mode
  and is loop-local is excluded from the hoist. A pure contracted intermediate is
  untouched -- its hoisting is the legitimate contracted-path win and any
  regression there is a bug.
- **OFF-path byte-identical at every task boundary** (`order_aware_recompute` OR
  `batch_spectator_indices` off): the hoist's existing `sl != kNoBatchScopeLevel`
  guard already excludes every node on the OFF path, so the new predicate is inert
  there.
- **Follow-on (b), NOT this plan:** external forest-merge across terms to recover
  the cross-term CSE that per-term external scatter drops (spec Sec. 8) -- its own
  spec/plan after this lands.
```
