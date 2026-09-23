# Lifetime-mask placement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the pre-slice-on-use placement apparatus (per-occurrence scalar `batch_scope_level`, the `consumed_upscope` bit, run-scope veto part (b), and `hoist_invariants`) with a per-canonical-node lifetime mask plus per-level placement, so cross-term CSE becomes operational under external batching.

**Architecture:** A forest-wide pass computes, per canonical eval node, the set of batch modes that slice it in *every* occurrence (a proto-aware cross-occurrence intersection). The runtime derives a node's placement level at eval time from that mask intersected with the live batch context; the run-scope builder admits only all-full nodes; each realized loop level materializes its own nodes locally, retiring the hoist pre-build pass.

**Tech Stack:** C++20, SeQuant `core/eval` + `core/optimize`. Tests: Catch2 (`tests/unit/test_eval_dryrun.cpp`, `test_optimize.cpp`). Build: `cmake-build-release` (factorizer-heavy dry-run must not run in Debug).

## Global Constraints

- **OFF-path byte-identity is mandatory.** On the order-blind path (`order_aware_recompute == false`) every node has an empty mask (all-full), the veto and per-level placement are inert, and the batched replay runs exactly as before. Any task that could touch the OFF path must assert byte-identity.
- **No `Co-Authored-By` trailers**; commit messages are plain.
- **ASCII only** in all committed text: no U+2013 en-dash, no U+00A0 nbsp (pre-commit rejects them).
- **Build gates use `__OPTIMIZE__`, never `NDEBUG`** (this project defines `NDEBUG` even in Debug builds; `SEQUANT_ASSERT` is elided only in release).
- **Cap builds at `-j6`.**
- **Runtime + emit only.** Do NOT change the DP cost/sizing (`szcell`, cell-restricted sizes). Cost-model lifetime-awareness is a separate later spec.
- **Acceptance gate (whole plan):** C60 `[.][dryrun-occ-veto]` (`SEQUANT_UT_DRYRUN_NTERMS=55`, `cmake-build-release`) replay peak `<= ~443 GB` (same-or-better); `avoidable_time` no worse; he10 CSV-CCSD energy correct in MPQC (Delta < 1e-7 vs `-0.33231474200227867`); OFF-path byte-identical.

## File Structure

- **Create** `SeQuant/core/eval/lifetime_mask.hpp` -- the forest-wide meet: given the eval forest, compute and stamp each canonical node's sliced-mode mask. One clear responsibility: mask computation. Header-only (matches the eval subsystem).
- **Modify** `SeQuant/core/eval/eval_expr.hpp` -- add the mask field + accessors on `EvalExpr`; remove `batch_consumed_upscope_` and `batch_scope_level_` (Task 4).
- **Modify** `SeQuant/core/eval/node_batch_annotation.hpp` -- remove `consumed_upscope` and `scope_level` (Task 4); keep `axes` (the meet's input) and `effective_count`.
- **Modify** `SeQuant/core/eval/cache_manager.hpp` -- veto part (b) reads the mask (Task 2).
- **Modify** `SeQuant/core/eval/eval.hpp` -- run the meet before evaluation; replace `hoist_invariants` with per-level placement (Task 3).
- **Modify** `SeQuant/core/optimize/cost_model.hpp` -- stop emitting `scope_level`/`consumed_upscope` (Task 4).
- **Test** `SeQuant/tests/unit/test_lifetime_mask.cpp` (new, added to `tests/unit/CMakeLists.txt`) -- unit tests for the meet on the survey fixtures.

---

## Task 1: Per-canonical sliced-mode mask + forest meet

Compute, for each canonical node, the set of batch modes (canonical `Index`es) that slice it in *every* occurrence, and stamp it on every occurrence of that node. Nothing consumes the mask yet, so behavior is unchanged; this task is purely additive and independently testable against the survey fixtures.

**Files:**
- Create: `SeQuant/core/eval/lifetime_mask.hpp`
- Modify: `SeQuant/core/eval/eval_expr.hpp` (add mask field + accessors, ~after line 294 `set_batched_here`)
- Create test: `SeQuant/tests/unit/test_lifetime_mask.cpp`
- Modify: `SeQuant/tests/unit/CMakeLists.txt`

**Interfaces:**
- Consumes: `EvalExpr::batched_here()` (the `svector<pair<Index,BatchModeType>>` External/Contracted stamps), `EvalExpr::canon_indices()`, `EvalExpr::hash_value()`, `Index::proto_indices()`.
- Produces:
  - `EvalExpr::sliced_modes() const noexcept -> container::svector<Index> const&` (canonical batch modes that slice this node under the meet; empty == all-full).
  - `EvalExpr::set_sliced_modes(container::svector<Index>) noexcept`.
  - `EvalExpr::mask_all_full() const noexcept -> bool` (== `sliced_modes().empty()`).
  - `template <meta::eval_node_range R> void stamp_lifetime_masks(R const& forest) noexcept;` in `namespace sequant` (header `lifetime_mask.hpp`).

- [ ] **Step 1: Write the failing test**

Create `SeQuant/tests/unit/test_lifetime_mask.cpp`. Build a two-occurrence synthetic forest where one canonical node is sliced by occ mode `i` in occurrence A and left full in occurrence B (meet => all-full), and a second canonical node is sliced by the same proto pair `(i,j)` in both occurrences (meet => sliced by `{i,j}`). Assert:

```cpp
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <catch2/catch_test_macros.hpp>
// ... build forest `nodes` (helpers as in test_eval_dryrun.cpp) ...
TEST_CASE("lifetime mask cross-occurrence meet", "[lifetime_mask]") {
  stamp_lifetime_masks(nodes);
  // node present sliced in A, full in B => all-full (meet demotes to full):
  CHECK(node_full.mask_all_full());
  // node sliced by (i,j) in every occurrence => sliced by exactly {i, j}:
  CHECK_FALSE(node_pair.mask_all_full());
  CHECK(as_set(node_pair.sliced_modes()) == index_set({i, j}));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build cmake-build-release -j6 --target unit_tests && ./cmake-build-release/tests/unit/unit_tests "[lifetime_mask]"`
Expected: FAIL to compile (`sliced_modes`/`stamp_lifetime_masks` undefined).

- [ ] **Step 3: Add the mask field + accessors to `EvalExpr`**

In `eval_expr.hpp`, after `set_batched_here` (~line 294), add:

```cpp
  /// \brief Canonical batch modes that slice this node in EVERY occurrence
  /// (the cross-occurrence meet; see stamp_lifetime_masks). Empty => all-full
  /// (block-agnostic, run-scope). Proto-aware: a composite slot contributes its
  /// proto indices. Set by stamp_lifetime_masks; empty by default (OFF path).
  [[nodiscard]] container::svector<Index> const& sliced_modes() const noexcept {
    return sliced_modes_;
  }
  void set_sliced_modes(container::svector<Index> m) noexcept {
    sliced_modes_ = std::move(m);
  }
  [[nodiscard]] bool mask_all_full() const noexcept {
    return sliced_modes_.empty();
  }
```

and the member (near `batch_axes_`, ~line 364):

```cpp
  /// See \c sliced_modes.
  container::svector<Index> sliced_modes_{};
```

- [ ] **Step 4: Implement the meet in `lifetime_mask.hpp`**

Create `SeQuant/core/eval/lifetime_mask.hpp`. Algorithm: (1) walk every tree top-down accumulating each node's occurrence-local sliced-mode set = ancestors' + self `batched_here()` External stamps, expanded proto-aware (a batched composite index contributes its proto indices); (2) group occurrences by `hash_value()`; (3) per canonical node, meet = set-intersection (by `Index` identity) of the per-occurrence sets; (4) stamp every occurrence with the meet.

```cpp
#ifndef SEQUANT_EVAL_LIFETIME_MASK_HPP
#define SEQUANT_EVAL_LIFETIME_MASK_HPP

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/meta.hpp>

namespace sequant {

/// Stamp each canonical eval node's cross-occurrence sliced-mode mask
/// (EvalExpr::sliced_modes). A mode slices a canonical node iff it slices EVERY
/// occurrence (proto-aware). Idempotent; no-op on the OFF path (no External
/// stamps => every set empty => every mask empty => all-full, byte-identical).
template <meta::eval_node_range R>
void stamp_lifetime_masks(R const& forest) noexcept {
  using Node = std::ranges::range_value_t<R>;
  // occ-local sliced modes for each occurrence, keyed by node identity.
  std::unordered_map<Node, container::svector<Index>,
                     TreeNodeHasher<Node>, TreeNodeEqualityComparator<Node>>
      meet;  // hash -> running intersection
  bool first_seen_init = false;  // per-key init handled inline
  auto ext_modes_of = [](Node const& n) {
    container::svector<Index> v;
    for (auto const& [ix, kind] : n->batched_here())
      if (kind == BatchModeType::External) {
        if (ix.has_proto_indices())
          for (auto const& p : ix.proto_indices()) v.push_back(p);
        else
          v.push_back(ix);
      }
    return v;
  };
  // Pass 1: per occurrence, accumulate ancestor+self external modes and meet.
  std::vector<std::pair<Node, container::svector<Index>>> occ;  // for pass 2
  auto walk = [&](auto&& self, Node const& n,
                  container::svector<Index> acc) -> void {
    if (!n.leaf()) {
      auto here = ext_modes_of(n);
      for (auto const& ix : here) acc.push_back(ix);
      occ.emplace_back(n, acc);
      auto it = meet.find(n);
      if (it == meet.end())
        meet.emplace(n, acc);  // first occurrence seeds the intersection
      else
        intersect_in_place(it->second, acc);  // by Index identity
    }
    if (!n.leaf()) {
      self(self, n.left(), acc);
      self(self, n.right(), acc);
    }
  };
  for (auto const& tree : forest) walk(walk, tree, {});
  // Pass 2: stamp every occurrence with its canonical meet.
  for (auto& [n, _] : occ) {
    if (auto it = meet.find(n); it != meet.end())
      n->set_sliced_modes(it->second);
  }
  (void)first_seen_init;
}

}  // namespace sequant
#endif
```

Provide `intersect_in_place(container::svector<Index>&, container::svector<Index> const&)` as a file-local helper: keep only elements present in both (by `Index` equality). Note in a comment: intersection by canonical `Index` identity is correct because canonicalization gives the shared proto pair consistent labels across occurrences, while a genuinely block-agnostic node (`s*C`) has occurrences with disjoint concrete modes and intersects to empty.

- [ ] **Step 5: Run the test to verify it passes**

Run: `cmake --build cmake-build-release -j6 --target unit_tests && ./cmake-build-release/tests/unit/unit_tests "[lifetime_mask]"`
Expected: PASS.

- [ ] **Step 6: Add the survey fixtures as a regression test**

Add to `test_lifetime_mask.cpp` a `[lifetime_mask][survey]` case driving the same `[.][dryrun-occ-veto]` term set (`SEQUANT_UT_DRYRUN_NTERMS` small, e.g. 12) and asserting the two ground-truth nodes from `.superpowers/sdd/lifetime-mask-ir-survey.md`: hash `1989507463377952644` (`s*C`) is `mask_all_full()`; hash `15545560759149115397` has `sliced_modes()` equal to its occ pair `{i_1, i_2}` (proto-aware, so the PNO slot does not add extra modes). Reuse the dry-run forest builder already in `test_eval_dryrun.cpp` (extract a shared helper if needed).

- [ ] **Step 7: Run + commit**

Run: `./cmake-build-release/tests/unit/unit_tests "[lifetime_mask]"` -> PASS.
```bash
git add SeQuant/core/eval/lifetime_mask.hpp SeQuant/core/eval/eval_expr.hpp \
        SeQuant/tests/unit/test_lifetime_mask.cpp SeQuant/tests/unit/CMakeLists.txt
git commit -m "eval: per-canonical sliced-mode lifetime mask + cross-occurrence meet"
```

---

## Task 2: Run-scope veto reads the mask

Switch the run-scope builder's batch-variant veto part (b) from the per-occurrence scalar `batch_scope_level() >= 0` to the mask (`!mask_all_full()`), and run the meet before the cache is built. Keep veto part (a) (contracted `sliced_batch_axis`) unchanged.

**Files:**
- Modify: `SeQuant/core/eval/cache_manager.hpp` (~line 668, the `batch_variant` computation)
- Modify: `SeQuant/core/eval/eval.hpp` (call `stamp_lifetime_masks` before `cache_manager(...)` is built for the evaluation)

**Interfaces:**
- Consumes: `EvalExpr::mask_all_full()` (Task 1); `stamp_lifetime_masks` (Task 1).
- Produces: no new symbols; a behavior change gated by the mask.

- [ ] **Step 1: Write the failing test**

Add to `test_lifetime_mask.cpp` a `[lifetime_mask][veto]` case: build the run-scope cache for a small forest via the same path `evaluate` uses (or call `cache_manager(...)` directly with the mask stamped), and assert that a node with a non-empty mask is NOT registered at run scope while an all-full repeated node IS. Assert an OFF-path forest (no External stamps) registers exactly the same set as before the change (capture the set with masks all-empty).

- [ ] **Step 2: Run to verify it fails**

Run: `./cmake-build-release/tests/unit/unit_tests "[lifetime_mask][veto]"` after building -> FAIL (veto still keys off `batch_scope_level`).

- [ ] **Step 3: Ensure the mask is stamped before the cache is built**

In `eval.hpp`, at the entry that builds the run-scope `cache_manager` for a batched evaluation (the caller of `cache_manager(...)` that feeds `make_batched_custom_evaluator`), call `stamp_lifetime_masks(nodes)` on the forest first. Guard: it is a no-op on the OFF path (empty masks), so it may run unconditionally.

- [ ] **Step 4: Switch the veto predicate**

In `cache_manager.hpp` replace (~line 668):

```cpp
    bool const batch_variant = sliced_batch_axis || n->batch_scope_level() >= 0;
```
with:
```cpp
    // (b): a node whose cross-occurrence mask is non-empty is sliced by some
    // enclosing external mode in every occurrence => batch-variant => refused
    // run-scope residence. all-full (empty mask; incl. the OFF path) is admitted.
    bool const batch_variant = sliced_batch_axis || !n->mask_all_full();
```
Update the surrounding doc comment (~lines 642-652) to describe the mask, not `batch_scope_level`.

- [ ] **Step 5: Run the test + OFF-path check**

Run: `./cmake-build-release/tests/unit/unit_tests "[lifetime_mask]"` -> PASS.
Run the existing OFF-path evaluator tests: `./cmake-build-release/tests/unit/unit_tests "[eval]"` -> PASS (byte-identical membership on the order-blind path).

- [ ] **Step 6: Commit**

```bash
git add SeQuant/core/eval/cache_manager.hpp SeQuant/core/eval/eval.hpp \
        SeQuant/tests/unit/test_lifetime_mask.cpp
git commit -m "eval: run-scope veto keys off the lifetime mask, not scalar scope_level"
```

---

## Task 3: Per-level placement replaces hoist_invariants

Replace the two `hoist_invariants(bs.cache, cache, ...)` calls (scatter branch ~line 1700, group/contracted branch ~line 1864 of `eval.hpp`) with per-level placement: at each realized loop level, materialize -- once per block -- the member-subtree nodes whose deepest sliced mode present in the live batch context is exactly this loop's mode, on this level's scratch. Deeper replays find them via the scope chain (fall-through + slice-on-use). This is THE behavioral task; the C60 gate is its spec.

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (remove `hoist_invariants` ~lines 1555-1643; add a per-level `place_at_this_level` helper; call it at the two branch sites)

**Interfaces:**
- Consumes: `EvalExpr::sliced_modes()` (Task 1); `CacheManager::batch_context()`, `access_at`, `ensure_hoist_slot`, `store`, `alive` (existing); `sliced_leaf`, `make_batched_scratch`, `slice_to_use` (existing).
- Produces: file-local `place_at_this_level(scratch_cache, parent_cache, member_roots, K)` (replaces `hoist_invariants`), invoked once per block after the block's `batch_context` is installed and before the deeper nest replays.

- [ ] **Step 1: Write the failing gate probe (wrong-block isolation)**

Add a `[.][dryrun-occ-veto]`-adjacent unit assertion (or extend the existing dry-run harness) that fails if a mask-scope-L node's stored value is served to a different L-block than the one it was built for. Concretely: run the small (`NTERMS=12`) occ-batched forest and assert the replay peak equals the DP model peak for that term set (the phase-1 fixture already records this equality). Before the change it will still hold via hoist; the test locks the invariant so Task 3 cannot regress it.

- [ ] **Step 2: Establish the current C60 baseline**

Run (record numbers, do not modify): `SEQUANT_UT_DRYRUN_NTERMS=55 ./cmake-build-release/tests/unit/unit_tests "[.][dryrun-occ-veto]"`. Expected: replay peak ~443.55 GB (the phase-1 post-slice-on-use value). This is the number Task 3 must hold (`<= ~443`).

- [ ] **Step 3: Add `place_at_this_level`**

In `eval.hpp`, add a helper that, given the current level's scratch (`bs.cache`, whose `batch_context` already has this block pushed), the member roots, and the picked mode `K`, walks each member subtree and for every non-leaf, non-volatile node whose sliced modes intersect the live `batch_context` with **deepest** element == `K` (i.e. this loop is the node's innermost enclosing sliced mode present in the nest), builds it once via `evaluate(node, sliced_leaf)`-equivalent (so it is sliced to the enclosing blocks) and stores it on `bs.cache` via `ensure_hoist_slot` + `store` (skip if `alive`). Nodes whose deepest live-context sliced mode is an OUTER level were placed by that outer level; nodes with no live-context sliced mode (all-full within this nest) belong to the root cache. Do NOT descend into a node once placed.

The predicate replaces hoist's `sl != kNoBatchScopeLevel && sl < depth && !ext_loop_local`: the mask makes `ext_loop_local` unnecessary (a node whose own result carries the External mode has `K` in its sliced set, so its deepest-live-sliced mode IS this level => placed here, loop-local, exactly right).

- [ ] **Step 4: Wire it at both branch sites; delete `hoist_invariants`**

Scatter branch (~1700): replace `hoist_invariants(bs.cache, cache, {&node});` with a `place_at_this_level` call inside the per-block loop, after the block's `batch_context` is installed (the reinstall path around ~1725-1734), before the deeper `evaluate(node, leaf_evaluator, bs.cache)`. Group branch (~1864): same, inside the per-member per-block replay (~1883-1900). Remove the `hoist_invariants` lambda (1555-1643). Keep `set_parent` wiring (the scope chain is still needed for fall-through).

- [ ] **Step 5: Run the small gate + OFF path**

Run: `./cmake-build-release/tests/unit/unit_tests "[.][dryrun-occ-veto]"` (small NTERMS via the Step-1 test) -> peak == model. Run `"[eval]"` -> OFF-path PASS (byte-identical: empty masks => `place_at_this_level` finds nothing => replay unchanged).

- [ ] **Step 6: Run the C60 gate**

Run: `SEQUANT_UT_DRYRUN_NTERMS=55 ./cmake-build-release/tests/unit/unit_tests "[.][dryrun-occ-veto]"`. Expected: replay peak `<= ~443 GB`; `avoidable_time` no worse than baseline (Step 2). If peak regresses, STOP and investigate placement (do not add compensating hacks) -- see systematic-debugging.

- [ ] **Step 7: Commit**

```bash
git add SeQuant/core/eval/eval.hpp SeQuant/tests/unit/test_eval_dryrun.cpp
git commit -m "eval: per-level placement replaces hoist_invariants (mask-driven)"
```

---

## Task 4: Retire consumed_upscope and scalar scope_level

Now that nothing consumes them, remove the `consumed_upscope` bit and the scalar `batch_scope_level` end to end: the emit (cost_model), the annotation struct, and the `EvalExpr` fields/accessors. Purely subtractive; compiles and passes only if Tasks 2-3 fully migrated the consumers.

**Files:**
- Modify: `SeQuant/core/optimize/cost_model.hpp` (~1931-1996: drop the `scope_level`/`consumed_upscope` emit; keep `axes`, `effective_count`)
- Modify: `SeQuant/core/eval/node_batch_annotation.hpp` (remove `scope_level`, `consumed_upscope`; keep `axes`, `effective_count`; keep `kNoBatchScopeLevel` only if still referenced -- else remove)
- Modify: `SeQuant/core/eval/eval_expr.hpp` (remove `batch_scope_level_`/`batch_consumed_upscope_` fields + accessors + setters)
- Modify: `SeQuant/core/eval/eval_expr.cpp` (~603-606: drop the two `set_batch_*` calls)

**Interfaces:**
- Consumes: nothing new.
- Produces: removal only. After this task `grep -rn "batch_scope_level\|consumed_upscope\|kNoBatchScopeLevel" SeQuant/core` returns nothing (or only an intentionally-kept `effective_count` neighbor).

- [ ] **Step 1: Remove the emit**

In `cost_model.hpp`, inside `if (order_aware_recompute) { ... }` (~1942-1996), delete the `placement`/`ann.scope_level` block, the `SEQUANT_NEST_PROBE` block, and the `ann.consumed_upscope` block. Keep the `effective_count` block. (If `effective_count` alone no longer needs the `order_aware_recompute` guard body, keep the guard for it.)

- [ ] **Step 2: Remove struct + EvalExpr members**

Delete `scope_level` and `consumed_upscope` from `NodeBatchAnnotation`. Delete `batch_scope_level_`, `batch_consumed_upscope_`, and their `batch_scope_level()`/`set_batch_scope_level()`/`batch_consumed_upscope()`/`set_batch_consumed_upscope()` from `eval_expr.hpp`. Delete the two `result.set_batch_*` lines in `eval_expr.cpp` (~604, 606). If `kNoBatchScopeLevel` is now unreferenced, remove it from `node_batch_annotation.hpp`.

- [ ] **Step 3: Verify no dangling references**

Run: `grep -rn "batch_scope_level\|consumed_upscope\|kNoBatchScopeLevel" SeQuant/core SeQuant/domain`
Expected: no output.

- [ ] **Step 4: Full build + suites**

Run: `cmake --build cmake-build-release -j6 --target unit_tests` -> builds clean (`-Werror` clean; no unused-variable fallout).
Run: `./cmake-build-release/tests/unit/unit_tests "[eval][optimize][lifetime_mask]"` -> PASS.
Run: `SEQUANT_UT_DRYRUN_NTERMS=55 ./cmake-build-release/tests/unit/unit_tests "[.][dryrun-occ-veto]"` -> peak still `<= ~443 GB`.

- [ ] **Step 5: Commit**

```bash
git add SeQuant/core/optimize/cost_model.hpp SeQuant/core/eval/node_batch_annotation.hpp \
        SeQuant/core/eval/eval_expr.hpp SeQuant/core/eval/eval_expr.cpp
git commit -m "eval: retire consumed_upscope and scalar batch_scope_level"
```

---

## Task 5: MPQC end-to-end energy validation

Confirm real CSV-CCSD energies are unchanged with the new placement, with aux+occ batching engaged. Not a code change; a validation gate. (Requires the MPQC tree repinned to this SeQuant branch -- do the repin+build only if not already current.)

**Files:** none (validation). Uses `he10-occtest.json`-style forced-batching input (occ_target_size small, aux_target_size small, `peak_threshold` tiny, `objective_function: dense_time_space`) -- recreate the temp input; do NOT commit it.

- [ ] **Step 1: Build the mpqc target** (NOT `MPQCmain`): `cmake --build cmake-build-release -j6 --target mpqc`.
- [ ] **Step 2: Run he10 CSV-CCSD forced-batched**, `MAD_NUM_THREADS=6`, serialized (no parallel MPQC runs). Confirm the trace shows external `BatchScatter` + aux `BatchGroup` events engaged.
- [ ] **Step 3: Check energy** == `-0.33231474200227867` within Delta < 1e-7 (the occ-batching accumulation-order Delta ~6.6e-8 is expected).
- [ ] **Step 4: Clean up** the temp input; record the result in the ledger.

---

## Self-Review

- **Spec coverage:** section 2 (mask/meet) -> Task 1; section 3 + 4 (placement rule + per-level mechanism) -> Task 3; section 5 retire list -> Tasks 2 (veto-b), 3 (hoist), 4 (consumed_upscope + scalar scope_level); section 6 (CSE payoff) -> observable as the C60 `avoidable_time` improvement in Task 3 Step 6; section 8 gate -> Global Constraints + Tasks 3/4/5. Section 7 scope boundary -> Global Constraint (no DP sizing change). Section 9 deferred -> out of plan by design.
- **Placeholder scan:** the meet code in Task 1 Step 4 is a working skeleton with a named helper (`intersect_in_place`) to implement; the per-level placement body in Task 3 is specified by predicate + hook sites + governing gate (the C60 peak and OFF-path byte-identity are the executable spec, as this is a research-grade orchestration change where the body is developed test-first against those gates). No "TBD"/"handle edge cases" left.
- **Type consistency:** `sliced_modes()`/`set_sliced_modes()`/`mask_all_full()` used identically in Tasks 1-3; `stamp_lifetime_masks` signature stable; removed names (`batch_scope_level`, `consumed_upscope`, `kNoBatchScopeLevel`) all gone by Task 4 Step 3's grep gate.
