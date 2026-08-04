# Phase 4b-2 -- hash threading + `remat_to_router` (additive) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the remat placement emittable as `PlacementRouter` overrides. Thread the value hash onto `ValueCell`, then add `remat_to_router` which, for MOVED cells only, emits one override per distinct occurrence key -> the cell's `HomeTarget`. Additive: ZERO production callers; nothing wired into runtime (that is 4b-3).

**Architecture:** A remat cell is one-per-VALUE (`hash_value()`, CSE identity). The router is keyed by the OCCURRENCE KEY (`occurrence_key(node, ctx_modes)` -> `SlotCanonicalizationMetadata`, batched-slot-aware). So one moved value emits MULTIPLE overrides (one per distinct occurrence key) all -> the same `HomeTarget{home_modes}` (the meet home). We thread the cheap value hash to link a `RematResult` cell back to its forest nodes; `remat_to_router` computes the expensive occurrence keys LAZILY, for moved cells only. Design: `doc/dev/specs/2026-08-04-phase4-threshold-o2-cutover-design.md` section 10.

**Tech Stack:** C++20 header, the dry-run backend, Catch2 (`tests/unit/`), CMake + Ninja.

## Global Constraints

- No en-dashes (U+2013) anywhere in source or comments.
- Format touched files with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`.
- No `Co-Authored-By` trailers in commits.
- RELEASE build in `cmake-build-release`; cap ninja at `-j6`.
- ADDITIVE / NON-REGRESSING: `[peak_profile]`, `[placement_remat]`, `[dryrun]`, `[occurrence_key]` stay
  green and their figures UNCHANGED (the new `hash` field + emitter do not perturb sizing / sweep).
  ZERO production callers of `remat_to_router`.

## File Structure

- MODIFY `SeQuant/core/eval/peak_profile.hpp` -- add `std::size_t hash` to `ValueCell`; populate it in
  `compute_dag_boulevard`'s grouping loop (from the already-computed `r.hash` / `hash_to_cell` key).
- MODIFY `SeQuant/core/eval/placement_remat.hpp` -- add `remat_to_router` (+ `#include
  <SeQuant/core/eval/occurrence_key.hpp>`, `placement_router.hpp`).
- MODIFY `tests/unit/test_placement_remat.cpp` -- the `remat_to_router` tests.

---

### Task T1 -- thread the value hash onto `ValueCell`

**Files:**
- Modify: `SeQuant/core/eval/peak_profile.hpp` (`ValueCell` struct ~line 197; the grouping loop ~350-385)
- Test: `tests/unit/test_peak_profile.cpp`

**Interfaces:**
- Consumes: `NodeRec::hash` (already carried in the walk), the `hash_to_cell` grouping.
- Produces: `ValueCell` gains `std::size_t hash;` -- the value's `EvalExpr::hash_value()`.

**Steps:**

- [ ] **Step 1: failing test** (`test_peak_profile.cpp`, `[peak_profile]`): on a small forest, run
  `compute_dag_boulevard`, and assert each `ValueCell::hash` equals the `hash_value()` of a node in that
  value's group (reuse an existing forest; for a leaf-shared value V, both occurrences' nodes have the
  same `hash_value()`, and the cell's `hash` must equal it). Compile fails: no `hash` member.

- [ ] **Step 2: add the member.** In `ValueCell` (`peak_profile.hpp`), add `std::size_t hash;` (document:
  the value's `EvalExpr::hash_value()`, the CSE identity that folds occurrences into this cell; distinct
  from the batched-slot-aware occurrence key the router uses).

- [ ] **Step 3: populate it.** In `compute_dag_boulevard`'s grouping loop, at the cell-CREATING branch
  (where `c.value_id = out.cells.size();` etc.), set `c.hash = r.hash;` (`r.hash` is the `NodeRec`'s
  `hash_value()`, already the `hash_to_cell` key). No change to subsequent-occurrence branch.

- [ ] **Step 4: run.** `ninja -C cmake-build-release unit_tests-sequant -j6`; run `[peak_profile]` +
  `[placement_remat]` -- green, and the anchor 481600 / end-to-end 16128->160 figures UNCHANGED (the new
  field does not affect sizing).

- [ ] **Step 5: commit.** `git commit -m "eval: thread value hash onto ValueCell (Phase 4b-2 T1)"`.

---

### Task T2 -- `remat_to_router` (moved cells -> overrides) + validation

**Files:**
- Modify: `SeQuant/core/eval/placement_remat.hpp`
- Test: `tests/unit/test_placement_remat.cpp`

**Interfaces:**
- Consumes: `ValueCell::{hash, home_modes, value_id}` (T1); `occurrence_key(node, ctx_modes)`
  (`occurrence_key.hpp`, `ctx_modes` is `container::svector<Index> const&`) -> `SlotCanonicalizationMetadata`;
  `PlacementRouter<TreeNode>::set_override(Key, HomeTarget)` with `HomeTarget{container::svector<Index>
  residency; std::size_t split_index;}` (`placement_router.hpp`); the post-order + ectx-accumulation walk
  shape from `compute_dag_boulevard`.
- Produces:
  ```cpp
  // Emit PlacementRouter overrides for the MOVED cells of a remat result: for each cell whose
  // home_modes differ from its seed home, one override per distinct occurrence key of that value ->
  // HomeTarget{final home_modes, 0}. Unmoved cells emit nothing.
  template <meta::eval_node_range R>
  PlacementRouter</*TreeNode of R*/> remat_to_router(
      container::svector<ValueCell> const& seed_cells,   // o2_cells(forest).cells (home = seed)
      container::svector<ValueCell> const& remat_cells,  // rematerialize_to_budget(...).cells (final)
      R const& forest);
  ```
  (Match the exact `PlacementRouter` template parameter to the forest's `TreeNode`, as `CacheManager` /
  the Phase-2 seam do; grep `PlacementRouter<` for the spelling.)

**Steps:**

- [ ] **Step 1: failing tests** (`test_placement_remat.cpp`, `[placement_remat]`):
  - **Moved cell emits its override.** Build the demoted-giant forest (reuse the T3 end-to-end forest),
    `o2_cells` -> `rematerialize_to_budget(threshold small enough to shrink the giant)` -> the giant cell
    moved. `remat_to_router(seed, result, forest)`: assert the router is NON-empty, and that looking up
    the giant's occurrence key (compute `occurrence_key(giant_node, ectx_modes)` for one of its
    occurrences) returns a `HomeTarget` whose `residency` == the shrunk `home_modes`.
  - **Unmoved cells emit nothing.** In the same result, assert a non-moved cell's occurrence key is
    absent from the router (no override).
  - **inf threshold => empty router.** `rematerialize_to_budget(..., +inf)` (no moves) -> `remat_to_router`
    returns an EMPTY router (`router.empty()`), i.e. the Phase-2 inert seed seam.
  - **Two keys, one home.** If practical on the constructed forest: a moved value with two occurrences of
    DIFFERENT occurrence keys emits TWO overrides -> the SAME `HomeTarget`. If the forest cannot produce
    two distinct keys for one value, assert the single-key case and note it.

- [ ] **Step 2: implement `remat_to_router`.**
  ```
  1. moved: unordered_map<size_t /*hash*/, svector<Index> /*final home_modes*/>.
     for i in [0, remat_cells.size()): if remat_cells[i].home_modes != seed_cells[i].home_modes
       moved[remat_cells[i].hash] = remat_cells[i].home_modes;   // (same value_id order in both)
  2. PlacementRouter<TreeNode> router;
  3. post-order walk the forest maintaining the ectx stack (push each node's batched_here loops,
     proto-expanded, on descent; pop on ascent) -- the SAME accumulation compute_dag_boulevard uses.
     For each node n:
       if moved.contains(n->hash_value()):
         ctx_modes = [e.first for e in ectx]   // svector<Index>
         key = occurrence_key(n, ctx_modes)
         router.set_override(key, HomeTarget{moved[n->hash_value()], 0})
  4. return router.
  ```
  (`set_override` is `insert_or_assign`, so repeat keys across occurrences coalesce harmlessly.)
  Add the `#include`s (`occurrence_key.hpp`, `placement_router.hpp`), project-headers-first.

- [ ] **Step 3: run.** `[placement_remat]` all green (incl. the new cases); `[peak_profile]`,
  `[occurrence_key]`, `[dryrun]` unchanged/green.

- [ ] **Step 4: commit.** `git commit -m "eval: remat_to_router emits overrides for moved cells (Phase 4b-2 T2)"`.

---

## Self-review checklist

- Coverage: design 10.2 step 1 -> T1; 10.2 step 2 -> T2; 10.3 validation -> T2 step 1 (moved/unmoved/
  inf-empty/two-keys). All covered.
- Types: `ValueCell::hash` (T1) consumed by `remat_to_router` (T2); `HomeTarget`/`occurrence_key`/
  `set_override` spellings verified against `placement_router.hpp` / `occurrence_key.hpp`.
- Additive: no runtime source touched; `remat_to_router` has zero production callers; existing figures
  unchanged.

## Deferred (NOT this plan)

- **Phase 4b-3** -- wire a dry-run remat pre-pass -> populate the router (via `remat_to_router`) ->
  `place_at_this_level` consults it; DELETE `has_demoted_external`; the default-threshold policy. Needs
  its own runtime-integration design pass. MPQC CCk (cross-repo, user-run) is its final energy gate.
