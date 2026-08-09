# CSE-aware rematerialization split -- Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the router-vs-hoist correctness hole by making one
slicing-signature criterion govern the hoist path, the router overlay, and
remat's pricing, so a CSE-folded value whose occurrences slice a relabeled mode
is materialized per occurrence (split) instead of sharing one wrong slice.

**Architecture:** Introduce a single shared slicing-signature helper and use it
in both runtime materialization paths (`make_batched_scratch`, already correct;
`place_at_this_level`, the hole). Make the occurrence key DAG-global (space +
batched-salt coloring) and the router home a label-free DAG-scope resolved to a
physical loop per use. remat emits one DAG-global overlay and prices the split.
All changes are behavior-neutral when the router is empty (every SeQuant unit
test path); the new behavior fires only under a populated remat router.

**Tech Stack:** C++20, SeQuant `core/eval` (`eval.hpp`, `cache_manager.hpp`,
`occurrence_key.hpp`, `placement_router.hpp`, `placement_remat.hpp`,
`peak_profile.hpp`, `lifetime_mask.hpp`); Catch2 unit tests under `tests/unit/`.

**Design spec:** `doc/dev/specs/2026-08-07-remat-cse-aware-split-design.md`
(read the Glossary and sections "The model", "The verified correctness hole",
and "Design" before starting).

## Global Constraints

- **Empty-router byte-identity:** every change must be a no-op when
  `cache.placement_router()` is null or `empty()`. All existing SeQuant unit
  tests run with no router and must stay byte-identical.
- **No en-dashes (U+2013) in committed files** -- the `forbid-en-dashes`
  pre-commit hook rejects them. Use ASCII `-`.
- **No AI-attribution trailers** in commit messages (`git log` on master has
  none; match it). No `Co-Authored-By`.
- **Build target:** `cmake --build cmake-build-release --target unit_tests-sequant -j6`
  (cap `-j` to avoid memory pressure). Run targeted filters, not the whole
  suite (an unrelated slow DP test in `test_eval_dryrun.cpp` hangs full runs).
- **New unit-test files** must be added to `tests/unit/CMakeLists.txt`.

---

### Task 1: Shared slicing-signature helper

Extract the "can these occurrences share one materialization?" criterion that
`make_batched_scratch` already implements (`eval.hpp` ~1443-1460) into a single
reusable helper, so the hoist path (Task 4), the router read (Task 6), and remat
(Task 5) all call one implementation instead of re-deriving it.

**Files:**
- Create: `SeQuant/core/eval/slicing_signature.hpp`
- Test: `tests/unit/test_slicing_signature.cpp`
- Modify: `tests/unit/CMakeLists.txt` (add the new test source)

**Interfaces:**
- Produces:
  - `sequant::eval::slicing_signature(TreeNode const& node, Index const& mode) -> std::optional<std::size_t>`
    -- the position of `mode` in `node->canon_indices()` (proto-expanded on the
    slot side, exactly as `stamp_lifetime_masks`/`make_batched_scratch` do), or
    `std::nullopt` if `mode` does not appear on this node's result.
  - `sequant::eval::signatures_consistent(range-of TreeNode occurrences, Index const& mode) -> bool`
    -- true iff `slicing_signature(occ, mode)` is equal across all `occurrences`
    (all present at the same position, or all absent).

- [ ] **Step 1: Read the current criterion.** Read `eval.hpp`
  `make_batched_scratch` (the doc comment ~1435-1460 and the signature-building
  loop below it) and note exactly how it derives "position of the batch mode in
  `canon_indices()`, or its absence." Copy that derivation faithfully; do not
  invent a new one.

- [ ] **Step 2: Write the failing test.**

```cpp
// tests/unit/test_slicing_signature.cpp
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <catch2/catch_test_macros.hpp>
// ... build two EvalExpr nodes that are canonically equal (same hash) but whose
// third-occ slot is physically i_3 in one and i_4 in the other (mirror the
// g.C legs; reuse the term-4 fixture helpers in test_eval_dryrun.cpp if handy).
TEST_CASE("slicing_signature: relabeled legs are inconsistent", "[slicing-signature]") {
  // A carries i_3 at canonical pos p, not i_4; B carries i_4 at pos p, not i_3.
  REQUIRE(sequant::eval::slicing_signature(A, i3).has_value());
  REQUIRE_FALSE(sequant::eval::slicing_signature(A, i4).has_value());
  REQUIRE_FALSE(sequant::eval::signatures_consistent(std::array{A, B}, i3));
  REQUIRE_FALSE(sequant::eval::signatures_consistent(std::array{A, B}, i4));
}
TEST_CASE("slicing_signature: shared-label mode is consistent", "[slicing-signature]") {
  // Both A and B carry i_1 (the amplitude pair index) at the same canonical pos.
  REQUIRE(sequant::eval::signatures_consistent(std::array{A, B}, i1));
}
```

- [ ] **Step 3: Run it to confirm it fails** (header does not exist yet):
  `cmake --build cmake-build-release --target unit_tests-sequant -j6` then
  `./cmake-build-release/tests/unit/unit_tests-sequant "[slicing-signature]"`.
  Expected: compile error / missing header.

- [ ] **Step 4: Implement `slicing_signature.hpp`** by lifting the exact
  derivation read in Step 1 into the two free functions above. Include only what
  it needs (`eval_expr.hpp`, `index.hpp`, `<optional>`). Keep it header-only
  (matches the `core/eval` style).

- [ ] **Step 5: Refactor `make_batched_scratch` to call the helper.** Replace
  the inline signature derivation in `eval.hpp` with a call to
  `slicing_signature`/`signatures_consistent`. This must be behavior-preserving:
  the registration decision is unchanged.

- [ ] **Step 6: Run tests.** `[slicing-signature]` passes; then run a batched
  eval smoke filter (e.g. `[eval]` or the specific batched-eval Catch2 tags) to
  confirm the `make_batched_scratch` refactor is byte-neutral.

- [ ] **Step 7: Commit.**

```bash
git add SeQuant/core/eval/slicing_signature.hpp tests/unit/test_slicing_signature.cpp tests/unit/CMakeLists.txt SeQuant/core/eval/eval.hpp
git commit -m "eval: extract shared slicing-signature helper"
```

---

### Task 2: DAG-global occurrence key (space + batched salt)

Change `occurrence_key` so batched slots are colored by space + a single
batched salt (renamable within the class) instead of space + label, so two
occurrences of one value bound to different externals map to the same key.

**Files:**
- Modify: `SeQuant/core/eval/occurrence_key.hpp` (`occurrence_key`, ~90-106)
- Test: `tests/unit/test_placement_router.cpp` (existing; add cases)

**Interfaces:**
- Produces: `occurrence_key(node, ctx_modes)` returning a
  `TensorNetwork::SlotCanonicalizationMetadata` that is **equal** for two
  occurrences differing only in the physical label bound to a batched slot.

- [ ] **Step 1: Read the coloring API.** Read `occurrence_key.hpp` 68-106 and the
  `TensorNetwork::canonicalize_slots` signature + how `named_indices` are colored
  (space + label) vs non-named (space alone) in `tensor_network.*`. Determine the
  supported way to request "color = space + a fixed salt, renamable within the
  class" for the batched slots -- either an existing coloring option or a minimal
  extension of `canonicalize_slots`'s named-index coloring. Record the chosen
  mechanism in the task's report before writing code.

- [ ] **Step 2: Write the failing test** (in `test_placement_router.cpp`):

```cpp
TEST_CASE("occurrence_key is DAG-global over batched labels", "[placement_router]") {
  // A: value with third-occ slot bound to i_3 (batched); B: same value, i_4.
  auto ka = sequant::eval::occurrence_key(A, /*ctx_modes with i_3 in scope*/);
  auto kb = sequant::eval::occurrence_key(B, /*ctx_modes with i_4 in scope*/);
  RouterKeyEqual eq;
  REQUIRE(eq(ka, kb));                 // same batched-slot structure -> same key
  // and a batched-occ slot differs from a non-batched-occ slot:
  auto kc = sequant::eval::occurrence_key(C_nonbatched_occ, /*no i in scope*/);
  REQUIRE_FALSE(eq(ka, kc));
}
```

- [ ] **Step 3: Run it to confirm it fails** (current space+label coloring makes
  `ka != kb`). `unit_tests-sequant "[placement_router]"`.

- [ ] **Step 4: Implement the coloring change** using the mechanism chosen in
  Step 1: color batched slots by space + a single shared batched salt (same salt
  for every space; the per-space color separates batched-occ from batched-aux),
  leaving non-batched slots colored by space alone. Update the doc comment
  (73-80) to describe space+salt.

- [ ] **Step 5: Run tests.** `[placement_router]` passes (new + existing). If an
  existing router test asserted per-label distinctness, update it to the DAG-
  global behavior and note why in the commit.

- [ ] **Step 6: Commit.**

```bash
git add SeQuant/core/eval/occurrence_key.hpp tests/unit/test_placement_router.cpp
git commit -m "eval: color batched occurrence-key slots by space+salt (DAG-global)"
```

---

### Task 3: DAG-scope home + per-use resolution in the router

Replace the router `HomeTarget`'s tree-local residency with a label-free
DAG-scope (ordered sequence of spaces), and make `home_depth` resolve it to a
physical loop via the use site's occurrence key.

**Files:**
- Modify: `SeQuant/core/eval/placement_router.hpp` (`HomeTarget`, `home_depth`)
- Modify: `SeQuant/core/eval/eval.hpp` (Enter-stage router consult, 732-762;
  `place_at_this_level` router branch, 1839-1858 -- update the `home_depth` call
  sites to pass the occurrence key / node)
- Test: `tests/unit/test_placement_router.cpp`

**Interfaces:**
- Produces:
  - `HomeTarget { container::svector<SpaceTag> dag_scope; }`
    where `SpaceTag` is the batch-loop space identifier (e.g. the `IndexSpace`
    base key / a small enum). `dag_scope` is an ordered nest prefix, and it
    alone identifies the home (split cells of one value differ in `dag_scope`).
  - `home_depth(HomeTarget const& home, BatchContext const& ctx, Key const& key) -> std::size_t`
    -- resolve `home.dag_scope` against `ctx` by mapping each DAG-scope position
    (a space) to the physical loop the occurrence binds there (via `key`'s slot
    canonicalization), returning `1 + deepest matched ctx level`, else 0.
- Consumes: `occurrence_key` (Task 2) at the resolution site.

- [ ] **Step 1: Read the slot mapping.** Read `SlotCanonicalizationMetadata`
  (in `tensor_network.*`) to find how it exposes the permutation from physical
  slots to canonical positions -- the data `home_depth` needs to map a DAG-scope
  position to the occurrence's physical index/loop. Record the exact accessor.

- [ ] **Step 2: Write the failing test.** Two occurrences A (`i_3`) and B (`i_4`)
  of one value, one overlay with `dag_scope = [occ]`, a `BatchContext` with both
  `i_3` and `i_4` live (nested):

```cpp
TEST_CASE("home_depth resolves one DAG-scope to per-use physical depth", "[placement_router]") {
  HomeTarget h; h.dag_scope = { SpaceTag::occ };
  // ctx = [ (i_3, ...), (i_4, ...) ]  (i_3 outer, i_4 inner)
  REQUIRE(router.home_depth(h, ctx, key_of_A) == 1);   // A's occ slot -> i_3 -> depth 1
  REQUIRE(router.home_depth(h, ctx, key_of_B) == 2);   // B's occ slot -> i_4 -> depth 2
}
```

- [ ] **Step 3: Run it to confirm it fails** (current `home_depth` takes no key
  and matches a canonical `Index` against `ctx` -- both A and B resolve the same).

- [ ] **Step 4: Implement `HomeTarget`/`home_depth`** per the interface, using the
  accessor from Step 1 to translate DAG-scope position -> physical loop -> ctx
  depth. Update the two `home_depth` call sites in `eval.hpp` (Enter-stage 740,
  hoist branch 1846) to pass the occurrence key already computed there.

- [ ] **Step 5: Run tests.** `[placement_router]` passes; empty-router eval smoke
  filter stays byte-identical (the branch is gated on a non-empty router).

- [ ] **Step 6: Commit.**

```bash
git add SeQuant/core/eval/placement_router.hpp SeQuant/core/eval/eval.hpp tests/unit/test_placement_router.cpp
git commit -m "eval: router home is a DAG-scope resolved per use to physical depth"
```

---

### Task 4: Per-occurrence hoist on signature divergence (the correctness fix)

Teach `place_at_this_level` to not dedup/share occurrences that are
signature-inconsistent at the home scope: register and build per occurrence
instead (the split), reusing the Task 1 helper.

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (`place_at_this_level`, 1796-1896: the
  `collect` dedup 1801-1802, `ensure_hoist_slot`/`alive` 1883-1884, build/store
  1892-1895)
- Test: `tests/unit/test_eval.cpp` (or the batched-eval test file) -- a numeric
  correctness test

**Interfaces:**
- Consumes: `signatures_consistent` (Task 1), the per-use DAG-scope resolution
  (Task 3). No cache re-keying: divergent occurrences resolve to DISTINCT
  scope-chain depths (distinct `cache_map_`s), so the canonical hash key already
  holds both materializations -- see spec Non-goals.

- [ ] **Step 1: Write the failing correctness test.** Build a minimal batched
  forest whose value has two occurrences slicing a relabeled mode (the g.C
  shape), with a router homing that value at the divergent sub-scope, and assert
  the numeric result equals the un-batched reference. Without the fix the two
  legs share one slice and the result is wrong.

```cpp
TEST_CASE("hoist splits a divergently-sliced CSE value", "[eval][split]") {
  // reference = evaluate(forest, /*no router*/)  -> correct
  // under a router homing g.C at the divergent sub-scope:
  auto got = evaluate(forest_with_router);
  REQUIRE(norm(got - reference) < 1e-12);   // per-occurrence build -> exact
}
```

- [ ] **Step 2: Run it to confirm it fails** (shared slice -> wrong number).

- [ ] **Step 3: Implement the split.** In `collect`, when a candidate occurrence
  is signature-inconsistent (`!signatures_consistent(...)`) with an existing
  target at the home scope, do not fold it via `eq`; register a distinct hoist
  target for it, built at its OWN router-resolved depth (a distinct scope-chain
  level, hence a distinct `cache_map_` entry -- no extra discriminator), and in
  the build loop skip the `alive(d)` early-out for that occurrence so it is built
  per occurrence. Keep the
  consistent case exactly as today.

- [ ] **Step 4: Run tests.** `[eval][split]` passes; the full batched-eval smoke
  filter stays green and byte-identical on the empty-router path.

- [ ] **Step 5: Commit.**

```bash
git add SeQuant/core/eval/eval.hpp SeQuant/core/eval/cache_manager.hpp tests/unit/test_eval.cpp
git commit -m "eval: split divergently-sliced CSE occurrences in the hoist path"
```

---

### Task 5: remat prices the split as two cells with a replication factor

Un-fold a divergent value into two real cells and price each cell's recompute by
its **replication factor** (the product of the block counts of every enclosing
loop it is homed-within but does not carry), replacing the flat 2x fudge that
`46b495eba` shipped as a placeholder. Make `shrink_candidates`/`apply_shrink`
CSE-aware (relabeled -> split-only), and make `remat_to_router` moved-detection
hash-aware so split cells (new `value_id`s, one shared hash) are still emitted as
one overlay.

> **HARD DEPENDENCY on Task 3.** The replication factor is a product over the
> levels a cell is *homed-within but does not carry* (leg B's `i_3`) -- exactly the
> levels a value-relative coordinate (`slot_positions`) cannot name, because they
> are not in that occurrence's `canon_indices`. This task **cannot be priced
> correctly until Task 3's DAG-scope home lands**; do Task 3 first. See spec
> Design section 3 ("remat prices the split") and "Scope must be a DAG-level
> coordinate".

**Files:**
- Modify: `SeQuant/core/eval/peak_profile.hpp` -- `ValueCell` (retain per-
  occurrence records `{carried, ectx, point, consumer_point, home}`, not just the
  first + union/min/max); `compute_dag_boulevard` grouping (stop collapsing them);
  `cell_footprint` (**delete** the `divergent ? 2*base : base` branch -- the sweep
  prices peak co-residency by liveness overlap once there are two real cells).
- Modify: `SeQuant/core/eval/placement_remat.hpp` -- `shrink_candidates`/
  `apply_shrink` (CSE-aware; `apply_shrink` grows the cell container instead of
  mutating one cell); `rematerialize_to_budget` (one->two growth in trial +
  commit, fresh `value_id`s, objective weighs DeltaPeak vs split recompute, restate
  termination); `remat_to_router` (hash-aware moved-detection).
- Test: `tests/unit/test_placement_remat.cpp`

**Interfaces:**
- Consumes: `slicing_signature` (Task 1); the DAG-scope `HomeTarget` and its
  per-use resolution (Task 3) -- **required**, see the dependency note above.
- Produces: an `apply_split` (replacing/extending `apply_shrink`) that replaces one
  divergent `ValueCell` with two same-hash cells; a recompute-priced `RematResult`;
  a hash-aware `remat_to_router` still emitting one overlay per value.

- [ ] **Step 1: Write the failing pricing test.** A `ValueCell` whose two
  occurrences bind a relabeled mode to `i_3` (outer, N3 blocks) and `i_4` (inner):
  assert `shrink_candidates` classifies that mode split-only; that applying the
  split yields two same-hash cells at distinct DAG-scope homes; and that the priced
  recompute matches the **replication factor**, not 2x -- cell A (homes `i_3`,
  carries it) costs 1x full build, cell B (homes `i_4`, invariant to `i_3`) costs
  **N3 x** full build, so total is `(1 + N3) x`, equal to 2x only at N3 = 1.

```cpp
TEST_CASE("remat splits along a relabeled mode, prices replication not 2x", "[placement_remat]") {
  auto cands = sequant::eval::shrink_candidates(cell);
  REQUIRE(is_split_only(cands, relabeled_mode));       // not an in-place shrink
  auto res = sequant::eval::apply_split(cells, ci, relabeled_mode);
  REQUIRE(res.cells_of_hash(cell.hash).size() == 2);   // two cells, one value
  // N3 = outer (i_3) block count; the inner cell is rebuilt per i_3 block:
  REQUIRE(res.modeled_recompute == Approx((1 + N3) * full_build));  // NOT 2 * full_build
}
```

- [ ] **Step 2: Run it to confirm it fails** (today: one cell, priced 2x).

- [ ] **Step 3: Retain per-occurrence records on the cell.** In
  `compute_dag_boulevard`, stop collapsing each occurrence's
  `{carried (canon_indices), ectx, point, consumer_point, home}` at grouping; keep
  a per-occurrence list on `ValueCell` sufficient to (a) classify a mode
  shared-label vs relabeled via `slicing_signature`, and (b) compute subset-local
  liveness and home for each split cell.

- [ ] **Step 4: Make `shrink_candidates`/`apply_shrink` CSE-aware + split.**
  Classify each candidate: shared-label (in `carried_isect`) -> in-place shrink
  (today's behavior). Relabeled (in `divergent_modes`) -> **split-only**: replace
  the one `ValueCell` in the working vector with **two same-hash cells**, one per
  occurrence subset, each with its subset-local `carried`/`home_modes`/liveness,
  `divergent_modes = {}`, and a fresh `value_id`, with the relabeled mode folded
  into each subset's home under its own physical label. `apply_shrink`'s signature
  changes from mutate-one-cell to grow-the-container.

- [ ] **Step 5: Price recompute by replication factor; delete the 2x.** Remove the
  `divergent ? 2*base : base` branch from `cell_footprint`. Add, to the cost each
  split cell contributes, a recompute term = `build_cost(cell) x PROD over every
  enclosing loop the cell is homed-within-but-does-not-carry of (that loop's block
  count)`. "Homed-within-but-not-carry" is read from the cell's DAG-scope home
  (Task 3) minus its `carried` set; the enclosing loops come from `enclosing_modes`
  / the nest. Peak co-residency is left to `peak_profile_sweep` (two liveness
  intervals). Confirm on the canonical fixture that the term is `(1 + N3) x`.

- [ ] **Step 6: Teach `rematerialize_to_budget` the one->two move.** A split erases
  one cell and appends two, so trial (`auto trial = cells; ...`) and commit must
  use `apply_split` and mint fresh `value_id`s. The split breaks the "all shrinks
  are zero-recompute, so the metric is pure `DeltaPeak`" invariant, so **record
  each split's recompute** (the replication factor) for the forecast. Keeping the
  selection objective **peak-only** is the current, correct default -- do NOT add
  recompute to the objective in this task; making recompute a *secondary* objective
  (prefer, among peak-feasible placements, the lowest-recompute one) is a separate,
  deliberate decision (see spec section 3, "Does this recompute term change the
  schedule?"). Restate termination: each split strictly reduces the forest's
  divergent-mode count (bounded); shrinks still grow `home_modes` toward `carried`
  (bounded).

- [ ] **Step 7: Make `remat_to_router` moved-detection hash-aware,** and still emit
  **one DAG-global overlay per value hash** (Task 2 key, Task 3 DAG-scope). Split
  cells carry fresh `value_id`s absent from the seed-home map, so `value_id`-keyed
  moved-detection would drop them -- key moved-detection by hash instead. Both split
  cells resolve to the *same* DAG-scope home, so assert they agree and emit once.

- [ ] **Step 8: Run tests.** `[placement_remat]` passes; the existing remat tests
  stay green (unmoved cells unchanged).

- [ ] **Step 9: Commit.**

```bash
git add SeQuant/core/eval/placement_remat.hpp SeQuant/core/eval/peak_profile.hpp tests/unit/test_placement_remat.cpp
git commit -m "remat: split along relabeled modes, priced by replication factor"
```

---

### Task 6: Defense-in-depth signature check on the router read

Have the router-directed read verify the fetched entry's slicing signature
matches the occurrence's before using it, recomputing on mismatch -- so no future
overlay can silently reintroduce the hole.

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (Enter-stage router consult, 732-762:
  after `access_at_hops` returns a hit, verify before `finalize`)
- Test: `tests/unit/test_eval.cpp`

**Interfaces:**
- Consumes: `slicing_signature` (Task 1).

- [ ] **Step 1: Write the failing test.** Construct a deliberately mis-emitted
  shared overlay (homes a divergent value at a single shared scope). Assert the
  router read detects the signature mismatch and recomputes (correct result)
  rather than serving the wrong slice.

```cpp
TEST_CASE("router read rejects a signature-mismatched share", "[eval][router-guard]") {
  auto got = evaluate(forest_with_bad_shared_overlay);
  REQUIRE(norm(got - reference) < 1e-12);   // guard forces recompute -> exact
}
```

- [ ] **Step 2: Run it to confirm it fails** (without the guard, the bad overlay
  corrupts).

- [ ] **Step 3: Implement the guard.** In the Enter-stage router hit branch
  (after line 743), compute the fetched entry's slicing signature for the
  in-scope batch mode and compare to the current occurrence's; if inconsistent,
  do not `finalize` the routed value -- fall through to the normal
  `access_at`/recompute path (set `routed = false`).

- [ ] **Step 4: Run tests.** `[eval][router-guard]` passes; empty-router path
  byte-identical.

- [ ] **Step 5: Commit.**

```bash
git add SeQuant/core/eval/eval.hpp tests/unit/test_eval.cpp
git commit -m "eval: router read verifies slicing signature before sharing"
```

---

### Task 7: Wet end-to-end validation (MPQC)

Prove the fix on a real system: drive w8 PNO-CCSD through a schedule that moves
the divergent value to a sub-scope and assert the CCk energy matches the
reference to lossless tolerance. This uses the MPQC `MPQC_SCHED_SIZE_OVERRIDE`
diagnostic already prototyped this session (uncommitted in `cck.ipp`/`sequant.h`).

**Files:**
- Reference: MPQC `src/mpqc/chemistry/qc/lcao/cc/cck.ipp` (schedule-size override),
  `repro/w8-occ8.json`, `repro/w8_rootseed.log`
- No new SeQuant source; this is an integration check.

**Interfaces:** consumes the full stack (Tasks 1-6) via a repinned SeQuant.

- [ ] **Step 1: Repin + rebuild.** Repin MPQC's SeQuant to the branch HEAD after
  Tasks 1-6; `cmake --build cmake-build-release --target mpqc -j6` (target is
  `mpqc`, not `MPQCmain`); confirm the binary mtime updated.

- [ ] **Step 2: Reference.** Confirm the trustworthy w8 PNO-CCSD reference:
  `CSV-CCk Energy = -2.0005382446100` (w8_rootseed.log:518; reproduce with a
  no-batch or own-schedule run if desired).

- [ ] **Step 3: Run the forced-split case.**

```bash
cd repro
MAD_NUM_THREADS=1 MPQC_SCHED_SIZE_OVERRIDE=c60 MPQC_SCHED_PEAK_THR_GB=500 \
  ../cmake-build-release/src/bin/mpqc/mpqc w8-occ8.json 2>&1 | tee w8-c60-split.out
```

  Confirm (via the remat diagnostic) that the schedule actually moves a divergent
  value to a sub-scope (`router_empty=0`, and the moved list includes a value
  whose occurrences slice a relabeled mode). If the g.C is not moved at these
  parameters, adjust `peak_threshold`/`occ_target_size` until a divergent value
  is moved -- the test is only meaningful when the split path fires.

- [ ] **Step 4: Assert lossless.** The converged `CSV-CCk Energy` matches
  `-2.0005382446100` to ~1e-8 (lossless batching tolerance). Pre-fix, this run
  diverges (the corruption); post-fix it must track.

- [ ] **Step 5: Record the result** in the spec's evidence section (or a short
  findings note) with the two energies and the moved-value confirmation. Do not
  commit the MPQC diagnostic scaffolding unless the team decides to keep it; it
  is a test harness, not a shipping feature.

---

## Self-review notes

- **Spec coverage:** occurrence-key coloring (spec Design #0) -> Task 2;
  hoist-path signature unification (#1) -> Tasks 1+4; one DAG-global overlay +
  per-use resolution (#2) -> Tasks 3+5; CSE-aware pricing/split (#3) -> Task 5;
  defense-in-depth router read (#4) -> Task 6; wet correctness + peak
  monotonicity (spec Testing 4-5) -> Task 7. Occurrence-key DAG-globality (spec
  Testing 1) -> Task 2; signature criterion (Testing 2) -> Tasks 1+4; router
  no-op invariance (Testing 3) -> the empty-router constraint checked in every
  task; defense-in-depth (Testing 6) -> Task 6.
- **Ordering:** Tasks 1-3 are behavior-neutral infrastructure (no-op under empty
  router); Task 4 is the correctness fix; Task 5 makes remat drive/price it;
  Task 6 is belt-and-suspenders; Task 7 validates end-to-end. **Task 5 has a HARD
  dependency on Task 3:** the split's recompute is a replication factor over the
  levels a cell is homed-within-but-does-not-carry, which only the Task 3 DAG-scope
  home can name -- Task 5 cannot be priced on Task 3's interim `slot_positions`.
- **Note on the current tree:** Task 3 shipped as `slot_positions`/`free_modes`
  (works intra-tree via a raw-Index `free_modes` patch) and Task 5 shipped only the
  2x pricing fudge (`46b495eba`) -- both are placeholders to be replaced by the
  DAG-scope home and the replication-factor split model described here.
- **Two API-read-first steps** (Task 2 Step 1: `canonicalize_slots` coloring;
  Task 3 Step 1: `SlotCanonicalizationMetadata` slot mapping) are load-bearing --
  the exact call shape must be read from `tensor_network.*` before coding, and
  the implementer records the chosen mechanism in the task report.
- **Peak monotonicity** (spec Testing 5) is exercised by Task 7's parameter
  sweep and by re-running the c60 dry-run `[c60-term-dump]` param scan after the
  fix; add an explicit assertion there if a committed regression guard is wanted.
