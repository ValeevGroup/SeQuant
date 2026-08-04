# Phase 4b-3 -- runtime cutover: dry-run remat pre-pass -> router -> place_at_this_level; delete the veto (cross-repo) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the remat placement drive runtime eval: an MPQC pre-pass runs remat on the real `enodes` (when `peak_threshold` is finite) and populates the `PlacementRouter` the runtime already consults; DELETE the demotion veto `has_demoted_external`. One coupled phase across SeQuant + MPQC, gated by MPQC CCk energy correctness.

**Architecture:** `place_at_this_level` ALREADY consults the router (`eval.hpp:1784-1790`, Phase 2 store seam) -- an empty router is inert (seed placement). 4b-3 (a) deletes the veto in SeQuant and (b) populates the router in MPQC's `build_cache_manager` via `remat_cells -> rematerialize_to_budget -> remat_to_router` on the real `EvalExprTA` forest with a dry-run `SizeRegime`. `threshold = inf` => no pre-pass, empty router, no batching (Section 5 of the design doc: the veto is inert without batching, so its deletion is safe by default). Design: `doc/dev/specs/2026-08-04-phase4-threshold-o2-cutover-design.md` section 11.

**Tech Stack:** C++20; SeQuant eval headers; MPQC lcao/expression + cc; Catch2 (SeQuant unit); MPQC ctest validation.

## Global Constraints

- No en-dashes (U+2013) in source or comments. clang-format touched files
  (`/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`). No `Co-Authored-By`. Project headers first.
- SeQuant: RELEASE build `cmake-build-release -j6`. MPQC: build the `mpqc` target (NOT `MPQCmain`) in
  `/Users/efv/code/mpqc4/cmake-build-release` (consumes local SeQuant via FETCHCONTENT_SOURCE_DIR_SEQUANT).
- **CORRECTNESS IS THE GATE.** Placement changes; CCk energies must NOT. The final gate is the cross-repo
  MPQC CCk + CSV-CCk validation (np1+np2), USER-run. If any energy moves it is a real regression -- STOP,
  do not re-baseline energies.
- Default `peak_threshold` stays `inf` (no memory-fraction default; arrays may be disk-backed). Batching
  requires an explicit finite budget (already validated by `validate_batch_config`).
- SeQuant commits land first (MPQC consumes them); repin note: this is a cross-repo change, coordinate the
  SeQuant/MPQC branches (both `evaleev/feature/multimode-batched-eval`).

## File Structure

- SeQuant `SeQuant/core/eval/eval.hpp` -- DELETE `has_demoted_external` (T1).
- SeQuant `tests/unit/test_eval_ta.cpp` -- re-baseline the veto cases (T1).
- SeQuant `tests/unit/test_placement_remat.cpp` (or a new test) -- the router-round-trip correctness
  invariant (T1).
- SeQuant `SeQuant/core/eval/cache_manager.hpp` -- OPTIONAL owned-router support (T2, if chosen for
  ownership).
- MPQC `src/mpqc/chemistry/qc/lcao/expression/sequant.h` (`build_cache_manager`) + `sequant.cpp` (shared
  dry-run `SizeRegime` construction) -- the pre-pass wiring (T2).

---

### Task T1 -- SeQuant: delete `has_demoted_external` + the router round-trip invariant test

**Files:**
- Modify: `SeQuant/core/eval/eval.hpp` (`has_demoted_external`: lambda ~1733, collect-gate conjunct ~1745, doc ~1780 -- grep-verify current lines)
- Test: `tests/unit/test_eval_ta.cpp` (re-baseline the veto cases), `tests/unit/test_placement_remat.cpp` (the invariant)

**Interfaces:**
- Consumes: `place_at_this_level`'s existing router store seam (`eval.hpp:1784-1790`), `residency_all_outer`, the `rl` walk.
- Produces: a veto-free `place_at_this_level` = router override if present, else deepest-`sliced_modes`-loop (seed).

**Steps:**

- [ ] **Step 1: grep-verify + read** the `has_demoted_external` lambda and its single use (the `collect`
  gate conjunct `&& !has_demoted_external(n)`), plus the store-seam router consult (`eval.hpp:1784-1790`).
  Confirm the veto is a SEPARATE conjunct removable without disturbing `residency_all_outer` / the rl walk.

- [ ] **Step 2: the correctness-invariant test FIRST** (`test_placement_remat.cpp`, `[placement_remat]`):
  build a forest with a moved value (reuse the CONSISTENT distinct-key fixture:
  `V{i_1;i_2}=A{i_1;x_1}*B{x_1;i_2}`, sliced `{i_1,i_2}` vs `{i_1}`), `remat_cells` ->
  `rematerialize_to_budget(threshold that shrinks V)` -> `remat_to_router`, then assert the runtime
  `place_at_this_level` would resolve V's home from the router: for a moved occurrence, look up
  `occurrence_key(V_node, ectx_modes)` in the router and assert `home_depth(HomeTarget, ectx)` equals the
  depth of V's remat home in that ectx. This pins that the pre-pass's keys == the runtime store seam's
  keys (design 11.4). (If the store-seam key path is not directly callable from a unit test, assert the
  equivalent via `router.route(occurrence_key(...))` + `router.home_depth(...)`, which is the same code
  the store seam runs.)

- [ ] **Step 3: delete `has_demoted_external`.** Remove the lambda and its `&& !has_demoted_external(n)`
  conjunct in the `collect` gate; remove/adjust the `:1780` doc reference. `place_at_this_level` now
  collects (hoists to the seed home) every node `residency_all_outer` accepts, unless the router overrides
  -- no veto.

- [ ] **Step 4: re-baseline the veto tests.** The `test_eval_ta.cpp` cases that asserted the veto homing a
  demoted giant LOCAL now assert the seed (hoisted) placement (with an empty router) -- update their
  expectations with a one-line justification per case (the veto is gone; remat, not the veto, now
  constrains the giant, and these unit tests have no remat pre-pass so they see the raw seed). Do NOT
  loosen a RESULT assertion -- only placement/scope expectations change.

- [ ] **Step 5: build + run.** `ninja -C cmake-build-release unit_tests-sequant -j6`; run `[eval]`,
  `[eval_ta]`, `[placement_remat]`, `[dryrun]`, `[cache_manager]` -> GREEN (veto cases re-baselined, the
  invariant test passes). Any `[eval_ta]` computed-RESULT change = STOP (real regression).

- [ ] **Step 6: clang-format + commit** `eval: delete the has_demoted_external veto; placement is
  router-or-seed (Phase 4b-3 T1)`.

---

### Task T2 -- MPQC: the remat pre-pass in `build_cache_manager` (populate + own the router)

**Files:**
- Modify (MPQC): `src/mpqc/chemistry/qc/lcao/expression/sequant.h` (`build_cache_manager`), `sequant.cpp`
  (factor the shared dry-run `SizeRegime` construction).
- Optional (SeQuant): `SeQuant/core/eval/cache_manager.hpp` -- an OWNED-router option, if chosen for
  ownership (see step 3).

**Interfaces:**
- Consumes: `sequant::eval::remat_cells`, `rematerialize_to_budget`, `remat_to_router` (`placement_remat.hpp`);
  `dryrun::CostModel` / `SizeRegime` (built from `ctx.idx_to_extent`); `CacheManager::set_placement_router`
  (`cache_manager.hpp:280`); `ctx.batch_policy.peak_threshold` + the batch target sizes for `block_of`.
- Produces: `build_cache_manager` returns a cache with a populated (non-empty) router when
  `peak_threshold` is finite; the router is OWNED with lifetime >= the cache.

**Steps:**

- [ ] **Step 1: factor the dry-run `SizeRegime` construction.** The existing dry-run cost-profile path in
  `sequant.cpp` (~440+) already builds a `SizeRegime` from `ctx.idx_to_extent` (+ inner_pow). Extract a
  helper `make_dryrun_size_regime(ctx)` (and a `block_of` from `ctx.batch_policy`'s per-mode batch target
  sizes) reused by BOTH the cost-profile path and the new pre-pass, so they cannot drift.

- [ ] **Step 2: the pre-pass in `build_cache_manager`.** After building the cache but before returning,
  when `std::isfinite(engine's peak_threshold)` (grep the engine for how it exposes the batch policy /
  threshold to `build_cache_manager` -- via `ctx.batch_policy.peak_threshold`):
  ```cpp
  if (std::isfinite(peak_threshold)) {
    auto const cm = dryrun::CostModel{make_dryrun_size_regime(ctx)};
    auto const in  = sequant::eval::remat_cells(nodes, cm, block_of);
    auto const res = sequant::eval::rematerialize_to_budget(in.cells, cm, block_of,
                                                            in.num_points, peak_threshold);
    auto router = sequant::eval::remat_to_router(in.cells, res.cells, nodes);
    // own `router` for the cache's lifetime, then:
    cache.set_placement_router(&owned_router);
  }
  ```
  `nodes` are the real `EvalExprTA` enodes (remat is backend-agnostic; design 11.1).

- [ ] **Step 3: router ownership (pick one; state the choice in the commit).** The cache holds a
  NON-owning `PlacementRouter const*`. Options: (a) add an OWNING `std::unique_ptr<PlacementRouter<...>>`
  member + `adopt_placement_router(...)` to SeQuant `CacheManager`, and move the router in so
  `build_cache_manager`'s return-by-value carries it (return type unchanged, lifetime tied to the cache --
  RECOMMENDED); or (b) return a `{CacheManager, PlacementRouter}` bundle from `build_cache_manager` and
  update the two callers (`cck.ipp:1443`, `eom_cck.ipp:493`). Do NOT leave the router as a local -> the
  cache's pointer would dangle.

- [ ] **Step 4: build.** Build the `mpqc` target in `/Users/efv/code/mpqc4/cmake-build-release`
  (`ninja -C cmake-build-release mpqc -j6`). Fix any glue drift (like the Phase-2 cache_manager/CacheConfig
  drift already fixed). Confirm it links.

- [ ] **Step 5: smoke + BatchGroup.** Run ONE batched CSV-CCk validation input (e.g.
  `he10-csv-cck-2-pao-batched-...`) with `eval:trace` on; confirm `BatchGroup` events still fire (batching
  engaged) and the energy matches its reference. commit `lcao/expression: remat pre-pass populates the
  PlacementRouter in build_cache_manager (Phase 4b-3 T2)`.

---

### Task T3 -- cross-repo validation (USER-run gate)

**Files:** none (validation only).

**Steps:**

- [ ] **Step 1: SeQuant full suite** (controller-runnable): `ctest`/the unit binary -> green.
- [ ] **Step 2: MPQC CCk + CSV-CCk (np1+np2)** -- the load-bearing gate. `ctest --test-dir
  cmake-build-release -R "cck.*-np1"` and `MAD_NUM_THREADS=6 ctest -R csv-cck`: ALL energies UNCHANGED vs
  reference (the batched variants included). A moved energy = real regression (STOP; do not re-baseline).
- [ ] **Step 3: `threshold = inf` byte-identity** -- a non-batched CCk input is byte-identical to pre-4b-3
  (no pre-pass, empty router).
- [ ] **Step 4: peak drop (informative)** -- a batched input with a finite `peak_threshold`: confirm the
  modeled/measured peak is lower than the seed (remat spilled), energies still exact. Note this is the
  first phase where the batched peak actually drops at runtime.

---

## Self-review checklist

- Coverage: design 11.2 -> T2; 11.3 -> T1; 11.4 -> T1 step 2 (invariant test); 11.5 -> T3. All covered.
- The veto deletion (T1) and the router population (T2) are the COUPLED unit; T1 alone (empty router) homes
  batched giants at the seed (peak-maximal) -- acceptable in SeQuant unit tests (no pre-pass) but the
  batched CCk regression only shows post-T2, hence T3 after both.
- Router lifetime > cache (T2 step 3) -- the one memory-safety invariant.

## Deferred (NOT this plan)

- **Phase 5 -- feedback:** consume remat's `RebatchNeeded` / `FactorizationInherent` to drive the
  factorizer (re-batch) + reconcile the static `PeakProfile` vs measured runtime peak + tune move ordering.
- Un-hoist / split evict moves (if shrink alone leaves witnesses over target); incremental re-sweep (O2a).
