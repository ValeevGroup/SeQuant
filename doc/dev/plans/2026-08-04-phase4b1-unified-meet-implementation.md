# Phase 4b-1 -- the unified meet (drop External filter, delete contracted_modes, retire seed_residency) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the runtime `sliced_modes` the cross-occurrence meet of ALL batched modes (drop the `External`-only filter), DELETE `contracted_modes` end-to-end, and retire the now-duplicate `seed_residency` -- a correctness-preserving unification. KEEP the veto (`has_demoted_external`).

**Architecture:** Today `place_at_this_level` homes a value at the deepest loop over `sliced_modes UNION contracted_modes`, where `sliced_modes` is External-only (`stamp_lifetime_masks` via `ext_modes_of`) and `contracted_modes` is a per-occurrence bolt-on. Dropping the `External` filter makes `sliced_modes` the all-batched-modes meet, which SUBSUMES `contracted_modes` (and corrects its latent per-occurrence aux under-meet), so `contracted_modes` is deleted wholesale and the union collapses to `sliced_modes` alone. `seed_residency` (Phase 3a's parallel all-modes field) then equals `sliced_modes` exactly and is retired. Design doc `doc/dev/specs/2026-08-04-phase4-threshold-o2-cutover-design.md` section 9; the CAT-1 site inventory is in `doc/dev/specs/2026-08-03-meet-based-home-scope-phase3-design.md` section 8.

**Tech Stack:** C++20 templated headers + `.cpp`, Catch2 (`tests/unit/`), CMake + Ninja, the dry-run backend.

## Global Constraints

- No en-dashes (U+2013) anywhere in source or comments.
- Format touched files with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`.
- No `Co-Authored-By` trailers in commits.
- RELEASE build in `cmake-build-release`; cap ninja at `-j6`.
- **CORRECTNESS IS THE GATE.** This phase CHANGES runtime placement (not byte-identical). The gate is
  that RESULTS are unchanged: SeQuant's `[eval]` / `[eval_ta]` correctness suites (which compute real
  contraction results) stay green. The final energy gate is the cross-repo MPQC CCk validation suite,
  which the USER runs -- FLAG in the report that MPQC CCk validation is required before this merges.
- **Witnesses may move, and that is EXPECTED.** `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]`
  peak figures MAY change (placement changed). RE-BASELINE them with a one-line justification in the
  commit (unified meet corrects the aux under-meet + drops the External/contracted asymmetry). Do NOT
  treat a moved witness figure as a regression -- but DO report the before/after numbers.
- The veto `has_demoted_external` STAYS in 4b-1 (its deletion is Phase 4b-3). Do NOT touch it.
- The CAT-1 audit line numbers PREDATE the Phase-4a refactor and have drifted. Treat the audit as a
  starting MAP; `grep` each symbol in the CURRENT tree and verify every site before editing.

## File Structure

- `SeQuant/core/eval/lifetime_mask.hpp` -- `stamp_lifetime_masks` selector (T1); delete `stamp_seed_residency` (T2); `home_scope` re-point (T2).
- `SeQuant/core/eval/eval_expr.{hpp,cpp}` -- delete `contracted_modes` field/accessors/threading (T1); delete `seed_residency_` field/accessors (T2).
- `SeQuant/core/eval/eval.hpp` -- `place_at_this_level` `in_union` / `residency_all_outer` contracted half (T1).
- `SeQuant/core/eval/node_batch_annotation.hpp`, `SeQuant/core/optimize/cost_model.hpp` -- `contracted_modes` annotation + emission (T1).
- `SeQuant/core/eval/schedule_dump.hpp` -- `contracted_modes` dump field (T1).
- `SeQuant/core/eval/peak_profile.hpp` -- `linearize_rich` calls `stamp_lifetime_masks`; `home_scope` reads `sliced_modes` (T2).
- `SeQuant/core/eval/occurrence_key.hpp` -- grep-verify (comment/helper only; likely no change).
- Tests: `test_eval_ta.cpp`, `test_eval_dryrun.cpp` (contracted_modes cases -> T1); `test_lifetime_mask.cpp`, `test_peak_profile.cpp`, `test_occurrence_key.cpp` (seed_residency cases -> T2).

---

### Task T1 -- unify `sliced_modes` (drop the External filter) + delete `contracted_modes`

**Files:**
- Modify: `SeQuant/core/eval/lifetime_mask.hpp` (the `ext_modes_of` selector)
- Modify: `SeQuant/core/eval/eval_expr.hpp`, `SeQuant/core/eval/eval_expr.cpp` (delete `contracted_modes`)
- Modify: `SeQuant/core/eval/node_batch_annotation.hpp`, `SeQuant/core/optimize/cost_model.hpp`
- Modify: `SeQuant/core/eval/eval.hpp` (`in_union` / `residency_all_outer`), `SeQuant/core/eval/schedule_dump.hpp`
- Test: `tests/unit/test_eval_ta.cpp`, `tests/unit/test_eval_dryrun.cpp`

**Interfaces:**
- Consumes: `EvalExpr::batched_here()`, `BatchModeType`, `EvalExpr::sliced_modes()`.
- Produces: `sliced_modes` = the all-batched-modes cross-occurrence meet (was External-only). `contracted_modes` no longer exists (`EvalExpr`, `NodeBatchAnnotation`, cost-model emission, `place_at_this_level`, `schedule_dump` all lose it).

**Steps:**

- [ ] **Step 1: grep-verify the current sites.** Run `grep -rn "contracted_modes" SeQuant/ tests/unit/` and `grep -rn "BatchModeType::External" SeQuant/core/eval/lifetime_mask.hpp`. Cross-check against the CAT-1 list (design 2026-08-03 section 8): `eval_expr.hpp` field/accessors + member; `node_batch_annotation.hpp`; `cost_model.hpp` emission; `eval_expr.cpp` threading; `eval.hpp` `residency_all_outer` contracted loop + `in_union` union half; `schedule_dump.hpp`; tests. Note the ACTUAL current line numbers.

- [ ] **Step 2: drop the External filter in `stamp_lifetime_masks`** (`lifetime_mask.hpp:148-153`): change the `ext_modes_of` lambda from
  ```cpp
  for (auto const& [ix, kind] : n->batched_here())
    if (kind == BatchModeType::External) detail::proto_expand_into(v, ix);
  ```
  to drop the `if` (proto-expand EVERY kind), matching `stamp_seed_residency`'s `all_batched_modes_of`. Rename the lambda to `all_batched_modes_of` and update the surrounding doc comment (it must no longer say "External-only" / "behavior unchanged / purely additive" -- it now IS the runtime residency = the all-modes meet).

- [ ] **Step 3: delete `contracted_modes` end-to-end.** Remove the `EvalExpr::contracted_modes()` get/set + the member (`eval_expr.hpp`), the `.cpp` threading (`eval_expr.cpp`), the `NodeBatchAnnotation::contracted_modes` field + its uses (`node_batch_annotation.hpp`), the cost-model emission (`cost_model.hpp`), and the `schedule_dump.hpp` field. In `eval.hpp` `place_at_this_level`: the `in_union` predicate becomes membership in `sliced_modes()` ALONE (drop the `contracted_modes` half); `residency_all_outer` drops its `contracted_modes` loop (it now checks only `sliced_modes`). Keep `residency_all_outer`'s structure and the veto conjunct (`has_demoted_external`) intact -- only the contracted input is removed.

- [ ] **Step 3a: verify the factorizer does not otherwise depend on the emission.** Read `cost_model.hpp:1885-1911` (the emission site): confirm it only WRITES the annotation (does not feed the DP's own peak/flops decisions). If it does feed a decision, STOP and report (that would make this more than a delete).

- [ ] **Step 4: update the tests that assert contracted-mode behavior.** `test_eval_ta.cpp` (the `contracted_modes` cases, audit ~2743-2760, 2837-2845, 2883, 4772-4773) and `test_eval_dryrun.cpp` (~4320): DELETE the assertions that read `contracted_modes` / expect the old `sliced UNION contracted` split; where a test asserted a placement that the unified meet changes, update the expectation to the unified-meet placement (justify in a comment). Do NOT touch the `has_demoted_external` veto cases.

- [ ] **Step 5: build + correctness run.** `ninja -C cmake-build-release unit_tests-sequant -j6`. Run `[eval]`, `[eval_ta]`, `[lifetime_mask]`, `[dryrun]` -- GREEN (results unchanged; this is the correctness gate). If an `[eval_ta]` RESULT (not a placement/scope expectation) breaks, STOP -- that is a real correctness regression, not a re-baseline.

- [ ] **Step 6: witness re-baseline.** Run the two slow witnesses FOREGROUND: `"[.][dryrun-occ-veto]"` and `"[.][dryrun-extmode-avoidable]"`. Record before (current: occ-veto 38897.4/6047.4/525.9, External-occ=244; extmode 6026.0/5999.7) and after. If any figure moved, UPDATE the witness expectations in `test_eval_dryrun.cpp` to the new values with a one-line justification comment. Paste both before/after in the report.

- [ ] **Step 7: clang-format + commit.**
  ```bash
  git add -A SeQuant/core/eval/ SeQuant/core/optimize/cost_model.hpp tests/unit/test_eval_ta.cpp tests/unit/test_eval_dryrun.cpp
  git commit -m "eval: unify sliced_modes to the all-batched-modes meet; delete contracted_modes (Phase 4b-1 T1)"
  ```

---

### Task T2 -- retire the now-duplicate `seed_residency`

**Files:**
- Modify: `SeQuant/core/eval/lifetime_mask.hpp` (delete `stamp_seed_residency`; `home_scope` -> `sliced_modes`)
- Modify: `SeQuant/core/eval/eval_expr.hpp` (delete `seed_residency_` member + accessors)
- Modify: `SeQuant/core/eval/peak_profile.hpp` (`linearize_rich` -> `stamp_lifetime_masks`)
- Test: `tests/unit/test_lifetime_mask.cpp`, `tests/unit/test_peak_profile.cpp`, `tests/unit/test_occurrence_key.cpp`

**Interfaces:**
- Consumes (post-T1): `sliced_modes()` is the all-modes meet.
- Produces: `home_scope(n)` returns `n->sliced_modes()`; `stamp_seed_residency` / `seed_residency_` deleted. `peak_profile`/remat behavior BYTE-IDENTICAL (they read the same set, now via `sliced_modes`).

**Steps:**

- [ ] **Step 1: grep-verify** `grep -rn "seed_residency\|stamp_seed_residency\|home_scope" SeQuant/ tests/unit/`. Expect: `lifetime_mask.hpp` (definition + `home_scope`), `eval_expr.hpp` (member/accessors), `peak_profile.hpp` (`linearize_rich` call), and tests.

- [ ] **Step 2: re-point `home_scope`** (`lifetime_mask.hpp`) to return `n->sliced_modes()` instead of `n->seed_residency()`; update its doc (the seed residency IS `sliced_modes` now). Change `peak_profile.hpp`'s `linearize_rich` to call `stamp_lifetime_masks(forest)` instead of `stamp_seed_residency(forest)`.

- [ ] **Step 3: delete** `stamp_seed_residency` (`lifetime_mask.hpp`) and `EvalExpr::seed_residency_` + its accessors (`eval_expr.hpp`).

- [ ] **Step 4: update tests.** In `test_lifetime_mask.cpp` / `test_peak_profile.cpp` / `test_occurrence_key.cpp`, replace `seed_residency()` reads with `sliced_modes()` and `stamp_seed_residency(...)` calls with `stamp_lifetime_masks(...)`. The `[lifetime_mask][seed]` equivalence tests (which asserted `seed_residency == sliced_modes` on External-only forests, and the Contracted-discriminating case) become tautological or move to assert the unified `sliced_modes` directly -- keep the Contracted-discriminating coverage (a Contracted mode on a slot IS now in `sliced_modes`), delete the now-vacuous "two fields agree" assertions.

- [ ] **Step 5: build + run.** `ninja -C cmake-build-release unit_tests-sequant -j6`; run `[lifetime_mask]`, `[peak_profile]`, `[placement_remat]`, `[occurrence_key]`, `[dryrun]`, `[eval_ta]` -- GREEN. Because `sliced_modes == ` the old `seed_residency` after T1, `[peak_profile]`/`[placement_remat]` figures (e.g. the anchor 481600, the end-to-end 16128->160) must be UNCHANGED from their T1-era values -- this byte-identity is T2's acceptance gate (the de-dup is inert). If any moved, STOP and report (the re-point is not equivalence-preserving).

- [ ] **Step 6: clang-format + commit.**
  ```bash
  git add -A SeQuant/core/eval/ tests/unit/test_lifetime_mask.cpp tests/unit/test_peak_profile.cpp tests/unit/test_occurrence_key.cpp
  git commit -m "eval: retire the duplicate seed_residency; home_scope reads sliced_modes (Phase 4b-1 T2)"
  ```

---

## Self-review checklist (run after writing, before dispatch)

- Coverage: design 9.1 step 1 -> T1 step 2; 9.1 step 2 -> T1 steps 3-4; 9.1 step 3 -> T2. Validation 9.3 -> T1 steps 5-6 (CCk + witness) and T2 step 5 (byte-identical de-dup). All covered.
- Types: `sliced_modes`, `home_scope`, `stamp_lifetime_masks` consistent across T1/T2. `contracted_modes`/`seed_residency`/`stamp_seed_residency` are fully removed by end of T2.
- The veto (`has_demoted_external`) is explicitly OUT (Phase 4b-3); no step touches it.

## Deferred (NOT this plan)

- **Phase 4b-2** -- thread the value hash / occurrence-key through `linearize_rich` -> `ValueCell`; `remat_to_router`.
- **Phase 4b-3** -- wire the dry-run remat pre-pass -> router -> `place_at_this_level`; DELETE `has_demoted_external`; the default-threshold policy. Needs its own runtime-integration design pass.
