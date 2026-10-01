# Batchability Role-Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make batchability role-aware end to end: a mode is batchable in the CONTRACTED role or the EXTERNAL role, decided by two separate caller-supplied building-block predicates, with "batchable in any role" derived from them. This removes the overloaded single predicate that let the factorizer batch modes in a role the runtime never realizes.

**Architecture:** Two building-block predicates everywhere they occur (`BatchPolicy`, `OptimizeOptions`, `CostParams`, `PeakBatchedModel`): `is_batchable_contracted_index` and `is_batchable_external_index`, both DEFAULT-DECLINE (return false), declared adjacently. "Batchable in any role" is DERIVED (a union), never a settable field: `BatchPolicy::is_batchable_index()` returns a `std::function` union (consumers store/pass it); `PeakBatchedModel::is_batchable(ix)` is a `bool` method (called internally). Consumers route by intent: the DP contracted-role filter and the cache veto take the contracted building block; the runtime accept and cost-profile/footprint mirror take the union; the DP external-role filter takes the external building block. Batching configuration is consolidated onto `CostParams` so the two role predicates sit together and `single_term_opt` stops carrying loose batching params. Callers are migrated to set both roles explicitly; only then is the historical "empty external => fall back to contracted" behavior removed.

**Tech Stack:** C++20, SeQuant core (`SeQuant/core/...`), MPQC CCk (`mpqc4/src/mpqc/chemistry/qc/lcao/cc/...`), Catch2.

## Background / why (VERIFIED root cause)
The single `is_batchable_index` predicate did triple duty: DP contracted-role predicate, DP external-role fallback (`cost_model.hpp:974`), and runtime accept (`eval.hpp:1883`). Callers put a mode into it to satisfy the runtime accept, which as a side effect declared it CONTRACTED-batchable to the factorizer. The DP then sliced that mode in its contracted cells (verified: `cell_slices` contains the contracted mode; the modeled `szcell` is ~225x below the true nominal `sz0`, which equals the replay footprint), a reduction the runtime never realizes -> the intermediate materializes whole. Full evidence: `.superpowers/sdd/contamination-role-predicate.md`, `.superpowers/sdd/d1-gate-rootcause.md`.

## Global Constraints
- Branch `evaleev/feature/suppress-heuristic-fallback`. Do NOT push, open PRs, or repin `external/versions.cmake`.
- No `Co-Authored-By` trailers. ASCII hyphens only (pre-commit rejects U+2013 / U+00A0). Run clang-format; re-add + re-commit if the hook reformats. Cap builds at `-j6`. One task per commit.
- SeQuant unit tests build in `SeQuant/cmake-build-debug` (`--target unit_tests-sequant`). MPQC builds in the user's MPQC build dir (consumes SeQuant via `FETCHCONTENT_SOURCE_DIR_SEQUANT` = `$HOME/code/SeQuant`, same working tree).
- Domain-neutral in all SeQuant code/docs: use "contracted" / "external", never occupied / aux / pao. Domain terms appear ONLY in the MPQC CCk layer (Task 3).
- The derived "any role" predicate is NEVER a settable field. The DP contracted-role filter and the cache veto NEVER receive the union (that re-creates the bug). Building blocks are default-decline; no `pred ? pred(x) : fallback` guards.
- No dependency-pin changes.

## Naming (uniform across BatchPolicy / OptimizeOptions / CostParams / PeakBatchedModel)
- contracted-role predicate: `is_batchable_contracted_index`
- external-role predicate: `is_batchable_external_index`
- derived "any role": `is_batchable_index()` (BatchPolicy, returns `std::function`) / `is_batchable(ix)` (PeakBatchedModel, `bool` method)

## Consumer routing map (load-bearing)
| Site | current | route to |
|---|---|---|
| `optimize.cpp:116-144` -> DP contracted role | `is_batchable_index` (loose arg) | **contracted** (via `CostParams`) |
| `eval.hpp:1883` runtime accept | `policy.is_batchable_index` | **union** `policy.is_batchable_index()` |
| `cost_profile.hpp:223` (feeds eval accept) | `policy.is_batchable_index` | **union** |
| `cost_profile.hpp` -> cache veto predicate | (same field, shared) | **contracted** (decouple; veto is Contracted-stamp-only) |
| `cache_manager.hpp:609` veto test | `is_batchable_index(ix)` | **contracted** `is_batchable_contracted_index(ix)` |
| MPQC `sequant.h:202` / `sequant.cpp:407` (feed the cache VETO) | `policy.is_batchable_index` | **contracted** `is_batchable_contracted_index` |
| MPQC `sequant_engine.cpp:203` (config-presence gate) | `policy.is_batchable_index` | either building block (`contracted \|\| external`) |
| DP role filter `cost_model.hpp:952/974/982`, seed probes `:1465/:2030` | `is_batchable` / `is_batchable_external` | **building blocks** |

---

### Task 1: Structural role-split API + config consolidation (BEHAVIOR-PRESERVING)

Rename to the two building blocks + derived accessors across all four structs; consolidate batching config onto `CostParams`; collapse `single_term_opt`; migrate the ~30 `PeakBatchedModel` construction sites; route every consumer per the map. This task KEEPS the empty-default + external->contracted fallback so behavior is byte-identical; Task 4 removes it. The building-block DEFAULTS stay empty-with-null-safe-derivation in THIS task (default-decline is adopted in Task 4, when the fallback is also removed) -- so the derived accessors are null-safe here.

**Files:**
- `SeQuant/core/batch_policy.hpp` (BatchPolicy), `SeQuant/core/optimize/options.hpp` (CostParams :129-178, OptimizeOptions), `SeQuant/core/optimize/cost_model.hpp` (PeakBatchedModel :551-662; role filter :952/974/982; seed probes :1465/2030), `SeQuant/core/optimize/single_term.hpp` (single_term_opt signature :69-77 + body; construction :135), `SeQuant/core/optimize/single_term_detail.hpp` (:311 guard), `SeQuant/core/optimize/optimize.cpp` (:96-147 CostParams build + call sites), `SeQuant/core/eval/eval.hpp:1883`, `SeQuant/core/eval/backends/dryrun/cost_profile.hpp:57/223`, `SeQuant/core/eval/cache_manager.hpp:533/609`.
- ~28 `PeakBatchedModel{...}` construction sites in `SeQuant/tests/unit/test_optimize.cpp` and `test_eval_dryrun.cpp` (grep `PeakBatchedModel.*{idxsz`).
- MPQC readers: `mpqc4/src/mpqc/chemistry/qc/lcao/expression/sequant.{h,cpp}`, `sequant_engine.cpp` (union accessor).
- Test: `SeQuant/tests/unit/test_optimize.cpp` new `[optimize][role-api]`.

**Interfaces / Produces:**
- `BatchPolicy`: `is_batchable_contracted_index`, `is_batchable_external_index` (fields), `std::function<bool(Index const&)> is_batchable_index() const` (union of the two, null-safe this task).
- `CostParams`: `is_batchable_contracted_index` (adjacent to renamed `is_batchable_external_index`), `batch_target_size`, `batch_persistent_only`, `inner_pow` added.
- `PeakBatchedModel`: `is_batchable_contracted_index` + `is_batchable_external_index` adjacent, OUT of the positional prefix; `bool is_batchable(Index const&) const` method.
- `single_term_opt(network, idxsz, subnet_cse, CostParams const& cost = {}, container::vector<NodeBatchAnnotation>* out_axes = nullptr)`.

- [ ] **Step 1: Failing test** `[optimize][role-api]`: a `BatchPolicy` with `is_batchable_contracted_index = pred_A`, `is_batchable_external_index = pred_B`; assert `p.is_batchable_index()(a)` and `(b)` both true; with only contracted set, union == contracted. A `PeakBatchedModel` with the two set: `m.is_batchable(x)` == OR. A `CostParams` carries both predicates + `batch_target_size`.
- [ ] **Step 2: Run, verify fail** (`"[role-api]"`) -> compile/assert fail.
- [ ] **Step 3: BatchPolicy + CostParams + OptimizeOptions.** Add building blocks + derived `is_batchable_index()` (null-safe OR) to `BatchPolicy`. In `CostParams` rename `is_batchable_external` -> `is_batchable_external_index`, add `is_batchable_contracted_index` adjacent, add `batch_target_size`/`batch_persistent_only`/`inner_pow`. Mirror on `OptimizeOptions`.
- [ ] **Step 4: PeakBatchedModel.** Remove `is_batchable` from the positional prefix; declare `is_batchable_contracted_index` + `is_batchable_external_index` adjacently outside the prefix (empty default this task); add `bool is_batchable(Index const&) const`. Role filter: contracted branch -> `is_batchable_contracted_index`, external branch KEEP fallback `is_batchable_external_index ? is_batchable_external_index : is_batchable_contracted_index` (removed in Task 4); `batchable_mode_list`/seed probes use the renamed building blocks. Migrate ALL ~30 `PeakBatchedModel{...}` sites: drop the 2nd positional arg, add `model.is_batchable_contracted_index = <that arg>;`.
- [ ] **Step 5: single_term_opt + optimize.cpp.** Collapse the signature to `(network, idxsz, subnet_cse, cost, out_axes)`; read `is_batchable_contracted_index`/`batch_target_size`/`inner_pow`/`batch_persistent_only` from `cost`; delete the `(void)` unpacking noise. In `optimize.cpp` build `CostParams` field-by-field from `BatchPolicy` (contracted from `is_batchable_contracted_index`, external from `is_batchable_external_index`, plus the moved fields) and drop the loose args at the 6 call sites.
- [ ] **Step 6: Route consumers.** `eval.hpp:1883` accept = `policy.is_batchable_index()`; `cost_profile.hpp:223` eval-accept = union, cache veto predicate = `policy.is_batchable_contracted_index`; `cache_manager.hpp:609` -> `is_batchable_contracted_index(ix)`; MPQC readers -> `policy.is_batchable_index()`.
- [ ] **Step 7: Build + behavior-preserving suites.** `cmake --build cmake-build-debug --target unit_tests-sequant -j6 && ./cmake-build-debug/tests/unit/unit_tests-sequant "[optimize],[eval],[cache_manager],[dryrun-objective],[role-api]"` -> `[role-api]` PASS; others assertion-count UNCHANGED (record). Build the MPQC CCk target to confirm the reader edits compile.
- [ ] **Step 8: Commit** (`-m "batchability: split into contracted/external building-block predicates + derived union; consolidate batching config on CostParams"`).

### Task 2: SeQuant witnesses + external test cases role-split (BEHAVIORAL, the fix in-tree)

Move every SeQuant caller that wants external batching to set `is_batchable_external_index` explicitly (occ-like modes -> external), so the DP no longer sees them as contracted-batchable. Fallback still present but now inert for these callers.

**Files:** `SeQuant/tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]` 4262-4280, `[.][dryrun-extmode-avoidable]` 4709-4717); sweep `test_optimize.cpp`/`test_eval_dryrun.cpp` for any `batch_spectator`/external case lacking an explicit external predicate. Test: new `[optimize][role-filter]`.

- [ ] **Step 1: Failing ground-truth test** `[optimize][role-filter]`: network with a CONTRACTED batchable mode `x` and an EXTERNAL batchable mode `e`; config `is_batchable_contracted_index = {x-space}`, `is_batchable_external_index = {e-space}`. Build the DP context; assert a node carrying contracted `x` has `x` NOT in its cell (`cell_union(B)` decode), i.e. contracted `x` is not sliced; and `e` is handled via the external path. (Reuse the `[optimize]` context idiom + the `cell_union`/`batchable_modes` decode.)
- [ ] **Step 2: Run, verify fail** (`"[role-filter]"`) -> FAIL if a contracted batchable mode is still put in the contracted predicate by a fixture; construct the fixture to demonstrate the corrected routing (contracted mode only in contracted predicate) already passes, and the FAILURE is a companion assertion that the witnesses' occ mode is external-only -- adjust per the actual fixture so RED is real, not vacuous.
- [ ] **Step 3: Role-split the witnesses.** In both `[.]dryrun` witnesses move the occ predicate to `is_batchable_external_index`; the contracted predicate keeps only genuinely-contracted batchable modes (aux). Sweep + fix any other external-reliant case.
- [ ] **Step 4: Build + run** `"[optimize],[eval],[dryrun-occ-veto],[role-filter]"` -> `[role-filter]` PASS; record any intentional shifts in external `[optimize]` cases (must be updated here, not left red).
- [ ] **Step 5: Commit** (`-m "tests: set the external batchability role explicitly for external-batching callers"`).

### Task 3: MPQC csv_batch_policy role-split + readers + test contract + NodeBatchAnnotation resync

**Files:** `mpqc4/src/mpqc/chemistry/qc/lcao/cc/csv_batch_policy.h` (81-113), `mpqc4/src/mpqc/chemistry/qc/lcao/expression/sequant.{h,cpp}`, `sequant_engine.cpp`, `mpqc4/tests/unit/cck_batch_policy_test.cpp` (92-116).

**Pre-req resync (pre-existing cross-repo debt, fix FIRST in this task):** SeQuant's `OptimizeOptions::term_batch_axes` is now `shared_ptr<unordered_map<Expr const*, container::vector<NodeBatchAnnotation>>>` (options.hpp:317). MPQC `sequant_engine.cpp:183-186` still builds the pre-`NodeBatchAnnotation` type (`vector<svector<pair<Index,BatchModeType>>>`), so MPQC does not compile on this branch. Scope is a SINGLE site: `sequant_engine.cpp:183-185` re-declares `BatchAxesVec = vector<svector<pair<Index,BatchModeType>>>` (old element type) then `make_shared<BatchAxesMap>()` into the SeQuant-typed field -> mismatch. Fix: make MPQC's element type `container::vector<NodeBatchAnnotation>` (match `OptimizeOptions::term_batch_axes`), or drop the local typedef and use the SeQuant type / `decltype` directly. The READER at `sequant.cpp:328-331` (`bopts.node_batch_axes = it->second`) already uses the SeQuant types (both `term_batch_axes` and `BinarizationOptions::node_batch_axes` are `vector<NodeBatchAnnotation>`), so it needs no change -- but VERIFY it compiles after the typedef fix. This is what makes a full MPQC build possible for the rest of Task 3.

**Veto-reader correction (Task-1 carryover):** `sequant.h:202` / `sequant.cpp:407` feed the cache veto; switch them from `is_batchable_index()` (union) to `is_batchable_contracted_index` per the corrected routing map. Byte-identical, correct intent.

- [ ] **Step 1: Update the test contract (RED).** In `cck_batch_policy_test.cpp` assert occ lands in `is_batchable_external_index` and NOT in `is_batchable_contracted_index`; aux/pao in `is_batchable_contracted_index`; the derived `is_batchable_index()(occ)` still true (runtime accept). `CHECK_FALSE(p_off.is_batchable_index()(occ))`.
- [ ] **Step 2: Run, verify fail** (MPQC build) `[cck_batch]` -> FAIL.
- [ ] **Step 3: Implement `csv_batch_policy.h`.** `is_batchable_contracted_index` = aux OR pao (targets>0); `is_batchable_external_index` = pure-occupied (occ_target>0). Rewrite the 82-90 comment: runtime accept is the derived union, so occ in the external role is still accepted at runtime; the optimizer no longer treats occ as contracted-batchable. Readers `sequant.{h,cpp}`/`sequant_engine.cpp` call the union accessor.
- [ ] **Step 4: Build + run** `[cck_batch]` PASS; build the CCk target.
- [ ] **Step 5: Commit** (`-m "cck: occ batchable in the external role only; aux/pao in the contracted role"`).

### Task 4: Enforce default-decline + remove the fallback (CLEANUP, byte-identical post-migration)

All callers now set roles explicitly, so the fallback is dead and defaults can decline.

**Files:** `SeQuant/core/batch_policy.hpp`, `SeQuant/core/optimize/options.hpp`, `SeQuant/core/optimize/cost_model.hpp` (:974 fallback), `single_term_detail.hpp` (:311 guard).

- [ ] **Step 1: Migrate the remaining fallback-reliant callers (Task-2 sweep gap).** Task 2's sweep keyed on `batch_spectator_indices` and MISSED callers that rely on the fallback WITHOUT setting that flag. Removing the fallback WILL break them; that is EXPECTED -- migrate each in THIS commit with the Task-2 pattern `model.is_batchable_external_index = <the contracted predicate>;` (byte-identical to the fallback). This is reproduce-driven: after removing the fallback, run `[optimize]` and for EACH failure (a `REQUIRE(ctx.m == N)` mismatch or an out-of-bounds/segfault from `ctx.m` shrinking) add the explicit external role to that caller, and REPEAT until `[optimize]` is fully green. Known starting sites in `tests/unit/test_optimize.cpp`: (a) the SECTION at ~line 753 (uses `allK = (1<<m)-1` from `batchable_mode_list(net, is_batchable)` -> `ctx.m` shrink -> out-of-bounds SEGFAULT), (b) the `REQUIRE(ctx.m == 2)` case at ~line 2735/2740. Do NOT migrate contracted-only cases (a model whose batchable mode never occurs externally is unaffected by fallback removal and needs no external role) -- only migrate a caller that actually breaks.
- [ ] **Step 2: Default-decline + remove fallback.** Give both building blocks (all four structs) a default `[](Index const&){ return false; }`; delete the null-guards in the derived accessors (now safe); delete the role-filter fallback (`cost_model.hpp:974` -> `is_batchable_external_index` only); simplify the `single_term_detail.hpp:311` guard.
- [ ] **Step 3: Build + full suites** -> byte-identical to the Task 3 state (`[optimize]`/`[eval]`/`[cache_manager]`/`[role-api]`/`[role-filter]` all green, counts unchanged vs Task 3). This task changes NO behavior (all callers explicit); it only removes dead fallback + tightens defaults.
- [ ] **Step 4: Commit** (`-m "batchability: default-decline building blocks; remove the external->contracted fallback"`).

### Task 5: C60 acceptance gate -- honest sizing, re-measure

**Files:** `SeQuant/tests/unit/test_eval_dryrun.cpp` `[.][dryrun-occ-veto]`.

- [ ] **Step 1: Add acceptance assertions.** Aux+occ leg: (a) NO node's DP cell slices a contracted occ index (promote the `cell_union`/`batchable_modes` decode to a permanent gated helper); (b) external-occ scatter STILL fires (external occ stamps present, occ-scatter op count > 0 -- the union accept did not drop external occ). Record the honest re-measured `cost_profile.peak_bytes` in the comment. The two aspirational gate CHECKs (`peak<100`, `avoidable<0.10`) stay RED and documented -- NOT forced here.
- [ ] **Step 2: Run** `"[dryrun-occ-veto]"`, capture numbers.
- [ ] **Step 3: Commit** (`-m "dryrun gate: assert contracted occ is never batch-sliced and external scatter still fires; record honest C60 sizing"`).

## Out of scope (follow-ups)
- Doc retractions of the contaminated tables/conclusions (code-first; do after Task 5 with clean numbers). Catalog in `.superpowers/sdd/contamination-role-predicate.md`.
- The D1 commits (`3ac0b6db5` propagate-External-to-descendants, `4a5c727d7` emit-External-outer-
  of-Contracted): DISPOSITION = **KEEP** (2026-07-27, post-role-split). Their session-level
  motivation ("drop the phantom C60 2947 GB peak") was contaminated, but the CODE is correct,
  independently tested (emit-level `[loop-tree]` RED->GREEN in test_optimize.cpp), and now
  load-bearing: they are general improvements to the EXTERNAL-mode emit path, which the role-split
  makes the active occ-batching path. They let the honest schedule stamp external occ on the
  descendant giants (`external_occ_stamps == 244` in the `[dryrun-occ-veto]` gate) and slice them on
  the external occ pair -- which LOWERS the honest peak versus materializing them whole (though it
  cannot reduce the `a<i_3,i_4>` leg, hence 5860 GB). Commit messages describe the mechanism
  accurately (no false peak claim), so no history rewrite; the contaminated peak-drop framing was
  corrected in the docs-retraction pass and the gitignored session notes. Reverting would remove
  tested-correct external-emit behavior and conflicts with the role-split's edits to the same
  `cost_model.hpp` emit region.
- The downstream design question (cross-pair giants exposed as genuinely over-budget need contracted-occ batching or a refactor).

## Self-review
- Coverage: API split + consolidation (T1), in-tree config fix (T2), MPQC config fix (T3), fallback removal + default-decline (T4), acceptance gate (T5).
- Safe ordering: migrate callers explicit (T1 keeps fallback, T2/T3 set roles) BEFORE removing the fallback (T4). Ground-truth gate (T5) last.
- Type consistency: `is_batchable_contracted_index` / `is_batchable_external_index` fields; `is_batchable_index()` (policy, std::function) / `is_batchable(ix)` (model, bool) derived.
- No union to the DP contracted filter or the veto (bug-reintroduction guard), stated in Global Constraints and the routing map.
