# Path A: measure the true C60 profile with the replay hold-gate off

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the dry-run C60 witnesses report the schedule's TRUE peak and avoidable recompute by measuring with the replay hold-gate off (`cfg.max_footprint = 0`, already the library default), retiring the hardcoded `1e11` (100 GB) gate that was distorting the measurement, and pinning the confound-free invariants.

**Architecture:** Path A of `doc/dev/specs/2026-08-09-remat-into-cost-profile-design.md`. NO library change -- `CacheConfig::max_footprint` already defaults to `0` (no gate). This is entirely test-side: the C60 witnesses hardcode `1e11`; make that an env override defaulting to gate-off, re-baseline the assertion bands to the true (gate-off) numbers, and add invariant tests. Batching itself is unaffected: it is driven by `policy.peak_threshold` (the DP), a SEPARATE knob.

**Tech Stack:** C++20, SeQuant `core/eval/backends/dryrun`; Catch2 tests in `tests/unit/test_eval_dryrun.cpp`.

## Global Constraints

- **No library change.** Only `tests/unit/test_eval_dryrun.cpp` is touched. `cost_profile` / `CacheConfig` are unchanged (default `max_footprint = 0` already means gate-off).
- **`policy.peak_threshold` is untouched** -- it is the DP's batching budget (per-node `subtree_peak` ceiling), a different knob from `cfg.max_footprint` (the replay hold-gate). Do not conflate them.
- **No en-dashes U+2013 or U+00A0** (pre-commit hooks reject them). Run `clang-format` on the touched file before committing (the commit hook also runs it; re-stage if it reformats).
- **No AI-attribution trailers** in commit messages.
- **Build:** `cmake --build cmake-build-release --target unit_tests-sequant -j6` (cap `-j`). The witnesses are `[.]` hidden (excluded from CI/default runs); run them BY NAME, not via the full suite (an unrelated slow DP test hangs full runs).
- **Measured numbers are captured, not invented.** Re-baseline bands come from running gate-off and reading the printed value; the plan gives the known ones and the exact command to capture the rest.

## Known gate-off numbers (occ-veto, 55 terms, `peak_threshold` = 100 GB)

Captured this session with `cfg.max_footprint` >= realized peak (gate inert):

| mode | avoidable | peak (GB) | total flops |
|---|---|---|---|
| unbatched | 0% | 50814 | 1.5798408944778614e16 |
| aux-only | 0% | 18585 | 1.586e16 |
| aux+occ | 75.80% | 563 | 2.0524411686101389e17 |

Contrast the gated (`1e11`) values (unbatched 60.8%, aux-only 19.2%, aux+occ 25.2%): the gate was inflating avoidable via evictions and hiding the true inherent-recompute tradeoff.

---

### Task 1: Env-gate the C60 witnesses' `max_footprint`, defaulting to gate-off

Make the replay hold-gate env-configurable in BOTH C60 witnesses, defaulting to `0` (gate off = measure the true profile). The 100 GB scenario becomes an opt-in override, mirroring the existing `SEQUANT_UT_DRYRUN_PEAK_THR_GB` pattern.

**Files:**
- Modify: `tests/unit/test_eval_dryrun.cpp` -- the `[.][dryrun-occ-veto]` run lambda's `cfg.max_footprint` (currently `1e11`, ~line 5239) and the `[.][dryrun-extmode-avoidable]` run lambda's `cfg.max_footprint` (currently `1e11`, ~line 5571).

**Interfaces:**
- Produces: both witnesses read `SEQUANT_UT_DRYRUN_MAXFP_GB` (GB, `* 1e9`) with default `0.` (gate off). A value of `100` reproduces the old gated behavior.

- [ ] **Step 1: Change the occ-veto gate to env-default-off.** Replace `cfg.max_footprint = 1e11;` (occ-veto run lambda) with:

```cpp
// Replay hold-gate. Default OFF (0) -> measure the schedule's TRUE peak and
// avoidable (see doc/dev/specs/2026-08-09-remat-into-cost-profile-design.md,
// path A). The gate is NOT the batching budget (that is policy.peak_threshold,
// a separate knob); it only decides what the replay HOLDS vs recomputes, and
// hardcoding 100 GB here was distorting the measurement. Set
// SEQUANT_UT_DRYRUN_MAXFP_GB=100 to reproduce the old gated numbers.
cfg.max_footprint = (std::getenv("SEQUANT_UT_DRYRUN_MAXFP_GB")
                         ? std::atof(std::getenv("SEQUANT_UT_DRYRUN_MAXFP_GB")) * 1e9
                         : 0.);
```

- [ ] **Step 2: Apply the identical change to the extmode-avoidable run lambda** (~line 5571), same comment.

- [ ] **Step 3: Build.** `cmake --build cmake-build-release --target unit_tests-sequant -j6`. Expect it to compile; the witnesses' OLD assertion bands will now FAIL (they were calibrated to the gated numbers) -- that is expected and fixed in Tasks 2-3.

- [ ] **Step 4: Commit.**

```bash
git add tests/unit/test_eval_dryrun.cpp
git commit -m "test: C60 witnesses measure gate-off by default (env-override for 100 GB)"
```

---

### Task 2: Re-baseline the occ-veto witness + pin the confound-free invariants

Update the three modes' avoidable/peak bands to the true gate-off numbers, and add explicit assertions that make the confound impossible to reintroduce.

**Files:**
- Modify: `tests/unit/test_eval_dryrun.cpp` -- the `[.][dryrun-occ-veto]` assertion block and the long narrative comment above it (add a short header note; do NOT rewrite the whole narrative -- append a pointer).

**Interfaces:**
- Consumes: the gate-off measurements from Task 1.

- [ ] **Step 1: Capture the current gate-off numbers** (confirm they match the table; they are seed-stable):

```
./cmake-build-release/tests/unit/unit_tests-sequant "dryrun occ batching wipes CSE (free-batchable-mode veto repro)" 2>&1 | grep 'occ-veto]'
```
Expected: unbatched `0%` / 50814 GB / total 1.5798...e16; aux-only `0%` / 18585 GB; aux+occ `75.8%` / 563 GB / total 2.05e17.

- [ ] **Step 2: Replace the stale avoidable bands** with the true invariants. The old `unbatched/aux avoidable < 0.05` (or whatever the current stale band is) becomes:

```cpp
// Gate-off = the schedule's TRUE profile (no eviction recompute). Perfect-CSE
// with enough memory: unbatched and aux-only hold everything and recompute
// NOTHING. aux+occ trades a low peak (563 GB) for real, INHERENT per-block
// re-form recompute (~76%) -- the honest cost of occ batching, not a gate
// artifact. See doc/dev/specs/2026-08-09-remat-into-cost-profile-design.md.
CHECK(unbatched.avoidable_time() == Catch::Approx(0.0).margin(1e-9));
CHECK(aux.avoidable_time() == Catch::Approx(0.0).margin(1e-9));
CHECK(both.avoidable_time() > 0.5);  // occ batching genuinely recomputes
```

- [ ] **Step 3: Pin the perfect-CSE floor** (the number the annotation-space fix landed on):

```cpp
// Unbatched, gate-off = the perfect-CSE floor: every value built once, full.
CHECK(unbatched.total_flops == Catch::Approx(1.5798408944778614e16).epsilon(1e-6));
```

- [ ] **Step 4: Pin gate-INDEPENDENCE (the blunt instrument is gone).** Run unbatched at two large gate values and assert identical avoidable -- proof the gate no longer governs once it clears the realized peak. Add a small helper run at an explicit `max_footprint` (do NOT rely on the env var inside the test; pass `cfg.max_footprint` directly for two values, e.g. `20e12` and `60e12`, and assert both give `avoidable_time() == 0` for unbatched). Keep it to the unbatched mode (fast).

- [ ] **Step 5: Append a header note** to the witness's narrative comment (one short paragraph): the witness now measures the gate-OFF true profile by default; `SEQUANT_UT_DRYRUN_MAXFP_GB=100` reproduces the historical gated numbers; the gate is the replay hold-gate, distinct from `policy.peak_threshold`. Point to the spec. Do not delete the existing role-stamp assertions (`Contracted-occ == 0`, `External-occ > 0`) -- those are batching-role guarantees, unaffected by the gate.

- [ ] **Step 6: Run + confirm green.**

```
./cmake-build-release/tests/unit/unit_tests-sequant "dryrun occ batching wipes CSE (free-batchable-mode veto repro)"
```

- [ ] **Step 7: Commit.**

```bash
git add tests/unit/test_eval_dryrun.cpp
git commit -m "test: re-baseline occ-veto to the true gate-off profile; pin CSE floor + gate independence"
```

---

### Task 3: Re-baseline the extmode-avoidable witness to gate-off

Same treatment for the contracted-vs-external witness: its `avoidable_time() < 0.05` bands were also gated. Re-baseline to the gate-off numbers.

**Files:**
- Modify: `tests/unit/test_eval_dryrun.cpp` -- the `[.][dryrun-extmode-avoidable]` assertion block + a header-note pointer.

- [ ] **Step 1: Capture the gate-off numbers.**

```
./cmake-build-release/tests/unit/unit_tests-sequant "dryrun external-mode occ batching matches contracted-mode avoidable" 2>&1 | grep 'extmode-avoidable]'
```
Read the printed contracted-occ / external-occ avoidable + peak.

- [ ] **Step 2: Replace the stale `avoidable_time() < 0.05` bands** with the captured gate-off values (as `Approx` or a bracketing range with a short comment: these are the schedule's inherent recompute, gate-off). Keep the structural checks (`replay_ops > 0`, the scatter/stamp markers, `model_flops == baseline`) unchanged.

- [ ] **Step 3: Append the same header note** (gate-off default; `SEQUANT_UT_DRYRUN_MAXFP_GB` override; spec pointer).

- [ ] **Step 4: Run + confirm green.**

```
./cmake-build-release/tests/unit/unit_tests-sequant "dryrun external-mode occ batching matches contracted-mode avoidable"
```

- [ ] **Step 5: Commit.**

```bash
git add tests/unit/test_eval_dryrun.cpp
git commit -m "test: re-baseline extmode-avoidable witness to the true gate-off profile"
```

---

### Task 4: Verify empty-budget regression (no other caller changed)

Confirm no non-witness `cost_profile` behavior moved -- the library default was already gate-off, so the visible suites must stay byte-identical.

**Files:** none (verification only).

- [ ] **Step 1: Run the real (non-hidden) suites.**

```
./cmake-build-release/tests/unit/unit_tests-sequant "[dryrun],[dryrun-costmodel],[composite-canon],[peak_profile],[placement_router],[placement_remat]"
```
Expected: all pass, counts unchanged from before Task 1 (the witnesses are `[.]` and not in these tags; nothing else consults `max_footprint`).

- [ ] **Step 2: Run `[eval]`** and confirm byte-identity (457 assertions / 28 cases).

- [ ] **Step 3: No commit** (verification only; if a non-witness test moved, STOP -- something other than the witnesses depended on the gated value, which contradicts the "test-side only" premise).

---

## Self-review notes

- **Spec coverage:** path A "confound gone" (spec Design A) -> Tasks 1-3; "perfect-CSE floor" (spec Testing 2) -> Task 2 Step 3; "no gate double-govern / independence" (spec Testing 1) -> Task 2 Step 4; "empty-budget regression" (spec Testing 5) -> Task 4. Path B (remat enforcement, swept B) is explicitly NOT in this plan -- a separate follow-up.
- **No placeholders:** the one genuinely runtime-derived value is the extmode gate-off band (Task 3), captured by the exact command in Step 1; every other assertion has its concrete number (0, 1.5798...e16, > 0.5).
- **Scope check:** single subsystem (dry-run witnesses), single file, no library change -- appropriately one plan.
- **Watch item:** the `both.avoidable_time() > 0.5` band (Task 2) is a floor, not an exact match -- aux+occ avoidable (75.8%) is schedule-dependent; assert a robust bound, not the exact percent, so the witness does not become brittle to unrelated DP changes. Same principle for the extmode bands.
