# Path B: remat-enforced budget in cost_profile (crux-first)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in budget-aware dry-run entry that shapes the C60 schedule to a single budget B by running remat ONCE on the final fused DAG (home-scope / split), then replays it with the hold-gate off, so the reported peak actually meets B -- and pins that the REPLAYED peak <= B.

**Architecture:** Path B of `doc/dev/specs/2026-08-09-remat-into-cost-profile-design.md`. The pieces exist (`remat_cells` -> `rematerialize_to_budget` -> `remat_to_router` -> `cost_profile(..., router, max_footprint=0)`); `cost_profile` already accepts a router. The load-bearing UNKNOWN is whether the remat MODELLED peak (`peak_profile_sweep`) equals the replay REALIZED peak (`cost_profile`'s high-watermark). This plan is CRUX-FIRST: Task 1 measures that equivalence and GATES the rest. Do not build the wiring before Task 1's result is known.

**Tech Stack:** C++20, SeQuant `core/eval` + `core/eval/backends/dryrun`; Catch2 in `tests/unit/`.

## Global Constraints

- **`policy.peak_threshold` (DP batching) is untouched.** Path B enforces B by PLACEMENT (remat), on top of whatever batching the DP chose. Both should take the same B, but they are different mechanisms (batch = which axes sliced; placement = where homed).
- **Remat runs ONCE, on the final fused DAG.** Never per-DP-candidate (combinatorial). See spec Cost section.
- **Opt-in.** The existing `cost_profile` (explicit router / gate) stays intact for the equivalence tests and the MPQC pre-pass; path B is a new entry point. Empty-budget callers stay byte-identical.
- **No en-dashes U+2013 / U+00A0; clang-format touched files; no AI-attribution trailers.** Build `unit_tests-sequant -j6`; run TARGETED filters (an unrelated slow DP test hangs full runs).

---

### Task 1: CRUX PROBE -- does modelled peak equal replayed peak? (decision gate)

Before any wiring, measure whether `rematerialize_to_budget`'s modelled peak agrees with the replay's realized peak on the real C60 forest. The answer decides whether path B is a hard bound (Task 2 as written) or needs reconciliation (Task 1b) first.

**Files:**
- Add: a `[.]` scratch/witness test in `tests/unit/test_eval_dryrun.cpp` (or a new `test_remat_cost_profile.cpp` added to `tests/unit/CMakeLists.txt`).

**Interfaces:**
- Consumes: the occ-veto C60 forest-loading + `policy` setup (reuse it); `remat_cells`, `rematerialize_to_budget`, `to_schedule`, `peak_profile_sweep`, `remat_to_router`, `cost_profile`.

- [ ] **Step 1: Build the aux+occ C60 forest** exactly as the occ-veto `run(true,true)` does (same `optimize` + `binarize` + `block_of = policy.batch_target_size`).
- [ ] **Step 2: Run remat at a budget B** (e.g. 1000 GB): `remat_cells` -> `rematerialize_to_budget(cells, cm, block_of, num_points, B)`. Record `peak_profile_sweep(to_schedule(res.cells,...)).peak_bytes` -- the MODELLED peak remat believes it achieved -- and `res.status` / `res.modeled_recompute`.
- [ ] **Step 3: Replay the placed schedule** with the derived router and gate off: `router = remat_to_router(seed, res.cells, forest)`; `cost_profile(forest, policy, cfg{max_footprint=0}, regime, nullptr, &router)`. Record the REALIZED `cp.peak_bytes`.
- [ ] **Step 4: Print and compare** MODELLED vs REALIZED peak (ratio + absolute) across a few B in {500, 1000, 2000} GB. `std::wcerr` a table. NO hard assertion yet -- this is a measurement.
- [ ] **Step 5: Decision.** If `realized ~= modelled` (within a small, understood factor) for all B -> the equivalence holds; SKIP Task 1b, proceed to Task 2. If they diverge materially, STOP and do Task 1b (diff the remat cell schedule against the replay `SCHEDULE_RUN_EVENT` dump; identify the divergence per spec "Likely divergence sources": tile granularity, interval-vs-scratch residency, replication count). Record the finding in the spec's evidence section before writing wiring that assumes a bound it does not have.

- [ ] **Step 6: Commit the probe** (as a `[.]` witness, so the measurement is reproducible).

```bash
git add tests/unit/... tests/unit/CMakeLists.txt
git commit -m "test: probe remat modelled-vs-replayed peak equivalence on C60 (path B gate)"
```

---

### Task 1b (CONDITIONAL -- only if Task 1 found divergence): reconcile the peak models

Close the gap the probe found, using the `2026-08-05-dryrun-wetrun-schedule-equivalence` instrument (the `SCHEDULE_RUN_EVENT` per-node build/fetch + ctx dump). Scope depends on the divergence; likely one of: align `cell_footprint`'s block extent with the replay's realized tile granularity; reconcile the interval-sweep co-residency with the batched-scratch high-watermark; or correct the replication count. This task is deliberately unscoped until Task 1 names the divergence -- do NOT pre-write it.

---

### Task 2: `cost_profile_to_budget` -- the opt-in wiring (GATED on Task 1/1b)

Add a budget-aware entry that runs the remat stage once and calls the existing `cost_profile` with the derived router and gate off.

**Files:**
- Modify: `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` -- add `cost_profile_to_budget(forest, policy, regime, double budget_bytes, ...)`.
- Test: `tests/unit/test_eval_dryrun.cpp` (or the new file).

**Interfaces:**
- Produces: `CostProfile cost_profile_to_budget(...)` = `remat_cells` -> `rematerialize_to_budget(B)` -> `remat_to_router` -> `cost_profile(..., cfg{max_footprint=0}, &router)`, returning the profile plus (via a field or out-param) `res.status` and `res.modeled_recompute`.
- Consumes: `block_of` from `policy.batch_target_size` (confirmed same partition in Task 1 Step 1).

- [ ] **Step 1: Write the failing test.** For the C60 aux+occ forest at budget B, assert `cost_profile_to_budget(...).peak_bytes <= B` (or the Task-1-documented slack). Without the wiring this does not compile.
- [ ] **Step 2: Run it to confirm it fails** (missing symbol).
- [ ] **Step 3: Implement `cost_profile_to_budget`** as the four-line pipeline above.
- [ ] **Step 4: Run tests.** The budget-enforcement assertion passes at the slack Task 1 established; `res.status` is `Feasible` (or the honest infeasible status if B is below what placement can reach).
- [ ] **Step 5: Commit.**

---

### Task 3: Budget-swept C60 witness + retire the last hardcoded budgets

Turn the path-A gate-off witnesses into a single-budget-B sweep where B drives BOTH `policy.peak_threshold` (batching) AND the remat budget (placement), and the reported REPLAYED peak <= B.

**Files:**
- Modify: `tests/unit/test_eval_dryrun.cpp` -- a new `[.]` witness (or extend occ-veto) that sweeps B and asserts `cost_profile_to_budget(...).peak_bytes <= B` per B, plus `res.modeled_recompute` tracks `avoidable_flops` within the Task-1 slack.

- [ ] **Step 1: Write the sweep** over B in {2000, 1000, 500} GB; per B, set `policy.peak_threshold = B`, run `cost_profile_to_budget(forest, policy, regime, B)`, assert `peak_bytes <= B` (+ slack).
- [ ] **Step 2: Run + confirm.** Document the realized peak / recompute per B (the corrected, budget-honest C60 scan).
- [ ] **Step 3: Commit.**

---

### Task 4: Empty-budget / regression

- [ ] **Step 1:** Run `[dryrun],[dryrun-costmodel],[composite-canon],[peak_profile],[placement_router],[placement_remat],[eval]`; confirm byte-identical (the new entry is additive; nothing else changed).
- [ ] **Step 2: No commit** (verification only).

---

## Self-review notes

- **Spec coverage:** crux equivalence (spec "The crux risk") -> Tasks 1/1b; the wiring (spec Design B) -> Task 2; budget enforcement (spec Testing 3) -> Task 3; regression (spec Testing 5) -> Task 4. Path A (gate-off measurement) is already shipped.
- **Crux-first, contingency honest:** Task 1 is a measurement with a decision, not a foregone conclusion; Task 1b is explicitly unscoped until Task 1 names the divergence; Tasks 2-3 are GATED on the equivalence holding (or being reconciled). No placeholder pretends the bound exists before it is measured.
- **Cost boundary respected:** remat runs once (Task 2), never per-candidate; `policy.peak_threshold` (the DP proxy) is untouched.
- **Watch item:** if Task 1 shows the peaks cannot be reconciled cheaply, the honest outcome (per spec) is to keep `max_footprint` as an OPTIONAL backstop and document the slack -- NOT to assert a hard `peak <= B` the replay does not honor. Task 3's assertion must use the slack Task 1 establishes.
