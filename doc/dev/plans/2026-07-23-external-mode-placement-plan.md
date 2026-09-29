# External-mode placement (stage 3, cheap variant) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make external modes first-class in the order-aware batching cost model by placing them on the fixed contracted schedule's over-budget nodes (the CHEAP variant: process the contracted root frontier before `select_root`), so external batching reaches the C60 peak-setting operand and drops its modeled peak.

**Architecture:** Two phases. Phase 1 = the existing contracted order-aware DP, with `build_cells` restricted to CONTRACTED modes. Phase 2 = a pass over the root frontier: for each candidate schedule, walk its tree, and at each node whose modeled peak exceeds `peak_threshold`, try injecting each carried external mode (append inner via `descend`), re-price the subtree, keep it if peak drops, iterate for cascade; `select_root` then chooses among the externally-augmented points. Reuses `escaped_outer` + resident-scan; replaces the root-level `seeded_forest_peak`.

**Tech Stack:** C++20, SeQuant, Catch2, the dry-run `cost_profile` backend. Build in `cmake-build-debug`.

## Global Constraints

- **Spec:** `doc/dev/specs/2026-07-23-external-mode-placement-design.md`. Scope = CHEAP variant only; the COMPLETE variant (external per DP candidate) is a documented TODO, OUT OF SCOPE.
- **Gated + byte-identical off:** phase 2 runs only when `order_aware_recompute && batch_spectator_indices`; both off (default) => byte-identical to today (`[optimize]`/`[eval]`/`[cache_manager]`/`[dryrun-objective]` unchanged).
- **`peak_threshold` = memory budget in bytes:** per-schedule in phase 1 (`select_root` ceiling, unchanged); per-node in phase 2 (a node's modeled peak > threshold triggers external injection).
- **Unified pricing:** an external mode in `B` is priced by the SAME `escaped_outer` (recompute for enclosed non-carriers) + resident-scan (hoist peak) as a contracted mode. No new cost formula.
- **External = `is_external_mode`:** open on the root, contracted at no node. `descend` appends it INNER (nest position determined by the injection node, not searched).
- **Annotation contract:** the DP emits placement; a well-behaved runtime respects it (evaluate the placed node's subtree one external-slice at a time). Not a DP guarantee.
- **Cross-repo:** validate on the C60 dry-run (DP-side modeled peak via `cost_profile`, NOT the replay `avoidable_time`) before any MPQC repin. No `Co-Authored-By`. ASCII hyphens only (pre-commit U+2013 detector).

---

## File structure

- `core/optimize/cost_model.hpp` -- ALL changes here. `build_context` / `Context::build_cells` (contracted-only enumeration), `PeakBatchedModel::is_external_mode` (unchanged, reused), `reconstruct_batched_modes` (the emit walk phase 2 extends / replaces `seeded_forest_peak` in), `seeded_forest_peak` (retired or repurposed), `select_root` (consumes the phase-2-augmented frontier).
- `tests/unit/test_optimize.cpp` -- DP-level unit gates (external placement drops a node's modeled peak; contracted-only cells).
- `tests/unit/test_eval_dryrun.cpp` -- the C60 giant DP-side measurement (reuse the `[.][ordered-key-c60-m]` fixture setup).

---

## Task S3.1: Investigate the phase-2 hook + re-price mechanism (INVESTIGATION)

The spec fixes the algorithm; the exact code shape depends on how the root frontier, the reconstruction walk, and subtree re-pricing are structured, and whether `seeded_forest_peak`'s sizing is reusable at a node. Settle this by reading before editing.

**Files (read-only):** `core/optimize/cost_model.hpp` -- `seeded_forest_peak` (the existing root seed + re-size + work-neutrality guard), `reconstruct_batched_modes` (the tree walk + the current external over_budget/adopt loop), `select_root`, the `State`/`BFrontPoint` frontier, `sliced_footprints` / `sz` (how a subset re-sizes under a sliced set).

**Deliverable** (write to `.superpowers/sdd/s3-1-note.md`): answer, with exact function+line references,
1. Where the root frontier lives after phase 1 and where `select_root` reads it -- the exact insertion point for a "process the frontier before selection" pass.
2. How to compute a NODE's modeled peak for the over-threshold trigger (the `BFrontPoint.peak` per subset/cell, or a walk), and how to re-price a subtree's peak with an external mode injected into `B` at that node -- specifically whether `seeded_forest_peak`'s re-size machinery (it already re-sizes the whole tree with a seed mode) can be applied at a subtree/node instead of the root, or whether a new local re-price is needed.
3. Whether the injection reuses `descend`/`escaped_outer`/resident-scan directly (it should: an external mode in `B` is a mode in `B`), and what changes to make `descend` accept an external mode (today `descend` is fed `aprime` = contracted-here only).
4. The minimal change to make `build_cells` enumerate CONTRACTED modes only (the external-mode bitmask must be computed in `build_context` via `is_external_mode` and passed to `build_cells`, which currently enumerates all `m`).

- [ ] **Step 1:** Read the functions; write the note answering (1)-(4) with exact change points and whether `seeded_forest_peak` is reusable node-level or must be generalized.
- [ ] **Step 2:** Confirm the note does not require changing the flag-off path (byte-identical). Report DONE with the note path.

## Task S3.2: `build_cells` enumerates contracted modes only

**Files:** `core/optimize/cost_model.hpp` (`Context` fields + `build_cells`, `build_context`); `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Produces: `Context` carries an `external_mask` (bitmask of `is_external_mode` bits); `build_cells` skips those bits when enumerating sequences, so ordered cells contain contracted modes only.

- [ ] **Step 1: Write the failing test** (`[optimize][ext-place]`): build a network with one contracted batchable mode `F` and one EXTERNAL batchable mode `x` (open on the root, e.g. a free index of the result). With `order_aware_recompute` on, assert that no ordered cell's sequence contains `x`'s bit -- i.e. `for all id: (cell_union(id) & (1<<x_bit)) == 0`. FAILS today (all `m` modes enumerated).

```cpp
// tests/unit/test_optimize.cpp  (sketch; reuse the [loop-tree] F-space idiom)
// network: g{r;F1} h{F1;x}  -> result carries x (external), F1 contracted.
// assert cell_union(id) never includes x's bit under order_aware.
```

- [ ] **Step 2: Run to confirm it fails.** `./tests/unit/unit_tests-sequant "[ext-place]"`.
- [ ] **Step 3: Implement:** in `build_context`, compute `ctx.external_mask` via `is_external_mode(ctx, k)` for each `k`; in `build_cells`, skip `k` when `external_mask & (1<<k)` during enumeration and descend. Keep the cell-count guard.
- [ ] **Step 4: Run to confirm it passes** + no regression: `./tests/unit/unit_tests-sequant "[optimize]"` unchanged (flag-off byte-identical; the `[ordered-key]` 1600->400 gate still passes since its modes are contracted).
- [ ] **Step 5: clang-format + Commit.**

## Task S3.3: Phase-2 external placement on the root frontier

**Files:** `core/optimize/cost_model.hpp` (a new `place_external_on_frontier` pass invoked before `select_root` in `reconstruct_batched_modes`, per S3.1); `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Consumes: S3.1's change note; S3.2's `external_mask`; `descend`, `escaped_outer`, resident-scan.
- Produces: for each root frontier point, an externally-augmented `(peak, flops)` variant added to the frontier when it cuts an over-threshold node's peak; `External` emitted for placed modes.

- [ ] **Step 1: Write the failing test** (`[optimize][ext-place]`): a network with an EXTERNAL mode `x` on a node whose modeled peak exceeds a set `peak_threshold`, where slicing `x` there drops the peak. With `order_aware_recompute + batch_spectator_indices` on and a finite `peak_threshold`, assert the chosen root peak is LOWER than with `batch_spectator_indices` off (external placement fired). State the expected pre/post peak from the `sz` of the node with/without `x` sliced (block/extent).
- [ ] **Step 2: Run to confirm it fails** (no external placement today at node level).
- [ ] **Step 3: Implement** per S3.1: before `select_root`, walk each frontier schedule's tree; at each node with modeled peak > `peak_threshold`, for each carried external mode inject it (`descend` append inner), re-price the subtree, `pareto_insert` if peak drops; iterate until no node over budget or no external helps. Gate on `order_aware_recompute && batch_spectator_indices`.
- [ ] **Step 4: Run to confirm it passes** + no regression (`[optimize]`, flag-off byte-identical).
- [ ] **Step 5: clang-format + Commit.**

## Task S3.4: Retire / subsume the root-level `seeded_forest_peak`

**Files:** `core/optimize/cost_model.hpp`; `tests/unit/test_optimize.cpp`.

- [ ] **Step 1:** Per S3.1, determine whether the phase-2 pass subsumes `seeded_forest_peak`'s root-level seeding (a node = the root is a special case of node-level placement). If so, route the root case through phase 2 and remove the separate `seeded_forest_peak` adopt-loop; if the work-neutrality guard is still needed anywhere, keep only that.
- [ ] **Step 2:** Run the existing external-mode tests (grep `[optimize]` for `External`/`spectator` cases) and confirm they still pass (or are updated to the node-level path with equal/better peak). No regression.
- [ ] **Step 3: clang-format + Commit.**

## Task S3.5: C60 giant DP-side measurement (GATE)

**Files:** `tests/unit/test_eval_dryrun.cpp` (extend the `[.][ordered-key-c60-m]` fixture setup into a `[.][ext-place-c60]` measurement).

- [ ] **Step 1: Write the measurement:** load the giant (fixture summand 38, DF spaces registered), build the batched context with `order_aware_recompute + batch_spectator_indices` on and a finite `peak_threshold` (e.g. 100e9), run `solve_single_term` + the phase-2 pass, and print the peak-setting node's modeled peak WITH external placement vs WITHOUT (`batch_spectator_indices` off). Assert the peak drops (external placement reached the `g.C` operand).
- [ ] **Step 2: Run** (hidden tag): `./tests/unit/unit_tests-sequant "[ext-place-c60]"`. Record the modeled peak drop.
- [ ] **Step 3: Decision gate:** if the giant's peak drops materially (external reached the operand), the cheap variant works -- record it. If the peak stays STUCK because the winning schedule was pruned pre-external, that is the trigger for the COMPLETE variant (the documented TODO), NOT more cheap-variant tuning -- escalate.
- [ ] **Step 4: Commit** the measurement + the recorded numbers.

---

## Notes for the executor

- **S3.1 gates S3.3.** The phase-2 re-price mechanism (reuse `seeded_forest_peak`'s sizing node-level, or a new local re-price) is the load-bearing unknown; do S3.1 first. If it finds the re-price is a much larger change than a node-level `seeded_forest_peak` call, escalate rather than force it.
- **Every batched change must keep flag-off byte-identical** and the `[ordered-key]` 1600->400 + `[batched-accum]` parity gates green.
- **The C60 metric is DP-side modeled peak (`cost_profile` / the DP frontier), never the replay `avoidable_time`** (untrusted per the session's trust boundary).
- **Do NOT implement the COMPLETE variant.** If S3.5 shows the cheap variant is insufficient, stop and escalate -- the complete variant is a separate, documented follow-on.
