# Phase 4a -- O2 shrink-first spill pass (standalone, non-regressing) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Naming:** the pass called **O2** here (the parent spec's ordinal) is named **`remat`** (rematerialization -- https://en.wikipedia.org/wiki/Rematerialization) in the code (`placement_remat.hpp`: `remat_cells`, `rematerialize_to_budget`, `RematResult`/`RematStatus`), because it reduces peak by REBUILDING values in smaller/shorter-lived pieces, not by spilling to storage. "O2" and "remat" are used interchangeably below.

**Goal:** Build the remat (spec "O2") pass, the whole-forest greedy placement-relaxation pass, as a STANDALONE static function `rematerialize_to_budget(cells, threshold) -> placement` that lowers the modeled peak by shrinking demoted giants -- validated on its own against the Phase-3b `PeakProfile`, with ZERO runtime wiring.

**Architecture:** O2 consumes the Phase-3b analysis (`peak_profile.hpp`: `home_scope` seed, `cell_footprint`, `peak_profile_sweep`). A cell is one value at a home mode-set. The seed homes each value at its meet `home_scope` (peak-maximal). O2 repeatedly, at the binding peak point, applies the cheapest peak-reducing move -- in this phase, only SHRINK: add a demoted carried batch-mode to a cell's home, block-sizing it at ZERO recompute -- until `peak <= threshold` or no shrink reduces the binding point (infeasible). Un-hoist/split (the evict moves) are DEFERRED (design section 3a: shrink first). Output = the final placement (per-value home mode-set) + its `PeakProfile`. Router-override emission + the runtime cutover are Phase 4b.

**Tech Stack:** C++20 header (`meta::eval_node_range`), the dry-run backend (`backends/dryrun/`), Catch2 (`tests/unit/`), CMake + Ninja.

## Global Constraints

- No en-dashes (U+2013) anywhere in source or comments.
- Format touched files with `/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`.
- No `Co-Authored-By` trailers in commits.
- RELEASE build in `cmake-build-release`; cap ninja at `-j6`.
- NON-REGRESSING: O2 is a NEW component with ZERO production callers; do NOT modify any runtime
  source (`eval.hpp`, `place_at_this_level`, `stamp_lifetime_masks`) or delete `contracted_modes` /
  `has_demoted_external` (all Phase 4b). The Phase-3b `peak_profile.hpp` flat `Schedule` output must
  stay BYTE-IDENTICAL where refactored (T1 guardrail). Existing `[peak_profile]`, `[dryrun]`,
  `[lifetime_mask]` suites + the `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]` witness
  figures stay green.
- SHRINK semantics (design section 3a): a shrink adds mode `m` to a cell's home where `m` is a
  DEMOTED carried batch mode -- `m` in the union of the cell's enclosing batch loops across
  occurrences, `m` in the cell's carried indices, `m` NOT already in the home. It block-sizes `m`
  (footprint via `cell_footprint`), leaves the liveness interval UNCHANGED, and adds ZERO recompute
  (result-mode slicing is free per the `W` model).

## File Structure

- CREATE `SeQuant/core/eval/placement_o2.hpp` -- the O2 pass: the rich per-value working record
  (`O2Cell`), the `Schedule` projection, the shrink move, the greedy loop, the placement result.
- MODIFY `SeQuant/core/eval/peak_profile.hpp` -- ONLY to expose the rich per-cell grouping that
  `linearize` already computes internally (carried / home_modes / enclosing-modes / interval), so O2
  reuses the ONE post-order walk instead of duplicating it; the flat `Schedule` stays a projection of
  it and byte-identical.
- CREATE `tests/unit/test_placement_o2.cpp` (add to `tests/unit/CMakeLists.txt`, same
  `eval_dryrun_test_sources` list `test_peak_profile.cpp` uses).

---

### Task T1 -- the O2 working representation (rich cells + Schedule projection + shrink candidates)

**Files:**
- Modify: `SeQuant/core/eval/peak_profile.hpp` (expose the rich per-cell records from `linearize`)
- Create: `SeQuant/core/eval/placement_o2.hpp` (`O2Cell`, `o2_cells`, `to_schedule`, shrink candidates)
- Create: `tests/unit/test_placement_o2.cpp` + register in `tests/unit/CMakeLists.txt`

**Interfaces:**
- Consumes: `detail::cell_footprint(carried, home_modes, cm, block_of)`,
  `struct Cell/Schedule`, `peak_profile_sweep` (all `peak_profile.hpp`); `home_scope(node)`,
  `stamp_seed_residency` (`lifetime_mask.hpp`); `dryrun::CostModel`.
- Produces:
  - `struct O2Cell { std::size_t value_id; container::svector<Index> carried;
    container::svector<Index> home_modes; container::svector<Index> enclosing_modes;
    std::size_t first_use; std::size_t last_use; };`
    (`home_modes` starts = the Phase-3b footprint home = `home_scope` MINUS the value's own-loop
    modes -- the SAME set `linearize` already feeds `cell_footprint`; `enclosing_modes` = the union
    over the value's occurrences of its enclosing `ectx` modes.)
  - `struct O2Input { container::svector<O2Cell> cells; std::size_t num_points = 0; };`
  - `template <meta::eval_node_range R> O2Input o2_cells(R const& forest,
    dryrun::CostModel const& cm, auto const& block_of);` -- `num_points` is the linearization's
    point count (same value `linearize` puts in `Schedule::num_points`), threaded through so the
    sweep can size its delta array.
  - `Schedule to_schedule(container::svector<O2Cell> const& cells, dryrun::CostModel const& cm,
    auto const& block_of, std::size_t num_points);` -- projects each `O2Cell` to a flat `Cell`
    (`footprint = cell_footprint(carried, home_modes, cm, block_of)`, interval copied).
  - `container::svector<Index> shrink_candidates(O2Cell const& c);` --
    `(enclosing_modes INTERSECT carried) MINUS home_modes`.

**Step 1: Write the failing test** (`test_placement_o2.cpp`, `[placement_o2]`)

```cpp
// A demoted value: carried {o_1,i_1}, meet-home {i_1} (o_1 dropped), o_1 enclosing in one occ.
// shrink_candidates == {o_1}; to_schedule footprint sizes o_1 FULL (matches the 3b demoted cell).
TEST_CASE("o2_cells exposes shrink candidates and projects to the 3b schedule",
          "[placement_o2]") {
  // build the SAME partial-demotion forest as test_peak_profile.cpp's partial-demotion case
  // (reuse its helper pattern); run o2_cells; find the V cell.
  // CHECK(shrink_candidates(V) == {o_1});
  // CHECK(to_schedule(cells,...).cells[V].footprint == full(o_1)*block(i_1)*8);  // == the 3b value
}
```

**Step 2: Run it to see it fail** -- `ninja -C cmake-build-release unit_tests-sequant -j6` then
`./cmake-build-release/tests/unit/unit_tests-sequant "[placement_o2]"` -- FAILS (no header).

**Step 3: Expose the rich records in `peak_profile.hpp`.** `linearize` already builds, per node, the
`hash / point / consumer_point / home(=home_scope) / carried / ectx / own_modes` records and groups by
hash. Refactor it so the grouping yields an intermediate rich per-value structure (value_id, carried,
`home_modes = home_scope MINUS own_modes_union`, `enclosing_modes = union of ectx.first over occ`,
first_use, last_use), and the existing flat `Schedule` becomes a projection of it. Keep `linearize`'s
public signature and its returned flat `Schedule` BYTE-IDENTICAL (guardrail below). Expose the rich
structure via a new `linearize_rich(...)` (or return both) that `placement_o2.hpp` consumes.

**Step 4: Implement `o2_cells` / `to_schedule` / `shrink_candidates`** in `placement_o2.hpp` on top of
the rich records. `o2_cells` = `linearize_rich` mapped to `O2Cell`. `to_schedule` re-projects (footprint
via `cell_footprint(carried, home_modes, cm, block_of)`). `shrink_candidates` = set difference/intersect.

**Step 5: Run the tests** -- `[placement_o2]` passes; and the guardrail:

```cpp
// GUARDRAIL in test_peak_profile.cpp (or _o2): linearize's flat Schedule is unchanged by the refactor.
TEST_CASE("linearize flat Schedule == to_schedule(o2_cells) on a forest", "[placement_o2]") {
  // same forest; assert peak_profile_sweep(linearize(f,...)) == peak_profile_sweep(to_schedule(o2_cells(f,...)))
  // (peak, binding_point, and per-cell footprints all equal) -- proves the rich split is faithful.
}
```
Also re-run `[peak_profile]` `[dryrun]` `[lifetime_mask]` -- green (refactor is behavior-preserving).

**Step 6: Commit** -- `git add SeQuant/core/eval/peak_profile.hpp SeQuant/core/eval/placement_o2.hpp
tests/unit/test_placement_o2.cpp tests/unit/CMakeLists.txt && git commit -m "eval: O2 working cells +
shrink candidates + Schedule projection (Phase 4a T1)"`.

---

### Task T2 -- the shrink move + the greedy spill loop + DeltaPeak + infeasibility

**Files:**
- Modify: `SeQuant/core/eval/placement_o2.hpp` (`apply_shrink`, `o2_spill`, the result struct)
- Modify: `tests/unit/test_placement_o2.cpp`

**Interfaces:**
- Consumes (T1): `O2Cell`, `o2_cells`, `to_schedule`, `shrink_candidates`; `peak_profile_sweep`.
- Produces:
  - `void apply_shrink(O2Cell& c, Index const& m);` -- appends `m` to `c.home_modes` (asserts `m` is a
    current shrink candidate). Footprint is re-derived by `to_schedule`; interval unchanged; recompute 0.
  - `enum class O2Status { Feasible, FactorizationInherent, RebatchNeeded };`
  - `struct O2Result { container::svector<O2Cell> cells; PeakProfile profile; O2Status status; };`
  - `O2Result o2_spill(container::svector<O2Cell> cells, dryrun::CostModel const& cm,
    auto const& block_of, std::size_t num_points, double threshold);`

**Step 1: Write the failing test** -- a forest whose seed peak exceeds a chosen threshold via one
demoted giant; assert `o2_spill(..., threshold)` returns `Feasible` with `profile.peak_bytes <=
threshold`, and that the binding giant's cell gained its demoted mode in `home_modes` (it was shrunk).
Add: `o2_spill(..., threshold = INF)` returns the seed unchanged (no move; `cells` home_modes
identical to `o2_cells`; `profile` == the seed sweep). Add: a forest whose binding cell has EMPTY
`shrink_candidates` (an irreducible giant) returns `RebatchNeeded` (not `Feasible`) and does NOT loop
forever.

**Step 2: Run to see it fail** (no `o2_spill`).

**Step 3: Implement the greedy** (design section 2/3b):

```cpp
inline O2Result o2_spill(container::svector<O2Cell> cells, dryrun::CostModel const& cm,
                         auto const& block_of, std::size_t num_points, double threshold) {
  for (;;) {
    Schedule s = to_schedule(cells, cm, block_of, num_points);
    PeakProfile pp = peak_profile_sweep(s);
    if (pp.peak_bytes <= threshold)
      return {std::move(cells), pp, O2Status::Feasible};
    // candidates: (cell alive at binding_point) x (its shrink_candidates). All shrinks are
    // zero-recompute, so pick the one maximizing DeltaPeak (re-sweep after the trial shrink).
    double best_drop = 0; std::size_t best_ci = 0; Index best_m;
    bool found = false;
    for (auto ci : pp.live_at_binding) {
      for (auto const& m : shrink_candidates(cells[ci])) {
        auto trial = cells; apply_shrink(trial[ci], m);
        double p2 = peak_profile_sweep(to_schedule(trial, cm, block_of, num_points)).peak_bytes;
        double drop = pp.peak_bytes - p2;
        if (drop > best_drop) { best_drop = drop; best_ci = ci; best_m = m; found = true; }
      }
    }
    if (!found || best_drop <= 0)
      // no placement move reduces the binding point: infeasible. A live cell with a nonempty
      // shrink set that still cannot drop the peak => factorization-inherent; none shrinkable at
      // all => rebatch-needed (the giant needs a batch loop the factorizer never added).
      return {std::move(cells), pp, classify_infeasible(cells, pp)};
    apply_shrink(cells[best_ci], best_m);  // full re-sweep next iteration (incremental = deferred)
  }
}
```
Implement `classify_infeasible`: `RebatchNeeded` if EVERY cell live at the binding point has empty
`shrink_candidates`; else `FactorizationInherent`. DeltaRecompute is 0 for every shrink, so the metric
is pure DeltaPeak here (design section 3b); the `avoidable_flops` term enters with the evict moves
(deferred).

**Step 4: Run the tests** -- feasible/inf/threshold-INF all pass; the loop terminates (each accepted
shrink strictly grows some `home_modes` toward the finite carried set, so the move supply is finite).

**Step 5: Commit** -- `git commit -m "eval: O2 shrink greedy spill loop + infeasibility classify
(Phase 4a T2)"`.

---

### Task T3 -- witness validation (O2 lowers the modeled peak; threshold=INF is a no-op)

**Files:**
- Modify: `tests/unit/test_placement_o2.cpp`

**Interfaces:** Consumes `o2_cells`, `o2_spill`, `peak_profile_sweep`; the witness forest builders in
`test_eval_dryrun.cpp` / the `[.][dryrun-occ-veto]` fixtures (read them for the exact forest + regime).

**Step 1: Write the validation tests** (`[placement_o2]`, and a hidden `[.][o2-witness]` for the slow
one):
- On the occ-veto witness forest + its faithful regime: `pp_seed = peak_profile_sweep(to_schedule(
  o2_cells(...)))`; pick a threshold BELOW `pp_seed.peak_bytes` (e.g. the documented target class) and
  assert `o2_spill(..., threshold).profile.peak_bytes <= threshold` with `status == Feasible`
  (O2 shrinks the giant under budget) -- OR, if it reports infeasible, assert the status and capture
  both numbers (an honest result, per design section 4; do NOT force feasibility by loosening the
  sizing). State which occurred in the test comment.
- `threshold = INF` is a strict no-op: `o2_spill(...).cells` home_modes are identical to `o2_cells(...)`
  and `profile.peak_bytes == pp_seed.peak_bytes`.
- A small CONSTRUCTED forest with a known demoted giant: hand-compute the post-shrink peak and assert
  `o2_spill` reaches exactly it.

**Step 2: Run** -- `[placement_o2]` green; the slow `[.][o2-witness]` foreground. Confirm the runtime
witnesses `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]` are UNCHANGED (O2 has no runtime
path; this is the non-regression proof).

**Step 3: Commit** -- `git commit -m "test: O2 lowers the witness modeled peak under threshold
(Phase 4a T3)"`.

---

## Self-review checklist (run after writing, before dispatch)

- Every task ends with a committed, independently testable deliverable. YES (T1 cells/projection,
  T2 the spill loop, T3 witness validation).
- Types consistent across tasks: `O2Cell`, `o2_cells`, `to_schedule`, `shrink_candidates`,
  `apply_shrink`, `o2_spill`, `O2Result`, `O2Status`. Used identically in T1/T2/T3.
- No placeholders: the greedy body + shrink semantics are shown in full; the infeasibility classifier
  is specified.
- Shrink correctness: only DEMOTED carried batch modes are candidates (T1 `shrink_candidates`);
  footprint via the tested `cell_footprint`; interval unchanged; recompute 0 -- matches design 3a and
  the Phase-3b own-loop-FULL fix (own modes are excluded from `home_modes` and from `enclosing_modes`).

## Deferred (NOT this plan)

- **Phase 4a-evict** -- the un-hoist + split moves (they change lifetimes / add recompute, so they need
  the `avoidable_flops` term in the metric and interval updates in the sweep). Design section 3a marks
  them lower-priority; add them if shrink alone does not bring the witnesses under target.
- **O2a incremental re-sweep** -- this plan re-sweeps fully after each move; optimize only if
  `moves x forest` is too slow (design section 3d).
- **Phase 4b** -- the CUTOVER: emit router overrides from the O2 placement, wire `place_at_this_level`
  to consult them, drop the `External` filter, DELETE `contracted_modes` + `has_demoted_external`, the
  default-threshold policy, and the CCk-correctness validation (design sections 2/5/7).
