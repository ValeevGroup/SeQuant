# Home-scope placement — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the placement-as-register-allocation design
(`doc/dev/specs/2026-08-02-home-scope-placement-design.md`), replacing the
caching-veto framing with correct cell placement under a peak budget.

**Architecture:** Five phases, each a working, testable increment; each later
phase gets its own detailed plan once the prior lands. Phase 1 (this plan,
detailed) corrects the peak *measurement* to the true co-resident sum — the
oracle everything else is validated against; Phases 2-5 are a roadmap here.

**Tech Stack:** C++20, SeQuant dry-run backend (`cost_profile()`), the real eval
loop (`eval.hpp`), `cache_manager.hpp`; Catch2 unit tests; Boost rational (for
`W` in a later phase).

## Global Constraints

- **No en-dashes** (U+2013) anywhere — the repo has a `forbid-en-dashes` commit
  hook; use `--` or `-`. Em-dashes are fine.
- **clang-format** touched C/C++ with `/opt/homebrew/opt/llvm/bin/clang-format
  --style=file -i` before committing.
- **No `Co-Authored-By` trailers**; keep dependency-pin bumps in their own commit.
- **Perfect-CSE default must stay byte-identical**: any keying/router default
  must reproduce today's behavior when no split/placement decision is made.
- **The replay is the ground-truth peak oracle**; the static model (Phase 3+)
  must be validated against it, never the reverse.

---

## Roadmap (phase = future plan unless detailed below)

- **Phase 1 (detailed here) — peak oracle.** Correct `cost_profile()`'s
  `peak_bytes` from `max(scratch_hwmark, cache_hwmark)` to the instant-resolved
  **sum** of co-resident residency across the live scope chain (spec §7c / O3b).
  Fixes the documented lower-bound under-count; nothing else depends on the
  router, so it stands alone.
- **Phase 2 — home_scope + router seed** (spec §7a, §7d / O1, O5). Compute
  `home_scope(value) = deepest scope over (sliced_modes ∪ demoted_external_modes)`
  and materialize the explicit `{value, use-site} -> (home, split=0)` router as
  the perfect-CSE seed; reads route via the router and reuse the existing
  `(use - home) INTERSECT carried` slicer instead of the `parent_` walk +
  `hops`. Byte-identical (seed == perfect CSE), but replaces the split authority.
- **Phase 3 — static peak profile** (spec §7c / O3a). The weighted-interval
  sweep (`[first-use, last-use]` per cell x footprint) over the router+schedule,
  validated to equal the Phase-1 replay oracle.
- **Phase 4 — the O2 greedy** (spec §7b). Spill loop: shrink/evict cells alive at
  the binding peak point by ΔPeak/ΔRecompute against the Phase-3 profile,
  emitting a router with non-trivial homes/splits.
- **Phase 5 — feedback** (spec §7b, O6). Detect-and-report the two infeasible
  modes (factorization-inherent vs. re-batch-needed); optional re-batch hint.

---

## Phase 1: correct the replay peak to the co-resident sum

**Why.** `cost_profile()`'s `peak_bytes` is
`max(profile.peak_bytes, peak.load(), cache.working_set_hwmark())`
(`cost_profile.hpp` replay loop) — the max of two *independent* high-watermarks
(the batched-scratch PeakSink and the outer cache), so when a persistent
cross-term cell co-resides with a batched-inner transient the two are additive
and the max under-reports (the field's own doc calls it a lower bound). The true
peak is `max over instants of (sum of all live residency across the scope
chain)`. Note `max(hwmark_a, hwmark_b)` is a lower bound and `hwmark_a + hwmark_b`
is an upper bound (the two peaks may occur at different instants); the correct
value is the sum taken *per instant*, then maxed — which requires summing the
scope chain at each `note_working_set` point, not combining two end-of-run
hwmarks.

**Files:**
- Modify: `SeQuant/core/eval/cache_manager.hpp` (add per-cache current-residency
  + scope-chain sum; feed the chain sum into the hwmark)
- Modify: `SeQuant/core/eval/eval.hpp` (the `note_working_set` call sites pass
  the chain-summed working set)
- Modify: `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` (`peak_bytes` doc;
  the fold no longer needs the separate `working_set_hwmark()` term once the
  chain sum subsumes it)
- Test: `tests/unit/test_cache_manager.cpp` (chain-residency unit), and
  `tests/unit/test_eval_dryrun.cpp` (the `[dryrun][peak]` / `[cost_profile]`
  co-resident assertion + witness re-baseline)

**Interfaces:**
- Produces: `size_t CacheManager::current_residency() const` — sum over *alive*
  entries of `entry.size_in_bytes()` in this cache only.
- Produces: `size_t CacheManager::chain_residency() const` — `current_residency()
  + (parent_ ? parent_->chain_residency() : 0)`.
- Consumes (Phase 3): the corrected `CostProfile::peak_bytes` as the oracle.

- [ ] **Step 0: Read the current peak path.** Read `cache_manager.hpp` lines
  ~185-305 (`working_set_hwmark_`, `note_working_set`, `reset`) and
  `eval.hpp:576-585, 1867-1874` (the `note_working_set(hwmark)` sites; `hwmark =
  log::bytes(cache, post).value`) and `cost_profile.hpp` replay fold
  (`profile.peak_bytes = std::max({..., peak.load(),
  cache.working_set_hwmark()})`). Confirm: each `CacheManager` tracks *its own*
  hwmark; the scratch PeakSink folds the scratch hwmark; the outer cache hwmark
  is separate; they are combined by `max` at the end.

- [ ] **Step 1: Write the failing test for `current_residency()`.**
  In `tests/unit/test_cache_manager.cpp`, in a new `TEST_CASE("cache_manager
  residency", "[cache_manager]")`: build a manager over two nodes, `store()` data
  into both, and assert
  `man.current_residency() == man.entry_size_in_bytes(k0) +
  man.entry_size_in_bytes(k1)` while both are alive, and that it drops as entries
  are drained (mirroring the existing `alive and entry_size_in_bytes` section).

- [ ] **Step 2: Run it, verify it fails** (`current_residency` undefined).
  Run: `ninja -C cmake-build-release unit_tests-sequant -j6 && \
  ./cmake-build-release/tests/unit/unit_tests-sequant "cache_manager residency"`.
  Expected: compile error / FAIL.

- [ ] **Step 3: Implement `current_residency()` + `chain_residency()`** in
  `cache_manager.hpp` (public):
  ```cpp
  [[nodiscard]] size_t current_residency() const noexcept {
    size_t s = 0;
    for (auto const& [k, e] : cache_map_)
      if (e.alive()) s += e.size_in_bytes();
    return s;
  }
  [[nodiscard]] size_t chain_residency() const noexcept {
    return current_residency() + (parent_ ? parent_->chain_residency() : 0);
  }
  ```

- [ ] **Step 4: Run the test, verify it passes.** Same command as Step 2.
  Expected: PASS.

- [ ] **Step 5: Commit.**
  ```bash
  git add SeQuant/core/eval/cache_manager.hpp tests/unit/test_cache_manager.cpp
  git commit -m "eval: add CacheManager current/chain residency accessors"
  ```

- [ ] **Step 6: Thread the chain sum into the hwmark.** In `eval.hpp`, at each
  `note_working_set(hwmark)` site (576-585 and the others found in Step 0),
  replace the local `hwmark` with the scope-chain total so the recorded peak is
  the co-resident sum. Concretely, where today `hwmark = log::bytes(cache,
  post).value` (this cache's working set), add the ancestors:
  `hwmark += (cache.parent() ? cache.parent()->chain_residency() : 0);`
  before the `note_working_set(hwmark)` call. (The local `log::bytes(cache,post)`
  already covers this cache; `chain_residency()` on the parent adds the outer
  live residency at the same instant.) In the scatter branch (1867-1874), fold
  `bs.cache.chain_residency()` instead of `bs.cache.working_set_hwmark()` into
  the PeakSink so the scratch peak already includes its ancestors.

- [ ] **Step 7: Simplify the `cost_profile()` fold.** In `cost_profile.hpp`,
  since the per-op hwmark now already sums the chain, drop the separate
  `cache.working_set_hwmark()` term from
  `profile.peak_bytes = std::max({profile.peak_bytes, peak.load(),
  double(cache.working_set_hwmark())})`, leaving
  `std::max(profile.peak_bytes, peak.load())` — `peak` (the PeakSink) now carries
  the summed instant maximum. Update the `peak_bytes` doc comment (lines ~189-201)
  to say it is the true co-resident peak, not a `max`-based lower bound.

- [ ] **Step 8: Write the failing co-resident test.** In `test_eval_dryrun.cpp`,
  adapt the `[dryrun][cost_profile]` fixture (or the `[dryrun][peak]` one) into a
  case where a persistent cross-term intermediate is alive while a batched-inner
  transient is alive: assert the new `cp.peak_bytes` equals the *sum* of the two
  live footprints at the co-resident instant, which is strictly greater than the
  old `max` of the two separate hwmarks. (Reuse the existing giant-term forest;
  the persistent gC-class intermediate co-resides with the batched W scratch.)
  Compute the expected sum from the two nodes' `memsize` at their realized
  extents, as the `[dryrun][cache]` test already sizes the giant.

- [ ] **Step 9: Run; verify PASS after Steps 6-7** (the summed peak is now
  reported). Run the `[dryrun]` filter.

- [ ] **Step 10: Re-baseline the witnesses.** The corrected (higher) peak shifts
  the documented `[.][dryrun-occ-veto]` / `[.][dryrun-extmode-avoidable]` peak
  figures (e.g. `~443.6 GB`, `~6026 GB`). Run each at nterms=55, record the new
  peaks, and update the doc figures + the `documented-RED research target (~X GB)`
  comments. These are diagnostics, not gates (the `peak < 100` CHECK stays RED);
  do **not** invent numbers — measure. Note in each that the peak rose because it
  is now the true co-resident sum, not the old lower bound.

- [ ] **Step 11: Full suite green.** `ctest --test-dir cmake-build-release` for
  the eval/cache/dryrun tests (`[cache_manager]`, `[lifetime_mask]`, `[dryrun]`,
  `[eval]`); the `[.]` witnesses' *acceptance* gates (role stamps, `abs < 0.02`)
  must stay green, only the documented-RED peak/avoidable targets move.

- [ ] **Step 12: clang-format + commit.**
  ```bash
  /opt/homebrew/opt/llvm/bin/clang-format --style=file -i \
    SeQuant/core/eval/eval.hpp SeQuant/core/eval/cache_manager.hpp \
    SeQuant/core/eval/backends/dryrun/cost_profile.hpp tests/unit/test_eval_dryrun.cpp
  git add -A && git commit -m "eval: peak_bytes = true co-resident sum, not max(scratch,cache)"
  ```

## Self-review checklist (run after Phase 1)

- The `max(a,b) -> per-instant sum` change is correct: the summed value is
  recorded *at each op* (an instant), then maxed, so it is neither the
  lower-bound `max(hwmarks)` nor the upper-bound `sum(hwmarks)`.
- No `[.]` acceptance gate (`contracted_occ == 0`, `external > 0`,
  `abs(ext-con) < 0.02`) moved; only documented-RED peak targets re-baselined
  from *measurement*.
- Perfect-CSE / OFF path unchanged: with no batching every scratch is absent, so
  `chain_residency()` == the single cache's residency == today's value.

## Notes for later phases (not tasks yet)

- Phase 2's router key is `{value, use-site}`; the use-site is the eval-tree node
  occurrence (pointer/position) — O5b / O1 residuals.
- Phase 3's sweep is validated to equal this phase's replay `peak_bytes` on the
  same forests.
- Phase 4 needs Phase 3's per-point live set (the sweep's argmax) for the spill
  candidates.
