# External-mode batching implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to execute this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax.

> **RETRACTED / CORRECTED (2026-07-27) -- contaminated "external-occ" data.**
> The D5 gate table below (and the "Gap 1 -- forest peak not bounded" conclusion
> built on it) was produced with occupied indices in the CONTRACTED-role
> batchability predicate (a single `is_batchable_index` was set;
> `is_batchable_external_index` was never supplied), so the DP ALSO batched
> CONTRACTED occ -- its cost cells sliced the contracted occ pair (`cell_slices`
> included `i_3 i_4`), a reduction the runtime never realizes. Every "external-occ"
> arm was thus "external-occ ON TOP OF contracted-occ", never "external-occ only".
> The tell is the IDENTICAL peak (2302.30 GB) reported for the "contracted-occ"
> and "external-occ" rows. Conclusions drawn from these numbers -- "forest peak
> not bounded", "external batching insufficient", "the peak floor is inevitable
> once external-occ pairs are batched", "needs forest-level / multi-mode
> co-batching" -- are unscoped and RETRACTED.
>
> The batchability role-split (SeQuant `79540c831`..`7f5240014`, MPQC
> `1685c7ac1e`..`5a3d94dc63`) confines occ to the EXTERNAL role. Honest re-measure
> of the C60 `[.][dryrun-occ-veto]` witness (nterms=55, post-role-split): aux+occ
> replay peak ~= **5860.9 GB** (NOT 2302 / 2947), avoidable_time ~= **44.6%**,
> `contracted_occ_stamps == 0`, `external_occ_stamps == 244`. The peak is HIGHER
> because the phantom contracted-occ reduction had HIDDEN the true cost of the
> cross-pair (two-PNO-leg) giants; external slicing cannot reduce them (the
> `a<i_3,i_4>` leg stays full). Full root cause:
> `.superpowers/sdd/contamination-role-predicate.md`. Historical content below is
> preserved as a record.

**Goal:** Make external-mode batching a first-class decision so the CSV/PNO-CCk
residual giants are bounded under a memory budget: the DP *sizes and selects* with
the external mode sliced, and the runtime *scatters* it — driving the C60 occ+aux
witness's modeled and replayed peak (and `avoidable_time`) down.

**Architecture (from `doc/dev/specs/2026-07-20-external-mode-batching-design.md`):**
The unified external-batch mechanism (emit + scatter) already exists (2026-07-11);
it is inert for these giants because (D1) the DP's sizing/selection ignore the
external slice (proven: modeled peak byte-identical flag-off/on), (D2) the runtime
cannot locate a proto-only external occ so the scatter no-ops, and (D3) the dry-run
backend cannot model the scatter. Fix order is **D1 first** — it is load-bearing
(the factorization must know true costs to schedule well; a runtime-only lever
cannot rescue a blind plan) and it is **dry-run-measurable alone** (predicted peaks
drop without the runtime). Then the **selection-policy experiment**, then **D2/D3**
(runtime), then **D4** (rename sweep) and **D5** (validation + MPQC).

**Tech stack:** C++20, SeQuant, Catch2, the dry-run cost-profile witness.

**Task types:** most tasks are TDD implementation with exact files/tests. Three are
explicitly **investigation** (D1.1, D2.1) or **experiment** (D1.3) — they carry a
*deliverable + verification* instead of forced TDD, because the design genuinely
depends on what they find. That is not a placeholder; it is the task.

## Global Constraints

- **Repo/branch:** `~/code/SeQuant`, branch `evaleev/feature/suppress-heuristic-fallback`.
- **Build/test (ONLY dir that compiles the eval tests):**
  `cd ~/code/SeQuant/cmake-build-debug && cmake --build . --target unit_tests-sequant -j6`,
  then `./tests/unit/unit_tests-sequant "<filter>"`. The dry-run witness
  `[.][dryrun-occ-veto]` is a hidden `[.]` test; select it explicitly. It is a
  *symbolic* dry-run, so Debug speed is fine.
- **Byte-identical fallback:** with `batch_external_modes` off (default; currently
  `batch_spectator_indices`), every path must stay byte-identical. Each task states
  its no-op fallback; do not regress the `[optimize]`/`[dryrun]`/eval suites.
- **Terminology:** the code still carries legacy `aux`/`spectator` identifiers;
  implement D1-D3 against the *current* names, then **D4** renames the touched files
  per `doc/dev/batching-mode-terminology.md` (do not interleave the rename into the
  functional diffs — it makes review harder). Spec/plan prose uses the target names;
  the terminology doc's map is the key.
- **Commit style:** plain messages, no `Co-Authored-By`/tool trailers. Pre-commit
  hook rejects en-dash U+2013 / U+00A0 (write ASCII; strip with
  `LC_ALL=C sed -i '' "s/$(printf '\342\200\223')/-/g" <files>`) and runs
  clang-format (re-`git add` + re-commit if it reformats). `git add` only the files
  a task changes; the repo has many untracked build dirs/`.DS_Store` — never
  `git add -A`.
- **No MPQC repin** until the SeQuant side validates (D5).
- Cap builds at `-j6`.

## Files

- `core/optimize/cost_model.hpp` — D1 (wire seeded external sizing into
  `select_root`; emit follows selection), D4 rename.
- `core/optimize/single_term_detail.hpp` — D1 plumbing (`seeded_root_peak_batched`
  / candidate selection), D4 rename.
- `core/eval/eval.hpp` — D2 proto-aware external locator, D4 rename.
- `core/eval/eval_expr.{hpp,cpp}` — D4 rename.
- `core/eval/backends/dryrun/result.hpp` — D3 `pre_sized_zeros_over_mode`.
- `tests/unit/test_optimize.cpp`, `tests/unit/test_eval_dryrun.cpp`,
  `tests/unit/test_eval_ta.cpp` — tests.

---

## Task D1.1 — Investigate + design the DP-selection wiring (INVESTIGATION)

**Files (read-only):** `core/optimize/cost_model.hpp` — `select_root` (~999-1087),
`seeded_root_peak_batched` (~1291-1341), the emit gate (~1136), `reconstruct_axes`
(the per-node emit walk, ~1122); `single_term_detail.hpp` (`sliced_footprints`,
`batchable_index_list`, `subset_open_aux`).

**Deliverable** (append a short "D1 wiring" note to the spec, or a plan-scratch
file): answer, with exact function/line references,
1. where `select_root` (and the perf-first ceiling) computes the **root peak** it
   uses for feasibility and reporting;
2. how to make that peak reflect **external seeding** for an over-budget term —
   i.e. how to route `seeded_root_peak_batched(seed_modes)`'s result into the
   feasibility/reported peak (it already rebuilds the sizing tables with the
   external mode's batchable predicate extended + a block-sized target, but is
   unused by `optimize()`/`select_root`);
3. how to make the **emit follow selection** — emit `External` for the modes the
   chosen feasible schedule actually relies on, replacing the current post-selection
   stamp gated on the *unsliced* `root_peak`.

- [ ] **Step 1:** Read the four functions and write the note (exact change points +
  the minimal wiring, plus the one initial policy to try in D1.2: "seed all the
  term's `is_external_mode`/`is_spectator` candidates").
- [ ] **Step 2:** Verify the note names concrete lines and does not require changing
  the contracted DP (external seeding is work-neutral — it must not perturb the
  min-flops factorization). Report DONE with the note path.

## Task D1.2 — Wire seeded external sizing into selection (simplest policy)

**Files:** `core/optimize/cost_model.hpp` (+ `single_term_detail.hpp` plumbing if
D1.1 finds it needed); `tests/unit/test_eval_dryrun.cpp`.

**Interface produced:** for a term whose unsliced root peak exceeds
`peak_threshold` with `batch_external_modes` on, the DP's *reported/feasibility*
peak is the **seeded** (external-sliced) peak, and `External` modes are emitted per
that selection. Off (default) => byte-identical.

**Initial policy (refined in D1.3):** seed *all* the term's external candidates
(the `is_spectator`/`is_external_mode` set), block size from
`batch_target_size(mode)`.

- [ ] **Step 1: Write the failing test.** In `test_eval_dryrun.cpp`, add a
  `[dryrun-extmode]` case (or extend the occ+aux witness path) that, with
  `batch_external_modes` on, asserts the C60 giant term's **DP-predicted** peak
  (from `cost_profile` / the DP's reported root peak) is strictly *less* than the
  flag-off value — and that with the flag off the predicted peak is byte-identical
  to baseline. (Baseline: the giant is modeled at ~3619 GB; flag-off must stay
  there, flag-on must drop.)
- [ ] **Step 2: Run to confirm it fails.**
  `./tests/unit/unit_tests-sequant "[dryrun-extmode]"` — fails (peak unchanged
  today, per the 2026-07-20 experiment).
- [ ] **Step 3: Implement the wiring** per D1.1's note: route
  `seeded_root_peak_batched(external candidates)` into the feasibility/reported peak
  for over-budget flag-on terms; make the emit follow selection. Keep the contracted
  DP untouched (work-neutral seed).
- [ ] **Step 4: Run to confirm it passes**, and confirm no regression:
  `./tests/unit/unit_tests-sequant "[optimize][dryrun-objective][dryrun-extmode]"`.
  Flag-off byte-identical.
- [ ] **Step 5: Commit** (against current identifier names; D4 renames later).

## Task D1.3 — Selection-policy experiment (EXPERIMENT)

**The real open problem (per the spec): which external modes to seed, how many, in
what order, at what block size.** Because D1.2 made this dry-run-measurable, settle
it on predicted peak/flops before any runtime work.

**Deliverable:** an evidence table (policy -> C60 giant predicted peak, total
predicted flops, #terms brought under budget) over these axes, and the chosen
policy wired into D1.2:
- candidate set: all external modes vs largest-footprint-first vs the term's
  external occ only;
- count: single outermost external seed vs jointly seeding both occ of a doubles
  residual;
- block size: sweep `batch_target_size`;
- order vs the contracted ceiling (external forest loop outermost).

- [ ] **Step 1:** Run the witness with each policy variant (env toggles or a small
  harness), collecting predicted peak/flops. `SEQUANT_UT_DRYRUN_*` knobs +
  `SEQUANT_DP_RECOMPUTE_DEBUG` as needed.
- [ ] **Step 2:** Pick the policy that best bounds predicted peak under budget with
  least added flops; wire it into D1.2 (update the code + keep `[dryrun-extmode]`
  green). Record the table in the plan-scratch / spec.
- [ ] **Step 3: Commit** the chosen policy.

## Task D2.1 — Confirm the runtime scatter no-op on a proto-only occ (INVESTIGATION)

**Deliverable:** direct evidence of the runtime behavior the experiment only
hinted at. With `External` emitted for the C60 proto occ, confirm whether the
scatter branch fires or no-ops, and why.

- [ ] **Step 1:** With D1 landed (external emitted for the giants), run the witness
  replay with tracing; check for `BatchScatter` markers and whether
  `pick_sliceable`/`find_leaf_carrying` locate the proto-only occ (add a temporary
  probe if needed). Confirm the hypothesis: `index_position` is not proto-aware, so
  the occ is skipped and the scatter no-ops.
- [ ] **Step 2:** Report the finding (fires / no-ops, with evidence). If it already
  fires (e.g. the occ is plain-locatable on leaves as in CSV), narrow D2.2's scope
  to only the genuinely proto-only case.

## Task D2.2 — Proto-aware external locator

**Files:** `core/eval/eval.hpp` (a proto-aware sibling to `index_position`, threaded
through `pick_sliceable`, `find_leaf_carrying`, the scatter branch's `dest_mode`,
the per-block leaf slicer, and `make_batched_scratch`'s external signature).

- [ ] **Step 1: Write the failing test.** In `test_eval_ta.cpp` (TA backend, where
  the scatter path is functional) or a dry-run case: a node/forest carrying a
  **proto-only** external occ (only as a subscript of a composite `a<i,j>`) — assert
  the scatter branch locates and slices it (result reconstructs the unbatched
  reference; peak scales ~block/extent). Include a **mixed leaf** (plain occ + proto
  legs) that does NOT carry the target, to guard the 2026-07-11 defect-1 over-throw.
- [ ] **Step 2: Run to confirm it fails** (proto-only occ un-locatable today).
- [ ] **Step 3: Implement** the proto-aware locator: position of the occ among a
  node's distinct composite outer-proto modes (reuse the canonical proto order; no
  `canonicalize_slots` rerun), feeding the existing ToT **outer-mode**
  `slice_mode`/`pre_sized_zeros_over_mode`/`write_into_slice`. Handle mixed leaves
  without over-throwing.
- [ ] **Step 4: Run to confirm it passes**, and the 2026-07-11 CSV / multi-mode /
  Hadamard external tests stay green:
  `./tests/unit/unit_tests-sequant "[eval]"` (or the specific tags).
- [ ] **Step 5: Commit.**

## Task D3.1 — Dry-run scatter modeling (`pre_sized_zeros_over_mode`)

**Files:** `core/eval/backends/dryrun/result.hpp` (`ResultDryRun` and
`ResultDryRunNested`); `tests/unit/test_eval_dryrun.cpp`.

- [ ] **Step 1: Write the failing test.** The witness occ+aux replay with external
  batching on currently throws (missing `pre_sized_zeros_over_mode`) or no-ops.
  Assert: the replay completes and the giant's per-op modeled size drops
  (external-sliced), and the measured replay peak / `avoidable_time` fall vs
  flag-off. (If the witness's avoidable-time parser keys only on
  `BatchGroup`/`BatchIter`, extend it to account for `BatchScatter`, or assert on
  peak directly — note which in the test.)
- [ ] **Step 2: Run to confirm it fails.**
- [ ] **Step 3: Implement** `pre_sized_zeros_over_mode` for both dry-run result
  classes: zero-data, widen the destination mode to the carrier's full (unsliced)
  extent (mirror the TA backend's outer-`TiledRange1` widening, in the dry-run's
  index/extent-override representation).
- [ ] **Step 4: Run to confirm it passes.** Full witness green; flag-off unchanged.
- [ ] **Step 5: Commit.**

## Task D4 — Terminology rename sweep (mechanical)

**Files:** the ones touched by D1-D3 (`cost_model.hpp`, `single_term_detail.hpp`,
`eval.hpp`, `eval_expr.{hpp,cpp}`, `backends/dryrun/result.hpp`) + their tests.

Rename per `doc/dev/batching-mode-terminology.md`: `ctx.aux`->`batchable_modes`,
`open_aux`->`open_modes`, `Acand`->`contracted_here`, `aprime`->`batched_here`,
`B`->`batched_enclosing`, `nbatch`->`nbatches`, `AxisKind`->`BatchModeType`,
`batch_axes()`/`set_batch_axes()`->`batched_here()`/`set_batched_here()`,
`is_spectator_axis`->`is_external_mode`,
`batch_spectator_indices`->`batch_external_modes`,
`reconstruct_axes`->`reconstruct_batched_modes`,
`batchable_index_list`->`batchable_mode_list`; "axis"->"mode", "spectator"->"external"
in comments/prose.

- [ ] **Step 1:** Mechanical rename (prefer `clang-rename`/IDE rename for
  identifiers with call sites across files; careful `sed` for the rest). Keep it a
  pure rename — no behavior change.
- [ ] **Step 2:** Build + run the full `[optimize]`/`[eval]`/`[dryrun]` suites;
  assert byte-identical behavior (identical assertion counts / outputs).
- [ ] **Step 3:** Add a one-line header pointer to `doc/dev/batching-mode-terminology.md`
  from `2026-07-07-*` and `2026-07-11-*` (legacy-term specs).
- [ ] **Step 4: Commit** as one reviewable rename. (The MPQC-side keyword rename
  `batch_spectator_indices`->`batch_external_modes` and the `cck` uses land in a
  paired MPQC commit at D5, with the repin.)

## Task D5 — Validation + MPQC enablement (final)

- [ ] **Step 1:** Full SeQuant suites green (targeted tags to avoid the known slow-DP
  hang): `[optimize]`, `[cache_manager]`, `[eval]`, `[dryrun-objective]`,
  `[dryrun-occ-veto]`, `[dryrun-extmode]`.
- [x] **Step 2 (GATE — RAN 2026-07-20; result: NECESSARY-BUT-NOT-SUFFICIENT).**
  Configured the external-batched replay (`batch_external_modes` on, occ block 8, D3
  scatter modeling) vs the contracted-occ baseline on the C60 residual forest (45
  terms, K@256, 100 GB, perf-first). The `[dryrun-extmode-avoidable]` witness in
  `tests/unit/test_eval_dryrun.cpp` records it. **The ~0 gate FAILED — and that RED
  reflects the physics, not a bug.** External batching *engages* (46 `BatchScatter`
  interceptions, 24 external stamps) and strictly improves both metrics, but does not
  eliminate the recompute or bound the peak:

  | config | replay ops | avoidable_time | scatter | modeled peak |
  |---|---|---|---|---|
  | contracted-occ | 92848 | 75.1% | 0 | 2302.30 GB |
  | external-occ | 66039 | **44.7%** | 46 | **2302.30 GB** |

  > **[RETRACTED 2026-07-27 -- see the top-of-doc note.]** The IDENTICAL
  > 2302.30 GB peak in both rows is the tell: the "external-occ" arm still batched
  > CONTRACTED occ (occ was in `is_batchable_index`, never
  > `is_batchable_external_index`), so this is "external-occ ON TOP OF
  > contracted-occ". Post-role-split honest re-measure: peak ~= 5860.9 GB,
  > avoidable ~= 44.6%, `contracted_occ_stamps == 0`, `external_occ_stamps == 244`.

  The test was **recast into a partial-win witness**: it asserts external strictly
  beats contracted (avoidable_time by >15 pts; fewer ops; scatter fires) and *records*
  the residual (`peak > 100 GB`, `avoidable > 0.10`) as the entry criterion for the
  deferred work — flip those two to `< 100`/`< 0.05` when it lands. TWO residual gaps,
  both deferred (see Follow-on): (1) the forest peak is set by a term slicing the
  proto-occ pair does not reach — needs **forest-level / multi-mode co-batching**; and
  (2) a **contracted middle-gap** node survives (the offender shifts from the gC
  intermediate to the 4-occ `I(..,i_3,K_2;a_3)*I(..,i_4,K_2;a_4)` contraction).
- [ ] **Step 3:** MPQC: enable the keyword on the CSV/PNO-CCk path (rename the
  keyword to `batch:external_...` if desired; keep back-compat), local CSV-CCSD
  sanity (converged energy unchanged vs unbatched), then **repin** SeQuant
  (`external/versions.cmake`, own commit).
- [ ] **Step 4:** Owl C60 occ+aux validation (separate session, SLURM only, refresh
  the TA pin): giants bounded, per-iteration time down, energy unchanged.

## Risks (carry into execution)

- **D1.1 may find the ceiling wiring is more than a peak swap** — e.g. the perf-first
  ceiling's fallback-to-global-min-flops path needs the seeded peak threaded into the
  feasibility test itself, not just the report. That is the load-bearing design
  moment; if it balloons, stop and escalate rather than force it.
- **Selection policy (D1.3)** is the real unknown; treat its table as the evidence
  that justifies the wired policy.
- **Proto-aware location on mixed leaves (D2.2)** is exactly the 2026-07-11 defect-1
  failure mode — the guard must not over-throw.
- **`avoidable_time` may not capture `BatchScatter`** — measure peak directly if so
  (D3.1).

## Follow-on after this plan: the residual middle gap (still needed)

External batching fixes the factorizer's pricing of EXTERNAL modes and bounds the
single dominant C60 giant W (D1's `[dryrun-extmode]` models it at ~8 GB). **The D5
gate proved that is not enough** (avoidable_time 75.1% -> 44.7%, but peak unchanged at
2302 GB). Two residual gaps remain — this is the fresh-planning-pass agenda:

> **[RETRACTED 2026-07-27 -- see the top-of-doc note.]** This "Gap 1" was inferred
> from the contaminated D5 table: the 2302 GB was an "external-occ ON TOP OF
> contracted-occ" figure, and the honest post-role-split peak is ~5860.9 GB. The
> true driver is the cross-pair two-PNO-leg giants that external slicing cannot
> reduce (the `a<i_3,i_4>` leg stays full) -- not an unaddressed "forest-level
> co-batching" gap. Re-derive any conclusion below against a role-separated config.

**Gap 1 — forest peak not bounded (NEW finding).** The 2302 GB peak is set by a term
that slicing the proto-occ pair (i_1,i_2) does not reach. Bounding it needs
**full forest-level / multi-mode co-batching** — co-batch i_3,i_4 and/or K_2, not just
one proto-occ pair — the "full forest-level co-batching" already flagged as required.
This gap is NOT addressed anywhere yet and needs its own design.

**Gap 2 — contracted middle-gap pricing + hoisting (already scoped).** The DP's
`charge_batch_recompute` (`esc = B & ~open_aux[n]`) is still ORDER-BLIND and prices a
node inside a *contracted* batch loop over an axis it does not carry at `rf = 1` when
the runtime rebuilds it per batch. Nor is the runtime multilevel hoisting built (the
Task-1 scope chain `parent_`/`set_parent`/fall-through `access` is committed on this
branch but INERT). The D5 witness's surviving offender — the 4-occ
`I(..,i_3,K_2;a_3)*I(..,i_4,K_2;a_4)` contraction — is exactly this case. Implement the
re-root-caused design in `mpqc4/doc/dev/specs/2026-07-17-nested-batch-group-join-design.md`:
(A) the ORDER-AWARE contracted cost (charge recompute by placement level), and
(B) MULTILEVEL scope-chain hoisting (store each intermediate at its placement level;
reuse the committed Task-1 fall-through primitive, do not rebuild). Supersedes the
stale depth-0 `2026-07-19-batched-hoisting-scope-chain.md` Tasks 3-6.
