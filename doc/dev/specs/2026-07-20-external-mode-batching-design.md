# External-mode batching as a first-class cost-model decision

**Status:** design
**Date:** 2026-07-20
**Terminology:** uses `doc/dev/batching-mode-terminology.md` (mode, batch vs
slice, `BatchModeType {Contracted, External}`, node-local vs forest-level).
**Identifier names below are the target (post-rename) names** — e.g.
`reconstruct_batched_modes`, `batch_external_modes`, `is_external_mode`,
`batchable_mode_list`, `contracted_here`. The current (pre-rename) code still
carries the legacy `aux`/`spectator` identifiers; that doc's map is the
translation key, and the rename is folded into this work (D4).
**Relation to prior work:** the unified external-batch mechanism
(`2026-07-11-unified-spectator-batch-axis-design.md`, *implemented* in commits
`d896b9841..874bfb43f`) landed the *emit* and the *scatter* branch. This spec
**corrects its central over-claim** — that "the peak model already accounts for
the sliced mode" — which the 2026-07-20 experiment disproves, and completes the
two holes that make external batching inert for the CSV/PNO-CCk giants.

---

## Problem

External modes — result-present modes (free or Hadamard), which for CSV/PNO-CCk
are the external occupied indices carried as protos of the PNO composite
`a<i,j>` — dominate the peak of the CC-residual giants (C60: individual summands
modeled at ~3.6 TB, far over a 100 GB budget). Bounding them requires batching
those external modes (batch the occ into blocks, evaluate per block, scatter into
the result). The unified mechanism emits `External` mode annotations and has a
runtime scatter branch, but **two holes make it inert** for these giants.

### Evidence (experiment, 2026-07-20; `experiment-enable-external-batching.md`)

Enabling `batch_external_modes` on the C60 occ+aux dry-run witness:

- **Emit fires** — `emit_external=1` for the 6 over-budget giant terms; 77
  `External` modes stamped, all on the external occ protos (40x `i_1`, 37x
  `i_2`). The candidate/emit path works.
- **But the DP is blind to the slice** — the modeled peak is *byte-identical*
  flag-off vs flag-on (giant `3619.8298417116248 GB -> 3619.8298417116248 GB`,
  full floating-point precision, three separate terms, flops identical).
  `select_root` runs on the unsliced sizing; the emit is computed *after*
  selection and only appends tags.
- **And the runtime does not scatter** — the replay raised *zero* scatter throws
  (the missing dry-run `pre_sized_zeros_over_mode` was never reached), consistent
  with the proto-only external occ being un-locatable at runtime (`index_position`
  is not proto-aware).

So the emit is a **stamp nothing acts on**: neither the DP's sizing/selection nor
the runtime scatter. This is precisely the 2026-07-11 over-claim, now measured.

## Current state: implemented vs. the holes

**Implemented (2026-07-11):**
- `batchable_modes` (legacy `ctx.aux`) admits external occ, including proto-only
  occ (`batchable_mode_list` proto + open-external passes).
- `reconstruct_batched_modes` emits `BatchModeType::External` per node, gated
  `batch_external_modes && perf_first && root_peak(unsliced) > peak_threshold`.
- The evaluator's scatter branch (slice carrying operands, `pre_sized_zeros_over_mode`,
  `write_into_slice`) and external modes in `make_batched_scratch`'s sharing signature.
- `seeded_root_peak_batched` (cost_model.hpp) *can* size an external slice — it
  rebuilds the sizing tables with the external mode's batchable predicate extended
  and a block-sized target — **but it is unused by `optimize()`/`select_root`.**

**The holes:**
1. **DP sizing/selection ignore external slicing.** The emit is post-selection;
   the modeled peak/feasibility never reflect the slice. (proven)
2. **Runtime cannot locate a proto-only external occ** — `index_position` /
   `find_leaf_carrying` are not proto-aware, so `pick_sliceable` skips the mode
   and the scatter no-ops. (strongly indicated; confirm as step 1 of D2)
3. **Dry-run cannot model scatter** (no `pre_sized_zeros_over_mode`), so the
   witness cannot measure external batching at all.

Both are needed, and **hole 1 (the DP) is fixed first**, for two reasons. First,
it is load-bearing: the factorization must know the true costs to drive the
schedule correctly — a runtime-only fix (D2) merely *executes* whatever schedule
the DP chose, so a DP blind to the external slice will keep choosing badly and no
runtime lever can recover a good schedule from a bad plan. Second, D1 is
**immediately measurable by dry-run alone**: it changes the *modeled/predicted*
peak (what `cost_profile` reports from the DP's decisions), so we can dry-run and
watch the predicted peaks drop **without** the runtime scatter (D2). D2 then makes
the runtime actually execute the schedule D1 chose (and bound the real replay
peak).

## Design

Make external-mode batching one coherent decision across the DP and the runtime,
per the two batching regimes: **node-local contracted** (already correct) and
**forest-level external** (the subject here).

### D1 — DP: size and select with the external slice (forest-level seed)

An external mode is batched **forest-level**, not node-local: it is carried onto
the result of every node up to the root, so batching it is a single outer
block-loop over the whole tree — work-neutral (no flops or recompute change),
scaling every carrying node's footprint by ~`block/extent`. It is therefore
represented as a **seed into the root batch context**, not as a per-node
`contracted_here` (legacy `Acand`) batch.

The sizing already exists — `seeded_root_peak_batched` yields the giant sized
sliced. Wire it into selection:

- **Feasibility under the perf-first ceiling.** For a term whose selected peak
  exceeds `peak_threshold` with external batching enabled, evaluate feasibility
  with the external mode(s) *seeded* (sliced). If the seeded peak fits the
  budget, the term is feasible via external batching.
- **Emit follows selection.** Emit `External` for a mode iff the chosen, feasible
  schedule relies on batching it — a selection outcome, not a fixed post-hoc
  stamp. The term's reported peak is the seeded (sliced) footprint.
- **No perturbation of the factorization.** External seeding is work-neutral, so
  it does not change the min-flops factorization choice — only the peak /
  feasibility. Sequencing: choose the min-flops factorization (contracted DP
  untouched); if over budget, apply the external seed to test feasibility, drive
  the emit, and report the sliced peak. Contracted batching (node-local) nests
  *inside* the external forest loop.

**External-mode selection is a genuine open problem — to be settled by
experiment, not assumed here.** We do not yet know *which* external modes to seed,
how many, in what order, or at what block size: candidates come from the
`is_external_mode` set, but which subset actually yields a good peak/flops
schedule (a single outermost external seed vs. jointly seeding both occ of a
doubles residual; interaction with the block size `batch_target_size`; ordering
against the contracted ceiling) is unknown. The virtue of D1-first is exactly that
it makes this cheap to explore: because D1 is dry-run-measurable, the plan drives
the selection choice from **predicted peak/flops on the witness**, iterating on
the selection policy before any runtime work. The spec fixes the *mechanism*
(seed external modes into the root batch context, feasibility gates the emit);
the *selection policy* is an experiment the plan carries.

### D2 — Runtime: scatter must fire on a proto-only external occ

The scatter branch resolves the external mode with `index_position` on
`canon_indices()`. A **proto-only** external occ (carried only as a subscript of a
composite `a<i,j>`) has no top-level slot, so `index_position` returns null,
`pick_sliceable` skips it, and the scatter never runs. Add a **proto-aware
locator**: the position of the occ among a node's distinct composite outer-proto
modes (reusing the canonical proto order), feeding the existing ToT **outer-mode**
`slice_mode` / `pre_sized_zeros_over_mode` / `write_into_slice`. Thread it through
`pick_sliceable`, `find_leaf_carrying`, the scatter branch's `dest_mode`, the
per-block leaf slicer, and `make_batched_scratch`'s external signature.

Step 1 of D2 is to **confirm the current no-op** on the C60 proto occ directly
(the experiment stopped before this). The 2026-07-11 forest path had
`spectator_outer_mode` / `outer_proto_position` for exactly this, retired on
unification; the unified per-node path needs the equivalent, now consistent with
the mainline mechanism — and it must handle **mixed leaves** (a plain occ leg plus
composite protos) without the old guard's over-throw (2026-07-11 defect 1).

### D3 — Dry-run: model scatter (so the witness can measure)

Add `pre_sized_zeros_over_mode` to the dry-run `Result` (flat `ResultDryRun` and
nested `ResultDryRunNested`), widening the destination mode to full (unsliced)
extent — zero-data, mirroring how the TA backend widens one outer `TiledRange1`.
Then the witness replay exercises the scatter branch. Ensure the witness's
avoidable-time / peak accounting handles the scatter markers (`BatchScatter`), or
measures the resident peak directly.

### D4 — Terminology rename (folded in)

As D1-D3 touch `cost_model.hpp`, `single_term_detail.hpp`, `eval_expr.*`,
`eval.hpp`, and the dry-run backend, apply the rename from
`doc/dev/batching-mode-terminology.md` (`aux`->`batchable_modes`,
`Acand`->`contracted_here`, `aprime`->`batched_here`, `B`->`batched_enclosing`,
`AxisKind`->`BatchModeType`, `nbatch`->`nbatches`, `batch_axes()`->`batched_here()`),
plus a one-line pointer from the old design docs to the terminology reference.

## Correctness

- **Work-neutral partition.** An external mode is contracted nowhere, so its
  blocks are disjoint slices of the result that reconstruct it exactly with no
  cross-block coupling and unchanged flops — free and Hadamard alike (Hadamard
  slices both operands to the same block; free slices one). Same identity the
  2026-07-11 mechanism already relies on.
- **Seeded sizing matches the runtime footprint.** The DP's seeded peak (giant
  sliced over the external mode) must equal the runtime scatter's per-block
  resident footprint; a validation check compares modeled-seeded vs measured
  peak.
- **Precision invariance.** External batching is a memory *schedule*, never an
  approximation; the converged energy must not change.

## Validation

1. **DP gate (fixes hole 1).** With external batching on, the modeled peak of the
   C60 giant *drops* under external seeding (was byte-identical) and drives
   `select_root` feasibility.
2. **Runtime gate (fixes holes 2+3).** The witness replay scatters the external
   occ (`BatchScatter` fires); the giants are sliced; measured peak and
   `avoidable_time` drop.
3. **Parity, external batching off** — byte-identical (no external emit), all existing
   `[optimize]`/`[dryrun]`/eval tests unchanged.
4. **CSV-CCSD (8-water)** — the 2026-07-11 validations still pass; the proto-aware
   locator does not regress the plain-occ path or the multi-mode / Hadamard tests.
5. **C60 occ+aux end-to-end** (real TA, separate session) — giants bounded under
   budget; converged energy unchanged.

## Risks / open questions

- **Selection policy (the real unknown).** Which external modes to seed, how many,
  in what order, at what block size — settled by dry-run experiment under D1, not
  by this spec (see D1). This, not the mechanism, is where the design risk lives.
- **Ceiling interaction.** Does the perf-first feasibility-with-external-seed
  compose cleanly with contracted (node-local) batching nesting inside the
  external forest loop?
- **Proto-aware location on mixed leaves** — the exact failure mode of 2026-07-11
  defect 1; must not over-throw.
- **Intermediate result pre-sizing** over the external mode's full outer extent —
  confirm that outer extent is not itself a bottleneck on any target system.
- **Residual middle-gap (out of scope here, but the interaction to watch).**
  External batching is work-neutral, so it introduces no recompute; it also
  *removes* the need for the contracted-occ batching that produced the original
  ~76% avoidable recomputation, by bounding the giants (e.g. the PPL `W`) so they
  need not be recomputed. If, after external batching, any avoidable_time survives
  (the D5 gate measures this), it is a **middle-gap** term — an intermediate inside
  a batch loop over an axis it does not carry — whose design (scope-chain fallback
  lookup + lifetime-scope + order-aware placement) is
  `mpqc4/doc/dev/specs/2026-07-17-nested-batch-group-join-design.md`, NOT this spec.

## Files touched (anticipated)

SeQuant:
- `core/optimize/cost_model.hpp` — wire seeded external sizing into `select_root`
  feasibility; emit follows selection; rename.
- `core/optimize/single_term_detail.hpp` — seeded-sizing plumbing; rename.
- `core/eval/eval.hpp` — proto-aware external locator across `pick_sliceable` /
  `find_leaf_carrying` / scatter branch / `make_batched_scratch`; rename.
- `core/eval/eval_expr.{hpp,cpp}` — rename (`batched_here()` etc.).
- `core/eval/backends/dryrun/result.hpp` — `pre_sized_zeros_over_mode`.
- tests: `test_optimize.cpp`, `test_eval_dryrun.cpp` (DP gate + runtime gate);
  keep the 2026-07-11 CSV/multi-mode/Hadamard cases green.

MPQC:
- enable `batch_external_modes` on the CSV/PNO-CCk occ-batching path (keyword
  already exists) once the SeQuant side is validated; repin.
