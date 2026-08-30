# Batched-DAG loop identity — implementation plan

> **For agentic workers:** Implement task-by-task. Each task ends with an
> independently-checkable validation. Tasks that touch internals begin with a
> scoped READ that produces a micro-design **before** any edit — do not fabricate
> diffs for code whose current shape you have not traced. Steps use checkbox
> (`- [ ]`) syntax.

**Goal:** Give batch loops a complete identity `(depth, loop_slot)` and track the
per-occurrence mode↔loop map through the pipeline, so a value carrying more than
one batched mode of one space is scheduled and executed correctly — fixing the
w20 aux+occ `is_range_set_congruent` crash without regressing the working cases.

**Architecture:** A loop's **identity** is `(depth, loop_slot)` (group + member),
its **layout** is `(altitude_ordinal, latitude_ordinal)` (nesting rank + PROCON
pass). `loop_slot` (the missing piece) is assigned at fusion by unifying
occurrences over producer→consumer slot connectivity; `build_ordered_schedule`
realizes one loop per slot (no collapse) and assigns layout; value-id and
occurrence-id color sliced modes by `space + (depth, loop_slot)`. Roles
(`per_axis`/`classify_axis`) are unchanged.

**Tech stack:** C++20, SeQuant `core/eval` + `core/optimize`. Validation:
`[ordered-schedule]` Catch2 tests; the MPQC `mpqc` binary on
`repro/w8-auxocc-ordered-*.json` (fast, local) and w20 aux+occ (remote z820); the
env-guarded dumps `SEQUANT_DUMP_PER_AXIS` / `SEQUANT_DUMP_SCHEDULE` already in the
tree.

**Spec:** `doc/dev/specs/2026-08-28-batched-dag-loop-identity-design.md` (read it
first; this plan argues from it).

## Global constraints

- **Frame purity.** No cross-frame index-**label** match. Cross-frame **space**
  comparison only as an abstract color (fusion group-matching); never a concrete
  space name in scheduling/identity logic.
- **Identity vs layout.** Seam, value-id, occurrence-id, atlas key on
  `(depth, loop_slot)` — never `altitude_ordinal`/`latitude_ordinal`.
- **Losslessness is the bar**: batched energy must match the unbatched reference
  to full precision, not merely "no crash". Reference energies:
  `w8-auxocc-fd-7b.json` (forest-descent) for the w8 configs.
- **No role changes** (`classify_axis`/`per_axis` kept as-is — measured correct).
- Build only `mpqc`, `-j6`; serialize mpqc runs. clang-format touched files
  (`/opt/homebrew/opt/llvm/bin/clang-format --style=file -i`; `.tpp` excluded).
  No en-dashes. No Co-Authored-By trailers. SeQuant→MPQC repin is its own commit.
- Local build dir consuming `~/code/SeQuant` = `cmake-build-release`. w20 on
  `ssh -i ~/.ssh/evaleev -p 54321 evaleev@204.111.133.10` (fetched SeQuant at
  `~/code/mpqc4/build/_deps/sequant-src`, sync edits via `ssh 'cat >'`).

## Files this change touches

- `SeQuant/core/eval/dag_scope.hpp` — `DagScopeLevel`, new `LoopKey`, the seam.
- `SeQuant/core/eval/peak_profile.hpp` — `compute_dag_boulevard` (fusion): assign
  per-occurrence `loop_slot`; `OccurrenceRec`/`ValueCell` carry it.
- `SeQuant/core/eval/ordered_schedule.hpp` — `build_ordered_schedule`
  (per-slot realization, layout, atlas), the level enumeration,
  `compute_sliced_mode_assignment`, `home_mode_depth`.
- `SeQuant/core/eval/value_id.hpp` — `value_id_coloring` / `occurrence_key`
  coloring on `(depth, loop_slot)`.
- `SeQuant/core/eval/eval.hpp` — `slice_to_use`, the completeness guard (uses of
  `ordinal` in diagnostics).
- `src/mpqc/.../cck.ipp`, `.../expression/sequant.h` — feed post-remat cells
  (Task 5).

---

### Task 1: Land the vocabulary (no behavior change)

**Why:** introduce `loop_slot`, `altitude_ordinal`, and the `latitude_ordinal`
rename with **zero** behavior change, so later tasks edit a clean base.

**Files:** `dag_scope.hpp`, `ordered_schedule.hpp`, `value_id.hpp`, `eval.hpp`,
and any other `.ordinal` user (grep-driven).

- [ ] **Step 1 — Rename `ordinal` → `latitude_ordinal`** across `core/eval`
  (grep `\bordinal\b` in `DagScopeLevel`/`ScopeBlock` contexts). Mechanical.
- [ ] **Step 2 — Add fields.** `DagScopeLevel`: add `int loop_slot = 0;` and
  `int altitude_ordinal = 0;`. Add `struct LoopKey { std::size_t depth; int
  loop_slot; };` with `==`/hash over `(depth, loop_slot)`.
- [ ] **Step 3 — Keep identity keying unchanged for now.** The level enumeration
  (`ordered_schedule.hpp:1745`) and every `DagScopeLevel::operator==` keep keying
  on the **full** tuple `(depth, space, loop_slot, latitude_ordinal)`. With
  `loop_slot == 0` everywhere and `latitude` as today, this is byte-identical.
  (The spec's "PROCON passes share a LoopId" is deferred to Task 6 — keeping them
  distinct here is harmless.)
- [ ] **Step 4 — Build + validate byte-identical.** `cmake --build
  cmake-build-release --target mpqc -j6`; `ctest -R ordered-schedule` green; run
  `repro/w8-auxocc-ordered-7b.json` and confirm energy **identical** to the
  pre-change run (`-1.6028511154362268`). Commit.

**Produces:** `LoopKey`; `DagScopeLevel{depth, space, loop_slot, altitude_ordinal,
latitude_ordinal}`; the vocabulary, behavior unchanged.

---

### Task 2: Assign per-occurrence `loop_slot` at fusion

> **Micro-design (2026-08-29), connectivity components (spec §4-§5).**
> Supersedes an earlier position-based note (deleted): position-in-own-frame was
> a mid-session shortcut, is exactly what e8bcee766's `slot_of` did, and it does
> NOT establish cross-frame member identity (spec §3: `loop_slot` "assigned by
> fusion").
>
> **Key decoupling (settled 2026-08-29):** `loop_slot` comes from PHYSICAL
> producer→consumer connectivity ONLY -- never the loop-colored value-id (which
> needs `loop_slot`; using it would be circular). The value graph is a DAG, so a
> single pass suffices -- no fixed point. **Symmetry stays OUT of this pass:** a
> symmetric tensor's two modes are still two distinct physical loops; folding
> `A(_,i)` with `A(i,_)` is a downstream VALUE-ID decision (Task 4) via the
> existing loop-colored `occurrence_key`/TNv3 graph path. So NO graph
> canonicalization here. Where connectivity is contradictory (§4 C/D) or symmetry
> leaves the mapping free, the **safe direction is to keep loops distinct**
> (reject the merge): an over-split loop is a missed fusion (still correct), a
> collapse is the crash. The design:
> - **Component nodes:** `(value_id, STRUCTURAL slot)` = a mode's position in
>   `canon_indices`. That order is invariant across a value's occurrences (same
>   structure => same canonical order; only the LABELS differ between terms), so
>   keying by position folds the occurrences across trees. (Keying by `cell.carried`
>   label lookup was tried and REJECTED -- it matches against the first
>   occurrence's labels, which come from a different term, so it strands every
>   relabeled mode; measured 4311 spurious skips on w8.)
> - **Edges:** per `OccurrenceRec` of value V with consumer C (`consumer_point`),
>   for each batched mode both carry, unite `(V, posV)`-`(C, posC)`. Mode
>   correspondence is a `Index ==` match in the **parent OCCURRENCE's own frame**
>   (child and parent share one tree, so labels agree there --
>   `contracted_indices` eval.hpp:549 relies on this). **Batched** = the mode is
>   sliced by an enclosing loop (in the occurrence's `ectx`); a carried-but-not-
>   enclosed mode is correctly left -1.
> - **Conflict = keep distinct (safe):** if a union would place two DISTINCT slots
>   of one value in one component, REJECT it (over-split, never collapse). This is
>   the §4 tracked choice; recording the loser's transposition so the transposed
>   occurrence reads the swapped slot is DEFERRED (improvement 2). w8 measures 36
>   such rejections (safe over-splits, e.g. `i` gets 3 global members not 2).
> - **Numbering:** each component -> one `loop_slot`, ranked per space (first-seen
>   over occurrences -- deterministic).
> - **Storage:** per-`carried`-position `loop_slot` (`svector<int>`, -1 where not
>   batched) on `OccurrenceRec`; additive, nothing consumes it yet.
>
> **Status (2026-08-29): position-based LOCKED, lossless on w8, within-value slots
> distinct.** Cross-occurrence consolidation (canonical-frame refinement) + the
> §4 transposition recording are deferred to be settled by Task 3's actual
> consumption (does w8 stay lossless / w20 stop crashing) rather than in isolation.
> The C/D/E transposition unit fixture (Step 3) is deferred with it.
> - The seam is already wired per-occurrence (`slice_to_use` ->
>   `loop_of_level(level) -> mode_of(hash, LoopId, consumer)`); the missing piece
>   is distinct `LoopId`s per slot, which `loop_slot` in the level provides.

**Why:** `loop_slot` is the missing identity. Fusion is where the correspondence
across frames is established (spec §4-§5).

**Files:** `peak_profile.hpp` (`compute_dag_boulevard`, `OccurrenceRec`,
`ValueCell`).

- [ ] **Step 1 — READ + micro-design (no edit).** Read how `compute_dag_boulevard`
  represents occurrences and their producer→consumer connectivity: `OccurrenceRec`
  (fields, `consumer_point`), how `ValueCell::occurrences` are folded, how
  `divergent_modes` is currently computed, and where a value's operand slots
  connect to its consumer's slots (the dep graph in
  `detail::ordered_schedule_dep_graph` may be the connectivity source). Write a
  half-page micro-design into this plan: the union-find domain (`(value, batched
  slot)` nodes), the edges (producer output slot ↔ consumer input slot, per the
  contraction), the conflict rule (when two edges disagree, keep the first-seen,
  record the loser's transposition), and where the resulting per-occurrence
  `position → loop_slot` is stored.
- [ ] **Step 2 — IMPLEMENT.** Add the per-occurrence `loop_slot` map to
  `OccurrenceRec` (and/or `ValueCell`), populated by the union-find of Step 1.
  Numbering: canonical (e.g. by first-seen order within each `(scope, space)`
  group). Replace nothing yet — this is additive data.
- [ ] **Step 3 — VALIDATE.** Extend `SEQUANT_DUMP_SCHEDULE` (or a new
  `SEQUANT_DUMP_LOOP_SLOT`) to print, per occurrence, `position → loop_slot`. On
  `w8-auxocc-ordered-7b`, a two-batched-mode value shows **distinct** slots `0,1`;
  a fixture matching the spec's C/D/E example reproduces the **transposition**
  (`C: 0→0,1→1` but `E: 0→1,1→0`). Add that fixture to `[ordered-schedule]`.
  Commit.

**Produces:** per-occurrence `position → loop_slot`, stable, order-free.
**Interfaces:** later tasks read `OccurrenceRec`'s `loop_slot` map.

---

### Task 3: Realize one loop per `loop_slot` in `build_ordered_schedule`

**Why:** stop the collapse (spec §5-§6). This is the change that fixes the w20
crash. Developed on the seed `rich` (fusion's `loop_slot` is already on it); Task 5
switches to post-remat.

**Files:** `ordered_schedule.hpp` (`types`, `depth_of_type`, the escape loop,
`make_block`, `compute_sliced_mode_assignment`).

- [ ] **Step 1 — READ + micro-design.** Re-read `types` (992-999),
  `depth_of_type` (1017-1023), the escape loop (1049-1068), `make_block`
  (1188-1221), and the assembly loop (1238+). Micro-design: `types` becomes one
  entry per `(space-group, loop_slot)`; `depth_of_type` becomes
  `level_of(space-group, loop_slot)`; the escape loop emits one escape per
  per-slot mode (no fold); `make_block` stamps `loop_slot` and `altitude_ordinal`
  (the canonical nesting order among a group's slots).
- [ ] **Step 2 — IMPLEMENT the nest.** One loop per `loop_slot`; assign
  `altitude_ordinal`; remove the escape-collapse. Keep the level enumeration
  keyed on the full tuple (Task 1).
- [ ] **Step 3 — IMPLEMENT the atlas.** In `compute_sliced_mode_assignment`,
  build the per-occurrence `(depth, loop_slot) → position` map from the Task 2
  `loop_slot` data + the occurrence's own frame — **no `ectx`/label/space match**.
- [ ] **Step 4 — VALIDATE.** `SEQUANT_DUMP_SCHEDULE` shows, for a two-member
  value, **two distinct `(depth, loop_slot)` loops**; `[ordered-schedule]` green;
  `w8-auxocc-ordered-7b` **lossless** vs `w8-auxocc-fd-7b` reference; **w20
  aux+occ** (remote) completes without `is_range_set_congruent` and matches
  reference. Commit.

**Produces:** a non-collapsed nest + the per-occurrence atlas.

---

### Task 4: Color value-id and occurrence-id on `(depth, loop_slot)`

**Why:** stable identity for cache/remat/router (spec §5.4); today depth-only.

**Files:** `value_id.hpp` (`value_id_coloring`, `occurrence_key`),
`ordered_schedule.hpp` (`home_mode_depth` must carry `loop_slot`, not just depth).

- [ ] **Step 1 — READ.** Confirm `value_id_coloring` (37-45), `occurrence_key`,
  and `home_mode_depth`'s current shape (`pair<Index,int depth>`).
- [ ] **Step 2 — IMPLEMENT.** `home_mode_depth` → carries `(Index, depth,
  loop_slot)`. `value_id_coloring` colors home-sliced modes by
  `space + (depth, loop_slot)`; `occurrence_key` colors **all** modes sliced at
  the occurrence the same way.
- [ ] **Step 3 — VALIDATE.** Two values of one node homed on different slots get
  **distinct** value-ids (dump); `w8-auxocc-ordered-7b` lossless; unit tests
  green. Commit.

---

### Task 5: Build the schedule from post-remat cells

**Why:** spec §6.5/§7.3 — the schedule must reflect remat's placement (splits,
homes), not the seed. Staged after Tasks 1-4 so the collapse fix is validated
first; this completes correctness where remat actually splits.

**Files:** `src/mpqc/.../cck.ipp` (schedule driver), `.../expression/sequant.h`
(remat pre-pass), possibly `ordered_schedule.hpp`.

- [ ] **Step 1 — READ + micro-design.** Re-read the driver (`cck.ipp:2222-2244`)
  and remat pre-pass (`sequant.h:269-276`): confirm `res.cells` is a
  `RichSchedule`-compatible cell set carrying `loop_slot`, and design feeding
  `build_ordered_schedule(RichSchedule{res.cells}, …)` instead of the seed.
  Determine whether the router overlay and the `remat_to_router` assert
  (`placement_remat.hpp:569`) are still needed or subsumed.
- [ ] **Step 2 — IMPLEMENT.** Feed the schedule from post-remat cells.
- [ ] **Step 3 — VALIDATE.** `w8-auxocc-ordered-7b` and a **remat-splitting**
  config (`w8-auxocc-ordered-p3` — after confirming Task 3 addressed its collapse
  aspect, or a config that actually splits) run lossless; the remat pre-pass
  assert no longer fires (or is intentionally removed). w20 lossless. Commit.

**Note:** if Step 1 finds the post-remat switch is entangled with the separate
`ordered_home_reads` bug (the w8-p3 PAO eviction), stop and report — that is a
different defect and out of this plan's scope.

---

### Task 6: Retire the consumer-disambiguation workaround (cleanup)

**Why:** `by_hash_consumer` / `mode_of`'s consumer arm existed only because
members shared a loop (spec §6.6). With `loop_slot`, each member is its own
`LoopId`. Optional; do only if it simplifies without regressing.

**Files:** `dag_scope.hpp` (seam), `ordered_schedule.hpp`, `eval.hpp`.

- [ ] **Step 1 — READ.** Identify every consumer-arm use of `mode_of` /
  `by_hash_consumer` and whether `loop_slot` now makes them single-mode-per-loop.
- [ ] **Step 2 — Optionally make PROCON passes share a `LoopId`** (spec §8):
  key the level enumeration on `(depth, loop_slot)` (drop `latitude` from the
  key), so the two passes fold. Validate byte-identical energy.
- [ ] **Step 3 — Retire the consumer arm** where `loop_slot` suffices; keep it
  only for genuinely-divergent occurrences if any remain.
- [ ] **Step 4 — VALIDATE.** w8 + w20 lossless; unit tests green. Commit.

---

### Task 7: End-to-end validation, cleanup, repin

- [ ] **Step 1 — w8 forced-multi-member** lossless (record both energies + delta).
- [ ] **Step 2 — w20 aux+occ** lossless on z820 (the original crash gone).
- [ ] **Step 3 — `[ordered-schedule]` suite green.**
- [ ] **Step 4 — Remove/guard instrumentation** (`SEQUANT_DUMP_*`) per user
  preference.
- [ ] **Step 5 — Repin.** Land SeQuant commits; separate MPQC repin commit
  (`external/versions.cmake`, `PREVIOUS_*` too). Push.

---

## Self-review

- **Spec coverage:** §3 identity/layout → Task 1; §4 atlas + §5 fusion timing →
  Task 2; §5-§6 collapse → Task 3; §5.4 coloring → Task 4; §6.5/§7.3 post-remat →
  Task 5; §6.6/§8 pass-sharing → Task 6; §9 validation → Tasks 3/7.
- **Ordering:** Task 1 (base) → 2 (data) → 3 (the fix, validated on w8+w20) → 4
  (coloring) → 5 (post-remat) → 6 (cleanup) → 7 (final). Tasks 3 and 7 are the
  losslessness gates.
- **Read-first tasks:** 2, 3, 5 each start with a scoped read producing a
  micro-design, because their internals must be traced before editing — the
  user's measurement-over-speculation constraint.
- **Risk flags:** Task 5 may collide with the separate `ordered_home_reads` bug
  (explicit stop-and-report). Task 6 is optional cleanup, not correctness.
