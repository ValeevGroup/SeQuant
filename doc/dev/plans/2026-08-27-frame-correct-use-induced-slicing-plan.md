# Frame-correct use-induced sliced-mode map — implementation plan

**Goal:** Make the ordered DAG evaluator's "which loop slices which mode of a
value occurrence" map frame-correct and per-instance, so aux+occ (external-occ)
CCSD wet evaluation runs crash-free and lossless (aux-only already works).

**Architecture:** The factorizer emits a *flat* per-instance set of batched
indices with roles (no same-space nest — settled by the hinge read). Stage 4
(`build_ordered_schedule`) currently **collapses** that set to one loop per
space; the fix is to **expand** it to one loop per batched index instance at
distinct `DagScopeLevel`s, then let every downstream consumer address loops by
`LoopId`/`{space,depth,ordinal}` and modes by the value's own-frame `per_axis`
instance — never by label, space, or `ectx` match. The runtime seam
(`LoopColoredSliceSeam` / `slice_to_use`) is unchanged in shape; it consults the
now-correct per-instance `occ_facts`.

**Tech stack:** C++20, SeQuant `core/eval` + `core/optimize`; validation via the
MPQC `mpqc` executable on `csv-cck-w8.json` (fast forced-occ repro) and the
w20 aux+occ wet run; SeQuant Catch2 `[ordered-schedule]` unit tests.

**Spec:** `SeQuant/doc/dev/specs/2026-08-27-frame-correct-use-induced-slicing-design.md`
(read it first — this plan argues from it and the two travel together).

## Global constraints

- **Frame purity (the binding invariant):** no code path outside a single
  value/occurrence frame may compare index **labels** or **space keys**
  (`base_key`), or walk `occ.ectx`, to bind a mode to a loop. Loops are
  addressed only by `LoopId` / `{space, depth, ordinal}`; modes only by the
  value's own-frame `per_axis` instance.
- **Losslessness is the acceptance bar**, not merely "no crash": batched energy
  must match the unbatched reference (`peak_threshold ≈ ∞`) to full precision.
- **Measurement over inference.** Every task that claims a code fact (a role, a
  level, a collapse) confirms it with an instrumented run on the w8 repro before
  building on it. Do not fabricate a diff for code whose current behavior has
  not been dumped. Where a task has a genuine design fork, it is called out as a
  **Decision point** to settle by measurement, not papered over.
- **No `node_slice_mask` / single-term / `binarize` / aux-path changes** (spec
  §10). The aux contracted path already works via `ectx`; leave it byte-identical.
- Build only the `mpqc` target, `-j6`. Serialize mpqc runs (each uses all cores).
- clang-format touched files (`/opt/homebrew/opt/llvm/bin/clang-format
  --style=file -i`); `.tpp` excluded. No en-dashes (pre-commit hook). No
  Co-Authored-By / AI trailers. A SeQuant→MPQC repin is its own commit.
- w20 wet runs execute on the remote box:
  `ssh -i ~/.ssh/evaleev -p 54321 evaleev@204.111.133.10`; log at
  `/home/evaleev/projects/pno/w20.aux+occ.log`. w8 runs locally.

## Files this change touches

- `SeQuant/core/eval/legality.hpp` — `classify_axis` (234-278): role producer.
  Verified/repaired in Task 1.
- `SeQuant/core/eval/ordered_schedule.hpp` — `build_ordered_schedule`: `types`
  (992-999), `depth_of_type` (1017-1023), escape-collapse (1049-1068),
  `make_block` (1188-1221), `home_mode_depth` (1395-1447),
  `compute_sliced_mode_assignment` (1680-1805). The core of the change.
- `SeQuant/core/eval/dag_scope.hpp` — `DagScopeLevel` / `SlicedModeAssignment`
  enumeration, if per-instance levels need a new field. Read before editing.
- `SeQuant/core/eval/value_id.hpp` — `home_mode_depth` consumer
  (`value_id_coloring`, 37-45): touched only if Task 4's decision needs it.
- `SeQuant/core/eval/ordered_executor.hpp` — `run_ordered_contracted_block`
  (244-565): scatter destination `dm[k]`.
- `SeQuant/core/eval/member_axis.hpp` — `member_external_axis` /
  `member_contracted_axis`: scatter axis resolution (Task 6).
- `SeQuant/core/eval/eval.hpp` — `slice_to_use` completeness guard (807-814),
  Task 7.

Runtime consumers already frame-correct and **not** edited except as noted:
`LoopColoredSliceSeam::mode_of`, `slice_to_use` ordered arm.

---

### Task 1: Confirm (and repair if needed) the per-instance role source

**Why first:** the whole design consumes `per_axis` roles (LoopLocal /
Reduction / LoopCarried) as the per-instance external-vs-contracted-vs-local
classification. The audit found the producer `classify_axis`
(`legality.hpp:234-278`) derives the role from `base_key` (space), `occ.ectx`,
and an across-occurrence label compare — the machinery this design forbids. Its
output may still be right within a non-divergent cell, or it may misclassify.
Settle this by measurement before anything else stands on it.

**Files:** instrument `SeQuant/core/eval/legality.hpp` (temporary dump in
`classify_axis` / `analyze_legality`); if repair is needed, edit `classify_axis`.

- [ ] **Step 1 — Instrument.** In `analyze_legality`, after `per_axis` is built
  for each value, dump per value: node hash, and for each axis its
  `{index-label, base_key, role}`. Guard behind an env flag
  (`SEQUANT_DUMP_PER_AXIS`) so it is off by default.
- [ ] **Step 2 — Measure on w8 forced-occ.** Run `csv-cck-w8.json` (small
  `occ_target_size`, `peak_threshold` large) with the flag set. Locate a doubles
  amplitude with two occ externals `{i_1, i_2}` and a contraction that sums an
  occ index (e.g. `g(i_2,i_1,Κ)·I(μ̃,i_2,Κ)→I(i_1,μ̃)`).
- [ ] **Step 3 — Check the invariant.** Confirm: both occ externals of the
  amplitude are **LoopCarried**; the summed occ index is **Reduction**; a
  home-consumed occ index is **LoopLocal**. Record the dump in the ledger.
- [ ] **Step 4 — Branch on the result.**
  - **If roles are correct per-instance:** `classify_axis` is *usable as-is*
    (mechanism suspect, output sound within the cell). Note this in the ledger
    and proceed; do not refactor it in this task.
  - **If any role is wrong:** repair `classify_axis` to derive the role from the
    value's **own build tree** only — LoopLocal = the axis is consumed within
    this value's build; Reduction = summed at or below this value; LoopCarried =
    survives into this value's own result — with **no** `base_key` compare, no
    `ectx` walk, no cross-occurrence label compare. Re-run Step 2-3 until the
    invariant holds.
- [ ] **Step 5 — Remove instrumentation** (or leave the env-guarded dump if it
  proves useful; decide with the user). Commit (repair only; the dump alone is
  not a commit).

**Decision point:** whether `classify_axis` is repaired or used as-is. This
gates the shape of every later task (they all read `per_axis`).

**Produces:** a trusted per-value, per-instance `role` for each build-site axis.

---

### Task 2: Expand Stage 4 from one-loop-per-space to one-loop-per-instance

**Why:** this is the core defect (spec §5, §6). `build_ordered_schedule` builds
its `DagScopeLevel` nest with **one representative index per space** (`types`,
992-999) and looks levels up by space (`depth_of_type`, 1017-1023). Two carried
occ externals `{i_1, i_2}` therefore share one depth and cannot be told apart
downstream. They must become **two distinct loops at distinct depths**, 1-to-1
with the batched index instances. The factorizer emits them flat (hinge, spec
§9-RESOLVED), so Stage 4 owns synthesizing the nest and choosing its order.

**Files:** `SeQuant/core/eval/ordered_schedule.hpp` (`types`, `depth_of_type`,
`make_block`); `SeQuant/core/eval/dag_scope.hpp` if a new level field is needed.

- [ ] **Step 1 — Read `dag_scope.hpp` fully** (the `DagScopeLevel` /
  `SlicedModeAssignment::levels` enumeration) and `make_block` (1188-1221) before
  editing, and dump the current w8 forced-occ schedule levels (env-guarded) so
  the "before" is on record: expect one occ depth even though two occ externals
  are batched.
- [ ] **Step 2 — Enumerate per-instance loops.** Replace the per-space
  representative in `types` with **one entry per batched index instance** carried
  by the scope (source: the flat per-instance batched set — the union over the
  scope's values of their LoopCarried/Contracted batched `per_axis` instances,
  from Task 1). Same-space instances get **distinct depths** (nested grid).
- [ ] **Step 3 — Fix the nesting order (Decision point).** Same-space instances
  are order-interchangeable (spec §6). Pick a **canonical, deterministic** order
  (e.g. by the instance's canonical index ordinal within its space) so the
  schedule is reproducible. Confirm the choice does not disturb the
  producer/consumer `ordinal` split (a *different* use of `ordinal`, over a
  single index — verify the two mechanisms compose by dumping a forced-split case
  that also has two occ externals).
- [ ] **Step 4 — Level lookup.** Turn `depth_of_type` (space→depth) into an
  **instance→level** resolution keyed on the batched index instance, since a
  space now spans several depths. Every caller of `depth_of_type` must pass the
  instance, not the space.
- [ ] **Step 5 — Validate.** `[ordered-schedule]` unit tests green
  (`ctest -R ordered-schedule` in the SeQuant build). Re-dump the w8 schedule:
  the two occ externals now occupy **two distinct occ levels**. Do NOT expect a
  correct wet run yet (escape-collapse and the seam still consume the old shape).
  Commit.

**Produces:** a schedule whose nest hosts one distinct `DagScopeLevel` per
batched index instance; `depth_of_type` replaced by an instance→level resolver.

---

### Task 3: Remove the escape-collapse

**Why:** even with per-instance levels, lines 1049-1068 fold "two carried occ
indices … to ONE escape at that depth, with LoopCarried dominating." That folds
`{i_1, i_2}` back into a single escape and destroys the 1-to-1 Task 2 just built.

**Files:** `SeQuant/core/eval/ordered_schedule.hpp` (escape-collapse block,
1049-1068).

- [ ] **Step 1 — Read the collapse block** and its consumers (what reads the
  escape set it produces) so removal does not orphan a downstream read.
- [ ] **Step 2 — Make each carried instance escape to its own loop.** Replace the
  same-depth fold with a per-instance escape: each LoopCarried batched instance
  escapes to *its own* level (the one Task 2 assigned it).
- [ ] **Step 3 — Validate.** `[ordered-schedule]` tests green; w8 schedule dump
  shows **two escapes** for the doubles amplitude, one per occ external. Commit.

**Produces:** per-instance escapes; the mode→loop binding survives Stage 4.

---

### Task 4: `home_mode_depth` resolves mode→carried via `per_axis`

**Why:** the home side (1395-1447) resolves mode→carried by matching each
level's **space** to the value's first-unconsumed carried mode of that space —
the same space anti-pattern, arbitrary when a space has two carried modes homed
at one scope (spec §1).

**Files:** `SeQuant/core/eval/ordered_schedule.hpp` (`home_mode_depth`,
1395-1447); `SeQuant/core/eval/value_id.hpp` only if the level-keying decision
requires it.

- [ ] **Step 1 — Replace the space match** with a `per_axis` LoopLocal lookup:
  the home-sliced modes are exactly the value's LoopLocal instances (Task 1), in
  the value's own frame — resolve mode→level per-instance from those, never by
  space.
- [ ] **Step 2 — Decision point: level-keying.** Determine whether
  `home_mode_depth` must key the full `{space,depth,ordinal}` (→ `value_id.hpp`
  change) or may keep `int depth` / slot-only. Per spec §9, per-scope scratches
  already separate cross-ordinal homes, so slot-only is the likely answer. Settle
  it by constructing (or dumping) two values of one node homed in sibling loops
  and checking they land in distinct scratches regardless of color. Record the
  finding; touch `value_id.hpp` only if the measurement demands it.
- [ ] **Step 3 — Validate.** `[ordered-schedule]` tests green; the home-slice
  path on a value with two LoopLocal occ modes resolves both correctly (dump).
  Commit.

**Produces:** frame-pure home-slice mode resolution.

---

### Task 5: `compute_sliced_mode_assignment` — read, don't match

**Why:** this is the seam that feeds `occ_facts` → `slice_to_use`. Today it does
a positional zip that reads `occ.ectx` (1775) — the ectx/label machinery this
design forbids, and the direct cause of every prior crash (spec §3).

**Files:** `SeQuant/core/eval/ordered_schedule.hpp`
(`compute_sliced_mode_assignment`, 1680-1805).

- [ ] **Step 1 — Rewrite the per-occurrence loop.** For each occurrence, for each
  of its **LoopCarried** `per_axis` instances (own frame, Task 1): the loop is
  that instance's realized `LoopId` from the Task 2 schedule, composed with the
  enclosing loop copy (build vs consumer pass) from
  `build_scope[occurrence.consumer]`, which already carries the full
  `DagScopeLevel` **including `ordinal`** (spec §4). Emit that as the occurrence's
  sliced-mode fact. **No `ectx` read, no `axis`/space/label compare anywhere.**
- [ ] **Step 2 — Exclude home-sliced.** LoopLocal modes stay with the value-id
  coloring (Task 4) and must **not** enter the use-induced seam (avoid
  double-slicing). Reduction modes are never sliced on fetch.
- [ ] **Step 3 — Validate.** w8 forced-occ: the fetch occurrence of the archetype
  top-homed leaf `g` (home_d=0, two hops) now reports **both** occ externals
  sliced to their two loops, and the contracted occ index reports **unsliced**
  (the "uniq over-sliced `i_2`" crash, spec §3, must not reappear). Dump
  `occ_facts` and confirm zero unresolved facts. Commit.

**Produces:** frame-correct per-occurrence use-induced sliced-mode facts.

---

### Task 6: Per-instance external-occ scatter

**Why:** with per-instance loops restored, a scatter output carrying two batched
occ externals must restrict **both** modes in its destination block; today
`dm[k]` / `member_axis` handle exactly one (spec §7.3) — the runtime face of the
same collapse.

**Files:** `SeQuant/core/eval/ordered_executor.hpp`
(`run_ordered_contracted_block`, scatter destination ~244-565);
`SeQuant/core/eval/member_axis.hpp` (`member_external_axis` /
`member_contracted_axis`).

- [ ] **Step 1 — Read `member_axis.hpp` and the `dm[k]` scatter path** and dump
  (env-guarded) the destination-block restriction it computes on w8 forced-occ —
  confirm it currently restricts one occ mode where two are batched.
- [ ] **Step 2 — Generalize `dm[k]` to a per-instance set.** The scatter must
  restrict every batched external instance the output carries, each to its own
  loop's current block, driven by the Task 5 per-instance facts (not a
  space-derived single axis). `member_axis.hpp`'s `base_key`-driven single-axis
  resolution (audit: all four functions space-tainted) is replaced by the
  per-instance facts for the batched modes.
- [ ] **Step 3 — Validate.** w8 forced-occ scatter writes the correct block for a
  two-occ-external output; energy sanity (a single-iteration residual matches the
  unbatched residual for that term to full precision). Commit.

**Produces:** scatter that restricts all batched external modes per-instance.

---

### Task 7: Re-express the completeness guard per-occurrence

**Why:** `slice_to_use`'s guard (`eval.hpp:807-814`) currently throws at
**hash-level** `participates()`. With per-instance facts, an occurrence that
legitimately carries an unsliced mode (invariant across a fetch) must not trip
merely because *another* occurrence of the same hash is sliced (spec §8).

**Files:** `SeQuant/core/eval/eval.hpp` (completeness guard, 807-814).

- [ ] **Step 1 — Rewrite the guard** to fire on "this occurrence carries a mode
  bound to the crossed loop but no fact resolved it," not on hash-level
  participation.
- [ ] **Step 2 — Validate.** w8 forced-occ runs with the guard active and does
  **not** throw on the correctly-unsliced case; it still throws (loudly) if a
  genuinely-required fact is missing (test by temporarily dropping one fact).
  Commit.

**Produces:** a guard that is per-occurrence, not hash-level.

---

### Task 8: End-to-end validation and repin

**Why:** the acceptance bar is lossless wet evaluation, not local dumps.

- [ ] **Step 1 — w8 forced-occ lossless.** `csv-cck-w8.json` with small
  `occ_target_size` and `peak_threshold ≈ ∞`: runs crash-free, zero unresolved
  facts, and batched CCSD energy matches the unbatched reference to full
  precision. Record the two energies and their delta in the ledger.
- [ ] **Step 2 — `[ordered-schedule]` unit suite green** in the SeQuant build.
- [ ] **Step 3 — w20 aux+occ (the original failure) lossless** on the remote box.
  Compare against `w20.aux+occ.log`'s prior TA-assert failure: it must now
  complete and match the aux-only / reference energy.
- [ ] **Step 4 — Repin.** Land the SeQuant commits, then a **separate** MPQC repin
  commit bumping `MPQC_TRACKED_SEQUANT_TAG` in `external/versions.cmake` (update
  `PREVIOUS_*` too). Push.

**Produces:** the shipped, validated fix.

---

## Self-review notes

- **Spec coverage:** §1 home-slicing → Task 4; §1/§4-caveat role source → Task 1;
  §5 collapse → Tasks 2-3; §6 per-instance model → Tasks 2-3-5; §7.1 Stage 4 →
  Tasks 2-3, §7.2 seam → Task 5, §7.3 scatter → Task 6; §8 guard → Task 7,
  validation → Task 8; §9-RESOLVED (flat factorizer) is the premise of Task 2.
- **Ordering rationale:** Task 1 gates all (everything reads `per_axis`); Tasks
  2-3 rebuild the schedule shape before any consumer (5, 6) reads it; Task 4
  (home) is independent of 5 (use-induced) but shares the `per_axis` resolver;
  the guard (7) is last before end-to-end because it can only be tuned once real
  facts flow.
- **Genuine open forks are Decision points, not fabricated code:** Task 1 (repair
  vs use-as-is), Task 2 Step 3 (nesting order canonicalization), Task 4 Step 2
  (level-keying). Each is settled by a named measurement, per the user's
  evidence-first constraint — deliberately not pre-committed here.
