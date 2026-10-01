# Order-aware multilevel batching -- runtime realization Implementation Plan

> **RETRACTED / CORRECTED (2026-07-27) -- contaminated "external-occ" C60 numbers.**
> This plan's OUTCOME and Tasks 7-8 rest on C60 measurements ("the C60 peak floor
> (2302 GB giant) is set by a separate, upstream defect", "the 4-occ recompute is
> inevitable once external-occ pairs are batched", "the 2302 GB transient is the
> external-forest lever") that were taken with occupied indices in the
> CONTRACTED-role batchability predicate (only `is_batchable_index` set;
> `is_batchable_external_index` never supplied). The DP ALSO batched contracted occ
> (cost cells sliced the contracted occ pair `i_3 i_4`), so every "external-occ"
> figure was "external-occ ON TOP OF contracted-occ" -- the IDENTICAL peak across
> the contracted/external arms is the tell. Post role-split (SeQuant
> `79540c831`..`7f5240014`, MPQC `1685c7ac1e`..`5a3d94dc63`) honest re-measure of
> `[.][dryrun-occ-veto]` (nterms=55): peak ~= **5860.9 GB** (NOT 2302 / 2947),
> avoidable_time ~= **44.6%**, `contracted_occ_stamps == 0`,
> `external_occ_stamps == 244`. The peak is HIGHER because the phantom
> contracted-occ reduction hid the true cost of the cross-pair two-PNO-leg giants;
> external slicing cannot reduce them (`a<i_3,i_4>` stays full). The "infrastructure
> landed" conclusion for Tasks 2-4 is unaffected; the C60 peak/lever attributions
> are RETRACTED. Root cause:
> `.superpowers/sdd/contamination-role-predicate.md`.

> **OUTCOME (2026-07-26).** Tasks 2, 3, 4, 6 are **implemented, reviewed clean,
> and landed** (`0c36ebe24` emit, `92c3c9001` scope-chain fall-through,
> `f42b641267` runtime hoist, `235309f01` veto reinterpret, `20dab38a9` sentinel
> test fix). The hoist mechanism is correct and **OFF-path byte-identical**
> (gated on `scope_level != kNoBatchScopeLevel`). It is landed as **infrastructure**.
>
> **Task 5 (B3 tight-count release) DEFERRED** -- premise was flawed: the emitted
> `effective_count` equals `rf` (the recompute factor, 1 for a fully-hoisted
> node), NOT the read-multiplicity a tight lifetime needs. Task 4 ships the
> spec-sanctioned `reset()` backstop instead (correct; the avoidable-time gate
> needs build-once, not tight release). A future tight-release pass must first fix
> the emit to carry the true read count (`prod nBatch` over loops enclosing the
> consumer but not the node, recoverable via the parent cell).
>
> **Tasks 7-8 (C60 gate + MPQC enablement) NOT the win here.** The `[dryrun-occ-veto]`
> gate does NOT reach `avoidable_time < 0.05`, and that is **expected**: two probes
> (2026-07-26) showed the C60 peak floor (2302 GB giant) is set by a *separate,
> upstream* defect -- the external-batching **DP selection never adopts/stamps the
> over-budget giants** (every >100 GB node emits `modes=(none)` or `K_2:con`, zero
> `External`), so they never enter external batching. The external *scatter* itself
> is sound (slices exactly `1/nblk`). The 4-occ recompute this plan's hoist targets
> is, per the user, **inevitable** once external-occ pairs are batched (batching a
> product over both free legs forces recompute) -- so it was never the peak lever.
> The peak-under-budget work lives in the external-mode-batching D1 hole
> (`2026-07-20-external-mode-batching-design.md`), a distinct workstream. This plan
> is closed at "infrastructure landed"; do not execute Tasks 5, 7, 8 as written.
>
> **PREREQUISITE before enabling order-aware batching through any runtime driver**
> (e.g. MPQC, the old Task 8): the runtime `scope_level >= 0` hoist/walk-up path
> (store an invariant at an INNER scratch) ships **untested** -- only the
> `scope_level == -1` path has a runtime test (F2). It is latent/unreachable today
> (`order_aware_recompute` is test-only in this branch, no shipped caller produces a
> nested 2+ mode order-aware schedule). Its release-mode null-deref risk (F3) is now
> guarded (unconditional throw, `13026b26f`). Land a runtime `scope_level >= 0` hoist
> test (a nested-nest TA network, built once, slice-exact) BEFORE any driver enables
> the ordered runtime path.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Supersedes** `2026-07-22-order-aware-multilevel-batching-plan.md`. That plan
> was written assuming the order-aware cost model (its Phase A) was unbuilt. A
> 2026-07-25 freshness pass against the branch found Phase A **already
> implemented, committed, and unit-tested**: the ordered-cell layer
> (`ac3f58bfd` stage 1, `388dcb7aa` stage 2, `e3198877d` depth cap, `049179059`
> external exclusion, `83028b0d4` node-level external placement) gives
> `PeakBatchedModel` an order-aware recompute charge (`escaped_outer` +
> `descend` + resident-scan), wired live into `relax`, gated by
> `order_aware_recompute` (default OFF), with `[loop-tree]` DP tests passing
> 7/7. The real remaining surface is the **annotation bridge + runtime
> hoisting**, which this plan covers. The prior plan's A1-A3 would re-derive
> committed code; do not execute them.

**Goal:** Make the batched evaluator realize the loop-invariant hoisting the
order-aware cost model already prices, driving the C60 PNO-CCSD residual's
contracted-occ middle-gap recompute (`[dryrun-occ-veto]` `both.avoidable_time()`)
from 0.437 to below 0.05, with the sliced peak unchanged.

**Architecture:** One annotation contract joining a DP that already chose the
schedule to a runtime that does not yet obey it. **The bridge (A4):** the emit
walk `reconstruct_batched_modes::build` (`cost_model.hpp:1722`) already descends
the *chosen* schedule tree with, at each node, the ordered enclosing cell `B` and
the contracted-here set `r.aprime`; extend it to also compute and emit each node's
**placement level** (its lifetime-scope) and **effective use count**
(`escaped_outer`-derived), threaded through the existing
`node_axes` -> `OptimizeOptions::term_batch_axes` ->
`BinarizationOptions::node_batch_axes` -> `EvalExpr::set_batched_here` channel.
**The runtime (B):** the batched evaluator (`eval.hpp`) consumes those fields via
the committed scope-chain primitive (`cache_manager.hpp` `parent_`/`set_parent`/
fall-through `access`, `75d240079`): set each nested scratch's parent, store each
loop-invariant descendant at its emitted scope level, release it by its emitted
effective count (count-zero, scope-reset backstop).

**Tech Stack:** C++20, SeQuant, Catch2 unit tests, the dry-run eval backend
(`core/eval/backends/dryrun/`). Build in `cmake-build-debug` (the only build dir
that compiles the eval unit tests -- it has the TiledArray deps; `build-test` /
`build-rwdi` have TA/BTAS OFF and skip those `.cpp`). The dry-run witness is a
symbolic zero-data replay, so Debug speed is fine.

## Global Constraints

- **Success gate (avoidable-time, in scope):** `[.][dryrun-occ-veto]` asserts
  `aux.avoidable_time() < 0.05` AND `both.avoidable_time() < 0.05` (fails today
  at 0.065 / 0.437). This is the Phase-B target and this plan's gate.
- **Peak is OUT of scope (external-forest follow-on).** The witness also asserts
  `both.peak_gb < 100`; today it is 2302 GB -- one transient
  `g(mu~,mu~,K_2)*C` full in the free PAO `mu~` and the proto-occ pair. Its lever
  is *external* occ batching (`2026-07-20-external-mode-batching-design.md`), not
  the contracted middle-gap cache this plan fixes. Task 7 **splits** the witness
  gate so the avoidable-time assertion is this plan's gate and the peak
  assertion is explicitly tagged as the external follow-on (not regressed, not
  this plan's to pass).
- **No regression / byte-identical OFF path.** `order_aware_recompute = false`
  and `charge_batch_recompute = false` and the unbatched / `peak_threshold=+inf`
  legacy heuristic path must stay byte-identical. Every new emit field is
  populated only under `order_aware_recompute`; every new runtime behavior is
  reached only when a node carries the new annotation, which is empty on the OFF
  path.
- **Peak stays sliced.** A hoist that restores sharing by caching a value *whole*
  across its own batched mode defeats batching and is caught by the witness's
  reported peak -- a hoisted node is held at its own footprint, sliced on the
  modes it carries.
- **`persistent_` stays one bit** (P = run-scope, the sole no-decay case). Inner
  scratches carry single-scope, count-released entries.
- **Determinism.** Co-contracted loop order among modes contracted at one node is
  the canonical ascending batched-slot ordinal already used by `descend`
  (`cost_model.hpp:854`); do not introduce a second ordering rule.
- **Cross-repo.** Any SeQuant behavior change needs an MPQC repin (own commit,
  update `PREVIOUS_*` too); validate on local `[.][dryrun-occ-veto]` before any
  Owl run. No `Co-Authored-By` trailers. ASCII hyphens only (pre-commit U+2013 /
  U+00A0 detectors). Cap builds at `-j6`. Work on branch
  `evaleev/feature/suppress-heuristic-fallback`; commit there; ONE task per
  commit; no push / PR / remote ops.

---

## File structure

- `core/optimize/cost_model.hpp` -- `reconstruct_batched_modes::build`
  (`:1722`, the emit walk); `escaped_outer` (`:869`), `descend` (`:857`),
  `cell_seq_`/`cell_union` (`:845-865`) already expose per-node placement. A4
  extends the emit here; no change to `relax`/selection.
- `core/optimize/options.hpp` -- `OptimizeOptions::term_batch_axes` payload;
  A4 widens the per-node element to carry the two new scalars.
- `core/eval/eval_expr.{hpp,cpp}` -- `batched_here()`/`set_batched_here`
  (`eval_expr.hpp:282-293`), `batch_axes_` (`:311`), and `BinarizationOptions::
  node_batch_axes` (`:330`). A4 adds the scope-level + effective-count fields and
  accessors; `binarize` stamps them.
- `core/eval/cache_manager.hpp` -- committed scope chain (`parent_`, `set_parent`
  `:205`, fall-through `access` `:257`), `entry` decay
  (`max_life`/`life_c`/`decay()` `:96-147`), the veto (`:555-620`). B1/B3/B6
  wire these.
- `core/eval/eval.hpp` -- `make_batched_custom_evaluator` (`:1361`),
  `make_batched_scratch` (`:1236`), the per-batch replay / reinstall
  (`:1535`, `:1696`). B2 sets scratch parents and stores invariants at scope
  level here.
- `tests/unit/test_optimize.cpp` -- `[loop-tree]` DP + emit tests.
- `tests/unit/test_cache_manager.cpp` -- scope-chain visibility + effective-count
  lifetime unit tests.
- `tests/unit/test_eval_ta.cpp` -- runtime build-count / reference-equality test.
- `tests/unit/test_eval_dryrun.cpp` -- `[.][dryrun-occ-veto]` witness (the gate).

---

## Task 1: Investigation -- placement recoverability + emit-vs-select charge (INVESTIGATION)

The freshness pass found a `2.24x` gap on the C60 witness between replayed
`total_exec` (1.02e17) and the DP cost profile `cp` (4.55e16), and the emit
walk's own debug block (`cost_model.hpp:1736`) computes the OLD order-blind
`esc = B & ~ctx.open_modes[n]`, not `ctx.escaped_outer`. Before emitting, settle
exactly what the chosen schedule carries and how to recover per-node placement.

**Files (read-only):** `core/optimize/cost_model.hpp` -- `relax` (`:445+`, the
`escaped_outer`/resident-scan charge, `:1156-1203`), `select_root` (the chosen
`best`), `reconstruct_batched_modes::build` (`:1722-1794`), `escaped_outer`
(`:869`), `descend` (`:857`), `cell_seq_`/`cell_union` (`:845-865`); the witness
`cp` computation in `tests/unit/test_eval_dryrun.cpp` (search `cost_profile` /
`cp `).

**Deliverable** (write to `.superpowers/sdd/oamb2-t1-note.md`), with exact
function+line references:

1. **Placement is recoverable in `build`.** Confirm that at each `build(n, B,
   idx)` node, `B` is the ordered enclosing cell id and `ctx.cell_seq_[B]` /
   `ctx.escaped_outer(B, ctx.open_modes[n])` yield (a) the node's placement level
   = the ordered position of its innermost carried mode, and (b) its escaped-outer
   set (enclosing loops it does not carry). State the exact expression for the
   **lifetime-scope level** (an integer nesting depth, or a canonical mode-slot
   identifying the scope) and for the per-node **recompute factor** `rf =
   prod(nbatches[k])` over `escaped_outer` bits.
2. **Effective use count.** State how to compute a node's effective use count
   (`Sigma over consumers Prod nBatch(L)` over loops enclosing the consumer but
   not the node) from the back-pointer tree in `build`. If it is not cheaply
   recoverable in the single post-order `build` walk, specify the minimal extra
   pass (e.g. a consumer-count accumulation keyed by node subset) -- do NOT
   guess; if it needs a second walk, say so and size it.
3. **Emit vs. select charge.** Determine whether the chosen `best` (from
   `select_root`) and the reported `cp` already reflect the order-aware
   `escaped_outer` charge, or whether `cp` is a `B=0`/set-based re-walk (the
   phantom). This decides whether Task 7's avoidable-time drop needs ONLY the
   runtime hoist (B) or also a `cp`/emit charge correction. Answer with the code
   path that produces `cp`.
4. **Threading the two scalars.** Confirm the channel
   `node_axes` (`:1720`) -> `OptimizeOptions::term_batch_axes` ->
   `BinarizationOptions::node_batch_axes` (`eval_expr.hpp:330`) -> `binarize`
   stamp -> `set_batched_here`, and state the minimal payload change (widen the
   per-node element from `svector<pair<Index,BatchModeType>>` to a struct adding
   `int scope_level; std::size_t effective_count;`, vs. parallel vectors). Pick
   one and justify.

- [ ] **Step 1:** Read the functions above; write the note answering (1)-(4) with
  exact line references and the chosen payload shape.
- [ ] **Step 2:** Confirm nothing in the note requires changing the
  `order_aware_recompute=false` path (must stay byte-identical). If (3) finds
  `cp` is phantom-priced, record the exact fix location for Task 7 to consume.
  Report DONE with the note path.

**Escalation:** if (2) finds the effective use count needs global cross-node
consumer analysis not available from the per-term back-pointer tree, STOP and
escalate -- the annotation is larger than a `build`-walk extension.

---

## Task 2: A4 -- emit `{scope_level, effective_count}` per node

Per Task 1's note, extend the emit so each contraction node carries its placement
level and effective use count alongside `batched_here()`.

**Files:** `core/optimize/cost_model.hpp` (`reconstruct_batched_modes::build`),
`core/optimize/options.hpp` (`term_batch_axes` payload),
`core/eval/eval_expr.{hpp,cpp}` (new fields + accessors +
`BinarizationOptions::node_batch_axes` stamp), `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Consumes: Task 1's payload shape; the ordered cell / `escaped_outer` per node.
- Produces: `EvalExpr::batch_scope_level() -> int` (the emitted lifetime-scope
  level; `-1` = run/term scope, i.e. carries no batched mode -> not hoisted into
  a batch scratch) and `EvalExpr::batch_effective_count() -> std::size_t`
  (emitted use count; `0` = unset/OFF path). Both default to the OFF-path values
  so `order_aware_recompute=false` emits nothing and is byte-identical.

- [ ] **Step 1: Write the failing test** (`test_optimize.cpp`, `[loop-tree]`):
  reuse the existing `[loop-tree]` `R{i,m}` network (batched `{k,i}`, `I2`-class
  node invariant to `i`). After `reconstruct_batched_modes`, assert the emitted
  node for the `I2`-class intermediate has `scope_level` equal to the `k`-level
  (above the `i` loop) and `effective_count == nBatch(i)`; assert a node
  contracting two batched modes emits them in ascending batched-slot order
  (matches `descend`). Follow the existing `[loop-tree]` test at
  `test_optimize.cpp:3394` for the network-build + `model.order_aware_recompute
  = true` idiom.
- [ ] **Step 2: Run to confirm it fails** (fields/accessors absent):
  `./cmake-build-debug/tests/unit/unit_tests-sequant "[loop-tree]"`.
- [ ] **Step 3: Implement.** Widen the `term_batch_axes` per-node payload per
  Task 1; in `build`, compute `scope_level` and `effective_count` from `B` /
  `escaped_outer` / the back-pointer tree; stamp them through
  `node_batch_axes` in `binarize`; add `batch_scope_level_` /
  `batch_effective_count_` members + accessors on `EvalExpr` next to
  `batch_axes_` (`eval_expr.hpp:311`). Populate only under
  `node_level_placement` / `order_aware_recompute`; leave OFF-path defaults.
- [ ] **Step 4: Run to confirm it passes.** Then no-regression:
  `./cmake-build-debug/tests/unit/unit_tests-sequant "[optimize]"` -- assertion
  count unchanged vs pre-change (OFF path byte-identical); record counts.
- [ ] **Step 5: clang-format** the changed files (`.tpp` excluded).
- [ ] **Step 6: Commit.**

---

## Task 3: B1 -- scope-chain visibility (read up on a local miss)

Test and harden the committed fall-through so a node stored at an ancestor cache
is found by an inner body, and confirm a level's reset does not clear ancestors.

**Files:** `core/eval/cache_manager.hpp` (`access` `:254`, `set_parent` `:205`,
`reset` `:231`), `tests/unit/test_cache_manager.cpp`.

**Interfaces:**
- Consumes: committed `parent_`/`set_parent`/fall-through primitive.
- Produces: verified fall-through read + level-local reset (no new API).

- [ ] **Step 1: Write the failing test** (`test_cache_manager.cpp`, follow the
  existing `set_parent` test at `:243`): build `outer` and `inner`
  `CacheManager`s registering the same node key; `inner.set_parent(&outer)`;
  `outer.store(node, v)`; assert `inner.access(node)` returns `v` (walk-up), then
  `inner.reset()` and assert `outer.access(node)` still returns `v` (ancestor
  untouched) while `inner`'s own local entries are cleared.
- [ ] **Step 2: Run to confirm** it passes if the primitive is complete, else
  fails: `./cmake-build-debug/tests/unit/unit_tests-sequant "[cache]"` (confirm
  the actual tag from the file header).
- [ ] **Step 3:** If any assertion fails, harden `access`/`reset` minimally (the
  fall-through is at `:257`; `reset` at `:231` already iterates only local
  `cache_map_`). Add cycle-guard note only if a real chain can loop (deferred
  minor from the prior ledger).
- [ ] **Step 4: Pass. Step 5: Commit.**

---

## Task 4: B2 -- store each invariant at its emitted scope level

Make the batched evaluator set nested scratch parents and store each
loop-invariant descendant once at its emitted scope level, instead of rebuilding
it on the per-batch scratch.

**Files:** `core/eval/eval.hpp` (`make_batched_custom_evaluator` `:1361`,
`make_batched_scratch` `:1236`, the per-batch reinstall `:1535`/`:1696`),
`tests/unit/test_eval_ta.cpp`.

**Interfaces:**
- Consumes: A4's `batch_scope_level()`; B1's scope chain.
- Produces: nested scratches whose `parent_` follows the loop nest; invariants
  built once at their scope level.

- [ ] **Step 1: Read** the per-batch replay + scratch construction; note the
  exact point a descendant is currently rebuilt per batch (the
  `evaluate(node, le_g, bs.cache)` reinstall loop `:1539`/`:1700`).
- [ ] **Step 2: Write the failing test** (`test_eval_ta.cpp`): an `{i,k}`-batched
  network whose `I2`-class node is invariant to `i` (mirror the `[loop-tree]`
  network with real TA tiles). Count builds of the `I2` node (dry-run op trace or
  a build counter) and assert it is built **once**, not `nBatch(i)` times; assert
  the batched result equals the unbatched reference exactly (slice-exact).
- [ ] **Step 3: Implement.** At a batching node: `bs.cache.set_parent(&cache)`
  before the per-batch loop (wire the nest); for each loop-invariant descendant
  (`batch_scope_level()` strictly outer to this loop), `store` it at its scope
  level (walk up the chain to that cache) with its emitted effective count,
  instead of on the per-batch scratch. Build a whole invariant node on a
  **fresh** cache via `evaluate(n, le)` (variadic -- `evaluate(n, le, cache)`
  re-enters the batched evaluator and SIGSEGVs; spec S5 gotcha). Keep the seeding
  guard (`carries_ext`) intact for external modes.
- [ ] **Step 4: Pass** + reference-equality. No-regression:
  `./cmake-build-debug/tests/unit/unit_tests-sequant "[eval]"` counts unchanged
  on the OFF path.
- [ ] **Step 5: clang-format. Step 6: Commit.**

---

## Task 5: B3 -- effective-count lifetime (release at zero, scope-reset backstop)

Register invariants with the emitted effective count so they release the instant
the count hits zero, backstopped by the scope reset.

**Files:** `core/eval/cache_manager.hpp` (entry registration with count),
`core/eval/eval.hpp` (pass the emitted count at store),
`tests/unit/test_cache_manager.cpp`.

**Interfaces:**
- Consumes: A4's `batch_effective_count()`; the existing `entry`
  `decay()`/`max_life`/`reset()`.

- [ ] **Step 1: Write the failing test** (`test_cache_manager.cpp`): register a
  node with `effective_use_count = n`; assert it survives `n-1` `access` calls
  and releases on the `n`-th (count-zero); assert an under-read node (fewer reads
  than `n`) is cleared by `reset()` (backstop) and never returns a wrong value.
- [ ] **Step 2: Run to confirm it fails.**
- [ ] **Step 3: Implement** wiring the emitted count into entry registration
  (the count is the `entry` constructor's `count`, `:101`). The accumulator
  read-modify-write MUST NOT decrement (route it outside the counted `access`;
  spec S4) -- add a test-visible assertion or a distinct non-decrementing path if
  the evaluator currently routes RMW through `access`. Keep `persistent_` one
  bit.
- [ ] **Step 4: Pass. Step 5: Commit.**

---

## Task 6: B4 -- finish the veto reinterpretation (context-level cacheable)

The committed veto (`2a52e063c`) already narrows to sliced-contracted nodes.
Complete it so a hoisted (batch-invariant, non-run-scope) node remains cacheable
at its context level, while a genuinely run-scope-illegal node (held whole across
its own batched mode) is still refused.

**Files:** `core/eval/cache_manager.hpp` (the veto, `:555-620`),
`tests/unit/test_cache_manager.cpp`.

**Interfaces:** Consumes B2/B3 (nodes now stored at context level).

- [ ] **Step 1: Write the failing test:** a hoisted batch-invariant node at its
  context level must NOT be vetoed from that cache; a node whose result carries a
  mode sliced AT this node (held-whole-across-its-own-batched-mode) is still
  refused. Assert both against the veto predicate directly.
- [ ] **Step 2: Run to confirm** the current here-only veto mis-handles the
  hoisted case (fails).
- [ ] **Step 3: Implement:** restate the veto as "cannot be run-scope if
  batched" -- permit caching at the node's emitted scope level (consult
  `batch_scope_level()`); refuse only run-scope residence of a batched node.
  Keep OFF-path byte-identical (no annotation -> old behavior).
- [ ] **Step 4: Pass. Step 5: Commit.**

---

## Task 7: B5 -- GATE: avoidable_time on the C60 dry run (and split the peak assertion)

**Files:** `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`, the gate at
`:4222`, `:4509`, `:4580`; the failing asserts `:4548` peak, `:4574` avoidable).

**Interfaces:** Consumes Tasks 2-6 (and Task 1's finding on whether `cp`/emit
also needs a charge correction).

- [ ] **Step 1: Split the gate.** Separate the avoidable-time assertion (this
  plan's gate) from the peak assertion (external follow-on). Keep
  `CHECK(both.avoidable_time() < 0.05)` and `aux.avoidable_time() < 0.05` as the
  binding gate; convert the `both.peak_gb < 100` CHECK into a `WARN` (or a
  separately-tagged `[.][dryrun-ext-peak]` expected-fail) with an inline comment
  pointing to `2026-07-20-external-mode-batching-design.md` -- documenting that
  the 2302 GB transient is the external-forest lever, not this plan's. Do NOT
  loosen avoidable-time.
- [ ] **Step 2: Run the witness:**
  `./cmake-build-debug/tests/unit/unit_tests-sequant "[dryrun-occ-veto]"`.
  Assert `both.avoidable_time() < 0.05` (was 0.437) with the reported peak still
  sliced (compare against the batched peak, not the unbatched 38897 GB). If it is
  still high, consult Task 1's finding: if `cp`/emit was phantom-priced, the
  runtime hoist alone may not register -- fix the emit charge per Task 1 before
  re-running (this is the one place the plan may loop back).
- [ ] **Step 3: No-regression:**
  `./cmake-build-debug/tests/unit/unit_tests-sequant "[optimize][eval][dryrun-objective]"`
  byte-identical; the `peak_threshold=+inf` heuristic path unchanged.
- [ ] **Step 4: Record** before/after `avoidable_time` (both legs), the sliced
  peak, and the per-iteration modeled time in the commit message. **Commit.**

---

## Task 8: C1 -- MPQC flag plumbing + repin + local validation

`order_aware_recompute` is currently reachable only from SeQuant test code; MPQC
never sets it. Plumb it, repin, validate locally.

**Files:** MPQC `src/mpqc/chemistry/qc/lcao/cc/csv_batch_policy.h` (and the
`BatchPolicy`/`CostParams` thread into `optimize`), `mpqc4/external/versions.cmake`
(own commit).

- [ ] **Step 1:** Thread `order_aware_recompute` from an MPQC CCk keyword through
  `BatchPolicy` -> `CostParams` (mirror how `batch_spectator_indices` /
  `peak_threshold` are threaded in `csv_batch_policy.h`). Default OFF (byte-
  identical). Build MPQC against local SeQuant (`FETCHCONTENT_SOURCE_DIR_SEQUANT`).
- [ ] **Step 2: Local validation:** run a CSV-CCSD sanity (e.g. an existing
  small CSV/PNO-CCk validation input) with the flag ON; assert converged energy
  is unchanged vs the flag OFF and vs unbatched (batching is a memory schedule,
  never an approximation).
- [ ] **Step 3: Repin** `MPQC_TRACKED_SEQUANT_TAG` in `external/versions.cmake`
  (update `PREVIOUS_*` too), its own commit paired with the MPQC code that
  requires it.
- [ ] **Step 4: Hand off** Owl C60 occ+aux validation (separate SLURM session,
  refresh the TA pin): per-iteration time down, energy unchanged. Peak-under-
  budget remains the external-forest follow-on, not this plan's gate.

---

## Notes for the executor

- **Task 1 is investigation-gated on purpose.** The emit's own debug block uses
  the OLD set `esc`, and the witness shows a `2.24x` replay/model gap -- settle
  whether the runtime hoist alone closes avoidable-time or whether `cp`/emit is
  also phantom-priced BEFORE building the runtime, so Task 7 does not chase a
  cost-model bug from the runtime side.
- **Every batched change must show the `[dryrun-occ-veto]` peak stays sliced** --
  a hoist that caches a value whole across its batched mode defeats the memory
  saving and is caught by the reported peak (even though peak-under-budget is not
  this plan's gate).
- **Phase A is done; do not touch `relax`/selection.** This plan only emits what
  the DP already computed and makes the runtime obey it.
- **The OFF path is the safety net.** `order_aware_recompute=false` must remain
  byte-identical at every task boundary -- it is how the non-batched and legacy
  heuristic paths stay untouched.
