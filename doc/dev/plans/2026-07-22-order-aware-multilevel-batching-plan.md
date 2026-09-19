# Order-aware multilevel batching Implementation Plan

> **SUPERSEDED (2026-07-25) by `2026-07-25-order-aware-multilevel-batching-runtime-plan.md`.**
> A freshness pass against the branch found this plan's **Phase A (order-aware
> cost model) already implemented, committed, and unit-tested** (ordered-cell
> layer `ac3f58bfd`/`388dcb7aa`/`e3198877d`/`049179059`/`83028b0d4`;
> `escaped_outer`+`descend`+resident-scan wired into `relax`, gated by
> `order_aware_recompute`; `[loop-tree]` DP tests 7/7). Its Phase A tasks
> (A1-A3) would re-derive committed code, and its "critical path" framing is
> inverted -- the real remaining surface is the **annotation bridge (A4) +
> runtime hoisting (Phase B)**. Do NOT execute this plan; use the 2026-07-25
> runtime plan. Kept for history and for the Phase-B mechanism narrative it
> shares with the spec.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the SeQuant factorizer price and realize batched contraction schedules by their true subtree-bound loop-tree recompute cost, so the C60 PNO-CCSD residual's contracted-occ middle-gap recompute drops from 76% to a few percent with peak still sliced.

**Architecture:** Two coupled halves joined by one annotation contract. **Phase A (cost model):** the DP (`PeakBatchedModel` in `core/optimize/cost_model.hpp`) charges each node's recompute from the loop tree induced by the contraction tree (a contracted index's batch loop wraps its contraction subtree; a node invariant to an enclosing loop is either hoisted -- resident -- or recomputed -- `flops x nBatch` -- whichever fits `peak <= budget`), replacing the order-blind `esc = B & ~open_modes[n]`. It emits per node `{lifetime-scope, effective-use-count, pinned co-contracted order}`. **Phase B (runtime):** the batched evaluator (`core/eval/eval.hpp`) obeys the annotation via the committed scope-chain primitive in `core/eval/cache_manager.hpp` (`parent_`/`set_parent`/fall-through `access`, `75d240079`): store each invariant at its emitted scope level, read up the chain, release by effective count.

**Tech Stack:** C++20, SeQuant, Catch2 unit tests, the dry-run eval backend (`core/eval/backends/dryrun/`). Build in `cmake-build-debug`; the dry-run instrument replays the real C60 forest in seconds.

## Global Constraints

- **Success gate (verbatim from spec S7):** `[.][dryrun-occ-veto]` asserts `aux.avoidable_time() < 0.05` AND `both.avoidable_time() < 0.05` (fails today at 0.75), **with peak still sliced, not materialized whole**.
- **No regression:** unbatched (~0%) and aux-only (~0.48%) stay clean; the `peak_threshold = +inf` / no-annotation legacy heuristic path stays **byte-identical**.
- **Contracted-index batching only** is in scope. External-index optimal loop placement is the forest-peak follow-on (Gap 1); this plan *prices* external modes but does not search their placement.
- **`charge_batch_recompute = false` and the unbatched path must remain byte-identical** (the new charge is active only when batching is engaged).
- **Co-contracted loop order** among indices contracted at one node is cost-neutral but MUST be pinned by a canonical deterministic rule (ascending batched-slot ordinal) -- for schedule determinism and bit-reproducibility.
- **`persistent_` stays one bit** (P = run-scope, the sole no-decay case).
- **Cross-repo:** any behavior change needs an MPQC repin; validate on local `[.][dryrun-occ-veto]` before any Owl run. No `Co-Authored-By` trailers. ASCII hyphens only (pre-commit U+2013 detector).

---

## File structure

- `core/optimize/cost_model.hpp` -- `PeakBatchedModel::relax` (recompute charge), the per-`B` frontier `State`, `select_root`, `reconstruct_batched_modes` (emit), `is_external_mode` (already distinguishes contracted vs external). All Phase A changes live here.
- `core/eval/eval_expr.{hpp,cpp}` -- the per-node stamp (`batched_here()` / `set_batched_here`); extended with the scope level, effective count, pinned co-contracted order.
- `core/eval/cache_manager.hpp` -- the committed scope-chain primitive (`parent_`, `set_parent`, entry decay); Phase B wires it as a consumer.
- `core/eval/eval.hpp` -- `make_batched_custom_evaluator` / `make_batched_scratch`; Phase B stores invariants at their scope level and sets scratch parents.
- `tests/unit/test_optimize.cpp` -- DP-level unit tests (loop-tree recompute charge on small hand-built networks).
- `tests/unit/test_eval_dryrun.cpp` -- the `[.][dryrun-occ-veto]` witness (Phase A checkpoint on modeled cost; Phase B final gate on replayed `avoidable_time`).
- `tests/unit/test_cache_manager.cpp` -- scope-chain visibility + effective-count lifetime unit tests.

---

## Phase A -- order-aware cost model

### Task A1: Investigate the loop-tree recompute charge (INVESTIGATION)

The spec (S3) fixes the *algorithm* (a node's `rf` = product of `nBatch` over enclosing loops it does not carry and is not hoisted above, read off the loop tree). The *exact* code change depends on how today's per-`B` frontier ties `B` to the tree's ancestor structure -- this task settles that by reading, before any edit.

**Files (read-only):** `core/optimize/cost_model.hpp` -- `relax` (~888-958, the `esc`/`rf` charge and the `C = B | Ap` descent), the `State`/`BFrontPoint` per-`B` frontier and `pareto_insert`, `sz`/`Lof`/`open_modes`/`nbatches` (~744-802), `is_external_mode` (~974-989), `select_root` (~1005+), `reconstruct_batched_modes` (the emit walk).

**Deliverable** (write to `.superpowers/sdd/oamb-a1-note.md`): answer, with exact function+line references,
1. Does the frontier index `B` at a subset already equal the set of batched modes contracted at-or-above that subset (i.e. is `B` at a node its true enclosing-loop set, via `C = B | Ap`), or can `B` float free of the tree? This decides whether the charge fix is "use the descent-propagated enclosing set instead of a free `B`" or "add enclosing-set tracking."
2. Where exactly is `rf` charged (`relax` line), and what is the minimal change so `rf` counts only `nBatch(x)` for enclosing loops `x` the node does not carry AND cannot be hoisted above under the peak budget -- i.e. how the hoist-vs-recompute choice enters (a per-node min over "resident peak add" vs "flops x nBatch", both already representable in the `(peak, flops)` frontier point).
3. How `is_external_mode` gates the charge: contracted modes get the subtree-bound charge; external modes are priced but their placement is not searched here.
4. The concrete shape of the emit: what `reconstruct_batched_modes` writes per node today, and where `{lifetime-scope, effective-count, co-contracted order}` attach.

- [ ] **Step 1:** Read the functions above; write the note answering (1)-(4) with exact line references and the minimal change points.
- [ ] **Step 2:** Confirm the note does not require changing the unbatched / `charge_batch_recompute=false` path (must stay byte-identical). Report DONE with the note path.

### Task A2: RED test -- the DP mis-prices the gC-class middle gap

Add a DP-level test on a *small hand-built* network isomorphic to the gC middle gap, so it runs in milliseconds and pins the charge without the full C60 forest.

**Files:** `tests/unit/test_optimize.cpp`.

**Interfaces:**
- Produces: a `[loop-tree]`-tagged Catch2 test asserting the modeled recompute factor of the invariant-carrying node.

The network: `R{i,m} = (A{i,k} B{k})(C{k,l} D{l,m})` reduced so one intermediate carries an aux-like batched mode `k` but not a batched contracted occ-like mode. Concretely reuse the spec's `R{i,m}` example with batched set `{k, i}` where `I2{k,m}` is invariant to `i`.

- [ ] **Step 1: Write the failing test.** Build the network via the existing `test_optimize.cpp` helpers (follow an existing `[optimize]` test for the network-construction idiom). Batch `{k, i}` with per-mode `nBatch` (k->1, i->8 via the batch-target lambda used by neighboring batched tests). Assert the DP's modeled cost for the schedule reflects `I2` being **recomputed `nBatch(i)` times** when it cannot be hoisted, i.e. the reported flops include the `I2`-rebuild term -- NOT the `rf=1` phantom. State the expected number from the roofline of `I2`'s contraction times `nBatch(i)`.

```cpp
// tests/unit/test_optimize.cpp  (sketch; use the file's actual network builder)
TEST_CASE("loop-tree recompute charge prices the middle gap", "[optimize][loop-tree]") {
  auto net = /* R{i,m} = (A{i,k} B{k})(C{k,l} D{l,m}) via existing builder */;
  auto const cost_phantom = model_cost(net, /*batch=*/{"k","i"}, /*charge=*/true);
  // gC/I2-class node is invariant to i; today esc drops i => rf=1 (phantom).
  // After the fix it is recomputed nBatch(i) times.
  CHECK(cost_phantom.flops >= expected_with_I2_rebuilt);  // FAILS today
}
```

- [ ] **Step 2: Run to confirm it fails.** `./tests/unit/unit_tests-sequant "[loop-tree]"` -- FAILS (today's `esc` prices `rf=1`, flops below `expected_with_I2_rebuilt`).

### Task A3: Implement the loop-tree recompute charge

Per Task A1's note, replace the order-blind `esc` charge in `relax` with the subtree-bound loop-tree charge + per-node hoist-vs-recompute, gated by `is_external_mode` (contracted modes only).

**Files:** `core/optimize/cost_model.hpp` (`relax`, and whatever A1 identifies as the enclosing-set source).

**Interfaces:**
- Consumes: A1's change note.
- Produces: a `relax` that, for a node inside enclosing contracted-mode loops it does not carry, either charges `flops x nBatch` (recompute) or a resident-peak add (hoist), whichever keeps the frontier point feasible/cheapest -- so the frontier already carries both realizations and `select_root` picks per budget.

- [ ] **Step 1:** Implement the charge per A1. Keep `charge_batch_recompute=false` and the unbatched path byte-identical (guard the new charge under the existing flag / batching-engaged condition).
- [ ] **Step 2: Run to confirm A2 passes.** `./tests/unit/unit_tests-sequant "[loop-tree]"` -- PASS.
- [ ] **Step 3: No-regression.** `./tests/unit/unit_tests-sequant "[optimize][dryrun-objective]"` -- assertion counts unchanged vs pre-change (the unbatched/legacy paths byte-identical). Record counts.
- [ ] **Step 4: clang-format** the changed file; **Step 5: Commit.**

### Task A4: Emit the annotation `{lifetime-scope, effective-use-count, co-contracted order}`

**Files:** `core/optimize/cost_model.hpp` (`reconstruct_batched_modes`), `core/eval/eval_expr.{hpp,cpp}` (the stamp).

**Interfaces:**
- Consumes: the chosen schedule + loop tree from A3.
- Produces: per-node stamp fields -- `lifetime_scope_level` (the outermost enclosing loop level containing all the node's carried batched modes), `effective_use_count` (`Sigma consumers Prod nBatch(L)` over loops L enclosing the consumer but not the node), and `co_contracted_order` (ascending batched-slot ordinal of the modes contracted at that node). Accessors alongside `batched_here()`.

- [ ] **Step 1: Write the failing test** (`tests/unit/test_optimize.cpp`, `[loop-tree]`): on the A2 network, assert `I2`'s emitted `lifetime_scope_level` is the `k`-level (above the `i` loop) and its `effective_use_count == nBatch(i)`; assert a node contracting two batched indices emits `co_contracted_order` in ascending slot ordinal.
- [ ] **Step 2: Run to confirm it fails** (fields absent).
- [ ] **Step 3: Implement** the stamp fields + populate them in `reconstruct_batched_modes` (extend the existing `set_batched_here` write; no new emit channel).
- [ ] **Step 4: Run to confirm it passes.** No-regression on `[optimize]`. **Step 5: clang-format + Commit.**

### Task A5: Phase-A checkpoint -- DP predicts the true cost on the C60 witness

**Files:** `tests/unit/test_eval_dryrun.cpp` (a `[.]`-hidden `[dryrun-loop-tree-cost]` case reusing the `[dryrun-occ-veto]` forest builder).

**Interfaces:** Consumes A3/A4.

- [ ] **Step 1: Write the test:** on the C60 occ+aux forest, with the loop-tree charge on, assert the **DP-modeled** recompute of the `gC`-class term is `~nBatch(i_3)` (not `rf=1`), and that `I2`-class invariants are modeled as hoisted (scope above their invariant loop). This is measurable from the DP's cost profile alone -- no runtime.
- [ ] **Step 2: Run** (hidden tag, explicit): `./tests/unit/unit_tests-sequant "[dryrun-loop-tree-cost]"`. Record the modeled `gC` rf and the hoisted-set.
- [ ] **Step 3: Commit.** This is the Phase-A gate: the cost model now prices the schedule correctly before any runtime work.

---

## Phase B -- runtime lifetime management

> Phase B consumes A4's emitted fields verbatim (Decision 0: the runtime never re-derives them). Each Phase B task begins by reading the specific `eval.hpp` region it touches (the batched evaluator internals are not restated here); the design in `2026-07-20-order-aware-multilevel-batching-design.md` S5 and `2026-07-17-nested-batch-group-join-design.md` gives the mechanism.

### Task B1: Scope-chain visibility -- read up on a local miss

**Files:** `core/eval/cache_manager.hpp` (the map-level `access(key)`), `tests/unit/test_cache_manager.cpp`.

**Interfaces:**
- Consumes: the committed `parent_`/`set_parent` primitive.
- Produces: `access(key)` returns a hit found in an ancestor cache when absent locally (fall-through already present at the map level per `75d240079`; this task *tests and hardens* it and confirms `store` still targets one level).

- [ ] **Step 1: Write the failing test:** two `CacheManager`s, inner `set_parent(&outer)`; `store` a node in `outer`; assert `inner.access(node)` returns it (walk-up), and that `inner`'s own `reset()` does not clear `outer`'s entry.
- [ ] **Step 2: Run to confirm** (fails if the fall-through is incomplete). **Step 3: Implement/harden** the fall-through in map-level `access`. **Step 4: Pass. Step 5: Commit.**

### Task B2: Store each invariant at its emitted scope level

**Files:** `core/eval/eval.hpp` (`make_batched_custom_evaluator` / `make_batched_scratch`), `tests/unit/test_eval_ta.cpp`.

**Interfaces:** Consumes A4's `lifetime_scope_level`; B1's scope chain.

- [ ] **Step 1: Read** the batched evaluator's per-batch replay + scratch construction; note where a descendant is currently rebuilt per batch.
- [ ] **Step 2: Write the failing test** (`test_eval_ta.cpp`): an `i,k`-batched network whose `I2`-class node is invariant to `i`; assert the node is built **once** (not `nBatch(i)` times) by counting builds (dry-run op trace or a build-counter), and the batched result equals the unbatched reference (slice-exact).
- [ ] **Step 3: Implement:** at a batching node, set each fresh scratch's `parent_` to the enclosing cache (B1), and `store` each loop-invariant descendant at its emitted scope level with its emitted effective count, instead of on the per-batch scratch. Build a whole invariant node on a **fresh** cache via `evaluate(n, le)` (variadic; `evaluate(n, le, cache)` re-enters the batched evaluator and SIGSEGVs -- spec S5 gotcha).
- [ ] **Step 4: Pass** + reference-equality. **Step 5: Commit.**

### Task B3: Effective-count lifetime (release at zero, scope-reset backstop)

**Files:** `core/eval/cache_manager.hpp` (register entries with the emitted count), `core/eval/eval.hpp` (pass the count), `tests/unit/test_cache_manager.cpp`.

**Interfaces:** Consumes A4's `effective_use_count`; the existing entry `decay()`/`max_life`/`reset()`.

- [ ] **Step 1: Write the failing test:** register a node with `effective_use_count = nBatch(i)`; assert it survives `nBatch(i)-1` reads and releases on the `nBatch(i)`-th (count-zero), and that an under-read node is cleared by scope `reset()` (backstop), never a wrong value.
- [ ] **Step 2: Run to confirm. Step 3: Implement** wiring the emitted count into entry registration (the accumulator RMW must NOT decrement -- route it outside the counted `access`, spec S4). Keep `persistent_` one bit.
- [ ] **Step 4: Pass. Step 5: Commit.**

### Task B4: Reinterpret the free-batchable veto

**Files:** `core/eval/cache_manager.hpp` (the veto, ~free-batchable eviction), `tests/unit/test_cache_manager.cpp`.

**Interfaces:** Consumes B2/B3.

- [ ] **Step 1: Write the failing test:** a hoisted batch-invariant (non-run-scope) node at its context level must remain cacheable -- assert the veto does NOT evict it, while a genuinely run-scope-illegal (held-whole-across-its-own-batched-mode) node is still refused.
- [ ] **Step 2: Run. Step 3: Implement:** narrow the veto to "cannot be run-scope if batched," permitting caching at the node's context level. **Step 4: Pass. Step 5: Commit.**

### Task B5: GATE -- avoidable_time on the C60 dry run

**Files:** `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`).

**Interfaces:** Consumes A3-B4.

- [ ] **Step 1: Run the witness:** `SEQUANT_UT_DRYRUN_NTERMS=<full> ./tests/unit/unit_tests-sequant "[dryrun-occ-veto]"`. Assert `aux.avoidable_time() < 0.05` AND `both.avoidable_time() < 0.05` (was 0.75), **and the reported peak is still sliced** (not materialized whole -- compare against the batched peak, not the unbatched footprint).
- [ ] **Step 2: No-regression:** `./tests/unit/unit_tests-sequant "[optimize][eval][dryrun-objective]"` byte-identical; the `peak_threshold=+inf` heuristic path unchanged.
- [ ] **Step 3: Record** the before/after `avoidable_time`, peak, and per-iteration modeled time in the commit message. **Commit.**

---

## Phase C -- cross-repo enablement

### Task C1: MPQC repin + local validation

**Files:** `mpqc4/external/versions.cmake` (own commit), local CSV/PNO-CCk path.

- [ ] **Step 1:** Build MPQC against the local SeQuant (`FETCHCONTENT_SOURCE_DIR_SEQUANT`); run a local CSV-CCSD sanity (converged energy unchanged vs unbatched).
- [ ] **Step 2:** Repin `MPQC_TRACKED_SEQUANT_TAG` in `external/versions.cmake` (update `PREVIOUS_*` too), its own commit paired with the code that requires it.
- [ ] **Step 3:** Hand off Owl C60 occ+aux validation (separate SLURM session, refresh the TA pin): per-iteration time down, energy unchanged. (Peak-under-budget is the forest-peak follow-on, not this plan's gate.)

---

## Notes for the executor

- **The DP-charge core (A3) is investigation-gated (A1) on purpose** -- the spec fixes the algorithm but the minimal code change depends on whether the per-`B` frontier already carries the enclosing-loop set. Do A1 before A3; if A1 finds `B` does NOT track the tree ancestors, escalate (the change is larger than a charge swap) rather than force it.
- **Phase A is the load-bearing, dry-run-measurable half.** A5 is a real gate: if the DP still prices `gC` at `rf=1`, stop -- Phase B cannot fix a blind plan.
- **Every batched change must show the `[dryrun-occ-veto]` peak stays sliced** -- a fix that restores sharing by caching a value whole across its batched mode defeats the memory saving and is caught by the reported peak.
