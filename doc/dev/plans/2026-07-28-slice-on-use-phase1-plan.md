# Slice-on-use (phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the batched evaluator slice a cached intermediate to the current
batch block when it is *used* in a deeper scope, so the C60 external-occ replay
peak stops diverging from the DP's cost model. **RESULT (achieved, SeQuant
`ac4f62377`): peak 5860.877 -> 443.55 GB (13.2x); the replay now MATCHES the model
(summand 12 replay 443.55 ~= model 445.46), eliminating the 167x gap.** (The "~35
GB" in earlier drafts was summand 46's per-term model, NOT the forest peak; the
forest modeled peak is ~445 GB, summand 12. Reducing 443 -> <100 is phase 2.)

**Architecture:** Implements Sec. 2.2c of
`doc/dev/specs/2026-07-27-external-loop-scope-semantics-design.md`. Today `le_g`
slices only *leaf* results; a cached intermediate fetched via the scope-chain
`access` fallback is served **whole**, so `s*C` (a full 8.7 GB top-cache entry)
contracts into the full 2930 GB giant. The fix: at the single value-return point
of `evaluate`, slice whatever was fetched to `(use scope) MINUS (lifetime scope)`
of the modes the node carries -- reading the *observed* lifetime scope from where
`access` resolved, so it is correct regardless of current placement.

**Tech Stack:** C++20, SeQuant, Catch2, the dry-run eval backend
(`core/eval/backends/dryrun/`). **Build/run in `cmake-build-release`** -- the
factorizer is too slow in Debug; the eval unit tests compile in Release.

## Global Constraints

- **Vocabulary** (per spec): `leaf_evaluator` (raw leaf compute, the `le` rename);
  `materialize(node) -> {value, lifetime_scope}` (unified cache-hit-or-leaf
  accessor, surfacing the OBSERVED resolved scope); `batch_context` (ordered
  enclosing `{axis, block}` stack); `slice_to_use(value, node, lifetime_scope)`
  (slice for the (use MINUS lifetime) modes the node carries, via `slice_mode`).
  `scope` = an axis-set ID; `scope_level` = a depth (int).
- **Slice from the OBSERVED lifetime**, not an emitted stamp: `materialize` reports
  where the value actually resolved on the scope chain (top/full for a freshly
  materialized leaf; the hit level for a cached intermediate). Phase 1 does NOT
  compute or change any placement -- it slices whatever the current machinery
  cached, wherever it cached it.
- **OFF-path byte-identical.** With `order_aware_recompute = false` OR
  `batch_spectator_indices = false` there is no external batch loop, so
  `batch_context` is empty and `slice_to_use` is a no-op -- byte-identical.
- **Slice-exact.** The batched result must equal the unbatched reference exactly;
  the C60 replay energy is unchanged.
- **Do NOT touch phase-2 concerns**: no mask computation, no placement changes, no
  removal of the veto / `hoist_invariants` / `consumed_upscope` / `scope_level`
  emit. Those stay in place and become redundant-but-harmless once slicing is
  correct; phase 2 retires them.
- **Process:** branch `evaleev/feature/suppress-heuristic-fallback`; commit there;
  ONE task per commit; NO push / PR / remote / repin. No `Co-Authored-By` trailers.
  ASCII hyphens only (pre-commit rejects U+2013 / U+00A0). clang-format is a
  pre-commit hook; re-add + re-commit if it reformats. Cap builds at `-j6`.

---

## File structure

- `core/eval/eval.hpp` -- `evaluate` Enter stage (hit `:619`, leaf `:659`),
  `make_batched_custom_evaluator` and the scatter/group branches where `le_g` is
  built and re-installed. Home of `materialize`, `slice_to_use`, the
  `batch_context` threading, and the `le -> leaf_evaluator` rename.
- `core/eval/cache_manager.hpp` -- `access` (`:276-287`): gains a form that
  surfaces the resolved chain level (the lifetime scope).
- `tests/unit/test_eval_ta.cpp` (`[eval]`) -- runtime slice-exactness test.
- `tests/unit/test_cache_manager.cpp` (`[cache_manager]`) -- access-reports-level test.
- `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`) -- C60 peak gate.

---

## Task 1: Rename `le -> leaf_evaluator`

Mechanical, behavior-neutral rename so the accessor's name is honest before the
substantive change lands (keeps the Task-3 diff about slicing, not renaming).

**Files:** `core/eval/eval.hpp` (the `evaluate` overloads' `F const& le` /
`F le` parameters at `:545,:810,:879,:937,:1362` and every in-body use +
Doxygen `\param le`), any other `.hpp`/`.cpp` in `core/eval/` that names the
parameter `le`.

**Interfaces:**
- Produces: the leaf-evaluator parameter is named `leaf_evaluator` everywhere; the
  concept `meta::leaf_node_evaluator` and the `DryRunLeafEvaluator` type are
  unchanged.

- [ ] **Step 1: Rename.** `grep -rn "\ble\b" core/eval/*.hpp core/eval/*.cpp` to
  find the parameter/local uses (it is the leaf-evaluator param and its
  invocations `le(...)`; do NOT touch unrelated identifiers). Rename the parameter
  and its uses to `leaf_evaluator`, and update the `\param le` Doxygen lines.
  Leave `le_g` alone (it is deleted in Task 3).
- [ ] **Step 2: Build.** `cmake --build cmake-build-release --target unit_tests-sequant -j6`.
- [ ] **Step 3: Run to confirm behavior-neutral.**
  `./cmake-build-release/tests/unit/unit_tests-sequant "[eval]~[tot]"` and
  `"[optimize]"` and `"[cache_manager]"` -- assertion counts unchanged from before
  the rename (record them).
- [ ] **Step 4: clang-format + Commit.**

---

## Task 2: INVESTIGATION -- pin the runtime edit points for slice-on-use

The Sec. 2.2c mechanism is designed; the *exact* edits depend on details of the
scope chain and the batched custom evaluator that must be read before touching
them. Settle them, do NOT change behavior.

**Files (read-only):** `core/eval/eval.hpp` -- the `evaluate` Enter stage
(`:614-672`: hit `:619`, custom-eval `:632`, leaf `:659`),
`make_batched_custom_evaluator` (`:1361` + its params) and the scatter (`:1564-1650`,
`le_g` at `:1635`) and group (`:1652+`) branches; `core/eval/cache_manager.hpp` --
`access` (`:276-287`), `parent_`/`set_parent`, `reset`.

**Deliverable** (write to `.superpowers/sdd/slice-on-use-t2-note.md`), with exact
file:line refs and the exact edit shape:

1. **How `access` surfaces the lifetime scope.** `access` currently returns
   `ResultPtr`, recursing up `parent_`. Specify how it reports WHERE it resolved
   (e.g. return `{ResultPtr, hops}` or the resolved `CacheManager*`, or a depth).
   The consumer needs enough to know which enclosing loops are already baked into
   the fetched value vs. not. State the representation and the exact `access`
   change (keep the old signature working for non-batched callers, or migrate them).
2. **`batch_context` threading.** The enclosing `{axis, block}` per realized loop
   lives implicitly in the composed `le_g` closures + `depth`. Specify how to make
   it explicit: a stack threaded into `make_batched_custom_evaluator` (a new param),
   pushed with `{K, [e_lo,e_hi]}` at the scatter (`:1602`) and the group
   (`BatchIter`) before the `depth+1` re-entry. Confirm it composes (each level
   appends its `{axis, block}`), and that on the OFF path it stays empty.
3. **`materialize` + `slice_to_use` + Enter-stage collapse.** Give the exact
   rewrite of the Enter stage: `materialize(node)` = cache-hit (`:619`) else
   leaf-eval (`:659`) else null, returning `{value, lifetime_scope}`; then
   `slice_to_use(value, node, lifetime_scope)` slices, for each `{axis, block}` in
   `batch_context` NOT already baked at `lifetime_scope` that `node` carries (via
   `index_position`), by `value->slice_mode(pos, lo, hi)`. Confirm order vs.
   `apply_phase`. Confirm a freshly materialized leaf's lifetime scope is
   "top/full" so `slice_to_use` reproduces exactly what `le_g` did for leaves
   (i.e. Task 3 is behavior-preserving for leaves).
4. **`le_g` deletion.** Confirm every `le_g` use (`:1616-1622` slicing, `:1627`
   re-install, `:1635`, `:1795`) becomes: pass the RAW `leaf_evaluator` to the
   re-installed evaluator, and rely on `slice_to_use` at the Enter stage for
   slicing. Confirm the leaf slicing that `le_g` did is fully covered by
   `slice_to_use` (leaf lifetime = top/full => sliced for all enclosing carried
   loops => identical to `le_g`).
5. **OFF-path proof.** Confirm that with `batch_context` empty, `materialize` +
   `slice_to_use` reduce to the current hit/leaf behavior byte-for-byte.

- [ ] **Step 1:** Read the functions; write the note answering (1)-(5) with exact
  edit points and the level/block-range representation decisions.
- [ ] **Step 2:** Report DONE with the note path. If (1) or (4) reveals the leaf
  slicing is NOT exactly reproduced by `slice_to_use` from a top/full lifetime
  (so deleting `le_g` would change leaf behavior), STOP and escalate -- that is a
  design gap, not an implementation detail.

**Escalation:** if `access` cannot surface a usable lifetime scope without a
larger refactor of the scope chain, STOP and report -- the representation is the
crux and a wrong choice cascades.

---

## Task 3: Implement slice-on-use

Per Task 2's note: `access` surfaces the lifetime scope; thread `batch_context`;
add `materialize` + `slice_to_use`; collapse the Enter stage; delete `le_g`.

**Files:** `core/eval/eval.hpp`, `core/eval/cache_manager.hpp`,
`tests/unit/test_eval_ta.cpp`, `tests/unit/test_cache_manager.cpp`.

**Interfaces:**
- Consumes: Task 2's note (exact edits), Task 1's `leaf_evaluator`.
- Produces: a runtime where a cached intermediate carrying a batched mode is sliced
  to the current block on use.

- [ ] **Step 1: Write the failing `[cache_manager]` test.** Build `outer -> mid ->
  inner` chained caches (`mid.set_parent(&outer)`, `inner.set_parent(&mid)`), store
  a key `X` in `outer`. Assert `inner.access(X)` (the level-surfacing form) reports
  the resolved level = `outer`'s level (the value came from 2 hops up), and
  `mid.access(X)` reports 1 hop. (Follow the existing `set_parent` test idiom.)
- [ ] **Step 2: Run to confirm it fails** (the level-surfacing form does not exist).
  `./cmake-build-release/tests/unit/unit_tests-sequant "[cache_manager]"`.
- [ ] **Step 3: Write the failing `[eval]` slice-exactness test** (`test_eval_ta.cpp`).
  Build a network `R{e} = M{e,c} * B{c}` where `M{e,c}` is a shared intermediate
  cached at an OUTER scope (full over `e`) and consumed inside the `e`-external
  scatter. Configure `order_aware_recompute = true`, `batch_spectator_indices =
  true`, `e` external-batched (block < full). Replay and assert: (a) `M` is
  realized at BLOCK-`e` extent at the consumer (no full-`e` contraction) -- via the
  dry-run op trace / a size probe; (b) the batched result equals the unbatched
  reference exactly (`norm2(ref - res) < 1e-12`). Use the `[eval]` replay idiom in
  `test_eval_ta.cpp` (copy an existing external-batched / `batch_spectator_indices`
  harness; if none, model the network on the `test_optimize.cpp` `[loop-tree]`
  builder and evaluate via the `[eval]` path).
- [ ] **Step 4: Run to confirm it fails** (today the cached `M` is served whole).
  `./cmake-build-release/tests/unit/unit_tests-sequant "[eval]"`.
- [ ] **Step 5: Implement per Task 2's note** (`.superpowers/sdd/slice-on-use-t2-note.md`
  is the authoritative edit spec -- follow it exactly). Summary: `access_at(key)
  -> {ptr, hops}` on `CacheManager` (keep `access` as a forwarder for the 3
  non-batched callers); `batch_context_` member on `CacheManager` (getter/setter
  mirroring `custom_evaluator_`), pushed `{K,{lo,hi}}` per block at scatter/group
  before the `depth+1` re-entry; `slice_to_use(value,node,hops)` slicing the `hops`
  innermost `batch_context` modes the node carries; Enter-stage hit -> `slice_to_use
  (apply_phase(ptr),node,hops)`, leaf -> store FULL then return
  `slice_to_use(stored,node,batch_context.size())`, custom-eval branch NOT sliced.
  **`le_g` is REPLACED, not deleted:** define `sliced_leaf = slice_to_use o
  leaf_evaluator` over `batch_context` and use it at the three internal sites
  (`pick_sliceable` probe, `carrier_full`, hoist build) that bypass the Enter stage
  -- passing the RAW evaluator there breaks the "K-not-re-picked" invariant (Q4).
  Keep the OFF path (empty `batch_context`) byte-identical.
- [ ] **Step 6: Run both new tests to confirm they pass** + reference-equality.
  No-regression: `"[eval]~[tot]"`, `"[cache_manager]"`, `"[optimize]"` unchanged on
  the OFF path -- record counts.
- [ ] **Step 7: clang-format + Commit** (note the `le_g`-deletion + slice-on-use in
  the message).

---

## Task 4: GATE -- C60 external-occ peak drops

**Files:** `tests/unit/test_eval_dryrun.cpp` (`[.][dryrun-occ-veto]`, external-occ
config: occ external, aux `K` contracted, `order_aware_recompute = true`,
`batch_spectator_indices = true`, `peak_threshold = 100 GB`).

**Interfaces:** Consumes Task 3.

- [ ] **Step 1: Run the witness** at full term count:
  `SEQUANT_UT_DRYRUN_NTERMS=55 ./cmake-build-release/tests/unit/unit_tests-sequant "[dryrun-occ-veto]"`.
  Assert the summand-46 `s*C` operand is served sliced (no `result=` near the
  2930 GB full giant / 8.7 GB full `s*C`), and the batched-scratch peak
  (`cost_profile(...).peak_bytes`) drops from 5860.877 GB to MATCH the DP model.
  **DONE: 5860.877 -> 443.55 GB; peak is now summand 12, replay 443.55 ~= model
  445.46 (the DP's honest cost) -- the 167x replay-vs-model gap is gone.** ("~35
  GB" was summand 46's per-term model, not the forest peak.) Record before/after.
- [ ] **Step 2: No-regression.**
  `./cmake-build-release/tests/unit/unit_tests-sequant "[optimize][eval][cache_manager][dryrun-objective]"`
  -- OFF-path byte-identical; the `peak_threshold = +inf` heuristic path unchanged.
- [ ] **Step 3: Record** before/after peak + `s*C`'s realized extent in the commit
  message. **Commit.** (If the peak does NOT drop to match the DP model, this is a
  MEASUREMENT gate -- record the honest number + probe evidence, report
  DONE_WITH_CONCERNS; do NOT fudge a green.)

---

## Phase-2 preview (NOT this plan -- its own spec+plan after phase 1 lands)

Recorded so the model is not lost (spec Sec. 10). Phase 2 = mask-based placement +
cost-model lifetime-awareness + retirement:

- **Per-slot lifetime mask.** Compute, per canonical node, a proto-aware per-slot
  SLICED/FULL bitmask by a meet over occurrences (proto-awareness is a *correctness*
  matter -- a mis-marked node lands in the wrong tier). The IR survey found this is
  dominated by within-occurrence structure, so it is per-occurrence-cheap.
- **Placement by rule.** all-full -> top; all-sliced -> deepest scratch; mixed ->
  deepest-sliced-mode scratch (full over the rest = hoisting). Tightens top-cache
  membership (only block-agnostic values), resolving the pair-specific-data hazard.
- **Cost model lifetime-aware.** The DP runs the SAME meet to size intermediates
  (residency = peak = cost); order choice must keep sliced modes outer to
  invariance loops.
- **Retire** the veto, `hoist_invariants`, `consumed_upscope`, and the
  `scope_level` emit (revert Tasks 1-3 of the prior effort) -- all subsumed by
  placement-by-rule + slice-on-use -- with a perf no-regression check.

---

## Notes for the executor

- **Phase 1 is a correctness fix bounded by a hard gate** (C60 peak drops,
  slice-exact). It does not touch placement; it slices from the observed lifetime.
- **`le_g` deletion is load-bearing** for the leaf path: `slice_to_use` from a
  top/full leaf lifetime must reproduce `le_g`'s leaf slicing exactly (Task 2 step
  4 / escalation) -- verify, do not assume.
- **OFF-path byte-identical at every task boundary** (`order_aware_recompute` OR
  `batch_spectator_indices` off => empty `batch_context` => `slice_to_use` inert).
