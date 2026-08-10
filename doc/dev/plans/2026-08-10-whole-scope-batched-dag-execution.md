# Whole-scope batched DAG execution -- implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a whole-scope batched DAG executor that walks the fused forest scope-by-scope (a value homed at a scope node is built once per block and reused by its subtree across all trees), replacing per-tree forest descent for batched schedules, so shared K-carrying intermediates are no longer rebuilt per tree.

**Architecture:** A scope-major SCHEDULER builds a scope tree from the existing DAG-scope placement (`compute_dag_boulevard`'s `ValueCell`s); a pure-realizer EXECUTOR walks that tree, reusing the existing `evaluate_impl` primitives (per-op eval, per-block scratch, accumulate/scatter, slice-on-use, lazy compute-at-home-on-first-hit). The executor makes NO placement decisions -- it honors the schedule exactly. It coexists with forest descent behind a flag until validated. See `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md`.

**Tech Stack:** C++20, SeQuant `core/eval`; Catch2 in `tests/unit/`.

## Global Constraints

- **Pure realizer.** The executor never makes a placement/spill decision; peak/recompute policy stays upstream in placement. Any "demote on pressure" in the executor is a bug.
- **Reuse primitives, replace the driver.** Do NOT rewrite per-op eval or the scratch/accumulate/scatter/slice-on-use mechanics in `evaluate_impl`; the executor is a new DRIVER over them.
- **Coexist, do not retire.** Forest descent (`evaluate(nodes,...)`, eval.hpp:1105) stays; whole-scope is selected by a flag. Retiring forest descent is out of scope.
- **Narrow placement is the target.** Consume the narrow canonical-chain scope tree the current machinery emits (one loop per index type, single order, set-based `home_modes`). General branching-tree placement (tree-LCA) is a later, separate project -- build the executor GENERAL but validate on the narrow input.
- **Correctness = small numerical noise, NOT byte-identical.** The schedule changes contraction order; validate whole-scope vs forest to a tight tolerance (energies/residual norms), not to the bit.
- **No en-dashes U+2013 / U+00A0.** `clang-format` touched files (re-stage after the hook). No AI-attribution trailers. Build `unit_tests-sequant -j6`; run TARGETED filters (an unrelated slow DP test hangs full runs). New test files go in `tests/unit/CMakeLists.txt`.

---

### Task 1: Scope-tree schedule -- data structure + builder

Build the scope tree and per-value home from the existing placement, with no execution yet.

**Files:**
- Create: `SeQuant/core/eval/scope_schedule.hpp`.
- Create test: `tests/unit/test_scope_schedule.cpp`; add to `tests/unit/CMakeLists.txt`.

**Interfaces:**
- Consumes: `sequant::eval::RichSchedule` / `ValueCell` (`peak_profile.hpp`) -- `home_modes`, `carried`, `value_id`, `first_use`/`last_use`, and the forest's `batched_here()` mode kinds.
- Produces:
  ```cpp
  namespace sequant::eval {
  struct ScopeNode {
    Index mode;                         // the batch loop mode (root: empty/sentinel)
    BatchModeType kind;                 // Contracted or External (root: n/a)
    container::svector<std::size_t> homed_values;  // value_ids homed AT this node
    container::svector<ScopeNode> children;         // deeper loops
  };
  struct ScopeSchedule {
    ScopeNode root;                     // top scope
    std::size_t num_values = 0;
  };
  // Build the narrow scope tree: one node per index type present, canonical
  // order; each value assigned to the node for its home_modes set (the deepest
  // level whose enclosing modes == home_modes). Root holds home_modes-empty
  // values.
  ScopeSchedule build_scope_schedule(RichSchedule const& rich,
                                     auto const& mode_order);
  }
  ```

- [ ] **Step 1: Write the failing test.** Two-term forest sharing one K-carrying intermediate, one K batch mode. Assert `build_scope_schedule` yields `root` with one child (the K loop, `kind==Contracted`), the shared intermediate's `value_id` in the K node's `homed_values`, and a K-independent value in `root.homed_values`.

```cpp
TEST_CASE("build_scope_schedule places values at their home scope", "[scope-schedule]") {
  // build a small RichSchedule via compute_dag_boulevard on a 2-term forest with
  // one Κ batch mode (reuse df_regime + a tiny hand-built forest, as the
  // peak_profile tests do)
  auto const sched = build_scope_schedule(rich, /*mode_order=*/{L"Κ"});
  REQUIRE(sched.root.children.size() == 1);
  CHECK(sched.root.children[0].mode.space().base_key() == L"Κ");
  CHECK(sched.root.children[0].kind == BatchModeType::Contracted);
  // the shared K-carrying intermediate homes in the K node; a K-independent
  // value homes at root
  CHECK(!sched.root.children[0].homed_values.empty());
}
```

- [ ] **Step 2: Run it to confirm it fails** (missing header/symbol).
- [ ] **Step 3: Implement `build_scope_schedule`** -- derive the mode set present, order canonically, build the chain of `ScopeNode`s, assign each `ValueCell` to the node matching its `home_modes` (root when empty). No execution.
- [ ] **Step 4: Run + confirm green.**
- [ ] **Step 5: Commit.**

---

### Task 2: Whole-scope executor -- top-scope-only walk (unbatched equivalence)

The executor skeleton: a scope-tree walk that, for a forest with NO batch loops (scope tree = root only), evaluates all root-homed values and sums -- provably equivalent to forest descent.

**Files:**
- Create: `SeQuant/core/eval/scope_executor.hpp`.
- Test: `tests/unit/test_scope_executor.cpp` (add to CMakeLists).

**Interfaces:**
- Consumes: `ScopeSchedule` (Task 1), the forest nodes, `leaf_evaluator`, `CacheManager`, layout -- mirroring `evaluate(nodes,...)`'s signature.
- Produces:
  ```cpp
  template <Trace EvalTrace = Trace::Default, ...>
  ResultPtr evaluate_whole_scope(Nodes const& forest, ScopeSchedule const& sched,
                                 auto const& layout, F const& leaf_evaluator,
                                 CacheManager<N,FHC>& cache);
  ```
  It walks `sched.root`; root-homed values are evaluated via `evaluate_impl` on the (unsliced) top scope and summed into the result. Reuses `evaluate_impl`; no new per-op code.

- [ ] **Step 1: Write the failing test.** An UNBATCHED forest (no `batched_here` stamps). Assert `evaluate_whole_scope(...)` returns a result equal to `evaluate(forest,...)` (forest descent) within `Approx` tolerance on the scalar/tensor norm.
- [ ] **Step 2: Run to confirm it fails** (missing symbol).
- [ ] **Step 3: Implement the root-only walk** -- for each root-homed value, `evaluate_impl` it and accumulate into the result (reusing the `add_inplace` + trace bookkeeping from `evaluate(nodes,...)`).
- [ ] **Step 4: Run + confirm green** (whole-scope == forest descent for unbatched).
- [ ] **Step 5: Commit.**

---

### Task 3: Single aux batch loop + lazy home-materialization -- the CSE win

Walk one batch loop (aux `Κ`) once for the whole forest: build each K-homed value once per block, share across trees; a value homed ABOVE (root) referenced mid-loop is materialized lazily at its home (reuse the existing compute-at-home-on-first-hit + slice-on-access), NOT rebuilt per block; accumulate contracted-K results on block exit.

**Files:**
- Modify: `SeQuant/core/eval/scope_executor.hpp`.
- Test: `tests/unit/test_scope_executor.cpp`.

**Interfaces:**
- Consumes: the lazy home-materialization already in `evaluate_impl` (a mid-loop reference to an above-homed value builds it full at its home and slices thereafter) -- confirm the seam and drive it from the walk.

- [ ] **Step 1: Write the failing test (equivalence).** The water-20 CSV doubles residual (reuse `df_regime(kWater20_pVDZF12)` + the aux-only batched forest from `test_eval_dryrun.cpp`). Assert `evaluate_whole_scope` and forest descent produce results agreeing to a tight tolerance (e.g. residual L2 norm within `1e-10` relative).
- [ ] **Step 2: Write the failing test (sharing).** Trace the build count of a shared K-carrying gC composite (e.g. `g{i₃,a₃<i₂,i₃>;K}`): under whole-scope it is built ONCE per K-block (`store == n_blocks`), NOT once per consumer group (`store == n_groups x n_blocks`). This is the regression the project targets.
- [ ] **Step 3: Run both to confirm they fail** (single-loop walk not implemented).
- [ ] **Step 4: Implement the single-loop walk.** `for each K-block b: bind context Κ=b; build each K-node-homed value for b (operands sliced-on-use from root-homed values, which materialize lazily at root on first hit); on block exit accumulate contracted-K results into the parent`. Free the K scratch after the loop.
- [ ] **Step 5: Run + confirm green** (equivalence + build-once sharing).
- [ ] **Step 6: Commit.**

---

### Task 4: Nested scopes (aux + occ)

Generalize the walk to a nested scope tree: recurse into child scope nodes, one batch loop per level, accumulate (contracted) / scatter (external) at each level's block exit.

**Files:**
- Modify: `SeQuant/core/eval/scope_executor.hpp`.
- Test: `tests/unit/test_scope_executor.cpp`.

**Interfaces:**
- Consumes: `ScopeNode::children` (the nested loops); the External-mode SCATTER primitive from `evaluate_impl` (block partials scattered into a pre-sized result), the Contracted-mode ACCUMULATE primitive.

- [ ] **Step 1: Write the failing test.** A forest with a nested (aux `Κ` outer, occ `i` inner) scope tree -- reuse the occ+aux batched forest construction from `test_eval_dryrun.cpp`. Assert nested `evaluate_whole_scope` vs forest descent agree to tolerance, and a value homed at the OUTER (K) level is built once per K-block and reused across all inner occ-blocks.
- [ ] **Step 2: Run to confirm it fails.**
- [ ] **Step 3: Implement the recursive walk** -- `walk(node, ctx)`: for each block of `node.mode`, build node-homed values, `for child in node.children: walk(child, ctx+block)`, then accumulate/scatter on exit. External modes scatter into a pre-sized result shared across the subtree; contracted modes accumulate.
- [ ] **Step 4: Run + confirm green.**
- [ ] **Step 5: Commit.**

---

### Task 5: Weighted use-count lifetimes -- remove `ensure_hoist_slot`

Replace the `ensure_hoist_slot` MAX-life hack (in BOTH eval.hpp's router-driven batched path and the whole-scope executor) with the correct WEIGHTED use count, so a value's cache life is deterministic and it is released as soon as its last in-block consumer is done (not when its home scope closes). This is the executor's core lifetime model and also fixes a long-standing misdesign in the forest-descent batched eval.

**Rationale (owner decision):** `ensure_hoist_slot(key)` (cache_manager.hpp:402) sets `life = MAX` because the emitted `effective_count` for a loop-invariant is 1 (it counts only the escaped-outer set), which as a life drains the entry after first use. The CORRECT life weights each consumer's access by the iteration counts of the nested loops between the value's home and the consumer:
```
life(V) = sum over consumers C of V:  product over loops L on the path (home(V), scope(C)]:  n_blocks(L)
```
With this life the use-counted cache holds V for exactly its in-block accesses, `reset()` at the home scope's block boundary rebuilds it for the next block, and V frees early. The correct store location is V's remat home (`node(home_meet)` + router).

**Files:**
- Investigate + modify: wherever `effective_count` / the escaped-outer set is computed and emitted for the batched eval (find it first -- likely the placement/schedule or eval.hpp's batched-annotation walk).
- Modify: `SeQuant/core/eval/eval.hpp` (the router hoist path ~1875-1930 -- drop `ensure_hoist_slot`, store at home with the weighted life), `SeQuant/core/eval/scope_executor.hpp` (the nested path -- drop `ensure_hoist_slot` + the fresh-per-level-cache workaround; store homed values at their home node's scratch with the weighted life).
- Modify: `SeQuant/core/eval/cache_manager.hpp` (remove `ensure_hoist_slot` once no callers remain).
- Test: `tests/unit/test_scope_executor.cpp`, and a `[eval]` regression run.

**Interfaces:**
- Produces: a `weighted_use_count(value, scope_tree, block_counts)` computation (Sum-over-consumers Product-over-intervening-loops), consumed as the cache life when a homed value is stored.

- [ ] **Step 1: Locate the current use-count/`effective_count` emission** for the batched eval (grep `effective_count`, `escaped`, the `ensure_hoist_slot` call sites) and document in the report where the life originates and why it is 1 for invariants.
- [ ] **Step 2: Write the failing tests.** (a) a unit test of `weighted_use_count` on a small nested scope tree (home at outer, consumer inside an n-block inner loop -> count == n, and a two-consumer case summing correctly); (b) the whole-scope build-once metric (Task 3/4: water-20 shared composite built == n_blocks) must still hold WITHOUT `ensure_hoist_slot`; (c) a `[eval]` numeric-equivalence pin (the forest-descent batched eval must stay BYTE-IDENTICAL in results -- a lifetime change must not alter numerics, only when entries free).
- [ ] **Step 3: Run to confirm (a) fails** (missing `weighted_use_count`) and note (b)/(c) as the invariants to preserve.
- [ ] **Step 4: Implement `weighted_use_count`; store homed values at their home with it; remove `ensure_hoist_slot` from eval.hpp and scope_executor.hpp; drop the fresh-per-level-cache workaround.**
- [ ] **Step 5: Run.** (a) green; (b) build-once still == n_blocks; (c) `[eval]` byte-identical numerics; and MEASURE the peak improvement from early release vs the old scope-close/MAX lifetime (report the number).
- [ ] **Step 6: Remove `ensure_hoist_slot` from cache_manager.hpp** (no callers remain) and confirm the build.
- [ ] **Step 7: Commit.**

---

### Task 6: Coexistence flag + cost-model peak-model selection

Select whole-scope vs forest descent by a flag at the driver entry, and make the dry-run cost model price the matching peak model (co-residency oracle for whole-scope).

**Files:**
- Modify: the driver dispatch (a `BatchPolicy` / eval-options flag), `SeQuant/core/eval/backends/dryrun/cost_profile.hpp` (peak-model selection).
- Test: `tests/unit/test_scope_executor.cpp` + a dry-run cost test.

**Interfaces:**
- Produces: a flag (e.g. `BatchPolicy::whole_scope_execution`, default false) routing `evaluate` to `evaluate_whole_scope`; `cost_profile` reads it to select the co-residency peak model (`peak_profile_sweep` over `home_modes`) rather than the batched-scratch replay high-watermark.

- [ ] **Step 1: Write the failing test.** With the flag ON, the driver routes to whole-scope (result matches Task 3/4); with it OFF, forest descent (byte-identical to today). And `cost_profile` under the flag reports the co-residency peak (which, per the spec, now MATCHES the realized whole-scope peak) rather than the forest replay peak.
- [ ] **Step 2: Run to confirm it fails.**
- [ ] **Step 3: Implement the flag + dispatch + cost-model peak selection.**
- [ ] **Step 4: Run + confirm green;** confirm the flag-OFF path is unchanged (regression on `[eval]`).
- [ ] **Step 5: Commit.**

---

### Task 7: Whole-branch validation (water-20 + C60)

End-to-end witnesses that whole-scope descent recovers sharing and that the cost model predicts it.

**Files:**
- Test: `[.]` witnesses in `tests/unit/test_scope_executor.cpp` (hidden; run by name).

- [ ] **Step 1: water-20 witness.** whole-scope vs forest agree to tolerance on the residual; the targeted shared composites build once (not per group); print the realized peak and confirm the co-residency oracle predicts it.
- [ ] **Step 2: C60 witness (`[.]`).** whole-scope runs on the C60 aux+occ residual; confirm correctness-to-tolerance and record the realized peak (feasibility of a given placement is the optimizer's concern, not the executor's -- just measure).
- [ ] **Step 3: Run + document** the recompute elimination vs forest descent.
- [ ] **Step 4: Commit** the witnesses (`[.]`, excluded from CI).

---

## Self-review notes

- **Spec coverage:** scheduler (spec "Schedule (new)") -> Task 1; pure-realizer executor + driver-swap (spec "load-bearing insight", "pure realizer") -> Tasks 2-4; lazy home-materialization (spec Open Questions) -> Task 3; nested general tree (spec "scope tree") -> Task 4; weighted use-count lifetimes / remove ensure_hoist_slot (owner decision during Task 4 review) -> Task 5; coexistence + cost-model binding (spec "paradox resolved", Open Questions) -> Task 6; validation-to-tolerance + recompute-elimination (spec Validation) -> Task 7.
- **Scope guardrails:** general branching-tree PLACEMENT and the intraloop-CSE-pricing OPTIMIZER are explicitly OUT (spec Scope boundary); no task touches them. C60 top-scope feasibility is measured, not enforced (executor is a pure realizer).
- **Greenfield honesty:** Tasks 3-4 depend on how `evaluate_impl`'s scratch/accumulate/scatter/lazy-home primitives compose under a whole-forest driver (a spec Open Question). Each task's FIRST deliverable is an equivalence test against forest descent, so the composition is validated empirically per increment rather than assumed. If a primitive does not compose (e.g. a scatter target cannot be shared across trees), that surfaces as a failing equivalence test at the smallest increment, not late.
- **Task 5 (use-count cleanup) invariant:** it changes only WHEN cache entries free, never numerics -- the `[eval]` byte-identical pin (Step 2c/5c) is the guard; it also removes `ensure_hoist_slot` from eval.hpp's existing forest-descent batched path, so `[eval]` regression coverage is load-bearing.
- **No silent regression:** Task 6 Step 4 pins the flag-OFF path byte-identical (forest descent unchanged) so whole-scope is purely additive until adopted.
- **Watch item:** the "build once" sharing assertions (Tasks 3/5/7) are the real success metric -- equivalence-to-tolerance alone would pass even if nothing shared. Both must hold.
