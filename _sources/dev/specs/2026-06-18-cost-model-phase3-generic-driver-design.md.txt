# Phase 3: CostModel generic single-term-DP driver (cost side)

_Status: design draft. Date: 2026-06-18._

Refines the cost side of the CostModel vision in
`2026-06-18-cost-model-batch-aware-design.md` (sections 3-4) into a concrete,
implementable, behavior-preserving refactor. The evaluator face
(model-is-evaluator) and mpqc wiring remain Phase 4.

## 1. Context and goal

After Phases 1-2, SeQuant has four single-term-optimization objectives -
`DenseFLOPs`, `DenseSize`, `DensePeakSize`, `DensePeakSizeBatched` - dispatched by
an `ObjectiveFunction` enum through an `if constexpr` chain in detail
`single_term_opt<Metric>` and a runtime `if`-chain in `opt_pure_product`. They are
two DP families: the additive DP (`single_term_opt_impl`, shared by FLOPs/Size via
a per-contraction `cost_fn`), and the peak family (`DensePeakSize`,
`DensePeakSizeBatched`, each its own structurally-different DP). Each DP
re-implements the same subset-lattice + bipartition orchestration.

**Goal:** factor that shared orchestration into one **generic, compile-time DP
driver** `run_single_term_opt<Model>`, and move each objective's recurrence and
reconstruction into a **`CostModel` type** (an associated `State` + `Context` and
a small set of methods). This:
- removes the duplicated lattice/bipartition driver code;
- encapsulates each objective's logic in one cohesive type;
- gives a clean compile-time extension point for custom objectives;
- positions the model to grow an evaluator face in Phase 4.

**Non-goals (Phase 4+):** the model-is-evaluator (`evaluate<Backend>`), replacing
the eval cache's `make_batched_custom_evaluator` wiring, mpqc changes, grouping
`OptimizeOptions` fields into per-model sub-structs, runtime-pluggable
(type-erased) custom models.

**Hard constraint:** behavior-preserving. Every built-in model driven through the
generic driver must reproduce its current DP's exact output. The existing
`[optimize]` tests and the Phase-1/2 oracles are the regression net.

## 2. The generic driver and the `CostModel` concept

The driver owns only what every single-term DP shares - the subset lattice
(subsets in increasing order) and the bipartition enumeration:

```cpp
template <class Model>          // Model satisfies the CostModel concept
EvalSequence run_single_term_opt(Model const& m, TensorNetwork const& network,
                                meta::range_of<Index> auto const& tidxs) {
  auto const nt = network.tensors().size();
  if (nt == 1) return EvalSequence{0};
  if (nt == 2) return EvalSequence{0, 1, -1};
  typename Model::Context ctx = m.build_context(network, tidxs);   // precompute
  container::vector<typename Model::State> st(size_t{1} << nt);
  for (size_t n = 1; n < st.size(); ++n) {
    if (std::popcount(n) == 1) { st[n] = m.leaf(ctx, n); continue; }
    typename Model::State acc = m.init(ctx, n);
    for (auto&& [lp, rp] : bits::bipartitions(n))
      if (lp && rp)
        m.relax(ctx, n, lp, rp, st[lp], st[rp], acc);   // fold one bipartition
    st[n] = std::move(acc);
    m.finalize(ctx, n, st);                              // post-bipartition hook
  }
  return m.reconstruct(ctx, st);
}
```

The **`CostModel` concept** (a built-in or custom objective satisfies it):

- associated `State` - the per-subset DP cell. The driver never inspects it.
  Scalar-ish for additive/peak; a `[B]`-vector for batched.
- associated `Context` - the model's precomputed tables AND mutable DP scratch
  (built once by `build_context`).
- `Context build_context(network, tidxs) const` - precompute open-index sets /
  footprint tables / subnet metadata / params, and validate (e.g. reject
  `subnet_cse` where unsupported).
- `State leaf(ctx, n) const` - singleton cell.
- `State init(ctx, n) const` - empty accumulator (worst).
- `void relax(ctx&, n, lp, rp, lp_state, rp_state, acc&) const` - fold one
  bipartition into `acc`, computing the candidate via the model's recurrence
  (order / `A'` / `[B]` all internal) and keeping the best, recording
  back-pointers.
- `void finalize(ctx&, n, states&) const` - per-subset post-processing (CSE
  bookkeeping); identity for most models.
- `EvalSequence reconstruct(ctx, states) const` - build the sequence from the
  recorded back-pointers.

`Context&` is passed mutable because `finalize` (and additive `relax` under CSE)
update DP scratch (`unique_meta_costs`). `value(State)` - the scalar objective
used for comparisons - is a model-internal helper used inside `relax`/`reconstruct`.

The `nt == 1`/`nt == 2` early returns match the current impls and keep trivial
products off the model path.

## 3. The four built-in models

`DenseFLOPs`/`DenseSize` collapse into one `AdditiveModel` parameterized by the
cost primitive. Each model's methods reproduce its current DP exactly.

### AdditiveModel (cost_fn = flops_counter | memsize_counter)
- `State`: `{double ops; size_t lp, rp; container::vector<size_t> subnets;}`.
- `Context`: open-index sets (`init_results`); `cost_fn`; `footprint_fn` +
  `footprint_weight`; `volatile_mask` + `volatile_weight`; subnet metadata
  (`meta_ids`, `unique_meta_costs`) when CSE.
- `relax`: `new = w*cost_fn(open[lp],open[rp],open[n]) + fp + (CSE ? sum_meta :
  lp.ops + rp.ops)`, with `w = (volatile_mask & n) ? volatile_weight : 1` and
  `fp = footprint_weight ? footprint_weight*footprint_fn(open[n]) : 0`; keep the
  min, record `lp,rp` (and combined `subnets` under CSE).
- `finalize`: CSE bookkeeping (update `unique_meta_costs[meta_ids[n]]`, insert the
  subnet id) - today's single_term_opt_impl lines ~356-371; identity when no CSE.
- `reconstruct`: the index-order tie-break emission (today's ~374-379).

### PeakModel
- `State`: `{double peak; size_t lp, rp; bool lp_first;}`.
- `Context`: `S[n]`, `L[n]` (`subset_footprints` + leaf-sums). `build_context`
  asserts `!subnet_cse`.
- `leaf`: `peak = S[n]`. `init`: `peak = max`.
- `relax`: pebble - `cand = min` over order of
  `max(L[other]+child.peak, S[child]+other.peak, S[lp]+S[rp]+S[n])`; keep min,
  record `lp_first`.
- `finalize`: identity. `reconstruct`: `lp_first`-ordered postorder
  (today's `single_term_opt_peak_impl`).

### PeakBatchedModel
- `State`: `container::vector<BatchedRes>` - the **`[B]`-vector**
  (`{peak,lp,rp,lp_first,aprime}` per `B`, `|B-space| = 2^m`).
- `Context`: `tables[B]` (`sliced_footprints`), per-`B` `L`, `open_aux[n]`
  (`subset_open_aux`), `volatile_mask`, `m`. `build_context` asserts
  `!subnet_cse`.
- `leaf`: per `B`, `peak = tables[B & open_aux[n]][n]`. `init`: per `B`, `max`.
- `relax`: per `B`, over `A' subset of A(lp,rp)` (persistence-gated), `C = B|A'`,
  pebble at `C`; update `acc[B]`, record `aprime`.
- `finalize`: identity. `reconstruct`: back-pointer walk from `(root, B=0)`
  (today's `single_term_opt_peak_batched_impl`).

The `[B]` dimension lives entirely inside `PeakBatchedModel` (its `State` is the
vector; its `relax` loops `B` and `A'`). The driver stores `State` per subset and
folds bipartitions, oblivious to it - this is what lets one generic driver serve
all three structures. `build_context` reuses the existing helpers
(`init_results`, `subset_footprints`, `sliced_footprints`, `subset_open_aux`,
`build_subnet_metadata`) verbatim, so no numbers move.

## 4. Dispatch transition

Models are compile-time types selected at runtime by the enum, so the dispatch
boundary stays where it is:

- `OptimizeOptions` keeps `objective_function` and its existing fields
  (`idx_to_extent`, `is_volatile_leaf`, `volatile_weight`, `footprint_weight`,
  `is_batchable_index`, `batch_target_size`, `CSE`). The refactor does NOT remove
  the enum (compile-time models need the runtime->type map). Fields stay flat,
  documented by consuming model.
- `opt_pure_product` becomes: `objective_function -> build Model from the opts
  fields it needs -> run_single_term_opt(model, network, tidxs)`. The detail
  `single_term_opt<Metric>` if-constexpr surface collapses into "construct
  `Model<Metric>`, call the driver."
- The public `single_term_opt(Product, ...)` wrapper (EvalSequence -> ExprPtr)
  stays; it builds the model and calls the driver.
- Each model's `build_context` reads only its own fields (Additive:
  cost-kind+weights+CSE; Peak: none extra; Batched:
  `is_batchable_index`+`batch_target_size`).
- **Custom objectives (compile-time):** `run_single_term_opt<Model>` is public, so
  a user-defined model type satisfying the concept can be driven directly. The
  four built-ins additionally route through the enum. Runtime-pluggable custom
  (an arbitrary model object via `OptimizeOptions`) is deferred (type erasure).

## 5. Scope boundaries

- **No evaluator face.** Models have only the cost/optimize side; the eval cache
  keeps its current separate `make_batched_custom_evaluator` wiring. (Phase 4.)
- **No mpqc changes.**
- **Behavior-preserving** by construction - see section 6.

## 6. Testing

- **Regression net (primary):** the existing `[optimize]` SECTIONs cover all four
  objectives plus the Phase-1/2 oracles and reconstruction checks. After the
  refactor they must all pass **unchanged** - that is the behavior-preservation
  proof. No test numbers change.
- **Concept conformance:** a compile-time check that each of `AdditiveModel`
  (both cost primitives), `PeakModel`, `PeakBatchedModel` satisfies the
  `CostModel` concept.
- **Extension point:** one tiny custom model (e.g. a trivial additive variant or
  a "size with a per-index weight" metric) driven directly via
  `run_single_term_opt<Custom>`, asserting it runs and returns a valid sequence -
  proving the public generic entry point works for a non-built-in.

## 7. File layout

Likely shape (final split decided in the plan):
- `core/optimize/cost_model.hpp` (new) - the `CostModel` concept, the four
  built-in model types, and `run_single_term_opt`.
- `core/optimize/single_term.hpp` - the existing DP helpers
  (`subset_footprints`, `peak_dp`, `sliced_footprints`, `subset_open_aux`,
  `build_subnet_metadata`, the cost/footprint counters) stay and are reused by
  the models' `build_context`/`relax`; the per-objective `*_impl` driver
  functions and the `single_term_opt<Metric>` if-constexpr are replaced by the
  generic driver + model dispatch.
- `core/optimize/optimize.cpp` - `opt_pure_product` builds the model and calls
  the driver.
- `tests/unit/test_optimize.cpp` - existing tests unchanged; add the
  concept-conformance and custom-model checks.

## 8. Risks and open questions

- **Rewrite risk / behavior preservation.** Moving three DPs behind one driver is
  the substantial-rewrite option (acknowledged when chosen). Mitigation: the
  models are mechanical extractions of existing loop bodies; the full existing
  test+oracle suite must stay green with no number changes; do it model-by-model
  (additive first, then peak, then batched) so a regression localizes.
- **`Context` mutability vs the mostly-const methods.** Only additive-CSE and
  `finalize` mutate; keep the mutation surface explicit (the driver passes
  `Context&`). Confirm the concept doesn't force needless non-const elsewhere.
- **`relax` signature (decided).** `relax` mutates `acc` in place (void return),
  as the driver sketch shows - this avoids copying the batched `[B]`-vector (`2^m`
  entries) on every bipartition. Settled, not open.
- **CSE generality.** Confirm the `finalize` hook fully reproduces the inlined CSE
  accounting (the `unique_meta_costs` update currently sits after the bipartition
  loop and also feeds later subsets); verify ordering is identical.
- **Whether the enum should eventually go.** It stays for Phase 3; if
  runtime-pluggable custom models arrive later, a type-erased path would
  supersede it. Out of scope now.
