# Phase 4: shared BatchPolicy - one batchability source for optimizer + eval

_Status: design draft. Date: 2026-06-19._

Completes the CostModel work (Phases 1-3) by removing the duplicated batchability
policy: today mpqc defines the same logical triple twice - once for the optimizer
(`OptimizeOptions` in `sequant_engine.cpp`) and once for the runtime batched
evaluator (`make_batched_custom_evaluator` in `cck.ipp`), even with the volatile
predicate in different types (Tensor- vs EvalNode-based). Phase 4 bundles the
triple into one `BatchPolicy` object both sides read.

## 1. Context and goal

The batchability policy is three things: which index spaces slice
(`is_batchable_index`), the per-index slice size (`batch_target_size`), and which
leaves are volatile/persistent (`is_volatile_leaf`). Phase 3 carried the first
two as loose `OptimizeOptions` fields (a scalar `batch_target_size`) and the
runtime evaluator took its own copies. **Goal:** define the policy once, in one
`BatchPolicy` object, consumed by both the optimizer and a thin eval-layer
adapter, so the two cannot drift; and generalize `batch_target_size` to be
per-index (keyed by `IndexSpace`).

**Non-goals:** changing the batched-eval algorithm (`make_batched_custom_evaluator`
stays; only its `target_batch_size` parameter generalizes to a function); the
old-standalone-DP removal (separate fast-follow); any objective/recurrence change.

**Hard constraint:** behavior-preserving. SeQuant: the full `[optimize]` suite
stays green. mpqc: a CSV-CCk reference run produces the same energy.

## 2. The `BatchPolicy` object

A small neutral header `SeQuant/core/batch_policy.hpp` (forward-declares `Index`,
`Tensor`; depends only on `<functional>`/`<cstddef>`), so both the optimize and
eval layers include it without coupling to each other:

```cpp
namespace sequant {
class Index;
class Tensor;

struct BatchPolicy {
  /// Which index spaces the runtime evaluator slices over (= eval accept_aux).
  std::function<bool(Index const&)> is_batchable_index = {};
  /// Per-index slice size (keyed by IndexSpace internally). Consulted only for
  /// batchable indices; slicing index Ki uses min(extent(Ki), batch_target_size(Ki)).
  std::function<std::size_t(Index const&)> batch_target_size = {};
  /// Canonical (Tensor-based) volatile-leaf predicate. The eval adapter lifts it
  /// to an EvalNode predicate (section 3).
  std::function<bool(Tensor const&)> is_volatile_leaf = {};
};
}  // namespace sequant
```

`is_batchable_index` and `batch_target_size` are kept parallel (a predicate plus
a size). `is_batchable_index` is used on its own in the optimizer's
internal/frontier gating where the size is irrelevant; `batch_target_size` is
read only when slicing. (They could collapse to one function - `batch_target_size>0`
implying batchable - but the explicit predicate reads clearer at the gating
sites; keep both.)

## 3. The eval-layer adapter

A free function in the eval layer (alongside `make_batched_custom_evaluator`),
templated on the backend leaf-evaluator `F`; this confines all backend-genericity
to the eval layer and keeps `cost_model.hpp`/`options.hpp` eval-agnostic:

```cpp
template <class F, class ScopeGuardFactory = /* default factory */>
[[nodiscard]] auto make_evaluator(BatchPolicy const& policy, F yielder,
                                  ScopeGuardFactory make_scope_guard = {}) {
  // Lift the canonical Tensor predicate to the EvalNode predicate the batched
  // evaluator needs. Reproduces cck.ipp's current is_volatile exactly when
  // is_volatile_leaf is the label test.
  auto is_volatile_node = [p = policy.is_volatile_leaf](auto const& n) {
    return n.leaf() && n->is_tensor() && p(n->as_tensor());
  };
  return make_batched_custom_evaluator(
      std::move(yielder), policy.batch_target_size, policy.is_batchable_index,
      std::move(make_scope_guard), std::move(is_volatile_node));
}
```

`make_batched_custom_evaluator` is updated so its `target_batch_size` parameter is
the per-index function `std::function<std::size_t(Index const&)>`: where it
currently calls `mode_batches(leaf_index, target_batch_size)` it calls
`mode_batches(leaf_index, target_batch_size(K))` for the batch axis `K` it
selects (the only behavioral substitution; everything else - `batch_axis`,
persistence gate, replay group, layering - is unchanged).

## 4. Optimizer side

`OptimizeOptions` embeds one `BatchPolicy batch_policy` in place of the three
loose Phase-3 fields. The model `build_context`s read from it:
- `PeakBatchedModel`: `batch_policy.is_batchable_index`, `batch_policy.batch_target_size`.
- `AdditiveModel`: `batch_policy.is_volatile_leaf` (for volatile weighting; `volatile_weight` stays a separate `OptimizeOptions` field).
- `PeakModel`: none.

The per-index `batch_target_size` threads through the batched optimizer path -
`batched_extent`/`sliced_footprints` apply `min(extent(ix), batch_target_size(ix))`;
the `batch` member of `peak_dp_batched`/`peak_cost_batched`/
`single_term_opt_peak_batched_impl`/`PeakBatchedModel` changes from `std::size_t`
to the function. Behavior is preserved when the function is the old constant.

## 5. mpqc wiring

Build one `BatchPolicy` in the CSV-CCk setup and feed both consumers:
- `is_batchable_index = [aux](Index const& ix){ return ix.space() == aux; }` (DF aux space).
- `batch_target_size = [n=csv_batch_aux_target_size_](Index const&){ return n; }` (constant for now; per-space variation is now expressible).
- `is_volatile_leaf = [](Tensor const& t){ return t.label() == L"t"; }`.

Threading: the policy is constructed once (in `CCk`/the CSV setup) and passed to
- the optimizer, via the `EvalContext` `SeQuantEngine::optimize` already takes (or a `CCk` member SeQuantEngine reads), setting `opts.batch_policy`;
- the eval cache: `cache.set_custom_evaluator(sequant::make_evaluator(policy, yielder, make_scope_guard))`.

Delete the now-duplicated `accept_aux` and `is_volatile` lambdas in `cck.ipp`
(and the separate volatile/batchable wiring in `sequant_engine.cpp`) - both now
derive from the single policy. The exact threading point (EvalContext vs member)
is a plan-level detail; the requirement is one construction site.

## 6. Scope: two stages

One spec, built in order (the SeQuant stage is self-contained and unit-testable;
the mpqc stage is the consumer rewire + the first mpqc validation in this work):

- **Stage A (SeQuant):** `batch_policy.hpp`, `make_evaluator`, the per-index
  `batch_target_size` ripple through the batched path and
  `make_batched_custom_evaluator`, and the `OptimizeOptions` rewire to embed
  `BatchPolicy`.
- **Stage B (mpqc):** construct one `BatchPolicy`, feed optimizer + eval cache,
  delete the duplicated definitions.

## 7. Testing

- **SeQuant (Stage A):**
  - The full `[optimize]` suite stays green; the batched equivalence tests now
    pass the policy/function `batch_target_size` form (a constant lambda),
    reproducing the prior scalar results - the behavior-preservation proof for
    the scalar->function change.
  - A new eval-layer test: `make_evaluator(policy, yielder)` produces an
    evaluator that yields the same result as a hand-built
    `make_batched_custom_evaluator(yielder, const_size, accept, guard, is_volatile_node)`
    on a small batched network (proves the adapter + the Tensor->EvalNode lift).
  - A test that a non-constant `batch_target_size` (different sizes per space) is
    honored per index (the new capability).
- **mpqc (Stage B):** a CSV-CCk reference calculation (e.g. an `he10` CSV-CCk
  validation input) produces the same energy as before the change - the
  bundling is behavior-preserving end to end. Requires an mpqc build + CSV-CCk
  run.

## 8. Risks and open questions

- **`mode_batches` per-index call site.** Confirm `make_batched_custom_evaluator`
  has exactly one `target_batch_size` use (the `mode_batches` probe on the
  selected axis `K`) so the scalar->function change is a one-line substitution;
  if the size is used elsewhere (e.g. the replay-group `mode_batches` equality
  probe at the candidate scan), each call must pass `batch_target_size(its own K)`
  consistently so the group-compatibility test still compares like with like.
- **`BatchPolicy` placement.** A dedicated `core/batch_policy.hpp` avoids an
  optimize<->eval dependency; confirm both `options.hpp` and the eval header can
  include it without cycles.
- **mpqc threading point.** Whether the policy rides on the existing
  `EvalContext` or a `CCk` member - pick the one that reaches both
  `SeQuantEngine::make_optimize_options` and the `cck.ipp` eval-cache setup with
  one construction.
- **Volatile lift generality.** The lift assumes a volatile leaf is a leaf
  tensor whose label the predicate tests; this matches all current callers.
  Document that `is_volatile_leaf` must be a leaf-tensor predicate.
- **Old-standalone-DP removal** remains a separate fast-follow (out of scope).
