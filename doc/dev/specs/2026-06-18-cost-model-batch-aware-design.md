# Batch-aware peak-memory cost model for SeQuant single-term optimization

_Status: design draft. Date: 2026-06-18._

## 1. Problem & motivation

SeQuant's single-term optimizer (STO, `single_term_opt_impl` in
`core/optimize/single_term.hpp`) chooses a contraction tree for one term by
minimizing a cost: either `DenseFLOPs` (flops) or `DenseSize` (summed dense
intermediate storage). For PNO-CCSD over DF integrals this objective does not
control peak memory, which is the binding constraint at cluster scale.

The empirical investigation (water 14/18/20-mer, c6h14; see
`project_pno_ccsd_water_trace_bottleneck`) established:

- **The binding resource is peak memory, not flops or total traffic.** water_20
  OOMs; water_14 fits.
- **The "giant"** `g(μ̃,μ̃,Κ)·C(i,i,μ̃;a<i,j>) → I(i,i,μ̃,Κ;a<i,j>)` (a
  half-transformed DF integral carrying a free aux index Κ and a free PAO index
  μ̃) is the dominant single intermediate (247 GB nominal at water_14). It is
  **batchable**: slicing over Κ bounds its footprint to `full/n_batches`.
- **The `n_replay` / `volatile_weight` knob is too blunt.** Cranking it slices
  the giant but defers t1 folding, which manufactures large *non-batchable*
  OSV-mode intermediates (79.7 GB, 108 GB at water_20). It trades one
  bottleneck for another.
- **`DenseSize` is worse, not better** (water_14: peak 290 GB vs 200 GB for
  flops+vw=100, and 96 s/iter). It is blind to batchability: it minimizes
  *gross* full-dense size, so it keeps the giant but in a factorization whose
  Κ-frontier is volatile — the runtime then materializes it whole.
- **Batchability is structural**, not per-index. The runtime
  (`make_batched_custom_evaluator`, `core/eval/eval.hpp`) slices a subtree only
  at a Κ-contraction whose subtree is persistent (the `subtree_any(is_volatile)`
  gate, eval.hpp:1128). A stateless "shrink Κ's extent everywhere" discount in
  the cost model is therefore wrong — it discounts intermediates the runtime
  would not batch.

Two further facts shape the design:

- **Gross footprint is useless; peak is what matters.** `footprint_weight`
  (options.hpp) adds `Σ(result size)` over intermediates — a total, not a peak.
- **The batchability control already exists, declaratively, at eval time.**
  `cck.ipp` (≈1166-1189) hands `make_batched_custom_evaluator` exactly: the
  batchable-space predicate `accept_aux` (`ix.space() == dfbs`), the slice size
  `batch:aux_target_size`, and the persistence predicate `is_volatile` (leaf
  labeled `t`). The factorization stage should consume the *same* triple so the
  cost model and the runtime agree by construction.

## 2. Goals and non-goals

**Goals.**

1. Make **peak memory** a first-class STO objective.
2. Make the peak objective **batchability-aware**: an intermediate carrying a
   free batchable index, under a *persistent* batch-index contraction, is costed
   at its sliced footprint, not its full one.
3. Make the objective **user-customizable** via injected behavior (lambdas /
   policy types), not a closed enum.
4. **Unify** the factorization cost policy and the runtime evaluation policy in
   one object holding the shared state, so they cannot drift.

**Non-goals.**

- **Cross-term CSE coupling.** STO optimizes one term at a time; cross-term
  common-subexpression sharing is a purely runtime concern and is *already*
  outside STO's scope. Intra-term CSE (`subnet_cse`) stays as is. We do not
  model how an intermediate's real peak changes when it is shared across terms.
- Changing the runtime batching mechanism itself (the `make_batched_custom_evaluator`
  algorithm). We change *who supplies its policy* and add a matching cost model.
- Sparse/CLR-aware costing. Extents remain dense (`approximate_size`).

## 3. Core abstraction: `CostModel`

One object owns the policy for both factorization and evaluation. The base holds
only the **problem description** that every model may consult; concrete models
add their own policy knobs.

```cpp
// backend-agnostic base — pure problem description
struct CostModelBase {
  index_to_extent_t       idx_to_extent;   // Index -> dense extent  (universal)
  std::function<bool(Tensor const&)> is_volatile;  // leaf changes between replays
};
```

A concrete cost model is a **policy type** (not a virtual base) satisfying a
concept with two faces:

**(A) STO cost algebra** — the incremental, DP-composable form of "score a
contraction tree". `State` is a model-defined associated type.

```cpp
State  leaf(idxs);                         // cost-state of a singleton subset
State  combine(State lp, State rp,         // cost-state of one binary contraction
               lhs_idxs, rhs_idxs, result_idxs);
double value(State);                       // scalar the DP minimizes
```

`cost(tree)` is the fold of `combine` over the tree; the DP never materializes a
whole tree, it folds bottom-up over subsets.

**(B) Evaluator** — the model *is* the evaluator (no separate object, no copied
state). Backend coupling is confined to one template member; eval-time
collaborators are arguments, not members:

```cpp
template <class Backend>
EvalResult evaluate(node, yielder, world, /* scope guard etc. */) const;
```

So mpqc constructs one `CostModel`, hands it to `optimize(expr, model)` and to
the eval cache; both faces read the same `is_batchable` / `batch_target_size` /
`is_volatile` members.

### Why a templated policy, not a virtual hierarchy

`State` differs per model (scalar for flops, a peak scalar for peak, a
`(peak_full, peak_slice)` pair for batched, arbitrary for custom). Type-erasing
`State` behind a virtual interface is awkward; instead the DP is **templated on
the CostModel concept**, matching today's style (`single_term_opt` is already
templated on the cost function and `ObjectiveFunction`). Built-ins are concrete
policy types; custom injection is a user policy type or a provided
`LambdaCostModel` adaptor.

## 4. Built-in models

| model | `State` | `combine` | eval |
|---|---|---|---|
| `DenseFLOPs` | scalar flops | `lp + rp + flops(contraction)` (order-independent) | plain |
| `DensePeakSize` | subtree peak | pebbling `max(...)` over both child orders | plain |
| `DensePeakSizeBatched` | `(peak_full, peak_slice)` | pebbling in both regimes + persistence-gated frontier swap | batched |

Policy knobs live on the derived models, not the base:

- `DenseFLOPs`: `volatile_flop_weight` (replay multiplier; today's
  `volatile_weight`).
- `DensePeakSizeBatched`: `is_batchable` (= eval's `accept_aux`),
  `batch_target_size`.

`DenseFLOPs` reproduces today's default behavior and is the migration target for
existing inputs.

## 5. Peak model

### 5.1 Peak is order-dependent

For a single binary contraction there is no order freedom. In a *tree*, a
computed subtree's result occupies memory while its sibling is computed, so peak
depends on the evaluation schedule. For `Z = L · R` with result sizes `S` and
subtree peaks `P`:

```
evaluate L first: peak = max( P_L, S_L + P_R, S_L + S_R + S_Z )
evaluate R first: peak = max( P_R, S_R + P_L, S_L + S_R + S_Z )
```

The middle terms differ (`S_L + P_R` vs `S_R + P_L`): hold the cheap-to-store
side, compute the peak-heavy side last. Flops and gross-sum are
order-*independent* (same intermediates, each counted once), which is why the
current DP fixes the eval order by index canonicalization (single_term.hpp
374-379) and never optimizes it. A peak objective makes evaluation order a real
lever; `combine` takes the `min` over the two child orders. (This sketch omits
the bystander-input term, added in 5.2.)

### 5.2 Recurrence and optimal substructure

**Memory model (all-co-resident).** Input tensors are materialized (DistArrays)
and resident; a tensor is live from its production (a leaf: from first use) until
its single consumption. The peak is the max over the schedule of the sum of all
live tensor sizes plus the result being formed. This is the realistic tensor
peak (the `opt_einsum`/`cotengra` "peak size"), not the register-allocation
(Sethi-Ullman) peak, which would omit not-yet-used input leaves.

Two build-independent tables per subset: `S[n]` = result footprint (product of
extents of `n`'s open indices) and `L[n]` = sum of leaf (singleton) sizes in
`n`. For subset `n = lp ⊕ rp`:

```
peak[n] = min over (bipartition, order) of
  lp-first: max( L[rp] + peak[lp], S[lp] + peak[rp], S[lp] + S[rp] + S[n] )
  rp-first: max( L[lp] + peak[rp], S[rp] + peak[lp], S[lp] + S[rp] + S[n] )
```

The `L[other]` term is the bystander cost: while one child evaluates, the other
child's inputs sit resident (constant throughout, since that subtree is
untouched - which is why the additive external context composes and the per-
subset `peak[n]` is well defined). `peak[singleton] = S[singleton]`.

`S[n]` and `L[n]` depend only on which tensors are in `n`, **not on how `n` was
built**, and the recurrence is monotone in the children's peaks, so minimizing
each child's peak minimizes the parent's. Hence **optimal substructure holds**:
a scalar minimization with per-subset back-pointers; no Pareto frontier needed.
(This assumes each child subtree is evaluated contiguously; the
brute-force-oracle regression test, which enumerates all schedules freely,
guards that assumption.)

The objective is `peak[root]` (the whole term; the final result is materialized,
so the root is not sliced).

## 6. Batchability in the peak model

### 6.1 Footprint as a monomial in batch extents

A tensor's footprint is `(∏ non-batch extents) · ∏ᵢ Kᵢ^{dᵢ}`, where `dᵢ` is the
number of times batchable index *i* appears in its free indices. Examples:

- giant `I(i,i,μ̃,Κ;a)` → `i²·μ̃·a · K¹` — degree 1 in K.
- μ̃⁴ `I(μ̃,μ̃,μ̃,μ̃)` → `μ̃⁴ · K⁰` — degree 0; nothing to slice.

Evaluating at `Kᵢ = full` gives the unbatched footprint; at `Kᵢ = batch_target_size`
gives the sliced footprint. This is the mechanism that makes "plug slice" a
one-line evaluation rather than a re-derivation, and it is what correctly keeps
μ̃⁴ at full size (it carries no free K) while discounting the giant.

### 6.2 Two-mode state and the persistence gate

`DensePeakSizeBatched::State = (peak_full, peak_slice)`:

- `peak_slice[n]` — the peak to evaluate subtree `n` **sliced** over the batch
  indices: the same recurrence with every `S` evaluated at `K = slice`.
- `peak_full[n]` — the peak to evaluate subtree `n` **whole**, except that a
  child which is a *batchable frontier* may contribute its `peak_slice` instead
  of its `peak_full`.

A node `n` is a **batchable frontier** when, in forming it, a batch index
becomes **internal** (it is in `lhs ∪ rhs` but not in `result` — i.e. contracted
here, so `n`'s result is K-free) **and** the subtree is **persistent**
(`is_volatile` false for every leaf below `n`; trackable as a per-subset
volatile bitmask exactly like today's `volatile_mask`). At such a node, the
runtime may slice the subtree, so `combine` for `peak_full[n]` may take
`peak_slice[child]` for that child. Both gates are local to the subset — no
ancestor lookups, dissolving the "DP knows below, not above" problem: we detect
the frontier *where the batch index is internalized*, looking down at the
children the DP has already solved.

This matches the runtime gate by construction: the cost model treats a subtree
as batchable iff a batch index is internal to it and it is persistent — the same
condition `make_batched_custom_evaluator` uses.

### 6.3 Worked check

- **giant under a persistent Κ-contraction:** frontier persistent → `peak_full`
  may take the giant subtree's `peak_slice`, which uses `i²·μ̃·a·K_slice` → small.
  The model keeps the giant and prices it sliced.
- **giant under a volatile Κ-contraction (vw=1 factorization):** frontier
  volatile → no swap → giant priced at `K_full` → the model knows it is
  expensive and will prefer a factorization that puts a persistent contraction
  above it (which is exactly where the runtime can batch it).
- **μ̃⁴:** degree 0 in K → `peak_slice == peak_full` → stays full, correctly
  costed as the non-batchable object it is.

## 7. Custom injection

A user supplies their own model two ways:

1. **Policy type** satisfying the `CostModel` concept (full control of `State`,
   `combine`, `value`, `evaluate`).
2. **`LambdaCostModel` adaptor** wiring user-supplied `leaf` / `combine` /
   `value` (and optionally `evaluate`) lambdas into the concept, for the common
   case where a user wants a bespoke additive or peak-like metric without
   defining a type.

Batchability is *not* something a custom lambda must re-implement: it is driven
by the shared declarative triple (`is_batchable`, `batch_target_size`,
`is_volatile`). A custom peak-style model reuses the framework's batch-plug; a
custom additive model simply ignores it.

## 8. Integration & migration

### 8.1 SeQuant

- Add `CostModelBase` + the `CostModel` concept in `core/optimize/`.
- Add the three built-in models. `DensePeakSizeBatched`'s `evaluate` wraps /
  replaces today's `make_batched_custom_evaluator` call path.
- Template `single_term_opt` / `single_term_opt_impl` on the `CostModel` concept;
  the DP folds `leaf`/`combine` and minimizes `value`. The peak models carry the
  two-mode state and the volatile bitmask (already computed today).
- `OptimizeOptions`: replace the `objective_function` enum + `footprint_weight`
  with a `CostModel` (held by value or `shared_ptr`). `idx_to_extent`,
  `is_volatile_leaf`, `CSE.subnet`, `reorder` stay. `volatile_weight` moves onto
  `DenseFLOPs`.
- Keep `subnet_cse` working with peak: the canonical-subnet cost folded into
  `unique_meta_costs` becomes a `State` rather than a scalar; the CSE accounting
  composes via `combine`. (Detail to validate during implementation — see
  Risks.)

### 8.2 mpqc

- `SeQuantEngine` constructs one `CostModel` from input
  (`sequant:optimize:...`) and the batchability triple it already has, and hands
  it to both `optimize()` and the eval cache. The `accept_aux` predicate and
  `batch:aux_target_size` currently built in `cck.ipp` for the evaluator now also
  parameterize the model — one source of truth.
- Keyword surface: `sequant:optimize:cost_model` ∈
  `{dense_flops, dense_peak_size, dense_peak_size_batched}` (default
  `dense_flops` to preserve current behavior). `dense_peak_size_batched` pulls
  `is_batchable`/`batch_target_size` from the existing `batch` group.

### 8.3 Backward compatibility

`DenseFLOPs` is the default and reproduces today's factorization and (plain)
evaluation. `dense_size` is dropped (superseded; it was never a good objective —
§1). Existing inputs that set `sequant:optimize:objective_function` get a
deprecation mapping (`dense_flops` → `DenseFLOPs`; `dense_size` → `DenseFLOPs`
with a warning).

## 9. Testing

- **Unit (SeQuant):** on small fixed networks, assert `combine`/`value` give the
  hand-computed peak; assert order-selection picks the lower-peak schedule;
  assert the monomial footprint plugs `full` vs `slice` correctly; assert the
  frontier gate fires only when a batch index is internal *and* the subtree is
  persistent (construct one volatile-frontier and one persistent-frontier case
  and check the giant is/ isn't discounted).
- **Optimal-substructure regression:** brute-force the optimal peak tree for
  n ≤ 6 tensors and compare to the DP.
- **Integration (mpqc):** c6h14 pVDZ trace under `dense_peak_size_batched` —
  confirm the giant is kept and priced sliced, OSV modes do not appear (folded
  early), and μ̃⁴ does not appear; compare peak `hw`/`rss` against
  `dense_flops` vw=1 and vw=100. Then water_14, expecting the giant sliced + no
  OSV family + lower peak than both 200 GB (vw=100) and 290 GB (`dense_size`).

## 10. Risks & open questions

- **subnet_cse × peak.** Intra-term CSE today counts a shared subnet's cost once
  (`unique_meta_costs`). Folding a *peak* `State` once (rather than a scalar)
  needs care: a CSE-shared subtree is live differently than a freshly recomputed
  one. Validate that the peak recurrence remains consistent when a subnet cost is
  shared; if not, gate `subnet_cse` off for peak models in v1.
- **Multiple batch indices / spaces.** The monomial form handles a vector of
  batch indices, but a frontier may internalize several at different nodes.
  v1 target: a single batchable space (DF aux), the production case; document the
  multi-space generalization as future work.
- **Evaluator type erasure.** If the eval cache slot is type-erased per backend,
  wiring the model's template `evaluate` needs a one-line `this`-bound adaptor at
  set time. Confirm against the cache's `set_custom_evaluator` signature.
- **Peak vs the gate fix.** A batch-aware peak objective steers factorization
  toward persistent batch frontiers, where the *current* runtime gate already
  batches — so it may subsume the separately-discussed eval.hpp:1128 gate fix.
  Decide during implementation whether the gate fix is still wanted as a
  belt-and-suspenders for factorizations that cannot avoid a volatile frontier.
- **Cost of two-mode + order in the DP.** Doubling state (full/slice) and adding
  the order `min` increases the per-subset constant factor; the subset count
  (`2ⁿ`, capped at n ≤ 24) is unchanged. Expected acceptable; measure.
