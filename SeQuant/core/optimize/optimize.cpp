#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/optimize/cost_model.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/optimize/sum.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {

namespace {

index_to_extent_t default_idx_to_size() {
  return [](Index const& ix) { return ix.space().approximate_size(); };
}

/// Diagnostic (env SEQUANT_FACTORIZER_DEBUG): for the chosen factorization
/// \p result of a single term, log each intermediate's result footprint AS THE
/// COST MODEL SIZES IT (idx_to_extent + inner_pow), plus the peak (the value
/// the DensePeakSize objective minimizes), so one can see why the factorizer
/// accepted a given intermediate -- e.g. an under-sized multi-composite tensor.
/// Footprints in mega-elements; outer{...} lists each free outer index extent,
/// inner{Np:e,...} lists each CSV/PNO composite's proto-index count N and the
/// extent e the model assigns it.
void log_chosen_factorization(ExprPtr const& result,
                              OptimizeOptions const& opts) {
  if (!opts.idx_to_extent || !result) return;
  auto describe = [&](ExprPtr const& e) -> std::pair<double, std::string> {
    auto g = get_unique_indices(e);
    std::vector<Index> all;
    all.insert(all.end(), g.bra.begin(), g.bra.end());
    all.insert(all.end(), g.ket.begin(), g.ket.end());
    all.insert(all.end(), g.aux.begin(), g.aux.end());
    auto const tot = tot_indices(all);
    double const fp = opt::detail::inner_aware_volume(tot, opts.idx_to_extent,
                                                      opts.inner_pow);
    std::string s = "outer{";
    for (auto const& o : tot.outer)
      s += std::to_string(opts.idx_to_extent(o)) + ",";
    s += "} inner{";
    for (auto const& i : tot.inner)
      s += std::to_string(i.proto_indices().size()) +
           "p:" + std::to_string(opts.idx_to_extent(i)) + ",";
    s += "}";
    return {fp, std::move(s)};
  };
  double peak = 0.;
  std::string peak_s;
  std::function<void(ExprPtr const&)> rec = [&](ExprPtr const& e) {
    if (!e->is<Product>()) return;
    auto const desc = describe(e);
    std::clog << "[FACTORIZE-NODE] fp=" << desc.first / 1e6 << "Me "
              << desc.second << "\n";
    if (desc.first > peak) {
      peak = desc.first;
      peak_s = desc.second;
    }
    for (auto const& f : e->as<Product>().factors()) rec(f);
  };
  std::clog << "[FACTORIZE] -------- chosen tree --------\n";
  rec(result);
  std::clog << "[FACTORIZE] PEAK fp=" << peak / 1e6 << "Me " << peak_s << "\n";
  std::clog.flush();
}

/// Optimize a Product that contains only Tensor and scalar factors.
ExprPtr opt_pure_product(Product const& prod, OptimizeOptions const& opts) {
  bool const subnet_cse = opts.CSE.subnet;
  // Build the cost knobs field-by-field from OptimizeOptions / its BatchPolicy.
  // Batching config (both role predicates, batch_target_size, inner_pow,
  // batch_persistent_only) now travels on CostParams rather than as loose args.
  CostParams cost;
  cost.is_volatile_leaf = opts.batch_policy.is_volatile_leaf;
  cost.volatile_weight = opts.volatile_weight;
  cost.footprint_weight = opts.footprint_weight;
  cost.peak_flops_tolerance = opts.peak_flops_tolerance;
  cost.roofline = opts.roofline;
  cost.accumulation_factor = opts.batch_policy.accumulation_factor;
  cost.peak_threshold = opts.batch_policy.peak_threshold;
  cost.prune_outer_products = opts.prune_outer_products;
  cost.batch_spectator_indices = opts.batch_policy.batch_spectator_indices;
  cost.order_aware_recompute = opts.batch_policy.order_aware_recompute;
  cost.is_batchable_contracted_index =
      opts.batch_policy.is_batchable_contracted_index;
  cost.is_batchable_external_index =
      opts.batch_policy.is_batchable_external_index;
  cost.batch_target_size = opts.batch_policy.batch_target_size;
  cost.inner_pow = opts.inner_pow;
  cost.batch_persistent_only = opts.batch_policy.persistent_only;
  // Filled only by the DensePeakSizeBatched arm below (via out_axes); every
  // other objective leaves it empty, so the term_batch_axes insertion at the
  // end is then a no-op-shaped empty-vector entry (harmless: Task 3.3 only
  // consumes entries for summands the batched objective actually annotated).
  container::vector<NodeBatchAnnotation> node_axes;
  auto run = [&]() -> ExprPtr {
    if (opts.objective_function == ObjectiveFunction::DenseFLOPs)
      return opt::single_term_opt<ObjectiveFunction::DenseFLOPs>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSize)
      return opt::single_term_opt<ObjectiveFunction::DenseSize>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTime)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTime>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseTimeSpace)
      return opt::single_term_opt<ObjectiveFunction::DenseTimeSpace>(
          prod, opts.idx_to_extent, subnet_cse, cost);
    if (opts.objective_function == ObjectiveFunction::DenseSpaceTimeBatched)
      return opt::single_term_opt<ObjectiveFunction::DenseSpaceTimeBatched>(
          prod, opts.idx_to_extent, subnet_cse, cost,
          opts.term_batch_axes ? &node_axes : nullptr);
    SEQUANT_ASSERT(opts.objective_function ==
                   ObjectiveFunction::DenseTimeSpaceBatched);
    return opt::single_term_opt<ObjectiveFunction::DenseTimeSpaceBatched>(
        prod, opts.idx_to_extent, subnet_cse, cost,
        opts.term_batch_axes ? &node_axes : nullptr);
  };
  ExprPtr result = run();
  if (opts.term_batch_axes)
    (*opts.term_batch_axes)[result.get()] = std::move(node_axes);
  if (std::getenv("SEQUANT_FACTORIZER_DEBUG"))
    log_chosen_factorization(result, opts);
  return result;
}

/// Deliberately non-identifier label prefix used to stand in for non-Tensor,
/// non-scalar factors during single-term optimization. Chosen so that no
/// user-defined tensor label can collide with it.
inline constexpr std::wstring_view placeholder_label_prefix = L"@__opt_";

/// Optimize a Product that contains some non-Tensor, non-scalar factors by
/// substituting placeholder tensors with target indices, optimizing the
/// resulting tensor-only product, then swapping the originals back in.
ExprPtr opt_mixed_product(Product const& prod, OptimizeOptions const& opts) {
  container::svector<ExprPtr> non_tensors(prod.size());
  container::svector<ExprPtr> new_factors;
  new_factors.reserve(prod.size());

  for (std::size_t i = 0; i < prod.size(); ++i) {
    auto&& f = prod.factor(i);
    if (f->is<Tensor>() || f->is_scalar()) {
      new_factors.emplace_back(f);
    } else {
      non_tensors[i] = f;
      auto target_idxs = get_unique_indices(f);
      new_factors.emplace_back(ex<Tensor>(
          std::wstring(placeholder_label_prefix) + std::to_wstring(i),
          bra(target_idxs.bra), ket(target_idxs.ket), aux(target_idxs.aux)));
    }
  }

  auto result = opt_pure_product(
      Product{prod.scalar(), new_factors, Product::Flatten::No}, opts);

  auto replacer = [&non_tensors](ExprPtr& out) {
    if (!out->is<Tensor>()) return;
    auto label = out->as<Tensor>().label();
    if (!label.starts_with(placeholder_label_prefix)) return;

    // The placeholder prefix is internal; anything carrying it must have been
    // emitted by this function, with a pure-decimal suffix indexing
    // non_tensors. Any deviation is a programming error.
    auto suffix_view = label.substr(placeholder_label_prefix.size());
    SEQUANT_ASSERT(!suffix_view.empty());
    std::size_t suffix = 0;
    for (wchar_t c : suffix_view) {
      SEQUANT_ASSERT(c >= L'0' && c <= L'9');
      suffix = suffix * 10 + static_cast<std::size_t>(c - L'0');
    }
    SEQUANT_ASSERT(suffix < non_tensors.size() && non_tensors[suffix]);
    out = non_tensors[suffix].clone();
  };

  result->visit(replacer, /* atoms_only = */ true);
  return result;
}

/// Recursive workhorse. \p parallel_outer controls whether the (single)
/// outermost Sum's summands are processed in parallel; nested recursive
/// calls always run sequentially to avoid `sequant::for_each` oversubscription.
ExprPtr optimize_impl(ExprPtr const& expr, OptimizeOptions const& opts,
                      bool reorder, bool parallel_outer) {
  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    bool pure = ranges::all_of(prod, [](auto&& x) {
      return x->template is<Tensor>() || x->is_scalar();
    });
    return pure ? opt_pure_product(prod, opts) : opt_mixed_product(prod, opts);
  }

  if (expr->is<Sum>()) {
    auto const& in_sum = expr->as<Sum>();
    Sum::summands_type new_smands(in_sum.size());

    auto do_term = [&](std::size_t i) {
      new_smands[i] = optimize_impl(in_sum.summand(i), opts,
                                    /*reorder=*/false,
                                    /*parallel_outer=*/false);
    };

    // Thread-safety of the parallel branch rests on two invariants; do NOT
    // break them without re-auditing:
    //   1. Each task writes a distinct, pre-allocated new_smands[i] slot, and
    //      the work below (single_term_opt ->
    //      TensorNetwork::canonicalize_slots) operates on per-task *clones* of
    //      the input tensors. The lazily populated `mutable` caches on
    //      Expr/Index (hash_value_, label_, ...) are unsynchronized, so
    //      concurrent work must never read/write them on a shared (non-cloned)
    //      node. Index comparison touches only immutable members, so building
    //      index sets over shared indices is safe.
    //   2. The binarize() pass below DOES read Index::label() (a lazy cache
    //      write) on the optimized summands, so it is run *sequentially, after*
    //      for_each() has joined -- never inside do_term().
    // The default Context and cardinal_tensor_labels must also be configured
    // before entering here (their writes are unsynchronized unless
    // SEQUANT_CONTEXT_MANIPULATION_THREADSAFE); optimize() only reads them.
    if (parallel_outer && in_sum.size() > 1) {
      auto indices = ranges::views::iota(std::size_t{0}, in_sum.size());
      sequant::for_each(indices, do_term);
    } else {
      for (std::size_t i = 0; i < in_sum.size(); ++i) do_term(i);
    }

    Sum new_sum(std::move(new_smands), Sum::move_only_tag{});
    if (!reorder) return ex<Sum>(std::move(new_sum));

    // Binarize once per optimized summand and hand the nodes to reorder()
    // so they aren't re-built inside clusters(). NOTE: this runs sequentially
    // by design -- see invariant (2) above.
    container::vector<FullBinaryNode<EvalExpr>> nodes;
    nodes.reserve(new_sum.size());
    // per-summand binarize for ordering only; positional head doesn't escape.
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    for (auto const& s : new_sum.summands()) nodes.push_back(binarize(s));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    return ex<Sum>(opt::reorder(new_sum, nodes));
  }

  return expr->clone();
}

}  // namespace

ExprPtr optimize(ExprPtr const& expr, OptimizeOptions opts) {
  if (!opts.idx_to_extent) opts.idx_to_extent = default_idx_to_size();
  return optimize_impl(expr, opts, opts.reorder == ReorderSum::Reorder,
                       /*parallel_outer=*/true);
}

OptimizeResult optimize_result(ExprPtr const& expr, OptimizeOptions opts) {
  if (!opts.idx_to_extent) opts.idx_to_extent = default_idx_to_size();
  OptimizeResult res;
  res.expr = optimize_impl(expr, opts, opts.reorder == ReorderSum::Reorder,
                           /*parallel_outer=*/true);
  return res;
}

ResultExpr& optimize(ResultExpr& expr, OptimizeOptions opts) {
  expr.expression() = optimize(expr.expression(), std::move(opts));
  return expr;
}

ResultExpr& optimize(ResultExpr&& expr, OptimizeOptions opts) {
  return optimize(expr, std::move(opts));
}

// backwards compatibility overloads

namespace {
inline OptimizeOptions compatibility_opts(bool reorder_sum) {
  return OptimizeOptions{.reorder = reorder_sum ? ReorderSum::Reorder
                                                : ReorderSum::NoReorder};
}
}  // namespace

ExprPtr optimize(ExprPtr const& expr, bool reorder_sum) {
  return optimize(expr, compatibility_opts(reorder_sum));
}

ResultExpr& optimize(ResultExpr& expr, bool reorder_sum) {
  return optimize(expr, compatibility_opts(reorder_sum));
}

ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum) {
  return optimize(std::move(expr), compatibility_opts(reorder_sum));
}

}  // namespace sequant
