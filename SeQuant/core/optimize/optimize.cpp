#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/optimize/sum.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

namespace sequant {

namespace {

index_to_extent_t default_idx_to_size() {
  return [](Index const& ix) { return ix.space().approximate_size(); };
}

/// Optimize a Product that contains only Tensor and scalar factors.
ExprPtr opt_pure_product(Product const& prod, index_to_extent_t const& idx2size,
                         OptFor opt_for) {
  if (opt_for == OptFor::Flops)
    return opt::single_term_opt<OptFor::Flops>(prod, idx2size);
  SEQUANT_ASSERT(opt_for == OptFor::Memsize);
  return opt::single_term_opt<OptFor::Memsize>(prod, idx2size);
}

/// Deliberately non-identifier label prefix used to stand in for non-Tensor,
/// non-scalar factors during single-term optimization. Chosen so that no
/// user-defined tensor label can collide with it.
inline constexpr std::wstring_view placeholder_label_prefix = L"@__opt_";

/// Optimize a Product that contains some non-Tensor, non-scalar factors by
/// substituting placeholder tensors with target indices, optimizing the
/// resulting tensor-only product, then swapping the originals back in.
ExprPtr opt_mixed_product(Product const& prod,
                          index_to_extent_t const& idx2size, OptFor opt_for) {
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
      Product{prod.scalar(), new_factors, Product::Flatten::No}, idx2size,
      opt_for);

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
ExprPtr optimize_impl(ExprPtr const& expr, index_to_extent_t const& idx2size,
                      OptFor opt_for, bool reorder, bool parallel_outer) {
  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    bool pure = ranges::all_of(prod, [](auto&& x) {
      return x->template is<Tensor>() || x->is_scalar();
    });
    return pure ? opt_pure_product(prod, idx2size, opt_for)
                : opt_mixed_product(prod, idx2size, opt_for);
  }

  if (expr->is<Sum>()) {
    auto const& in_sum = expr->as<Sum>();
    Sum::summands_type new_smands(in_sum.size());

    auto do_term = [&](std::size_t i) {
      new_smands[i] =
          optimize_impl(in_sum.summand(i), idx2size, opt_for,
                        /*reorder=*/false, /*parallel_outer=*/false);
    };

    if (parallel_outer && in_sum.size() > 1) {
      auto indices = ranges::views::iota(std::size_t{0}, in_sum.size());
      sequant::for_each(indices, do_term);
    } else {
      for (std::size_t i = 0; i < in_sum.size(); ++i) do_term(i);
    }

    Sum new_sum(std::move(new_smands), Sum::move_only_tag{});
    if (!reorder) return ex<Sum>(std::move(new_sum));

    // Binarize once per optimized summand and hand the nodes to reorder()
    // so they aren't re-built inside clusters().
    container::vector<FullBinaryNode<EvalExpr>> nodes;
    nodes.reserve(new_sum.size());
    for (auto const& s : new_sum.summands()) nodes.push_back(binarize(s));
    return ex<Sum>(opt::reorder(new_sum, nodes));
  }

  return expr->clone();
}

ExprPtr optimize_dispatch(ExprPtr const& expr, index_to_extent_t idx2size,
                          OptFor opt_for, ReorderSum reorder) {
  if (!idx2size) idx2size = default_idx_to_size();
  return optimize_impl(expr, idx2size, opt_for, reorder == ReorderSum::Reorder,
                       /*parallel_outer=*/true);
}

}  // namespace

ExprPtr optimize(ExprPtr const& expr, OptFor opt_for, ReorderSum reorder) {
  return optimize_dispatch(expr, default_idx_to_size(), opt_for, reorder);
}

ExprPtr optimize(ExprPtr const& expr, index_to_extent_t idx2size,
                 OptFor opt_for, ReorderSum reorder) {
  return optimize_dispatch(expr, std::move(idx2size), opt_for, reorder);
}

ResultExpr& optimize(ResultExpr& expr, OptFor opt_for, ReorderSum reorder) {
  expr.expression() = optimize(expr.expression(), opt_for, reorder);
  return expr;
}

ResultExpr& optimize(ResultExpr& expr, index_to_extent_t idx2size,
                     OptFor opt_for, ReorderSum reorder) {
  expr.expression() =
      optimize(expr.expression(), std::move(idx2size), opt_for, reorder);
  return expr;
}

ResultExpr& optimize(ResultExpr&& expr, OptFor opt_for, ReorderSum reorder) {
  return optimize(expr, opt_for, reorder);
}

ResultExpr& optimize(ResultExpr&& expr, index_to_extent_t idx2size,
                     OptFor opt_for, ReorderSum reorder) {
  return optimize(expr, std::move(idx2size), opt_for, reorder);
}

}  // namespace sequant
