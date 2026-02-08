#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/single_term.hpp>
#include <SeQuant/core/optimize/sum.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/range/access.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/view.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>

namespace sequant {

namespace opt {

///
/// \param expr  Expression to be optimized.
/// \param idxsz An invocable object that maps an Index object to size.
/// \param reorder_sum If true, the summands are reordered so that terms with
///                    common sub-expressions appear closer to each other.
/// \return Optimized expression for lower evaluation cost.
template <typename IdxToSize, typename = std::enable_if_t<std::is_invocable_r_v<
                                  size_t, IdxToSize, const Index&>>>
ExprPtr optimize(ExprPtr const& expr, IdxToSize const& idx2size,
                 bool reorder_sum) {
  using ranges::views::transform;
  if (expr->is<Product>()) {
    if (ranges::all_of(*expr, [](auto&& x) {
          return x->template is<Tensor>() || x->template is<Variable>();
        }))
      return opt::single_term_opt(expr->as<Product>(), idx2size);
    else {
      auto const& prod = expr->as<Product>();

      container::svector<ExprPtr> non_tensors(prod.size());
      container::svector<ExprPtr> new_factors;

      for (auto i = 0; i < prod.size(); ++i) {
        auto&& f = prod.factor(i);
        if (f.is<Tensor>() || f.is<Variable>())
          new_factors.emplace_back(f);
        else {
          non_tensors[i] = f;
          auto target_idxs = get_unique_indices(f);
          new_factors.emplace_back(
              ex<Tensor>(L"I_" + std::to_wstring(i), bra(target_idxs.bra),
                         ket(target_idxs.ket), aux(target_idxs.aux)));
        }
      }

      auto result = opt::single_term_opt(
          Product(prod.scalar(), new_factors, Product::Flatten::No), idx2size);

      auto replacer = [&non_tensors](ExprPtr& out) {
        if (!out->is<Tensor>()) return;
        auto const& tnsr = out->as<Tensor>();
        auto&& label = tnsr.label();
        if (label.at(0) == L'I' && label.at(1) == L'_') {
          size_t suffix = std::stoi(std::wstring(label.data() + 2));
          out = non_tensors[suffix].clone();
        }
      };

      result->visit(replacer, /* atoms_only = */ true);
      return result;
    }
  } else if (expr->is<Sum>()) {
    auto smands = *expr | transform([&idx2size](auto&& s) {
      return optimize(s, idx2size, /* reorder_sum= */ false);
    }) | ranges::to_vector;
    auto sum = Sum{smands.begin(), smands.end()};
    return reorder_sum ? ex<Sum>(opt::reorder(sum)) : ex<Sum>(std::move(sum));
  } else
    return expr->clone();
}

}  // namespace opt

ExprPtr optimize(ExprPtr const& expr, bool reorder_sum) {
  return opt::optimize(
      expr, [](Index const& ix) { return ix.space().approximate_size(); },
      reorder_sum);
}

ResultExpr& optimize(ResultExpr& expr, bool reorder_sum) {
  expr.expression() = optimize(expr.expression(), reorder_sum);

  return expr;
}

ResultExpr& optimize(ResultExpr&& expr, bool reorder_sum) {
  return optimize(expr, reorder_sum);
}

}  // namespace sequant
