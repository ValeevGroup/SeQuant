//
// Created by Eduard Valeyev on 3/31/18.
//

#ifndef SEQUANT2_WICK_IMPL_HPP
#define SEQUANT2_WICK_IMPL_HPP

namespace sequant2 {

namespace detail {

inline std::map<Index, Index> compute_index_replacement_rules(const std::shared_ptr<Product> &product,
                                                              const container::vector<Index> &external_indices) {
  std::map<Index, Index> result;
  for(const auto& factor : product->factors()) {
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      const auto& tensor = static_cast<const Tensor&>(*factor);
      if (tensor.label() == L"S") {
        assert(tensor.bra().size() == 1);
        assert(tensor.ket().size() == 1);
      }
    }
  }
  return result;
}

inline void apply_index_replacement_rules(std::shared_ptr<Product> &expr, const std::map<Index, Index> &replacement_rules) {
  // not yet implemented
}

/// resolves Kronecker deltas (=overlaps between indices in orthonormal spaces) in sums
inline void reduce_wick_impl(std::shared_ptr<Product> &expr, const container::vector<Index> &external_indices) {
  const auto replacement_rules = compute_index_replacement_rules(expr, external_indices);
  apply_index_replacement_rules(expr, replacement_rules);
}

}  // namespace sequant2::detail

template<Statistics S>
void
WickTheorem<S>::reduce(ExprPtr &expr) const {
  // there are 2 possibilities: expr is a single Product, or it's a Sum of Products
  if (expr->type_id() == Expr::get_type_id<Product>()) {
    auto expr_cast = std::static_pointer_cast<Product>(expr);
    detail::reduce_wick_impl(expr_cast, external_indices_);
    expr = expr_cast;
  } else {
    assert(expr->type_id() == Expr::get_type_id<Sum>());
    for (auto &subexpr: *expr) {
      assert(subexpr->type_id() == Expr::get_type_id<Product>());
      auto subexpr_cast = std::static_pointer_cast<Product>(subexpr);
      detail::reduce_wick_impl(subexpr_cast, external_indices_);
      subexpr = subexpr_cast;
    }
  }
}

}  // namespace sequant2

#endif //SEQUANT2_WICK_IMPL_HPP
