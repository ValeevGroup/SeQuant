//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT2_EXPR_ALGORITHM_HPP
#define SEQUANT2_EXPR_ALGORITHM_HPP

namespace sequant2 {

inline ExprPtr
operator*(const ExprPtr& left, const ExprPtr& right) {
  // naive version is to just make a Product
  // TODO why is ExprPtrList needed?
  auto result = std::make_shared<Product>(ExprPtrList{left,right});
  return result;
}

inline ExprPtr
operator+(const ExprPtr& left, const ExprPtr& right) {
  // naive version is to just make a Sum
  // TODO why is ExprPtrList needed?
  auto result = std::make_shared<Sum>(ExprPtrList{left,right});
  return result;
}

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; will be replaced if canonicalization is impure
inline void canonicalize(ExprPtr& expr) {
  const auto biproduct = expr->canonicalize();
  if (biproduct && biproduct->type_id() == Expr::get_type_id<Constant>()) {
    const auto constant_ptr = std::static_pointer_cast<Constant>(biproduct);
    expr = biproduct * expr;
  }
}

namespace detail {
struct expand_visitor {
  void operator()(ExprPtr& expr) {
    if (debug) std::wcout << "expand_visitor received " << to_latex(expr) << std::endl;
    // apply expand() iteratively until done
    while(expand(expr)) {
      if (debug)  std::wcout << "after 1 round of expansion have " << to_latex(expr) << std::endl;
    }
    if (debug) std::wcout << "expansion result = " << to_latex(expr) << std::endl;
    // now need to flatten!
    // and simplify?! e.g. extract constants out of products, etc.
  }

  /// expands a Product
  bool expand_product(ExprPtr& expr) {
    auto& expr_ref = *expr;
    std::shared_ptr<Sum> result;
    const auto nsubexpr = ranges::size(*expr);
    for(std::size_t i=0; i != nsubexpr; ++i) {
      if (expr_ref[i]->type_id() == Expr::get_type_id<Sum>()) {
        // allocate the result, if not done yet
        if (!result)
          result = std::make_shared<Sum>();
        ExprPtr subexpr_to_expand = expr_ref[i];
        for(auto& subsubexpr: *subexpr_to_expand) {
          auto expr_clone = expr->clone();
          (*expr_clone)[i] = subsubexpr;
          result->append(std::move(expr_clone));
        }
        expr = std::static_pointer_cast<Expr>(result); // expanded one Sum, return
        return true;
      }
    }
    return false;
  }

  /// expands a Sum
  bool expand_sum(ExprPtr& expr) {
    bool expanded = false;
    auto& expr_ref = *expr;
    const auto nsubexpr = ranges::size(*expr);
    if (debug) std::wcout << "in expand_sum: expr = " << to_latex(expr) << std::endl;
    for(std::size_t i=0; i != nsubexpr; ++i) {
      if (expr_ref[i]->type_id() == Expr::get_type_id<Product>()) {
        expanded |= expand_product(expr_ref[i]);
        if (debug) std::wcout << "in expand_sum: after expand_product(" << (expanded ? "true)" : "false)") << " expr = " << to_latex(expr) << std::endl;
      }
    }
    return expanded;
  }

  // @return true if expanded Product of Sum into Sum of Product
  bool expand(ExprPtr& expr) {
    const auto type_id = expr->type_id();
    if (type_id == Expr::get_type_id<Product>()) {
      return expand_product(expr);
    }
    else if (type_id == Expr::get_type_id<Sum>()) {
      return expand_sum(expr);
    } else
      return false;
  }

  bool debug = false;
};
};  // namespace detail

/// Recursively expands products of sums
inline void expand(ExprPtr& expr) {
  detail::expand_visitor expander{};
  expr->visit(expander);
  expander(expr);
}

}  // namespace sequant2

#endif //SEQUANT2_EXPR_ALGORITHM_HPP
