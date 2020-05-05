#ifndef SEQUANT_FACTORIZER_HPP
#define SEQUANT_FACTORIZER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor_network.hpp>

// IndexSpace type based hashing of tensors
#include <SeQuant/domain/evaluate/eval_tree.hpp>

#include <tuple>

///
/// Find common sub-networks between a pair of tensor networks.
/// @author Bimal Gaudel
/// @version May 2020
///

namespace sequant::factorize {
///
/// Find the largest common sub network between a pair of tensor networks.
///
/// @param exprA Expression to find common subnetwork of. All of the
/// subexpressions of @exprA must be ExprPtr to sequant Tensors'.
///
/// @param exprB Expression to find common subnetwork of. All of the
/// subexpressions of @exprB must be ExprPtr to sequant Tensors'.
///
/// @return Tuple of position vectors of the common subexpressions.
///
std::tuple<container::svector<size_t>, container::svector<size_t>>
largest_common_subnet(const ExprPtr& exprA, const ExprPtr& exprB);

/// Generate a subexpression of an expr made of those present at @c indices
/// positions in the @c expr.
/// @tparam Container A container type that allows range-based for looping.
/// @param expr The expression to extract subexpression from.
/// @param indices Vector of indices of subexpressions in @c expr.
template <typename Container>
ExprPtr _getSubExpr(const ExprPtr& expr, const Container& indices) {
  if (expr->is<Product>()) {
    auto& prod = expr->as<Product>();
    auto result = std::make_shared<Product>();
    for (auto idx : indices) {
      result->append(1, prod.factor(idx));
    }
    return result;
  } else if (expr->is<Sum>()) {
    auto& sum = expr->as<Sum>();
    auto result = std::make_shared<Sum>();
    for (auto idx : indices) {
      result->append(sum.summand(idx));
    }
    return result;
  } else if (expr->is<Tensor>() && *indices.begin() == 0) {
    return expr;
  } else {
    throw std::logic_error(
        "Only Tensor, Sum or Product type expression expected!");
  }
}  // _getSubExpr

}  // namespace sequant::factorize

#endif /* #ifndef SEQUANT_FACTORIZER_HPP */
