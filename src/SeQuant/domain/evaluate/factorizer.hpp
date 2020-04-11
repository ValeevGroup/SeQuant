#ifndef SEQUANT_EVALUATE_FACTORIZER_HPP
#define SEQUANT_EVALUATE_FACTORIZER_HPP

///
/// Find optimal evaluation sequences of EvalTensor evaluation tree that result
/// into lower ops count.
///
/// @author Bimal Gaudel
/// @version March 2020
///

#include "eval_tensor_builder.hpp"
#include "eval_tensor_fwd.hpp"

#include <SeQuant/core/expr.hpp>
#include <algorithm>
#include <tuple>

namespace sequant::evaluate {
namespace detail {
/// Generate a subexpression of an expr made of those present in indices
/// positions @c idx_vec.
/// @param expr The expression to extract subexpression from.
/// @param idx_vec Vector of indices of subexpressions in @c expr.
ExprPtr getSubExpr(const ExprPtr& expr,
                   const container::vector<size_t>& idx_vec) {
  if (expr->is<Product>()) {
    auto& prod = expr->as<Product>();
    auto result = std::make_shared<Product>();
    for (auto idx : idx_vec) {
      result->append(1, prod.factor(idx));
    }
    return result;
  } else if (expr->is<Sum>()) {
    auto& sum = expr->as<Sum>();
    auto result = std::make_shared<Sum>();
    for (auto idx : idx_vec) {
      result->append(sum.summand(idx));
    }
    return result;
  } else if (expr->is<Tensor>() && *idx_vec.begin() == 0) {
    return expr;
  } else {
    throw std::domain_error(
        "Only Tensor, Sum or Product type expression expected!");
  }
}  // getSubExpr

}  // namespace detail

///
/// Find the largest common sub network between a pair of tensor networks.
///
/// @param exprA Expression to find common subnetwork of.
/// @param exprB Expression to find common subnetwork of.
/// @param builder EvalTensorBuilder object to give a context while interpreting
/// tensors.
///
/// @return Tuple of ExprPtr that correspond to equivalent data tensors.
template <typename T>
std::tuple<ExprPtr, ExprPtr> largest_common_subnet(
    const ExprPtr& exprA, const ExprPtr& exprB,
    const EvalTensorBuilder<T>& builder) {
  // finding the positions of the common sub-expressions
  container::svector<size_t> commonIdxA;
  container::svector<size_t> commonIdxB;
  //
  auto i = 0;
  for (const auto& subA : *exprA) {
    auto j = 0;
    for (const auto& subB : *exprB) {
      if (j >= i) {  // skip redundant computations
        if (builder.build_tree(subA)->get_hash_value() ==
            builder.build_tree(subB)->get_hash_value()) {
          commonIdxA.push_back(i);
          commonIdxB.push_back(j);
        }
      }
      ++j;
    }
    ++i;
  }

  // finding the sub-networks
  //
  // hold the best subnet found so far
  auto optimal_subnetA = ExprPtr(nullptr), optimal_subnetB = ExprPtr(nullptr);
  // the operation counts for the optimal_subnets
  auto optimal_ops_count = 0;
  //
  // choose determines the size of a combination
  auto choose = commonIdxA.size();  // or commonIdxB.size()
  while (choose > 0) {
    // generating unique combinations, of decreasing size,
    // from the vector of
    //
    // https://stackoverflow.com/questions/9430568/generating-combinations-in-c
    //           -answer by mitchnull

    // the combination will be generated in this vector
    container::vector<size_t> combination;

    container::vector<size_t> grid(commonIdxA.size());
    std::fill(grid.begin(), grid.begin() + choose, true);

    do {
      // clear previous combination
      combination.clear();
      //
      for (auto k = 0; k < grid.size(); ++k) {
        if (grid[k]) {
          combination.push_back(k);
        }
      }
      // generated a combination, do something with it
      //
      // get the two subnets
      container::vector<size_t> index_selectorA, index_selectorB;
      for (auto& idx : combination) {
        index_selectorA.push_back(commonIdxA[idx]);
        index_selectorB.push_back(commonIdxB[idx]);
      }
      auto subNetA = detail::getSubExpr(exprA, index_selectorA);
      auto subNetB = detail::getSubExpr(exprB, index_selectorB);
      //
      // now check if the two subnets are equivalent
      // using eval_tensor for now
      auto evalTreeA = builder.build_tree(subNetA);
      auto evalTreeB = builder.build_tree(subNetB);

      if (evalTreeA->get_hash_value() == evalTreeB->get_hash_value()) {
        // found identical subnets
        return std::make_tuple(subNetA, subNetB);
      }
    } while (std::prev_permutation(grid.begin(), grid.end()));
    // if nothing found so far proceed to find a smaller subnet
    --choose;
  }
  // if not even a single subnet is found in common
  return std::make_tuple(nullptr, nullptr);
}

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_FACTORIZER_HPP */
