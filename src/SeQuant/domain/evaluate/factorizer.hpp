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
  } else {
    throw std::domain_error("Only Sum or Product type expression expected!");
  }
}  // getSubExpr

}  // namespace detail

///
/// Find the largest common sub network between a pair of tensor networks.
/// If multiple tensor networks of the largest size are found, return the
/// one that results into minimal operations count while evaluating.
///
/// @param exprA Expression to find common subnetwork of.
/// @param exprB Expression to find common subnetwork of.
/// @param space_size_map Map from index space to its size.
///
/// @return ExprPtr to the largest common subnetwork with minimum ops count.
template <typename T>
ExprPtr largest_common_subnet(
    const ExprPtr& exprA, const ExprPtr& exprB,
    const container::map<IndexSpace::TypeAttr, size_t>& space_size_map,
    const EvalTensorBuilder<T>& builder) {
  // finding the positions of the common sub-expressions
  container::svector<size_t> commonIdxA;
  container::svector<size_t> commonIdxB;
  //
  auto i = 0;
  for (const auto& subA : *exprA) {
    auto j = 0;
    for (const auto& subB : *exprB) {
      if (subA == subB) {
        commonIdxA.push_back(i);
        commonIdxB.push_back(j);
      }
      ++j;
    }
    ++i;
  }

  // finding the sub-networks
  // optimal_subnet holds the best subnet found so far
  auto optimal_subnet = ExprPtr(nullptr);
  // the operation counts for the optimal_subnet
  auto optimal_ops_count = 0;
  // builds eval_tensor from sequant expressions assuming
  // the tensors refer to a real valued data tensor

  // choose determines the size of a combination
  auto choose = commonIdxA.size();  // or commonIdxB.size()
  while (choose > 0) {
    // generating unique combinations, of decreasing size,
    // from the common indices
    //
    // https://stackoverflow.com/questions/9430568/generating-combinations-in-c
    //           -- answer by mitchnull

    // the combination will be generated in this vector
    container::vector<size_t> combination;

    container::vector<size_t> grid(commonIdxA.size());
    std::fill(grid.begin(), grid.begin() + choose, true);

    do {
      // clear previous combination
      combination.clear();
      //
      for (auto i = 0; i < grid.size(); ++i) {
        if (grid[i]) {
          combination.push_back(i);
        }
      }
      // generated a combination, do something with it
      //
      // get the two subnets
      auto subNetA = detail::getSubExpr(exprA, combination);
      auto subNetB = detail::getSubExpr(exprB, combination);
      //
      // now check if the two subnets are equivalent
      // using eval_tensor for now
      auto evalTreeA = builder.build_tree(subNetA);
      auto evalTreeB = builder.build_tree(subNetB);

      if (evalTreeA->get_hash_value() != evalTreeB->get_hash_value()) {
        continue;
      }
      // found identical subnets
      // update the running subnet
      if (optimal_subnet) {
        // getting the operations count
        evalTreeA->set_ops_count(space_size_map);
        if (evalTreeA->get_ops_count() < optimal_ops_count) {
          // this one is better than the previous one, update
          optimal_subnet = subNetA;
          optimal_ops_count = evalTreeA->get_ops_count();
        }  // else nothing needs to be done, proceed for next subnet
      } else {
        // this is the first time we have encountered the common subnet
        optimal_subnet = subNetA;
        evalTreeA->set_ops_count(space_size_map);
        optimal_ops_count = evalTreeA->get_ops_count();
      }
    } while (std::prev_permutation(grid.begin(), grid.end()));
    // before finding a smaller subnet, if we have found one so far return it
    if (optimal_subnet) {
      return optimal_subnet;
    }
    // if nothing found so far proceed to find a smaller subnet
    --choose;
  }
  // if not even a single subnet is found in common
  return nullptr;
}

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_FACTORIZER_HPP */
