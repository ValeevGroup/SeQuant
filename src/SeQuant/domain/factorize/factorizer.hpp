//
// Created by Bimal Gaudel on 12/20/19.
//

#ifndef SEQUANT_FACTORIZER_HPP
#define SEQUANT_FACTORIZER_HPP

#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <limits>

namespace sequant::factorize {
namespace detail {

using FlopsType = unsigned long long;

/// @brief uncontracted @c indices and @c flops holder
struct ContractionCostResult {
  container::svector<Index> indices;
  FlopsType flops = 0;

  ContractionCostResult() = default;
  explicit ContractionCostResult(const ExprPtr&);
};

/// @brief stores sizes of IndexSpace::Type's and
/// defines the parenthesis operator( )
struct ContractionCostCounter {
  const std::shared_ptr<container::map<IndexSpace::Type, size_t>> map_ptr = nullptr;
  /// @return ContractionCostResult
  ContractionCostResult operator()(const ContractionCostResult&,
                                   const ContractionCostResult&) const;
};

/// @brief a PathTree @c path and its corresponding @c flops
struct PathCostResult {
  PathCostResult() = default;
  std::shared_ptr<PathTree> path = std::make_shared<PathTree>(0);
  FlopsType flops = std::numeric_limits<FlopsType>::max();
};

/// @brief computes the flops for a path
/// @return ContractionCostResult
ContractionCostResult compute_path_cost(const std::shared_ptr<PathTree>&,
                                        const Product&,
                                        const ContractionCostCounter&);

/// @brief finds the optimal path for a given TensorNetwork
/// Note: a PathCostResult object has to be initialized and passed
/// as the last parameter and it will hold the optimal path result
void optimal_path(container::svector<std::shared_ptr<PathTree>>&,
                  const Product&, const ContractionCostCounter&,
                  std::shared_ptr<PathCostResult>&);

/// @brief translate a PathTree to a TensorNetwork object
/// @return ExprPtr
ExprPtr path_to_product(const std::shared_ptr<PathTree>&, const Product&);

}  // namespace detail

/// @brief optimize a TensorNetwork
ExprPtr factorize_product(
    const Product&,
    const std::shared_ptr<container::map<IndexSpace::Type, size_t>>&);

ExprPtr factorize_expr(const ExprPtr&,
                       const std::shared_ptr<container::map<IndexSpace::Type, size_t>>&,
                       bool factorize=true);

}  // namespace sequant::factorize

#endif  // SEQUANT_FACTORIZER_HPP
