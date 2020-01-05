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

namespace sequant {
namespace factorize {
namespace detail {

/// @brief uncontracted @c indices and @c flops holder
struct ContractionCostResult {
  sequant::container::svector<sequant::Index> indices;
  unsigned long long flops = 0;

  ContractionCostResult() = default;
  ContractionCostResult(const ExprPtr&);
};

/// @brief stores sizes of IndexSpace::Type's and
/// defines the parenthesis operator( )
struct ContractionCostCounter {
  size_t nocc, nvirt;
  /// @return ContractionCostResult
  ContractionCostResult operator()(const ContractionCostResult&, const ContractionCostResult&);
};

/// @brief a PathTree @c path and its corresponding @c flops
struct PathCostResult {
  PathCostResult() = default;
  std::shared_ptr<PathTree> path  = std::make_shared<PathTree>(0);
  unsigned long long flops = std::numeric_limits<unsigned long long>::max();
};

/// @brief computes the flops for a path
/// @return ContractionCostResult
ContractionCostResult compute_path_cost(const std::shared_ptr<PathTree>&,
                        const sequant::Product&, const ContractionCostCounter&);

/// @brief finds the optimal path for a given TensorNetwork
/// Note: a PathCostResult object has to be intialized and passed
/// as the last parameter and it will hold the optimal path result
void optimal_path(sequant::container::svector<std::shared_ptr<PathTree>>&,
                  const sequant::Product&, const ContractionCostCounter&,
                  std::shared_ptr<PathCostResult>&);

/// @brief translate a PathTree to a TensorNetwork object
/// @return ExprPtr
sequant::ExprPtr path_to_product(const std::shared_ptr<PathTree>&,
    const sequant::Product&);

}  // namespace detail

/// @brief optimize a TensorNetwork
sequant::ExprPtr factorize_product(const sequant::Product&,
    const detail::ContractionCostCounter&);
}  // namespace factorize
}  // namespace sequant

#endif  // SEQUANT_FACTORIZER_HPP
