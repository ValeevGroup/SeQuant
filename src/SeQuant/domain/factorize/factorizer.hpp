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

namespace sequant {
namespace factorize {
namespace detail {

struct CostResult {
  sequant::container::svector<sequant::Index> indices;
  unsigned long flops = 0;

  CostResult() = default;
  CostResult(const ExprPtr&);
};

struct CostCounter {
  size_t nocc, nvirt;
  CostResult operator()(const CostResult&, const CostResult&);
};

CostResult compute_cost(const std::shared_ptr<PathTree>&,
                        const sequant::Product&, const CostCounter&);

void path_scanner(sequant::container::svector<std::shared_ptr<PathTree>>&,
                  const sequant::Product&, const CostCounter&);

}  // namespace detail
}  // namespace factorize
}  // namespace sequant

#endif  // SEQUANT_FACTORIZER_HPP
