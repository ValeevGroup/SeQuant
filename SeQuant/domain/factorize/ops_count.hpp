#ifndef SEQUANT_FACTORIZE_OPS_COUNT_HPP
#define SEQUANT_FACTORIZE_OPS_COUNT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/domain/factorize/eval_sequence.hpp>

namespace sequant::factorize {

/**
 * Store flops and the remaining indices after a binary contraction.
 */
struct OpsCalcResult {
  using ops_type = unsigned long long;
  using idx_container_type = container::set<Index, Index::LabelCompare>;

  ops_type flops;

  idx_container_type indices;
};

/**
 * Compute the number of operations/flops when a product is evaluated.
 *
 * @param prod Pointer to a Product.
 *
 * @param tree Rooted tree that dictates the sequence of evaluation.
 *
 * @param nocc Number of occupied orbitals.
 *
 * @param nvirt Number of virtual orbitals.
 *
 * @return OpsCalcResult.
 */
OpsCalcResult ops_count(const ExprPtr& prod, const rooted_tree& tree,
                        size_t nocc, size_t nvirt);

}  // namespace sequant::factorize

#endif  // SEQUANT_FACTORIZE_OPS_COUNT_HPP
