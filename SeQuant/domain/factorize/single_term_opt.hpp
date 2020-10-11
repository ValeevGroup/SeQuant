#ifndef SEQUANT_FACTORIZE_SINGLE_TERM_OPT_HPP
#define SEQUANT_FACTORIZE_SINGLE_TERM_OPT_HPP

#include <SeQuant/core/expr.hpp>
#include "eval_sequence.hpp"
#include "ops_count.hpp"

namespace sequant::factorize {

/**
 * Exhaustively scan all possible evaluation sequences for the given product to
 * find the optimal sequence.
 *
 * @param expr Expression pointer to product.
 *
 * @param nocc Number of occupied orbitals.
 *
 * @param nvirt Number of virtual orbitals.
 *
 * @return Product made of the copies of factors in @c expr.
 */
ExprPtr sto_exhaustive_scan(const ExprPtr& expr, size_t nocc = 2,
                            size_t nvirt = 3);

/**
 * Repack a product in accordance with a evaluation sequence.
 *
 * @param prod Reference product to be repacked.
 * @param tree Dictates the way of repacking.
 * @param scale Scale the result by the scalar of @c prod. True by default.
 *
 * @return A new expression.
 */
ExprPtr repack_prod(const ExprPtr& prod, const factorize::rooted_tree& tree,
                    bool scale = true);

/**
 * To be used as a function object to find the optimal evaluation sequence for a
 * given prodcut.
 */
struct OptimalRootedTree {
  /** Number of occupied orbitals. */
  size_t nocc;

  /** Number of virtual orbitals. */
  size_t nvirt;

  /** Reference to a product. */
  const ExprPtr& prod;

  /** Optimal evaluation sequence. */
  std::optional<rooted_tree> tree_;

  /** Flops for the optimal evaluation. */
  OpsCalcResult::ops_type flops;

  /**
   * Ctor.
   *
   * @param prod Reference to a product.
   * @param nocc The number of occupied orbitals.
   * @param nvirt The number of virtual orbitals.
   */
  OptimalRootedTree(const ExprPtr& prod, size_t nocc, size_t nvirt);

  /**
   * Compute flops for a given evaluation sequence.
   *
   * @param tree Evaluation sequence for @c prod.
   */
  void operator()(const rooted_tree& tree);

  /**
   * @return The optimal evaluation sequence.
   */
  rooted_tree tree();

};  // struct OptimalRootedTree

}  // namespace sequant::factorize
#endif  // SEQUANT_FACTORIZE_SINGLE_TERM_OPT_HPP
