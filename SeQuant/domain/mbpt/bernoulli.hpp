#ifndef SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
#define SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP

#include <SeQuant/core/expr.hpp>
#include <cstddef>

namespace sequant::mbpt::bernoulli {

namespace detail {

/// Operator-valued Wick reduction: applies Wick's theorem to @p expr retaining
/// PARTIAL contractions, reducing a product of normal-ordered operators to a
/// sum of normal-ordered operators (each = coefficient tensor × one residual
/// NormalOperator). Unlike the expectation-value path it keeps operators rather
/// than collapsing to a scalar VEV.
ExprPtr wick_reduce(ExprPtr expr);

/// Normal-ordered commutator [A, B] = wick_reduce(A·B − B·A). NOT the bare
/// algebraic commutator: the operator product is Wick-reduced, so contractions
/// between A and B generate the lower-rank terms the Bernoulli expansion relies
/// on. B's summed indices are reindexed to fresh temporaries first, making them
/// disjoint from A's.
ExprPtr wick_commutator(const ExprPtr& A, const ExprPtr& B);

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli

#endif  // SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
