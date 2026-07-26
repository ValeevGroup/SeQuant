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

/// Rewrites every general (non-base) index of the residual NormalOperator as
/// the sum over the hole/particle base spaces it spans (occupied/virtual), an
/// identity in the single-reference setting where the other base spaces are
/// empty. After expansion every residual index is definite so the
/// N/R classifier can act on it. Idempotent on block-resolved input.
ExprPtr expand_to_blocks(const ExprPtr& expr);

/// Block-resolved N part (O_N of 10.1063/1.5030344, defined above Eq. (43)):
/// the terms whose single residual
/// NormalOperator is a pure excitation or pure de-excitation of rank ≤
/// @p cutoff. Applies expand_to_blocks first.
ExprPtr N_part(const ExprPtr& expr, std::size_t cutoff);

/// Block-resolved R (rank-preserving remainder) part:
/// expand_to_blocks(expr) minus N_part(expr, cutoff).
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff);

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli

#endif  // SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
