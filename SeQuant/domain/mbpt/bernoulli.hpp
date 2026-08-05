#ifndef SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
#define SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP

#include <SeQuant/core/expr.hpp>
#include <cstddef>

namespace sequant::mbpt::bernoulli {

/// Tensor-level H̄ = Σ_{k=0..rank} H̄^k in the Bernoulli expansion, for
/// σ = T−T† of rank N.
///
/// The Bernoulli expansion rewrites the non-terminating UCC
/// similarity-transform series so that Bernoulli numbers appear as the
/// expansion coefficients; the rank-by-rank operators H̄⁰..H̄⁴ are Eqs. (46)-(50)
/// of 10.1063/1.5030344.
///
/// @warning Single-reference only. The N/R split expands general indices over
/// the hole and particle spaces alone (see detail::expand_to_blocks), dropping
/// any other base space the registry defines. That is harmless only because the
/// single-reference projection manifolds annihilate the dropped terms. Under a
/// multireference registry they contribute, and both the N and the R part come
/// out wrong. Nothing checks for this.
///
/// @param N cluster/excitation rank (also the N/R rank cutoff)
/// @param rank highest Bernoulli order H̄^k to include (0..4)
/// @param skip1 exclude singles from T
/// @throw Exception if @p rank > 4
ExprPtr hbar(std::size_t N, std::size_t rank, bool skip1);

namespace detail {

/// Applies Wick's theorem to @p expr retaining PARTIAL contractions,
/// reducing a product of normal-ordered operators to a sum of normal-ordered
/// operators (each = coefficient tensor × at most one residual NormalOperator;
/// fully-contracted terms carry none). Unlike the expectation-value path it
/// keeps operators rather than collapsing to a scalar VEV.
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

/// Block-resolved N part (O_N of 10.1063/1.5030344): the terms whose single
/// residual NormalOperator is a pure excitation or pure de-excitation of rank ≤
/// @p cutoff. Applies expand_to_blocks first.
ExprPtr N_part(const ExprPtr& expr, std::size_t cutoff);

/// R part (O_R of 10.1063/1.5030344: expr minus its N part). Unlike N_part
/// the result is NOT block-resolved; it stays in compact general-index form.
/// This is exact because expand_to_blocks is an identity, and it is much
/// cheaper for the nested commutators that consume R.
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff);

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli

#endif  // SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
