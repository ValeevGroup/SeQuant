#ifndef SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
#define SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP

#include <SeQuant/core/expr.hpp>
#include <cstddef>

namespace sequant::mbpt::bernoulli {

/// Tensor-level H̄ = Σ_{k=0..rank} H̄^k in the Bernoulli expansion, for
/// σ = T−T† of rank N. H̄⁰..H̄⁴ are Eqs. (46)-(50) of 10.1063/1.5030344.
///
/// The result is coefficient tensors times normal-ordered operators, not
/// `mbpt::op` operators, and nothing is screened out of it, so the caller
/// projects every term.
///
/// @warning Single-reference only, and nothing checks for it. The N/R split
/// expands general indices over the hole and particle spaces alone (see
/// detail::expand_to_blocks) and classifies each one as wholly occupied or
/// wholly unoccupied. An active space, as in make_mr_spaces(), is neither, so
/// no term touching an active index can ever be classified N.
///
/// @pre an HF reference: F is taken to have no occupied-virtual block, which
/// keeps F out of H̄² and higher (see the F-cancellation in bernoulli.cpp)
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
/// @note @p expr is left untouched; the reduction runs on a clone.
ExprPtr wick_reduce(const ExprPtr& expr);

/// Normal-ordered commutator [A, B] = wick_reduce(A·B − B·A). NOT the bare
/// algebraic commutator: the operator product is Wick-reduced, so contractions
/// between A and B generate the lower-rank terms the Bernoulli expansion relies
/// on. Every index of B is reindexed to a fresh temporary first, making them
/// disjoint from A's.
ExprPtr wick_commutator(const ExprPtr& A, const ExprPtr& B);

/// Rewrites every general (non-base) index of the residual NormalOperator as
/// the sum over the hole/particle base spaces it spans, so that every residual
/// index is definite and the N/R classifier can act on it. The registry's other
/// base spaces are dropped (see the @warning on hbar). Idempotent on
/// block-resolved input.
/// @pre the registry specifies both a hole and a particle space
ExprPtr expand_to_blocks(const ExprPtr& expr);

/// Block-resolved N part (O_N of 10.1063/1.5030344): the terms whose single
/// residual NormalOperator is a pure excitation or pure de-excitation of rank ≤
/// @p cutoff. Applies expand_to_blocks first.
ExprPtr N_part(const ExprPtr& expr, std::size_t cutoff);

/// R part (O_R of 10.1063/1.5030344: expr minus its N part). Unlike N_part
/// the result is NOT block-resolved; it stays in compact general-index form,
/// which is much cheaper for the nested commutators that consume R.
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff);

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli

#endif  // SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
