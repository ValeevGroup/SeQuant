#ifndef SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
#define SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>

#include <cstddef>
#include <string>

namespace sequant::mbpt::bernoulli {

/// Tensor-level @f$\bar{H} = \sum_{k=0}^{\mathrm{rank}} \bar{H}^{k}@f$ in
/// the Bernoulli expansion, for @f$\sigma = T-T^\dagger@f$ of rank @p N.
/// @f$\bar{H}^{0},\ldots,\bar{H}^{4}@f$ are Eqs. (46)-(50), and higher orders
/// are generated from Eq. (44), of 10.1063/1.5030344.
///
/// The result is coefficient tensors times normal-ordered operators, not
/// `mbpt::op` operators, and nothing is screened out of it, so the caller
/// projects every term.
///
/// @warning Single-reference only, and nothing checks for it. The N/R split
/// expands general indices over the hole and particle spaces alone (see
/// detail::expand_to_blocks) and classifies each one as wholly occupied or
/// wholly unoccupied relative to the single-product vacuum. That classification
/// does not represent multireference excitation semantics; for example, the
/// active space in make_mr_spaces() is vacuum-unoccupied.
///
/// @pre a HF reference: @f$F@f$ is taken to have no occupied-virtual block,
/// which keeps @f$F@f$ out of @f$\bar{H}^{2}@f$ and higher (see the
/// F-cancellation in bernoulli.cpp)
///
/// @param N cluster/excitation rank (also the N/R rank cutoff)
/// @param rank highest Bernoulli order @f$\bar{H}^{k}@f$ to include
/// @param skip1 exclude singles from T
/// @throw Exception if CSV is enabled
ExprPtr hbar(std::size_t N, std::size_t rank, bool skip1);

namespace detail {

/// Exact coefficients keyed by a nested-commutator partition path.
///
/// A path begins with the partition applied to V and contains one additional
/// character for the partition after each commutator with @f$\sigma@f$: `A`
/// keeps all terms, `N` keeps pure excitation/de-excitation terms, and `R`
/// keeps the remainder. Thus an order-@f$k@f$ path has @f$k+1@f$ characters
/// and ends in `A`.
using PartitionPathCoefficients = container::map<std::string, rational>;

/// Generates the partition paths and exact rational coefficients of every
/// homogeneous Bernoulli contribution through @p max_order.
///
/// For @f$b_n@f$ defined by
///
/// @f[
///   \frac{x}{\exp(x)-1} = \sum_{n=0}^{\infty} b_n x^n,
/// @f]
/// this applies Eq. (44) of 10.1063/1.5030344 in the form
/// @f[
///   V_m = (-1)^m b_m \operatorname{ad}_{\sigma}^{m}(V)
///         - \sum_{j=1}^{m} b_j
///           \operatorname{ad}_{\sigma}^{j}((V_{m-j})_R).
/// @f]
///
/// Duplicate paths are combined exactly, zero coefficients are removed, and
/// pairs related by @f$A-R=N@f$ are represented by the corresponding `N` path.
/// @return one coefficient map per order, indexed from 0 through @p max_order
container::svector<PartitionPathCoefficients> nested_commutator_coefficients(
    std::size_t max_order);

/// Applies Wick's theorem to @p expr retaining PARTIAL contractions,
/// reducing a product of normal-ordered operators to a sum of normal-ordered
/// operators, each consisting of a coefficient tensor and at most one residual
/// NormalOperator; fully-contracted terms carry none. Unlike the
/// expectation-value path it keeps operators rather than collapsing to a
/// scalar VEV.
/// @note @p expr is left untouched; the reduction runs on a clone.
ExprPtr wick_reduce(const ExprPtr& expr);

/// Normal-ordered commutator
/// @f$[A,B] = \operatorname{wick\_reduce}(AB-BA)@f$. NOT the bare algebraic
/// commutator: the operator product is Wick-reduced, so contractions between
/// @f$A@f$ and @f$B@f$ generate the lower-rank terms the Bernoulli expansion
/// relies on. Every index of @f$B@f$ is reindexed to a fresh temporary first,
/// making them disjoint from @f$A@f$'s.
ExprPtr wick_commutator(const ExprPtr& A, const ExprPtr& B);

/// Rewrites every general (non-base) index of the residual NormalOperator as
/// the sum over the hole/particle base spaces it spans, so that every residual
/// index is definite and the N/R classifier can act on it. The registry's other
/// base spaces are dropped (see the @warning on hbar). Idempotent on
/// block-resolved input.
/// @pre the registry specifies both a hole and a particle space
ExprPtr expand_to_blocks(const ExprPtr& expr);

/// Block-resolved @f$N@f$ part (@f$O_N@f$ of 10.1063/1.5030344): the terms
/// whose single residual NormalOperator is a pure excitation or pure
/// de-excitation of rank in
/// @f$[\mathtt{min\_rank},\mathtt{cutoff}]@f$. Applies expand_to_blocks first.
/// @param expr expression to partition
/// @param cutoff highest excitation/de-excitation rank in @f$O_N@f$
/// @param min_rank lowest rank @f$\sigma@f$ carries, so 2 when singles are
/// skipped. Below it there is no amplitude and the condition
/// @f$\bar{V}_N=0@f$ that justifies calling a term @f$N@f$ does not hold.
ExprPtr N_part(const ExprPtr& expr, std::size_t cutoff,
               std::size_t min_rank = 1);

/// @f$R@f$ part (@f$O_R@f$ of 10.1063/1.5030344: @p expr minus its @f$N@f$
/// part). Unlike N_part the result is NOT block-resolved; it stays in compact
/// general-index form, which is much cheaper for the nested commutators that
/// consume @f$R@f$.
/// @param expr expression to partition
/// @param cutoff highest excitation/de-excitation rank in @f$O_N@f$
/// @param min_rank lowest rank @f$\sigma@f$ carries
ExprPtr R_part(const ExprPtr& expr, std::size_t cutoff,
               std::size_t min_rank = 1);

}  // namespace detail

}  // namespace sequant::mbpt::bernoulli

#endif  // SEQUANT_DOMAIN_MBPT_BERNOULLI_HPP
