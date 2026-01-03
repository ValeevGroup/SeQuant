#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <Eigen/Core>

namespace sequant {

namespace {
static constexpr double default_biorth_threshold = 1e-12;
}  // namespace

[[nodiscard]] ResultExpr biorthogonal_transform_copy(
    const ResultExpr& expr, double threshold = default_biorth_threshold);

[[nodiscard]] container::svector<ResultExpr> biorthogonal_transform_copy(
    const container::svector<ResultExpr>& exprs,
    double threshold = default_biorth_threshold);

void biorthogonal_transform(ResultExpr& expr,
                            double threshold = default_biorth_threshold);

void biorthogonal_transform(container::svector<ResultExpr>& exprs,
                            double threshold = default_biorth_threshold);

/// performs symbolic biorthogonal transform of CC-like equation using
///(for rank-3 and higher
/// Wang-Knizia biorthogonalization (https://arxiv.org/abs/1805.00565) is used
/// uses hardcoded coefficients for ranks 1-6, computed coefficients for higher
/// ranks
[[nodiscard]] ExprPtr biorthogonal_transform(
    const ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups = {},
    double threshold = default_biorth_threshold);

/// \brief Computes the non-null space (NNS) projection matrix
///
/// \param n_particles The rank of external index pairs
/// \param threshold The threshold to compute the pseudoinverse matrix
///        (set to default_biorth_threshold)
/// \return The computed NNS projection matrix
[[nodiscard]] Eigen::MatrixXd compute_nns_p_matrix(
    std::size_t n_particles, double threshold = default_biorth_threshold);

/// \brief Provides permuted indices using libperm unrank function
///
/// \param indices The indices to permute
/// \param perm_rank The rank of the permutation
/// \param n_particles The rank of external index pairs
///
/// \return The permuted indices
container::svector<size_t> compute_permuted_indices(
    const container::svector<size_t>& indices, size_t perm_rank,
    size_t n_particles);

}  // namespace sequant

#endif
