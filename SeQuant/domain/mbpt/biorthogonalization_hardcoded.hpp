#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP

#include <SeQuant/core/math.hpp>

#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace sequant {

/// \brief Provides the first row of the biorthogonal coefficients matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The Myrvold-Ruskey unrank1 algorithm (doi.org/10.1016/S0020-0190(01)00141-7)
/// is used to order permutations, then the permutational overlap matrix M is
/// constructed with elements (-2)^{c} × (-1)^{n_particles}, where c is the
/// number of cycles in the relative permutation.
///
/// The biorthogonal coefficients are obtained from the normalized pseudoinverse
/// of M: first compute M_pinv (the pseudoinverse), then normalize it by the
/// factor ((n_particles)!/rank(M)).
/// Finally, biorthogonal coefficients = normalized_M_pinv · e_1,
/// where e_1 is the first unit vector.
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of rational coefficients representing the first row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
std::vector<sequant::rational> hardcoded_biorth_coeffs_first_row(
    std::size_t n_particles);

/// \brief Constructs the full biorthogonal coefficient matrix from the first
/// row
///
/// \param first_row The first row of rational biorthogonal coefficients matrix
/// \param n_particles The rank of external index pairs
///
/// \return The complete biorthogonal coefficient matrix
Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
make_hardcoded_biorth_coeffs_matrix(
    const std::vector<sequant::rational>& first_row, std::size_t n_particles);

/// \brief Provides the hardcoded biorthogonal coefficient matrix
///
/// \param n_particles The rank of external index pairs
///
/// \return The biorthogonal coefficient matrix as a matrix of rational numbers
Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
hardcoded_biorth_coeffs_matrix(std::size_t n_particles);

#include "biorthogonalization_hardcoded.ipp"

}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
