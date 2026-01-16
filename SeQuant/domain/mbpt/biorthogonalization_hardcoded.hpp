#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP

#include <SeQuant/core/math.hpp>

#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace sequant {

// clang-format off
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
/// See [DOI 10.48550/ARXIV.1805.00565](https://doi.org/10.48550/ARXIV.1805.00565)
/// for more details.
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of rational coefficients representing the first row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
// clang-format on
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
/// \return Optional vector of NNS projector weights (std::nullopt if
/// n_particles not in [1,5])
Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
hardcoded_biorth_coeffs_matrix(std::size_t n_particles);

/// \brief Provides one row of the NNS projector matrix,
/// hardcoded from Mathematica to avoid numerical precision loss.
///
/// The NNS projector weights are obtained from the normalized pseudoinverse
/// of M: first compute M_pinv (the pseudoinverse), then normalize it by the
/// factor ((n_particles)!/rank(M)).
//// Finally, NNS projector = normalized_M_pinv · M.
///
/// \param n_particles The rank of external index pairs
///
/// \return Vector of NNS projector weights representing the last row
///
/// \throw std::runtime_error if n_particles is not in the range [1,5]
template <typename T>
  requires(std::floating_point<T> || meta::is_complex_v<T>)
std::optional<std::vector<T>> hardcoded_nns_projector(std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<T>{T(1) / T(1)};

    case 2:
      return std::vector<T>{T(1) / T(1), T(1) / T(1)};

    case 3:
      return std::vector<T>{T(-1) / T(5), T(-1) / T(5), T(-1) / T(5),
                            T(-1) / T(5), T(-1) / T(5), T(1) / T(1)};

    case 4:
      return std::vector<T>{
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(1) / T(7),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(1) / T(7),   T(-2) / T(7),
          T(-1) / T(14), T(1) / T(7),   T(-1) / T(14), T(-2) / T(7),
          T(1) / T(7),   T(-1) / T(14), T(-1) / T(14), T(-2) / T(7),
          T(-2) / T(7),  T(-2) / T(7),  T(-2) / T(7),  T(1) / T(1)};

    case 5:
      return std::vector<T>{
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(-1) / T(14), T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(2) / T(21),  T(-1) / T(14), T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(2) / T(21),  T(-1) / T(14), T(-1) / T(14),
          T(-1) / T(14), T(-1) / T(14), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(0) / T(1),
          T(2) / T(21),  T(2) / T(21),  T(-1) / T(21), T(2) / T(21),
          T(0) / T(1),   T(2) / T(21),  T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(-1) / T(21), T(1) / T(7),   T(0) / T(1),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(-1) / T(3),
          T(2) / T(21),  T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(0) / T(1),   T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(-1) / T(21), T(1) / T(7),   T(2) / T(21),  T(-1) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(0) / T(1),
          T(1) / T(7),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(-1) / T(21), T(-1) / T(21), T(-1) / T(21), T(-1) / T(21),
          T(1) / T(7),   T(-1) / T(21), T(2) / T(21),  T(2) / T(21),
          T(2) / T(21),  T(0) / T(1),   T(-1) / T(21), T(2) / T(21),
          T(2) / T(21),  T(2) / T(21),  T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(1) / T(7),   T(-1) / T(3),
          T(0) / T(1),   T(1) / T(7),   T(1) / T(7),   T(0) / T(1),
          T(-1) / T(3),  T(1) / T(7),   T(0) / T(1),   T(1) / T(7),
          T(0) / T(1),   T(-1) / T(3),  T(1) / T(7),   T(1) / T(7),
          T(0) / T(1),   T(0) / T(1),   T(-1) / T(3),  T(-1) / T(3),
          T(-1) / T(3),  T(-1) / T(3),  T(-1) / T(3),  T(1) / T(1)};

    default:
      return std::nullopt;
  }
}

}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
