#ifndef SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
#define SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP

#include <SeQuant/core/math.hpp>

#include <Eigen/Core>
#include <cstddef>
#include <vector>

namespace sequant {

std::vector<sequant::rational> get_first_row_biorth_coeffs_rational(
    std::size_t n_particles);

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
biorth_coeffs_from_first_row_rational(
    const std::vector<sequant::rational>& first_row, std::size_t n_particles);

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
get_hardcoded_biorth_coeffs_rational(std::size_t n_particles);

#include "biorthogonalization_hardcoded.ipp"

}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_BIORTHOGONALIZATION_HARDCODED_HPP
