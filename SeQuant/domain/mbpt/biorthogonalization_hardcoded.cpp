#include <SeQuant/domain/mbpt/biorthogonalization_hardcoded.hpp>

#include <libperm/Permutation.hpp>
#include <libperm/Rank.hpp>
#include <libperm/Utils.hpp>

namespace sequant {

std::vector<sequant::rational> hardcoded_biorth_coeffs_first_row(
    std::size_t n_particles) {
  switch (n_particles) {
    case 1:
      return std::vector<sequant::rational>{ratio(1, 2)};

    case 2:
      return std::vector<sequant::rational>{ratio(1, 3), ratio(1, 6)};

    case 3:
      return std::vector<sequant::rational>{ratio(17, 120), ratio(-7, 120),
                                            ratio(-1, 120), ratio(-1, 120),
                                            ratio(-1, 120), ratio(-7, 120)};

    case 4:
      return std::vector<sequant::rational>{
          ratio(43, 840), ratio(-19, 1680), ratio(-19, 1680),
          ratio(-1, 105), ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(1, 120),    ratio(-1, 105),
          ratio(1, 120),  ratio(-1, 105),   ratio(-19, 1680),
          ratio(-1, 105), ratio(1, 120),    ratio(1, 120),
          ratio(13, 840), ratio(-1, 105),   ratio(-1, 105),
          ratio(1, 120),  ratio(-19, 1680), ratio(-19, 1680),
          ratio(13, 840), ratio(-19, 1680), ratio(1, 120)};

    case 5:
      return std::vector<sequant::rational>{
          ratio(59, 3780),   ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-31, 7560),  ratio(-5, 3024),
          ratio(-5, 3024),   ratio(-23, 30240), ratio(19, 7560),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-5, 3024),   ratio(19, 7560),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-31, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-23, 30240),
          ratio(-13, 7560),  ratio(-5, 3024),   ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(37, 15120),
          ratio(19, 7560),   ratio(-23, 30240), ratio(19, 7560),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-23, 30240), ratio(-31, 7560),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(19, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-31, 7560),
          ratio(37, 15120),  ratio(-31, 7560),  ratio(37, 15120),
          ratio(-5, 3024),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(37, 15120),
          ratio(-31, 7560),  ratio(37, 15120),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(19, 7560),   ratio(37, 15120),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-13, 7560),
          ratio(-23, 30240), ratio(37, 15120),  ratio(-13, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(-23, 30240),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-23, 30240),
          ratio(19, 7560),   ratio(-13, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(19, 7560),   ratio(-31, 7560),  ratio(-31, 7560),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-5, 3024),
          ratio(37, 15120),  ratio(-13, 7560),  ratio(37, 15120),
          ratio(-13, 7560),  ratio(-23, 30240), ratio(-5, 3024),
          ratio(19, 7560),   ratio(-23, 30240), ratio(-5, 3024),
          ratio(37, 15120),  ratio(-5, 3024),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-23, 30240), ratio(-13, 7560),
          ratio(19, 7560),   ratio(19, 7560),   ratio(-23, 30240),
          ratio(-23, 30240), ratio(-13, 7560),  ratio(-5, 3024),
          ratio(19, 7560),   ratio(-5, 3024),   ratio(-23, 30240),
          ratio(37, 15120),  ratio(37, 15120),  ratio(-13, 7560),
          ratio(-13, 7560),  ratio(37, 15120),  ratio(-23, 30240)};

    default:
      throw std::runtime_error(
          "hardcoded biorthogonal coefficients only available for ranks 1-5, "
          "requested rank is : " +
          std::to_string(n_particles));
  }
}

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
make_hardcoded_biorth_coeffs_matrix(
    const std::vector<sequant::rational>& first_row, std::size_t n_particles) {
  const auto n = first_row.size();
  Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic> M(n, n);

  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < n; ++col) {
      perm::Permutation row_perm = perm::unrank(n - 1 - row, n_particles);
      perm::Permutation col_perm = perm::unrank(col, n_particles);

      col_perm->preMultiply(row_perm);

      std::size_t source_idx = perm::rank(col_perm, n_particles);
      M(row, col) = first_row[source_idx];
    }
  }
  return M;
}

Eigen::Matrix<sequant::rational, Eigen::Dynamic, Eigen::Dynamic>
hardcoded_biorth_coeffs_matrix(std::size_t n_particles) {
  auto first_row = hardcoded_biorth_coeffs_first_row(n_particles);
  return make_hardcoded_biorth_coeffs_matrix(first_row, n_particles);
}

}  // namespace sequant
