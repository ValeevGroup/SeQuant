#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/transform_expr.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>

namespace sequant {

Eigen::MatrixXd permutational_overlap_matrix(std::size_t n_particles) {
  const auto n = boost::numeric_cast<Eigen::Index>(factorial(n_particles));

  // The matrix only contains integer entries but all operations we want to do
  // with the matrix will (in general) require non-integer scalars which in
  // Eigen only works if you start from a non-integer matrix.
  Eigen::MatrixXd M(n, n);
  M.setZero();

  std::size_t n_row = 0;
  container::svector<int> v(n_particles), v1(n_particles);
  std::iota(v.begin(), v.end(), 0);
  std::iota(v1.begin(), v1.end(), 0);

  container::svector<double> permutation_vector;
  permutation_vector.reserve(n);
  do {
    permutation_vector.clear();
    do {
      permutation_vector.push_back(std::pow(-2, sequant::count_cycles(v1, v)));
    } while (std::next_permutation(v.begin(), v.end()));

    // TODO: M is symmetric -> we could make use of that in its construction
    M.row(n_row) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
        permutation_vector.data(), permutation_vector.size());

    ++n_row;
  } while (std::next_permutation(v1.begin(), v1.end()));

  if (n_particles % 2 != 0) {
    M *= -1;
  }

  return M;
}

container::svector<rational> compute_biorth_coeffs(std::size_t n_particles,
                                                   double threshold) {
  auto perm_ovlp_mat = permutational_overlap_matrix(n_particles);
  assert(perm_ovlp_mat.rows() == perm_ovlp_mat.cols());
  assert(perm_ovlp_mat == perm_ovlp_mat.transpose());

  // Find Pseudo Inverse, get 1st row only
  auto decomp =
      Eigen::CompleteOrthogonalDecomposition<decltype(perm_ovlp_mat)>();
  decomp.setThreshold(threshold);
  decomp.compute(perm_ovlp_mat);
  Eigen::MatrixXd pinv = decomp.pseudoInverse();

  // We need to normalize to the amount of non-zero eigenvalues via
  // normalization = #eigenvalues / #non-zero eigenvalues
  // Since perm_ovlp_mat is symmetric, it is diagonalizable and for every
  // diagonalizable matrix, its rank equals the amount of non-zero eigenvalues.
  double normalization =
      static_cast<double>(perm_ovlp_mat.rows()) / decomp.rank();

  container::svector<double> bt_coeff_dvec;
  bt_coeff_dvec.resize(pinv.rows());
  Eigen::VectorXd::Map(&bt_coeff_dvec[0], bt_coeff_dvec.size()) =
      pinv.row(0) * normalization;

  return bt_coeff_dvec | ranges::views::transform([&](double d) {
           return to_rational(d, threshold);
         }) |
         ranges::to<container::svector<rational>>();
}

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups,
    const double threshold) {
  assert(!ext_index_groups.empty());
  const auto n_particles = ext_index_groups.size();

  using sequant::container::svector;

  container::svector<rational> coeffs =
      compute_biorth_coeffs(n_particles, threshold);

  // Transformation maps
  container::svector<container::map<Index, Index>> bt_maps;
  {
    container::svector<Index> idx_list(ext_index_groups.size());

    for (std::size_t i = 0; i != ext_index_groups.size(); ++i) {
      idx_list[i] = *ext_index_groups[i].begin();
    }

    const container::svector<Index> const_idx_list = idx_list;

    do {
      container::map<Index, Index> map;
      auto const_list_ptr = const_idx_list.begin();
      for (auto& i : idx_list) {
        map.emplace(*const_list_ptr, i);
        const_list_ptr++;
      }
      bt_maps.push_back(map);
    } while (std::next_permutation(idx_list.begin(), idx_list.end()));
  }

  // If this assertion fails, change the threshold parameter
  assert(coeffs.size() == bt_maps.size());

  // Checks if the replacement map is a canonical sequence
  auto is_canonical = [](const container::map<Index, Index>& idx_map) {
    bool canonical = true;
    for (auto&& pair : idx_map)
      if (pair.first != pair.second) return false;
    return canonical;
  };

  // Scale transformed expressions and append
  Sum bt_expr{};
  auto coeff_it = coeffs.begin();
  for (auto&& map : bt_maps) {
    const auto v = *coeff_it;
    if (is_canonical(map))
      bt_expr.append(ex<Constant>(v) * expr->clone());
    else
      bt_expr.append(ex<Constant>(v) *
                     sequant::transform_expr(expr->clone(), map));
    coeff_it++;
  }
  ExprPtr result = std::make_shared<Sum>(bt_expr);
  return result;
}

}  // namespace sequant
