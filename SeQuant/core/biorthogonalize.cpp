#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/permutation.hpp>
#include <SeQuant/core/utility/transform_expr.hpp>

#ifdef SEQUANT_HAS_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <algorithm>

namespace sequant {

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr,
    const container::svector<container::svector<sequant::Index>>&
        ext_index_groups,
    const double threshold) {
  assert(!ext_index_groups.empty());
  const auto n_particles = ext_index_groups.size();

  using sequant::container::svector;

  // Coefficients
  container::svector<rational> bt_coeff_vec;
  {
#ifdef SEQUANT_HAS_EIGEN
    using namespace Eigen;
    // Dimension of permutation matrix is n_particles!
    const auto n = boost::numeric_cast<Eigen::Index>(factorial(n_particles));

    // Permutation matrix
    Eigen::Matrix<double, Dynamic, Dynamic> M(n, n);
    {
      M.setZero();
      size_t n_row = 0;
      svector<int> v(n_particles), v1(n_particles);
      std::iota(v.begin(), v.end(), 0);
      std::iota(v1.begin(), v1.end(), 0);
      do {
        container::svector<double> permutation_vector;
        do {
          permutation_vector.push_back(
              std::pow(-2, sequant::count_cycles(v1, v)));
        } while (std::next_permutation(v.begin(), v.end()));
        Eigen::VectorXd pv_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
            permutation_vector.data(), permutation_vector.size());
        M.row(n_row) = pv_eig;
        ++n_row;
      } while (std::next_permutation(v1.begin(), v1.end()));
      M *= std::pow(-1, n_particles);
    }

    // Normalization constant
    double scalar;
    {
      auto nonZero = [&threshold](const double& d) {
        using std::abs;
        return abs(d) > threshold;
      };

      // Solve system of equations
      SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
      container::svector<double> eig_vals(eig_solver.eigenvalues().size());
      VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
          eig_solver.eigenvalues();

      double non0count =
          std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
      scalar = eig_vals.size() / non0count;
    }

    // Find Pseudo Inverse, get 1st row only
    MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
    container::svector<double> bt_coeff_dvec;
    bt_coeff_dvec.resize(pinv.rows());
    VectorXd::Map(&bt_coeff_dvec[0], bt_coeff_dvec.size()) =
        pinv.row(0) * scalar;
    bt_coeff_vec.reserve(bt_coeff_dvec.size());
    ranges::for_each(bt_coeff_dvec, [&bt_coeff_vec, threshold](double c) {
      bt_coeff_vec.emplace_back(to_rational(c, threshold));
    });

//    std::cout << "n_particles = " << n_particles << "\n bt_coeff_vec = ";
//    std::copy(bt_coeff_vec.begin(), bt_coeff_vec.end(),
//              std::ostream_iterator<rational>(std::cout, " "));
//    std::cout << "\n";
#else
    // hardwire coefficients for n_particles = 1, 2, 3
    switch (n_particles) {
      case 1:
        bt_coeff_vec = {ratio(1, 2)};
        break;
      case 2:
        bt_coeff_vec = {ratio(1, 3), ratio(1, 6)};
        break;
      case 3:
        bt_coeff_vec = {ratio(17, 120), ratio(-1, 120), ratio(-1, 120),
                        ratio(-7, 120), ratio(-7, 120), ratio(-1, 120)};
        break;
      default:
        throw std::runtime_error(
            "biorthogonal_transform requires Eigen library for n_particles > "
            "3.");
    }
#endif
  }

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
  assert(bt_coeff_vec.size() == bt_maps.size());

  // Checks if the replacement map is a canonical sequence
  auto is_canonical = [](const container::map<Index, Index>& idx_map) {
    bool canonical = true;
    for (auto&& pair : idx_map)
      if (pair.first != pair.second) return false;
    return canonical;
  };

  // Scale transformed expressions and append
  Sum bt_expr{};
  auto coeff_it = bt_coeff_vec.begin();
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
