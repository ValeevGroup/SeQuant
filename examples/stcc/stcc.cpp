#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/formalism.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <clocale>

#include <Eigen/Eigenvalues>

using namespace sequant;

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr, int n_particles,
    const std::vector<std::vector<sequant::Index>>& ext_index_groups = {{}},
    double threshold = 1.e-12);
ExprPtr symmetrize_expr(
    ExprPtr& expr,
    const container::vector<container::vector<Index>>& ext_index_groups = {{}});

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__                    \
        << " in closed-shell spin-traced coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());                   \
  }

int main(int argc, char* argv[]) {
  std::setlocale(LC_ALL, "en_US.UTF-8");
  std::wcout.precision(std::numeric_limits<double>::max_digits10);
  std::wcerr.precision(std::numeric_limits<double>::max_digits10);
  std::wcout.sync_with_stdio(false);
  std::wcerr.sync_with_stdio(false);
  std::wcout.imbue(std::locale("en_US.UTF-8"));
  std::wcerr.imbue(std::locale("en_US.UTF-8"));
  std::wcout.sync_with_stdio(true);
  std::wcerr.sync_with_stdio(true);

  sequant::set_default_context(
      SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;

  /// Make external index
  auto ext_idx_list = [](const int i_max) {
    container::vector<container::vector<Index>> ext_idx_list;

    for (size_t i = 1; i <= i_max; ++i) {
      auto label = std::to_wstring(i);
      auto occ_i = Index::make_label_index(
          IndexSpace::instance(IndexSpace::active_occupied), label);
      auto virt_i = Index::make_label_index(
          IndexSpace::instance(IndexSpace::active_unoccupied), label);
      container::vector<Index> pair = {occ_i, virt_i};
      ext_idx_list.push_back(pair);
    }
    return ext_idx_list;
  };

  // Spin-orbital coupled cluster
  auto cc_r = sequant::mbpt::sr::cceqs{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Closed-shell spintrace (fast)
  //
  std::cout << "\nClosed-shell coupled cluster:\n";
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    auto ext_idx = ext_idx_list(i);
    cc_st_r[i] = sequant::closed_shell_CC_spintrace(cc_r[i]);
    canonicalize(cc_st_r[i]);

    // Remove S operator
    for (auto& term : *cc_st_r[i]) {
      if (term->is<Product>()) term = remove_tensor(term->as<Product>(), L"S");
    }

    // Biorthogonal transformation
    cc_st_r[i] = biorthogonal_transform(cc_st_r[i], i, ext_idx);

    // The symmetrizer operator is required for canonicalizer to give the
    // correct result
    if (i != 1) cc_st_r[i] = symmetrize_expr(cc_st_r[i], ext_idx);
    simplify(cc_st_r[i]);

    // Remove S operator
    //    for (auto& term : *cc_st_r[i]) {
    //      if (term->is<Product>()) term = remove_tensor(term->as<Product>(),
    //      L"S");
    //    }

    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu time: %5.3f sec.\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }

  if (NMAX == 4) {
    runtime_assert(cc_st_r.size() == 5)
        runtime_assert(cc_st_r.at(1)->size() == 30)    // T1
        runtime_assert(cc_st_r.at(2)->size() == 78)    // T2
        runtime_assert(cc_st_r.at(3)->size() == 567)   // T3
        runtime_assert(cc_st_r.at(4)->size() == 2150)  // T4
  } else if (NMAX == 3) {
    runtime_assert(cc_st_r.size() == 4)
        runtime_assert(cc_st_r.at(1)->size() == 30)   // T1
        runtime_assert(cc_st_r.at(2)->size() == 73)   // T2
        runtime_assert(cc_st_r.at(3)->size() == 490)  // T3
  }
}

ExprPtr biorthogonal_transform(
    const sequant::ExprPtr& expr, const int n_particles,
    const std::vector<std::vector<sequant::Index>>& ext_index_groups,
    const double threshold) {
  assert(n_particles != 0);
  assert(!ext_index_groups.empty());

  using sequant::container::svector;

  // Coefficients
  std::vector<double> bt_coeff_vec;
  {
    auto factorial = [](auto x) {
      if (x > 20)
        throw std::invalid_argument(
            "20! cannot be represented by std::intmax_t");
      return boost::numeric_cast<std::intmax_t>(sequant::factorial(x));
    };

    using namespace Eigen;
    // Dimension of permutation matrix is n_particles!
    const auto n = factorial(n_particles);

    // Permutation matrix
    Eigen::Matrix<double, Dynamic, Dynamic> M(n, n);
    {
      M.setZero();
      size_t n_row = 0;
      svector<int, 6> v(n_particles), v1(n_particles);
      std::iota(v.begin(), v.end(), 0);
      std::iota(v1.begin(), v1.end(), 0);
      do {
        std::vector<double> permutation_vector;
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
        return abs(d) > threshold;
      };

      // Solve system of equations
      SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
      std::vector<double> eig_vals(eig_solver.eigenvalues().size());
      VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
          eig_solver.eigenvalues();

      double non0count =
          std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
      scalar = eig_vals.size() / non0count;
    }

    // Find Pseudo Inverse, get 1st row only
    MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
    bt_coeff_vec.resize(pinv.rows());
    VectorXd::Map(&bt_coeff_vec[0], bt_coeff_vec.size()) = pinv.row(0) * scalar;
  }

  // Transformation maps
  std::vector<std::map<Index, Index>> bt_maps;
  {
    std::vector<Index> idx_list(ext_index_groups.size());

    for (auto i = 0; i != ext_index_groups.size(); ++i) {
      idx_list[i] = *ext_index_groups[i].begin();
    }

    const std::vector<Index> const_idx_list = idx_list;

    do {
      std::map<Index, Index> map;
      auto const_list_ptr = const_idx_list.begin();
      for (auto& i : idx_list) {
        map.emplace(std::make_pair(*const_list_ptr, i));
        const_list_ptr++;
      }
      bt_maps.push_back(map);
    } while (std::next_permutation(idx_list.begin(), idx_list.end()));
  }

  // If this assertion fails, change the threshold parameter
  assert(bt_coeff_vec.size() == bt_maps.size());

  // Checks if the replacement map is a canonical sequence
  auto is_canonical = [](const std::map<Index, Index>& idx_map) {
    bool canonical = true;
    for (auto&& pair : idx_map)
      if (pair.first != pair.second) return false;
    return canonical;
  };

  // Scale transformed expressions and append
  Sum bt_expr{};
  auto coeff_it = bt_coeff_vec.begin();
  for (auto&& map : bt_maps) {
    const auto v = to_rational(*coeff_it, threshold);
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

// Generate S operator from external index list
ExprPtr symmetrize_expr(
    ExprPtr& expr,
    const container::vector<container::vector<Index>>& ext_index_groups) {
  container::vector<Index> bra_list, ket_list;
  for (auto&& idx_group : ext_index_groups) {
    bra_list.push_back(*idx_group.begin());
    ket_list.push_back(*(idx_group.begin() + 1));
  }

  assert(bra_list.size() == ket_list.size());
  auto S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  return ex<Tensor>(S) * expr;
}
