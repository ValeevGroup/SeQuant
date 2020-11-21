#include "../sequant_setup.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace sequant;

#define CLOSED_SHELL_SPINTRACE 1
#if CLOSED_SHELL_SPINTRACE
container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold);
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const container::vector<container::vector<Index>> ext_index_groups);
ExprPtr symmetrize_expr(ExprPtr& expr, const container::vector<container::vector<Index>> ext_index_groups = {{}});
#endif

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
  sequant::detail::OpIdRegistrar op_id_registrar;

  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // set_num_threads(1);

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 3;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  // change to true to print out the resulting equations
  constexpr bool print = true;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

    ranges::for_each(std::array<bool, 2>{false, true}, [=](const bool screen) {
      ranges::for_each(
          std::array<bool, 2>{false, true}, [=](const bool use_topology) {
            ranges::for_each(std::array<bool, 2>{false, true},
                             [=](const bool canonical_only) {
                               tpool.clear();
                               // comment out to run all possible combinations
                               if (screen && use_topology && canonical_only)
                                 compute_all{NMAX}(print, screen, use_topology,
                                                   true, canonical_only);
                             });
          });
    });

#if CLOSED_SHELL_SPINTRACE
  auto cc_r = cceqvec{ 4, 4}(true, true, true, true, true);

  /// Make external index
  auto ext_idx_list = [] (const int i_max){
    container::vector<container::vector<Index>> ext_idx_list;

    for(size_t i = 1; i <= i_max; ++i) {
      auto label = std::to_wstring(i);
      auto occ_i = Index::make_label_index(IndexSpace::instance(IndexSpace::active_occupied), label);
      auto virt_i = Index::make_label_index(IndexSpace::instance(IndexSpace::active_unoccupied), label);
      container::vector<Index> pair = {occ_i, virt_i};
      ext_idx_list.push_back(pair);
    }
    return ext_idx_list;
  };

  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (int i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    const auto list = ext_idx_list(i);
    cc_st_r[i] = closed_shell_spintrace(cc_r[i], list);
    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu spintrace: %5.6f sec.\n", i, cc_st_r[i]->size(), time_elapsed.count());
    canonicalize(cc_st_r[i]);
    tstop = std::chrono::high_resolution_clock::now();
    time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu canonicalize: %5.6f sec.\n", i, cc_st_r[i]->size(), time_elapsed.count());

    // Remove S operator
    for(auto&& term: *cc_st_r[i]){
      if(term->is<Product>()) term = remove_tensor_from_product(term->as<Product>(), L"S");
    }

    // Checks if the replacement map is a canonical sequence
    auto is_canonical = [] (const std::map<Index, Index>& idx_map){
      bool canonical = true;
      for(auto&& pair: idx_map) if(pair.first != pair.second) return false;
      return canonical;
    };

    // Get coefficients and replacement maps
    auto btc = biorthogonal_tran_coeff(list.size(), 1.e-12);

    auto idx_map = biorthogonal_tran_idx_map(list);
    assert(btc.size() == idx_map.size());

    // Append scale and append all transformations
    Sum bt_expr{};
    auto btc_ptr = btc.begin();
    for(auto&& map: idx_map){
      ExprPtr transformed_expr{};
      if(is_canonical(map))
        transformed_expr = ex<Constant>(*btc_ptr) * cc_st_r[i];
      else
        transformed_expr = ex<Constant>(*btc_ptr) * transform_expression(cc_st_r[i], map);
      btc_ptr++;
      bt_expr.append(transformed_expr);
    }
    cc_st_r[i] = std::make_shared<Sum>(bt_expr);

    tstop = std::chrono::high_resolution_clock::now();
    time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu biorthogonal transform time: %5.6f sec.\n", i, cc_st_r[i]->size(), time_elapsed.count());

    if(i != 1)
      cc_st_r[i] = symmetrize_expr(cc_st_r[i], list);

    tstop = std::chrono::high_resolution_clock::now();
    time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu S time: %5.6f sec.\n", i, cc_st_r[i]->size(), time_elapsed.count());

    expand(cc_st_r[i]);
    canonicalize(cc_st_r[i]);
    rapid_simplify(cc_st_r[i]);

    tstop = std::chrono::high_resolution_clock::now();
    time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu simplify time: %5.6f sec.\n\n", i, cc_st_r[i]->size(), time_elapsed.count());
  }
#endif

  return 0;
}

#if CLOSED_SHELL_SPINTRACE
// Generate S operator from external index list
ExprPtr symmetrize_expr(ExprPtr& expr, const container::vector<container::vector<Index>> ext_index_groups){

  container::vector<Index> bra_list, ket_list;
  for(auto&& idx_group : ext_index_groups) {
    bra_list.push_back(*idx_group.begin());
    ket_list.push_back(*(idx_group.begin() + 1));
  }

  assert(bra_list.size() == ket_list.size());
  auto S = Tensor(L"S", bra_list, ket_list, Symmetry::nonsymm);
  return ex<Tensor>(S) * expr;
}

container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold){
  using namespace Eigen;

  int n = std::tgamma(n_particles + 1); // <- Dimension of permutation matrix is n_particles!
  // Permutation matrix
  Eigen::MatrixXd M(n,n);
  {
    M.setZero();
    size_t n_row = 0;
    container::svector<int, 6> v(n_particles), v1(n_particles);
    std::iota(v.begin(), v.end(), 0);
    std::iota(v1.begin(), v1.end(), 0);
    do {
      container::vector<double> permutation_vector;
      do {
        auto cycles = count_cycles(v1, v);
        permutation_vector.push_back(std::pow(-2, cycles));
      } while (std::next_permutation(v.begin(), v.end()));
      Eigen::VectorXd pv_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(permutation_vector.data(),
                                                                             permutation_vector.size());
      M.row(n_row) = pv_eig;
      ++n_row;
    } while (std::next_permutation(v1.begin(), v1.end()));
    M *= std::pow(-1, n_particles);
    // std::cout << "permutation_matrix:\n" << M << std::endl;
  }

  // Normalization constant
  double scalar;
  {
    // inline bool nonZero(double d) { return abs(d) > threshold ? true : false; }
    auto nonZero = [&threshold] (const double d) { return abs(d) > threshold ? true : false; };

    // Solve system of equations
    SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
    container::vector<double> eig_vals(eig_solver.eigenvalues().size());
    VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
        eig_solver.eigenvalues();

    double non0count = std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
    scalar = eig_vals.size() / non0count;
  }
  // std::cout << "scalar: " << scalar << std::endl;
  // Find Pseudo Inverse, get 1st row only
  MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
  container::vector<double> result(pinv.rows());
  VectorXd::Map(&result[0], result.size()) = pinv.row(0) * scalar;
  return result;
}

/// @brief Biorthogonal transformation map
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const container::vector<container::vector<Index>> ext_index_groups){
  //Check size of external index group; catch exception otherwise
  if(ext_index_groups.size() == 0) throw( "Cannot compute index map since " && "ext_index_groups.size() == 0");
  assert(ext_index_groups.size() > 0);

  std::vector<Index> idx_list;
  for(auto&& idx_group : ext_index_groups) idx_list.push_back(*idx_group.begin());

  const container::vector<Index> const_idx_list = idx_list;
  // Do permutations and append to map
  std::vector<std::map<Index, Index>> result;
  do{
    std::map<Index, Index> map;
    auto const_list_ptr = const_idx_list.begin();
    for(auto&& i : idx_list){
      map.emplace(std::make_pair(*const_list_ptr, i));
      const_list_ptr++;
    }
    result.push_back(map);
  } while(std::next_permutation(idx_list.begin(), idx_list.end()));
  return result;
}
#endif

