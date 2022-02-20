#include <SeQuant/core/op.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/eqs/cceqs.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <clocale>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace sequant;

container::vector<double> biorthogonal_tran_coeff(const int n_particles, const double& threshold);
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(const container::vector<container::vector<Index>> ext_index_groups);
ExprPtr symmetrize_expr(ExprPtr& expr, const container::vector<container::vector<Index>> ext_index_groups = {{}});

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

  using sequant::eqs::compute_all;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // set_num_threads(1);

#ifndef NDEBUG
  const size_t DEFAULT_NMAX = 4;
#else
  const size_t DEFAULT_NMAX = 4;
#endif
  const size_t NMAX = argc > 1 ? std::atoi(argv[1]) : DEFAULT_NMAX;
  // change to true to print out the resulting equations
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  ranges::for_each(std::array<bool, 2>{false, true}, [=](const bool screen) {
    ranges::for_each(
        std::array<bool, 2>{false, true}, [=](const bool use_topology) {
          ranges::for_each(std::array<bool, 2>{false, true},
                           [=](const bool canonical_only) {
                             // TODO tpool was in scope before
                             // separting header and source files
                             // tpool.clear();
                             // comment out to run all possible combinations
                             if (screen && use_topology && canonical_only)
                               compute_all{NMAX}(print, screen, use_topology,
                                                 true, canonical_only);
                           });
        });
  });

  // reset index tags
  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };

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
  auto cc_r = sequant::eqs::cceqvec{3, 3}(true, true, true, true, true);

  auto A3 = ex<Tensor>(Tensor(L"A",
                              WstrList{L"i_1", L"i_2", L"i_3"},
                              WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm));
  // Bra or ket for P operators don't matter
  auto P_ab = ex<Tensor>(Tensor(L"P",WstrList{L"a_1", L"a_2"},{}));
  auto P_ac = ex<Tensor>(Tensor(L"P",WstrList{L"a_1", L"a_3"},{}));
  auto P_bc = ex<Tensor>(Tensor(L"P",WstrList{L"a_2", L"a_3"},{}));
  auto P_ij = ex<Tensor>(Tensor(L"P",WstrList{L"i_1", L"i_2"},{}));
  auto P_ik = ex<Tensor>(Tensor(L"P",WstrList{L"i_1", L"i_3"},{}));
  auto P_jk = ex<Tensor>(Tensor(L"P",WstrList{L"i_2", L"i_3"},{}));

  // G and t are kept consistent with the SeQuant CC notation
  auto G_oovv = ex<Tensor>(Tensor(L"g",
                                  WstrList{L"i_4", L"i_5"}, WstrList{L"a_4", L"a_5"},
                                  Symmetry::antisymm));

  auto t2_a1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_4", L"a_1"},
                                 WstrList{L"i_4", L"i_1"}, Symmetry::antisymm));
  auto t3_a1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_2", L"a_3"},
                                 WstrList{L"i_5", L"i_2", L"i_3"}, Symmetry::antisymm));

  auto t2_a2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_2"},
                                 WstrList{L"i_1", L"i_4"}, Symmetry::antisymm));
  auto t3_a2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_4", L"a_5", L"a_3"},
                                 WstrList{L"i_2", L"i_5", L"i_3"}, Symmetry::antisymm));

  auto t2_b1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_4", L"a_5"},
                                 WstrList{L"i_1", L"i_2"}, Symmetry::antisymm));
  auto t3_b1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_2", L"a_3"},
                                 WstrList{L"i_4", L"i_5", L"i_3"}, Symmetry::antisymm));

  auto t2_b2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_2"},
                                 WstrList{L"i_4", L"i_5"}, Symmetry::antisymm));
  auto t3_b2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_4", L"a_5", L"a_3"},
                                 WstrList{L"i_2", L"i_5", L"i_3"}, Symmetry::antisymm));

  auto t2_c1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_4", L"a_1"},
                                 WstrList{L"i_4", L"i_5"}, Symmetry::antisymm));
  auto t3_c1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_5", L"a_2", L"a_3"},
                                 WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm));

  auto t2_c2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_4"},
                                 WstrList{L"i_1", L"i_2"}, Symmetry::antisymm));
  auto t3_c2 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_2", L"a_5", L"a_3"},
                                 WstrList{L"i_4", L"i_5", L"i_3"}, Symmetry::antisymm));

  auto t2_d1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_1", L"a_4"},
                                 WstrList{L"i_1", L"i_4"}, Symmetry::antisymm));
  auto t3_d1 = ex<Tensor>(Tensor(L"t",
                                 WstrList{L"a_5", L"a_2", L"a_3"},
                                 WstrList{L"i_5", L"i_2", L"i_3"}, Symmetry::antisymm));

  // Terms
  auto a1 = ex<Constant> (-0.5) * (ex<Constant>(1.0) - P_ij - P_ik) *
           G_oovv * t2_a1 * t3_a1;
  auto a2 = ex<Constant> (-0.5) * (ex<Constant>(1.0) - P_ij - P_ik) *
           (ex<Constant>(1.0) - P_ac - P_bc) *
           G_oovv * t2_a2 * t3_a2;
  auto b1 = ex<Constant> (0.25) * (ex<Constant>(1.0) - P_ij - P_ik) *
           G_oovv * t2_b1 * t3_b1;
  auto b2 = ex<Constant> (0.25) * (ex<Constant>(1.0) - P_bc - P_ac) *
           G_oovv * t2_b2 * t3_b2;
  auto c1 = ex<Constant> (-0.5) * (ex<Constant>(1.0) - P_ab - P_ac) *
            G_oovv * t2_c1 * t3_c1;
  auto c2 = ex<Constant> (-0.5) * (ex<Constant>(1.0) - P_ik - P_jk)  *
            (ex<Constant>(1.0) - P_ab - P_ac)  *
             G_oovv * t2_c2 * t3_c2;
  auto d1 = (ex<Constant>(1.0) - P_ij - P_ik)  *
            (ex<Constant>(1.0) - P_ab - P_ac)  *
            G_oovv * t2_d1 * t3_d1;

  for(auto &t : {a1,a2,b1,b2,c1,c2,d1}){
    std::wcout << to_latex(t) << std::endl;
  }

  // Full term that needs to be subtracted from SeQuant derived CCSDT R3
  auto w_t2_t3 = a1 + a2 + b1 + b2 + c1 + c2 + d1;

  // p_CCSDT correction
  auto pCCSDT_correction = [&](const std::vector<int>& params){
    assert(params.size() == 3);
    auto correction = ex<Constant>(0.5) * a1 +
                      ex<Constant>(0.5) * a2 +
                      ex<Constant>(params[0]) * (ex<Constant>(0.5) * a1 + b1) +
                      ex<Constant>(params[1]) * (ex<Constant>(0.5) * a2 + b2) +
                      ex<Constant>(params[2]) * (c1 + c2 + d1);
    simplify(correction);
    return correction;
  };

  // Operations to simplify the correction term
#if 0
  expand(c1);
  c1 = expand_P_op(c1);
  c1 = A3 * c1;
  expand(c1);
  canonicalize(c1);
  c1->visit(reset_idx_tags);
  c1 = remove_tensor(c1, L"A");
#endif

  // pCCSDT R3 expression
  auto r3 = cc_r[3] - w_t2_t3 + pCCSDT_correction({1, 1, 0});


#if 0
  //
  // Closed-shell spintrace (fast)
  //
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (size_t i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    auto ext_idx = ext_idx_list(i);
    cc_st_r[i] = sequant::closed_shell_CC_spintrace(cc_r[i]);
    canonicalize(cc_st_r[i]);

    // Remove S operator
    for (auto& term : *cc_st_r[i]) {
      if (term->is<Product>())
        term = sequant::remove_tensor_from_product(term->as<Product>(), L"S");
    }

    // Biorthogonal transformation
    cc_st_r[i] = biorthogonal_transform(cc_st_r[i], i, ext_idx);

    // The symmetrizer operator is required for canonicalizer to give the
    // correct result
    if (i != 1) cc_st_r[i] = symmetrize_expr(cc_st_r[i], ext_idx);
    simplify(cc_st_r[i]);

    // Remove S operator
    for (auto& term : *cc_st_r[i]) {
      if (term->is<Product>())
        term = remove_tensor_from_product(term->as<Product>(), L"S");
    }

    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%lu size: %lu time: %5.3f sec.\n", i,
                                  cc_st_r[i]->size(), time_elapsed.count());
  }

  //
  // Open-shell spintrace
  //
  std::cout << "Open-shell coupled cluster: nterms per spin blocks: " << std::endl;
  std::vector<std::vector<ExprPtr>> os_cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    Tensor A = cc_r[i]->as<Sum>().summand(0)->as<Product>().factors()[0]->as<Tensor>();
    assert(A.label() == L"A");
    auto P_vec = open_shell_P_op_vector(A);
    auto A_vec = open_shell_A_op(A);
    assert(P_vec.size() == i+1);

    std::vector<Sum> concat_terms(i + 1);
    size_t n_spin_orbital_term = 0;
    for (auto& product_term : *cc_r[i]) {
      auto term = remove_tensor(product_term->as<Product>(), L"A");
      std::vector<ExprPtr> os_st(i+1);

      // Apply the P operators on the product term without the A,
      // Expand the P operators and spin-trace the expression
      // Then apply A operator, canonicalize and remove A operator
      for (int s = 0; s != os_st.size(); ++s) {
        os_st.at(s) = P_vec.at(s) * term;
        expand(os_st.at(s));
        os_st.at(s) = expand_P_op(os_st.at(s), false, true);
        os_st.at(s) = open_shell_spintrace(os_st.at(s),
                                           ext_idx_list(i), s).at(0);
        if (i > 2) {
          os_st.at(s) = A_vec.at(s) * os_st.at(s);
          simplify(os_st.at(s));
          os_st.at(s) = remove_tensor(os_st.at(s), L"A");
        }
      }

      for (size_t j = 0; j != os_st.size(); ++j) {
        concat_terms.at(j).append(os_st.at(j));
      }
      ++n_spin_orbital_term;
    }

    // Combine spin-traced terms for the current residual
    std::vector<ExprPtr> expr_vec;
    std::cout << "CC R" << i << ": ";
    for (auto& spin_case : concat_terms) {
      auto ptr = sequant::ex<Sum>(spin_case);
      expr_vec.push_back(ptr);
      std::cout << ptr->size() << " ";
    }

    os_cc_st_r.at(i) = std::move(expr_vec);
    std::cout << "\n";
  }

#define runtime_assert(tf)                                         \
  if (!(tf)) {                                                     \
    std::ostringstream oss;                                        \
    oss << "failed assert at line " << __LINE__ << " in SRCC example"; \
    throw std::runtime_error(oss.str().c_str());                   \
  }

  if (NMAX == 4) {
    runtime_assert(os_cc_st_r.size() == 5);
    runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30);
    runtime_assert(os_cc_st_r.at(2).at(1)->size() == 130);
    runtime_assert(os_cc_st_r.at(2).at(2)->size() == 74);
    runtime_assert(os_cc_st_r.at(3).at(1)->size() == 249);
    runtime_assert(os_cc_st_r.at(3).at(3)->size() == 124);
    runtime_assert(os_cc_st_r.at(4).at(1)->size() == 356);
    runtime_assert(os_cc_st_r.at(4).at(2)->size() == 386);
    runtime_assert(os_cc_st_r.at(4).at(4)->size() == 156);
  } else if (NMAX == 3) {
    runtime_assert(os_cc_st_r.size() == 4);
    runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30);
    runtime_assert(os_cc_st_r.at(2).at(0)->size() == 65);
    runtime_assert(os_cc_st_r.at(2).at(1)->size() == 122);
    runtime_assert(os_cc_st_r.at(3).at(2)->size() == 209);
    runtime_assert(os_cc_st_r.at(3).at(3)->size() == 75);
  }
#endif
}

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

