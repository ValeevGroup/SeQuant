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

  using sequant::eqs::compute_all;

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
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;
#if !CLOSED_SHELL_SPINTRACE
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
#else
  auto cc_r = sequant::eqs::cceqvec{3, 3}(true, true, true, true, true);

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


#define PERMUTATION_TEST 1

#if PERMUTATION_TEST
  auto P13_b = ex<Tensor>(L"P", WstrList{L"i_98", L"i_99"},
                          WstrList{L"a_1", L"a_3"},
                          Symmetry::nonsymm);
  auto P13_k = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                          WstrList{L"a_98", L"a_99"},
                          Symmetry::nonsymm);
  auto P12_b = ex<Tensor>(L"P", WstrList{L"i_98", L"i_99"},
                          WstrList{L"a_1", L"a_2"},
                          Symmetry::nonsymm);
  auto P12_k = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                          WstrList{L"a_98", L"a_99"},
                          Symmetry::nonsymm);

  auto P23_b = ex<Tensor>(L"P", WstrList{L"i_98", L"i_99"},
                          WstrList{L"a_2", L"a_3"},
                          Symmetry::nonsymm);
  auto P23_k = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                          WstrList{L"a_98", L"a_99"},
                          Symmetry::nonsymm);

  auto P4_1313 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                            WstrList{L"a_1", L"a_3"},
                            Symmetry::nonsymm);
  auto P4_1323 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                            WstrList{L"a_2", L"a_3"},
                            Symmetry::nonsymm);
  auto P4_2313 = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                            WstrList{L"a_1", L"a_3"},
                            Symmetry::nonsymm);
  auto P4_2323 = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                            WstrList{L"a_2", L"a_3"},
                            Symmetry::nonsymm);

  auto P4_1212 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                            WstrList{L"a_1", L"a_2"},
                            Symmetry::nonsymm);
  auto P4_1213 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                            WstrList{L"a_1", L"a_3"},
                            Symmetry::nonsymm);
  auto P4_1312 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                            WstrList{L"a_1", L"a_2"},
                            Symmetry::nonsymm);

  auto p_aab = ex<Constant>(1.) - P13_b - P23_b - P13_k - P23_k +
      P4_1313 + P4_1323 + P4_2313 + P4_2323;

  auto p_abb = ex<Constant>(1.) - P13_b - P12_b - P13_k - P12_k +
      P4_1212 + P4_1213 + P4_1312 + P4_1313;

  auto A_12 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
  auto A_23 = ex<Tensor>(L"A", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"}, Symmetry::antisymm);
  auto A3 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm);


  auto occA = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha);
  auto virA = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::alpha);
  auto occB = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta);
  auto virB = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::beta);
  const auto i1A = Index(L"i⁺_1", occA);
  const auto i2A = Index(L"i⁺_2", occA);
  const auto a1A = Index(L"a⁺_1", virA);
  const auto a2A = Index(L"a⁺_2", virA);
  const auto a2B = Index(L"a⁻_2", virB);
  const auto a3B = Index(L"a⁻_3", virB);
  const auto i2B = Index(L"i⁻_2", occB);
  const auto i3B = Index(L"i⁻_3", occB);
  auto A2_aa = Tensor(L"A", {i1A, i2A}, {a1A, a2A}, Symmetry::antisymm);
  auto A2_bb = Tensor(L"A", {i2B, i3B}, {a2B, a3B}, Symmetry::antisymm);

  size_t n_aab = 0;
  size_t n_abb = 0;
  for (auto& product_term : *cc_r[3]) {
    if(product_term->is<Product>()){
      auto input = remove_tensor_from_product(product_term->as<Product>(), L"A");
      std::wcout << "input: " << to_latex(input) << std::endl;
      auto aab_input = p_aab * input;
      expand(aab_input);
      aab_input = expand_P_operator(aab_input, false);
      aab_input->visit(reset_idx_tags);
      simplify(aab_input);

      auto aab_result = open_shell_spintrace(aab_input,ext_idx_list(3));
      expand(aab_result[1]);
      aab_result[1] = expand_A_operator(aab_result[1]);
      aab_result[1]->visit(reset_idx_tags);
      simplify(aab_result[1]);
      n_aab += aab_result[1]->size();
      // std::wcout << "aab_result:\n" << to_latex(aab_result[1]) << "\n";

      auto aab_final = ex<Tensor>(A2_aa) * aab_result[1];
      // std::wcout << "aab_final:\n" << to_latex(aab_final) << "\n";
      expand(aab_final);
      aab_final = expand_A_operator(aab_final);
      aab_final->visit(reset_idx_tags);
      simplify(aab_final);
      // std::wcout << "aab_final:\n" << to_latex(aab_final) << "\n";

      auto abb_input = p_abb * input;
      expand(abb_input);
      abb_input = expand_P_operator(abb_input, false);
      abb_input->visit(reset_idx_tags);
      simplify(abb_input);

      auto abb_result = open_shell_spintrace(abb_input,ext_idx_list(3));
      expand(abb_result[2]);
      abb_result[2] = expand_A_operator(abb_result[2]);
      abb_result[2]->visit(reset_idx_tags);
      simplify(abb_result[2]);
      n_abb += abb_result[2]->size();
      // std::wcout << "abb_result:\n" << to_latex(abb_result[2]) << "\n";

      auto abb_final = ex<Tensor>(A2_bb) * abb_result[2];
      // std::wcout << "abb_final:\n" << to_latex(abb_final) << "\n";
      expand(abb_final);
      abb_final = expand_A_operator(abb_final);
      abb_final->visit(reset_idx_tags);
      simplify(abb_final);
      // std::wcout << "abb_final:\n" << to_latex(abb_final) << "\n";

//       auto st_with_A = open_shell_spintrace(product_term, ext_idx_list(3));

//      if(to_latex(aab_final) != to_latex(st_with_A[1])){
//        std::wcout << "aab_final:    " << to_latex(aab_final) << "\n";
//        std::wcout << "st_with_A[1]: " << to_latex(st_with_A[1]) << "\n\n";
//      }
//      if(to_latex(abb_final) != to_latex(st_with_A[2])){
//        std::wcout << "abb_final:    " << to_latex(abb_final) << "\n";
//        std::wcout << "st_with_A[2]: " << to_latex(st_with_A[2]) << "\n\n";
//      }

    }

  }
  std::cout << n_aab << " " << n_abb << std::endl;
  return 0;
#endif


#define PRINT_PRODUCT_TERMS 0
#if PRINT_PRODUCT_TERMS
  for (int i = 3; i < cc_r.size(); ++i){
    size_t term_counter = 0;
    for(auto term : *cc_r[i]) {
      std::wcout /*<< term_counter << ": " */ << "$$" << to_latex(term) << "$$ \n";
      ++term_counter;
    }
    std::cout << "\n";
  }
  // return 0;
#endif


#define OPEN_SHELL_TEST 0
#if OPEN_SHELL_TEST
  auto occA = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha);
  auto virA = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::alpha);
  auto occB = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta);
  auto virB = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::beta);
  const auto i1A = Index(L"i⁺_1", occA);
  const auto i2A = Index(L"i⁺_2", occA);
  const auto i3A = Index(L"i⁺_3", occA);
  const auto a1A = Index(L"a⁺_1", virA);
  const auto a2A = Index(L"a⁺_2", virA);
  const auto a3A = Index(L"a⁺_3", virA);

  const auto i1B = Index(L"i⁻_1", occB);
  const auto i2B = Index(L"i⁻_2", occB);
  const auto i3B = Index(L"i⁻_3", occB);
  const auto a1B = Index(L"a⁻_1", virB);
  const auto a2B = Index(L"a⁻_2", virB);
  const auto a3B = Index(L"a⁻_3", virB);

  auto A2_alpha = Tensor(L"A", {i1A, i2A}, {a1A, a2A}, Symmetry::antisymm);
  auto A2_beta = Tensor(L"A", {i1B, i2B}, {a1B, a2B}, Symmetry::antisymm);

  auto A3_alpha = Tensor(L"A", {i1A, i2A, i3A}, {a1A, a2A, a3A}, Symmetry::antisymm);
  auto A3_beta = Tensor(L"A", {i1B, i2B, i3B}, {a1B, a2B, a3B}, Symmetry::antisymm);

  auto A3_aab = Tensor(L"A", {i1A, i2A, i3B}, {a1A, a2A, a3B}, Symmetry::antisymm);
  auto A3_abb = Tensor(L"A", {i1A, i2B, i3B}, {a1A, a2B, a3B}, Symmetry::antisymm);


  for (int i = 3; i < 4 /*cc_r.size()*/; ++i) {
  size_t counter = 0;
  std::vector<size_t> n_st_terms(i+1,0);
  std::vector<size_t> n_st_with_A_terms(i+1,0);
  const auto list = ext_idx_list(i);

  for (auto& product_term : *cc_r[i]) {
    ExprPtr term;
    if(product_term->is<Product>())
      term = remove_tensor_from_product(product_term->as<Product>(), L"A");

    // auto term = product_term;
    // std::wcout << "Input " << counter << ": " << to_latex(term) << std::endl;

    auto os_st = open_shell_spintrace(term, list);
    auto st_with_A = open_shell_spintrace(product_term, list);

    // for (size_t j = 0; j != os_st.size(); ++j) {
    for (size_t j : {0,3} /* all same spin blocks only*/) {
      // std::wcout << __LINE__ << " " << os_st[j]->size() << ": " << to_latex(os_st[j]) << "\n";

      ExprPtr new_expr{};
      if (j==0)
        new_expr = ex<Tensor>(A3_alpha) * os_st[j];
      else if (j==1)
        new_expr = ex<Tensor>(A3_aab) * os_st[j];
      else if (j==2)
        new_expr = ex<Tensor>(A3_abb) * os_st[j];
      else if (j==3)
        new_expr = ex<Tensor>(A3_beta) * os_st[j];

      expand(new_expr);
      new_expr = expand_A_operator(new_expr);
      new_expr->visit(reset_idx_tags);
      canonicalize(new_expr);
      rapid_simplify(new_expr);

      assert(new_expr->size() == st_with_A[j]->size());
      assert(to_latex(new_expr) == to_latex(st_with_A[j]));
      n_st_terms[j] += os_st[j]->size();
      n_st_with_A_terms[j] += st_with_A[j]->size();
    }
      ++counter;
  }
  std::cout << "Full expansion: " << n_st_with_A_terms[0] << " " << n_st_with_A_terms[3] << "\n";
  std::cout << "Preserve A2: " << n_st_terms[0] << " " << n_st_terms[3] << "\n";

  return 0;

  std::cout << "CC R" << i << ": ";
  size_t n_terms_R = 0;
  for (auto& n_terms : n_st_terms) {
    n_terms_R += n_terms;
    std::cout << n_terms << " ";
  }
  std::cout << ": " << n_terms_R << std::endl;
}
  return 0;
#endif

#if 0 // Open-shell
  for (int i = 1; i < cc_r.size(); ++i) {
    const auto list = ext_idx_list(i);
    auto temp = open_shell_spintrace(cc_r[i], list);
    std::cout << "R" << i << ": ";
    for(auto& t : temp){
      std::cout << t->size() << " ";
      // std::wcout << to_latex(t) << "\n";
    }
    std::cout << "\n";
  }
  return 0;
#endif

  // Closed-shell coupled cluster spin-trace
  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (int i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    const auto list = ext_idx_list(i);
    cc_st_r[i] = closed_shell_spintrace(cc_r[i], list);
    // cc_st_r[i] = closed_shell_cc_spintrace(cc_r[i]);
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
    simplify(cc_st_r[i]);

    tstop = std::chrono::high_resolution_clock::now();
    time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu simplify time: %5.6f sec.\n\n", i, cc_st_r[i]->size(), time_elapsed.count());
  }
  std::wcout << "CCSDT: " << to_latex(*cc_st_r[3]) << std::endl;
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

