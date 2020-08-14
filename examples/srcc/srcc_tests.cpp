//
// Created by Nakul Teke on 8/4/20.
//

#include "../sequant_setup.hpp"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "SeQuant/domain/mbpt/spin.hpp"

using namespace sequant;

/// @brief Find coefficients for biorthogonal transformation
/// @detailed Given the number of external indices, this function calculates the
/// permutation matrix, counts it's eigenvalues to find normalization constant
/// and calculates a pseudoinverse to find the coefficients corresponding to the
/// n! permutations. For details, see: http://arxiv.org/abs/1805.00565
/// @param n_particles Number of external index group
/// @param threshold Cut-off for counting number of non-zero eigen values
/// @return A vector<double> of biorthogonal transformation coefficients
container::vector<double> biorthogonal_tran_coeff(const int n_particles,
                                                  const double& threshold) {
  using namespace Eigen;

  int n = std::tgamma(n_particles +
                      1);  // <- Dimension of permutation matrix is n_particles!
  // Permutation matrix
  Eigen::MatrixXd M(n, n);
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
      Eigen::VectorXd pv_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
          permutation_vector.data(), permutation_vector.size());
      M.row(n_row) = pv_eig;
      ++n_row;
    } while (std::next_permutation(v1.begin(), v1.end()));
    M *= std::pow(-1, n_particles);
    // std::cout << "permutation_matrix:\n" << M << "\n";
  }

  // Normalization constant
  double scalar;
  {
    // inline bool nonZero(double d) { return abs(d) > threshold ? true : false;
    // }
    auto nonZero = [&](const double d) {
      return abs(d) > threshold ? true : false;
    };

    // Solve system of equations
    SelfAdjointEigenSolver<MatrixXd> eig_solver(M);
    container::vector<double> eig_vals(eig_solver.eigenvalues().size());
    VectorXd::Map(&eig_vals[0], eig_solver.eigenvalues().size()) =
        eig_solver.eigenvalues();

    double non0count = std::count_if(eig_vals.begin(), eig_vals.end(), nonZero);
    scalar = eig_vals.size() / non0count;
  }

  // Find Pseudo Inverse, get 1st row only
  MatrixXd pinv = M.completeOrthogonalDecomposition().pseudoInverse();
  container::vector<double> result(pinv.rows());
  VectorXd::Map(&result[0], result.size()) = pinv.row(0) * scalar;
  return result;
}

/// @brief Biorthogonal transformation map
std::vector<std::map<Index, Index>> biorthogonal_tran_idx_map(
    const std::initializer_list<IndexList> ext_index_groups = {{}}) {
  // Check size of external index group; catch exception otherwise
  if (ext_index_groups.size() == 0)
    throw("Cannot compute index map since " && "ext_index_groups.size() == 0");
  assert(ext_index_groups.size() > 0);

  container::vector<Index> idx_list;
  for (auto&& idx_group : ext_index_groups)
    idx_list.push_back(*idx_group.begin());

  const container::vector<Index> const_idx_list = idx_list;
  // Do permutations and append to map
  std::vector<std::map<Index, Index>> result;
  do {
    std::map<Index, Index> map;
    auto const_list_ptr = const_idx_list.begin();
    for (auto&& i : idx_list) {
      map.emplace(std::make_pair(*const_list_ptr, i));
      const_list_ptr++;
    }
    result.push_back(map);
  } while (std::next_permutation(idx_list.begin(), idx_list.end()));

  return result;
}

/// @brief R2 equation coded manually
ExprPtr r2();

int main(int argc, char* argv[]) {
  // SEQUANT SETUP
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

  const int n = 2;
  const int p = 2;

  if (n == 2)
    std::cout << "CCSD";
  else if (n == 3)
    std::cout << "CCSDT";
  else if (n == 4)
    std::cout << "CCSDTQ";
  else if (n == 5)
    std::cout << "CCSDTQP";
  std::cout << " expressions:\n";

  auto cc_r = cceqvec{n, p}(true, true, true, true);

  std::vector<ExprPtr> cc_st_r(cc_r.size());
  for (int i = 1; i < cc_r.size(); ++i) {
    std::wcout << "Spin-Orbit R" << i << "\n" << to_latex(cc_r[i]) << "\n";
    const auto tstart = std::chrono::high_resolution_clock::now();
    std::initializer_list<IndexList> external_indices = {{}};
    if (i == 1)
      external_indices = {{L"i_1", L"a_1"}};
    else if (i == 2)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
    else if (i == 3)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
    else if (i == 4)
      external_indices = {{L"i_1", L"a_1"},
                          {L"i_2", L"a_2"},
                          {L"i_3", L"a_3"},
                          {L"i_4", L"a_4"}};

    // cc_st_r[i] = closed_shell_spintrace(cc_r[i], external_indices);
    cc_st_r[i] = spintrace(cc_r[i], external_indices);
    canonicalize(cc_st_r[i]);
    printf("R%d Spin-orbit: %lu terms;\nSPINTRACED: With S operator: %lu;", i,
           cc_r[i]->size(), cc_st_r[i]->size());
    cc_st_r[i] = expand_S_operator(cc_st_r[i]);
    rapid_simplify(cc_st_r[i]);
    canonicalize(cc_st_r[i]);
    printf("S expanded: %lu\n", cc_st_r[i]->size());

    // Checks if the replacement map is a canonical sequence
    auto is_canonical = [&](const std::map<Index, Index>& idx_map) {
      bool canonical = true;
      for (auto&& pair : idx_map)
        if (pair.first != pair.second) return false;
      return canonical;
    };

    // Get coefficients and replacement maps
#if 0
    auto btc = biorthogonal_tran_coeff(external_indices.size(), 1.e-12);
    // for(auto&& c: btc) std::cout << c << " ";
    // std::cout << "\n";
#else
    container::vector<double> btc;
    if(i== 1) btc = {0.5};
    else if(i == 2) btc = {1./3, 1./6};
    else if(i == 3) btc = {17./120, -1./120, -1./120, -7./120, -7./120, -1./120};
#endif
    // TODO: external_indices are getting modified before.
    if (i == 1)
      external_indices = {{L"i_1", L"a_1"}};
    else if (i == 2)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
    else if (i == 3)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
    else if (i == 4)
      external_indices = {{L"i_1", L"a_1"},
                          {L"i_2", L"a_2"},
                          {L"i_3", L"a_3"},
                          {L"i_4", L"a_4"}};

    auto idx_map = biorthogonal_tran_idx_map(external_indices);
    assert(btc.size() == idx_map.size());

    // Append scale and append all transformations
    Sum bt_expr{};
    auto btc_ptr = btc.begin();
    for (auto&& map : idx_map) {
      ExprPtr transformed_expr{};
      if (is_canonical(map))
        transformed_expr = ex<Constant>(*btc_ptr) * cc_st_r[i];
      else
        transformed_expr =
            ex<Constant>(*btc_ptr) * transform_expression(cc_st_r[i], map);
      btc_ptr++;
      bt_expr.append(transformed_expr);
    }
    ExprPtr bt_expr_p = std::make_shared<Sum>(bt_expr);
    expand(bt_expr_p);
    canonicalize(bt_expr_p);
    rapid_simplify(bt_expr_p);
    cc_st_r[i] = bt_expr_p;
    printf("Biorthogonal transform: %lu\n", cc_st_r[i]->size());

    if (i == 1)
      external_indices = {{L"i_1", L"a_1"}};
    else if (i == 2)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}};
    else if (i == 3)
      external_indices = {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}};
    else if (i == 4)
      external_indices = {{L"i_1", L"a_1"},
                          {L"i_2", L"a_2"},
                          {L"i_3", L"a_3"},
                          {L"i_4", L"a_4"}};

    cc_st_r[i] = factorize_S_operator(cc_st_r[i], external_indices, true);
    printf("Factorize S: %lu\n", cc_st_r[i]->size());
    canonicalize(cc_st_r[i]);
    rapid_simplify(cc_st_r[i]);
    const auto tstop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu time: %5.3f sec.\n\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }

  // Testing CCSD R2 expression
  auto r2_expr = r2();
  std::wcout << __LINE__ << " " << r2_expr->size() << " " << to_latex(r2_expr) << "\n";

#if 0
  for(auto&& term: *cc_st_r[2]){
    std::wcout << term->hash_value() << " " << term->to_latex() << "\n";
    for(auto&& tensor: term->as<Product>()) {
      if(tensor->as<Tensor>().bra_rank() == 1) {
        std::wcout << tensor->to_latex() << " "
                   << tensor->hash_value() << " "
                   << (tensor->as<Tensor>().symmetry() == Symmetry::nonsymm) << " "
                   << (tensor->as<Tensor>().braket_symmetry() == BraKetSymmetry::conjugate) << " "
                   << (tensor->as<Tensor>().particle_symmetry() == ParticleSymmetry::symm) << " "
                   << "\n";
      }
    }
  }
  std::cout << "\n";

  for(auto&& term: *r2_expr){
    std::wcout << term->hash_value() << " " << term->to_latex() << "\n";
    for(auto&& tensor: term->as<Product>()) {
      if(tensor->as<Tensor>().bra_rank() == 1) {
        std::wcout << tensor->to_latex() << " "
                   << tensor->hash_value() << " "
                   << (tensor->as<Tensor>().symmetry() == Symmetry::nonsymm) << " "
                   << (tensor->as<Tensor>().braket_symmetry() == BraKetSymmetry::conjugate) << " "
                   << (tensor->as<Tensor>().particle_symmetry() == ParticleSymmetry::symm) << " "
                   << "\n";
      }
    }
  }

  container::set<size_t> hash_list;
  for(auto&& term : *cc_st_r[2]) hash_list.insert(term->hash_value());
  // for(auto&& term : *r2_expr) hash_list.insert(term->hash_value());

  std::cout << "hash_list size: " << hash_list.size() << "\n";
  for(auto i = r2_expr->begin(); i != r2_expr->end(); ++i){
    auto h0 = (*i)->hash_value();
    auto it = std::find(hash_list.begin(), hash_list.end(), h0);
    if(it == hash_list.end()) std::wcout << (*i)->to_latex() << "\n";
  }
  std::cout << "\n";
//
   return 0;
#endif

  cc_st_r[2] = ex<Constant>(-1.0) * cc_st_r[2];
  expand(cc_st_r[2]);
  r2_expr->as<Sum>().append(cc_st_r[2]->clone());
  std::wcout << __LINE__ << " " << r2_expr->size() << " " << to_latex(r2_expr) << "\n";
  simplify(r2_expr);
  std::wcout << __LINE__ << " " << r2_expr->size() << " " << to_latex(r2_expr) << "\n";

  return 0;

/*
  auto result = std::make_shared<Sum>();
  result->append(ex<Constant>(0.0));
  for(auto i = r2_expr->begin(); i != r2_expr->end(); ++i){
    auto h0 = (*i)->hash_value();
    auto scalar = (*i)->as<Product>().scalar();
    bool match = false;
    if((*i)->is<Product>())
      for(auto j = i + 1; j != r2_expr->end(); ++j){
        if((*j)->hash_value() == h0){
          match = true;
          scalar += (*j)->as<Product>().scalar();
          *j = ex<Constant>(0.0);
        }
      }
    // Change scalar of (*i
    // std::wcout << scalar.real() << " " << (*i)->to_latex() << "\n";

    if(match){
      Product new_product{};
      new_product.scale(scalar);
      for(auto&& tensor: (*i)->as<Product>()) new_product.append(tensor);
      result->append(ex<Product>(new_product));
    } else {
      result->append(*i);
    }
  }
  ExprPtr result_ptr = result->clone();
  rapid_simplify(result_ptr);
  std::wcout << "\n" << __LINE__ << " " << result_ptr->size() << " " << to_latex(result_ptr) << "\n\n";


  auto t1 = *cc_st_r[2]->begin();
  auto T1 = *r2_expr->begin();
  auto g1 = t1->as<Product>().factor(1);
  auto G1 = T1->as<Product>().factor(1);
  std::wcout << g1->to_latex() << " " << G1->to_latex() << "\n";
  std::cout << (g1->as<Tensor>().braket_symmetry() == G1->as<Tensor>().braket_symmetry()) << "\n";
  std::cout << (g1->as<Tensor>().symmetry() == G1->as<Tensor>().symmetry()) << "\n";
  std::cout << (g1->as<Tensor>().particle_symmetry() == G1->as<Tensor>().particle_symmetry()) << "\n";

  {
    Sum test_sum{};
    test_sum.append(g1);
    test_sum -= g1->as<Expr>();
    std::wcout << __LINE__ << " " << to_latex(test_sum) << "\n";
    auto ptr = ex<Sum>(test_sum);
    expand(ptr);
    canonicalize(ptr);
    rapid_simplify(ptr);
    std::wcout << __LINE__ << " " << to_latex(ptr) << "\n";
  }

  auto g_diff = g1->clone() + ex<Constant>(-1.0) * g1->clone() + ex<Constant>(1.0);
  simplify(g_diff);
  std::wcout << __LINE__ << " " << g_diff->size() << " " << to_latex(g_diff) << "\n";

  auto g_diff2 = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm);
  simplify(g_diff2);
  std::wcout << __LINE__ << " " << to_latex(g_diff2) << "\n";

  g_diff2 = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) -
      ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm);
  simplify(g_diff2);
  std::wcout << __LINE__ << " " << to_latex(g_diff2) << "\n";

  {
    // Check hash values
    auto g = Tensor(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm);
    auto expr1 = ex<Constant>(1.0) * ex<Tensor>(g);
    auto expr2 = ex<Constant>(-1.0) * ex<Tensor>(g);
    auto expr3 = ex<Constant>(-1.0/12.0) * ex<Tensor>(g);
    std::wcout << g.hash_value() << " " << expr1->hash_value() << " " << expr2->hash_value() << " " << expr3->hash_value() <<"\n";

  }
*/
//  for(auto i = cc_st_r[2]->begin(); i != cc_st_r[2]->end(); ++i)
//    std::cout << (*i)->hash_value() << "\n";
//  std::cout << "\n";
//
//  for(auto i = r2_expr->begin(); i != r2_expr->end(); ++i)
//    std::cout << (*i)->hash_value() << "\n";
//  std::cout << "\n";

  return 0;
}

ExprPtr r2(){
  //﻿http://dx.doi.org/10.1063/1.4907278

  // 39
  auto W_ijam = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
  // std::wcout << "39 W_ijam: " << to_latex(W_ijam) << "\n\n";

  // 36
  auto F_em = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"i_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "36 F_em: " << to_latex(F_em) << "\n\n";

  // 35
  auto F_ea = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"a_1"}) -
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) -
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_5"}, WstrList{L"a_6", L"a_1"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_6"});
  // std::wcout << "35 F_ea: " << to_latex(F_ea) << "\n\n";

  // 34
  auto F_im = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
      F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_5"}) +
      (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
  // std::wcout << "34 F_im: " << to_latex(F_im) << "\n\n";

  // 37
  auto W_ejmn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
  // std::wcout << "37 W_ejmn: " << to_latex(W_ejmn) << "\n\n";

  // 37'
  auto W_eimn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "W_eimn " << to_latex(W_eimn) << "\n\n";

  // 37"
//    auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) +
//        ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
  // std::wcout << "W_iemn " << to_latex(W_iemn) << "\n\n";

  // 43 // Antisymmetrized
  auto W_eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
      W_eimn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
      ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) +
      ex<Constant>(0.25) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
          (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
      ex<Constant>(0.25) *  ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);
  // std::wcout << "W_eima " << to_latex(W_eima) << "\n\n";

  // 44
  auto W_iema = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
      W_iemn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
      ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) -
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"});
  // std::wcout << "W_iema " << to_latex(W_iema) << "\n\n";

  // 49
  auto W_ijmn_temp = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
          (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));
  auto W_ijmn_temp_c = W_ijmn_temp->clone();
  expand(W_ijmn_temp_c);
  // std::wcout << "W_ijmn_temp_c: " << to_latex(W_ijmn_temp_c) << "\n\n";
  std::map<Index, Index> P_ij_mn;
  {
    Index i(L"i_1");
    Index j(L"i_2");
    Index m(L"i_5");
    Index n(L"i_6");
    P_ij_mn.emplace(std::make_pair(i, j));
    P_ij_mn.emplace(std::make_pair(j, i));
    P_ij_mn.emplace(std::make_pair(m, n));
    P_ij_mn.emplace(std::make_pair(n, m));
  }
  // std::wcout << "W_ijmn_temp transformed: " << to_latex(transform_expression(W_ijmn_temp_c->clone(), P_ij_mn)) << "\n\n";

  auto W_ijmn = W_ijmn_temp_c->clone() + transform_expression(W_ijmn_temp_c->clone(), P_ij_mn);
  // std::wcout << "W_ijmn " << to_latex(W_ijmn) << "\n\n";

  //
  auto temp_ab_ = W_iema->clone() * ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm);
  auto temp_ab_c = temp_ab_->clone();
  // std::wcout << "temp_ab_c " << to_latex(temp_ab_c) << "\n\n";
  expand(temp_ab_c);

  std::map<Index, Index> P_ab;
  {
    Index a(L"a_1");
    Index b(L"a_2");
    P_ab.emplace(std::make_pair(a, b));
    P_ab.emplace(std::make_pair(b, a));
  }
  auto temp_ab = ex<Constant>(0.5) * temp_ab_c->clone() + transform_expression(temp_ab_c->clone(), P_ab);
  // std::wcout << "temp_ab " << to_latex(temp_ab) << "\n\n";

  // CCSD Z2
  auto temp1 = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) +
      ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) -
      W_ijam->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}) +
      F_ea->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2"}, WstrList{L"a_5", L"a_2"}) -
      F_im->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_2"}, WstrList{L"a_1", L"a_2"}) +
      ex<Constant>(0.5) * (ex<Constant>(2.0) * W_eima->clone() - W_iema->clone()) *
          (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_2", L"a_5"}, Symmetry::nonsymm)) -
      temp_ab->clone() +
      ex<Constant>(0.5) * W_ijmn->clone() * (ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"}, Symmetry::nonsymm)) +
      ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) *
          (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}, Symmetry::nonsymm));

  auto temp1_c = temp1->clone();
  expand(temp1_c);


  auto ccsd_z2 = ex<Tensor>(L"S", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}) * temp1_c->clone(); // + transform_expression(temp1, P_ab_ij);
  expand(ccsd_z2);
  // std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";

  ccsd_z2 = swap_bra_ket(ccsd_z2);
  // std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";

  canonicalize(ccsd_z2);
  rapid_simplify(ccsd_z2);
  std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";
  return ccsd_z2;
}
