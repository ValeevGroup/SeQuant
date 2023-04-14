#include <SeQuant/core/op.hpp>
#include <SeQuant/core/timer.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/formalism.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <clocale>

using namespace sequant;

#define runtime_assert(tf)                                       \
  if (!(tf)) {                                                   \
    std::ostringstream oss;                                      \
    oss << "failed assert at line " << __LINE__                  \
        << " in open-shell spin-traced coupled cluster example"; \
    throw std::runtime_error(oss.str().c_str());                 \
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
  mbpt::set_default_formalism();
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
  auto cc_r = sequant::mbpt::sr::so::cceqs{NMAX}.t();
  for (auto i = 1; i < cc_r.size(); ++i) {
    std::cout << "Spin-orbital CC R" << i << " size: " << cc_r[i]->size()
              << "\n";
  }

  //
  // Open-shell spintrace
  //
  std::cout << "\nOpen-shell coupled cluster: nterms per spin blocks: "
            << std::endl;
  std::vector<std::vector<ExprPtr>> os_cc_st_r(cc_r.size());
  for (auto i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    Tensor A =
        cc_r[i]->as<Sum>().summand(0)->as<Product>().factors()[0]->as<Tensor>();
    assert(A.label() == L"A");
    auto P_vec = open_shell_P_op_vector(A);
    auto A_vec = open_shell_A_op(A);
    assert(P_vec.size() == i + 1);

    std::vector<Sum> concat_terms(i + 1);
    size_t n_spin_orbital_term = 0;
    for (auto& product_term : *cc_r[i]) {
      auto term = remove_tensor(product_term->as<Product>(), L"A");
      std::vector<ExprPtr> os_st(i + 1);

      // Apply the P operators on the product term without the A,
      // Expand the P operators and spin-trace the expression
      // Then apply A operator, canonicalize and remove A operator
      for (int s = 0; s != os_st.size(); ++s) {
        os_st.at(s) = P_vec.at(s) * term;
        expand(os_st.at(s));
        os_st.at(s) = expand_P_op(os_st.at(s), false, true);
        os_st.at(s) =
            open_shell_spintrace(os_st.at(s), ext_idx_list(i), s).at(0);
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
    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf(" time:  %5.3f sec.\n", time_elapsed.count());
  }

  if (NMAX == 4) {
    runtime_assert(os_cc_st_r.size() == 5)
        runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30)   // T1a
        runtime_assert(os_cc_st_r.at(2).at(1)->size() == 130)  // T2ab
        runtime_assert(os_cc_st_r.at(2).at(2)->size() == 74)   // T2bb
        runtime_assert(os_cc_st_r.at(3).at(1)->size() == 249)  // T3aab
        runtime_assert(os_cc_st_r.at(3).at(3)->size() == 124)  // T3bbb
        runtime_assert(os_cc_st_r.at(4).at(1)->size() == 356)  // T4aaab
        runtime_assert(os_cc_st_r.at(4).at(2)->size() == 384)  // T4aabb
        runtime_assert(os_cc_st_r.at(4).at(4)->size() == 156)  // T4bbbb
  } else if (NMAX == 3) {
    runtime_assert(os_cc_st_r.size() == 4)
        runtime_assert(os_cc_st_r.at(1).at(0)->size() == 30)   // T1a
        runtime_assert(os_cc_st_r.at(2).at(0)->size() == 65)   // T2aa
        runtime_assert(os_cc_st_r.at(2).at(1)->size() == 122)  // T2ab
        runtime_assert(os_cc_st_r.at(3).at(2)->size() == 209)  // T3abb
        runtime_assert(os_cc_st_r.at(3).at(3)->size() == 75)   // T3bbb
  }
}
