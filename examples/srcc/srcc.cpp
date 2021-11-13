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

#if 0
  // First 4 terms only
    cc_r[1] = cc_r[1]->as<Sum>().take_n(5);
    cc_r[2] = cc_r[2]->as<Sum>().take_n(5);
  //  cc_r[3] = cc_r[3]->as<Sum>().take_n(4);

  for (int i = 1; i < cc_r.size(); ++i) {
  // size_t counter = 1;
  std::vector<size_t> n_st_terms(i+1,0);
  for (auto& product_term : *cc_r[i]) {
    // auto term = remove_tensor_from_product(product_term->as<Product>(), L"A");
    auto term = product_term;
     std::wcout << "Input: " << to_latex(term) << "\n";
    const auto list = ext_idx_list(i);
    auto os_st = open_shell_spintrace(term, list);
    for (size_t j = 0; j != os_st.size(); ++j) {
       std::wcout<< "st: " << to_latex(os_st[j]) << "\n";
      n_st_terms[j] += os_st[j]->size();
    }
      std::wcout << "\n";
  }
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
  std::vector<ExprPtr> cc_st_r(cc_r.size(),nullptr);
  for (int i = 1; i < cc_r.size(); ++i) {
    const auto tstart = std::chrono::high_resolution_clock::now();
    cc_st_r[i] = closed_shell_CC_spintrace(cc_r[i]);
    auto tstop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = tstop - tstart;
    printf("CC R%d size: %lu closed-shell spin-trace with biorthogonal transformation time: %5.6f sec.\n\n", i, cc_st_r[i]->size(),
           time_elapsed.count());
  }
  std::wcout << "CCSDT: " << to_latex(*cc_st_r[3]) << std::endl;
#endif

  return 0;
}
