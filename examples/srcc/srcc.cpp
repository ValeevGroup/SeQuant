
#include "../sequant_setup.hpp"

#define CCSD_r1 1
#define CCSD_r2 1
#define CCSDT 1
#define CCSDTQ 0

using namespace sequant;

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
  constexpr bool print = false;
  // change to true to print stats
  Logger::get_instance().wick_stats = false;

  auto ccsd = cceqvec{2, 2}(true, true, true, true);
# if CCSD_r1
  expand(ccsd[1]);
  rapid_simplify(ccsd[1]);
  // std::wcout << to_latex(ccsd[1]) << "\n\n";

  {
    auto n_term = 0;
    auto result = std::make_shared<Sum>();
    for (auto& term : *ccsd[1]) {
      ++n_term;
      auto traced = ex<Constant>(0.5) * spintrace(term, {{L"i_1", L"a_1"}});
      std::wcout << n_term << ": " << to_latex(term) << "\n";
      expand(traced);
      rapid_simplify(traced);
      canonicalize(traced);
      std::wcout << n_term << ": " << to_latex(traced) << "\n";
      result->append(traced);
    }
    // std::wcout << "R1 traced: " << to_latex(result) << "\n\n";
  }
#endif

# if CCSD_r2
  expand(ccsd[2]);
    rapid_simplify(ccsd[2]);
    // std::wcout << to_latex(ccsd[2]) << "\n";
    {
      auto n_term = 0;
      auto result = std::make_shared<Sum>();
      for (auto &term : *ccsd[2]) {
        ++n_term;
          // std::wcout << n_term << ": " << to_latex(term) << "\n";
          auto traced = ex<Constant>(1) * spintrace(term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}} );
          expand(traced);
          rapid_simplify(traced);
          canonicalize(traced);
          // std::wcout << n_term << "R:  " << to_latex(traced) << "\n\n";
          result->append(traced);
      }
      // std::wcout << "R2 traced: " << to_latex_align(result) << "\n";
      // std::wcout << "R2 traced: " << to_latex(result) << "\n";
      std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                       {Index{L"i_2"}, Index{L"i_1"}}};

      //  auto transformed_result = transform_expression(result, idxmap);
      // 1/3 R + 1/6 R' for simpler equations
      auto temp_expr =  transform_expression(result, idxmap);
      // std::wcout << "temp_expr: " << to_latex(temp_expr) << "\n";
      auto simpler_R2 = ex<Constant>(1.0/3.0) * result +  ex<Constant>(1.0/6.0) * temp_expr;
      // std::wcout << "simpler_R2:\n" << to_latex(simpler_R2) << "\n";
      expand(simpler_R2);
      rapid_simplify(simpler_R2);
        if(simpler_R2->is<Tensor>())
          ranges::for_each(simpler_R2->as<Tensor>().const_braket(),
                           [&](const Index& idx) { idx.reset_tag(); });
      canonicalize(simpler_R2);
      simpler_R2->rapid_canonicalize();
      rapid_simplify(simpler_R2);
      // std::wcout << "simpler_R2:\n" << to_latex_align(simpler_R2,10,3) << "\n";
      // std::wcout << "simpler_R2:\n" << to_latex(simpler_R2) << "\n";
    }
#endif

#if CCSDT
  {
    auto ccsdt = cceqvec{3, 3}(true, true, true, true);
    auto n_term = 0;
    auto result = std::make_shared<Sum>();
    for (auto &term : *ccsdt[1]) {
      ++n_term;
      std::wcout << n_term << ": " << to_latex(term) << "\n";
//       auto traced = ex<Constant>(1) * spintrace(term, {{L"i_1", L"a_1"}});
//       auto traced = ex<Constant>(1) * spintrace(term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
       auto traced = ex<Constant>(1) * spintrace(term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      // std::wcout << n_term << ": " << to_latex(traced) << "\n";
      // if(n_term > 1) {
        expand(traced);
        rapid_simplify(traced);
        canonicalize(traced);
        std::wcout << n_term << "R:  " << to_latex(traced) << "\n\n";
        result->append(traced);
      // }
    }
  }
#endif

#if CCSDTQ
  {
    auto ccsdtq = cceqvec{4, 4}(true, true, true, true);
    auto n_term = 0;
    auto result = std::make_shared<Sum>();

    for (auto &term : *ccsdtq[4]) {
      ++n_term;
      std::wcout << n_term << ": " << to_latex(term) << "\n";
      auto traced = ex<Constant>(1) * spintrace(term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
//      expand(traced);
//      rapid_simplify(traced);
//      canonicalize(traced);
//      std::wcout << n_term << "R:  " << to_latex(traced) << "\n\n";
//      result->append(traced);
    }
  }
#endif



  /*

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
  */
  return 0;
}

