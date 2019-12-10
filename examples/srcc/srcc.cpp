#include "../contract/sequant_setup.hpp"

using namespace  sequant;

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
  //set_num_threads(1);

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

  auto ccsd = cceqvec{2, 2}(true, true, true, true);
  // std::wcout << "CCSD R1 =" << to_latex_align(ccsd[1], 20) << std::endl;
  // auto function_check = spintrace(ccsd[2]);

  Tensor g = Tensor(L"g", {Index{L"a_1"}, Index{L"a_2"}},
                    {Index{L"i_1"}, Index{L"i_2"}}, Symmetry::antisymm);
  Tensor T2 = Tensor(L"t",{Index{L"i_4"}, Index{L"i_5"}}, {Index{L"a_4"}, Index{L"a_5"}}, Symmetry::antisymm);
  Tensor T3 = Tensor(L"t",{Index{L"i_1"}, Index{L"i_2"}, Index{L"i_3"}}, {Index{L"a_1"}, Index{L"a_2"}, Index{L"a_3"}}, Symmetry::antisymm);
  Tensor T4 = Tensor(L"t",{Index{L"i_6"}, Index{L"i_7"}, Index{L"i_8"}, Index{L"i_9"}}, {Index{L"a_6"}, Index{L"a_7"}, Index{L"a_8"}, Index{L"a_9"}}, Symmetry::antisymm);
  // std::wcout << g.to_latex() << "\t" << T3.to_latex() << std::endl;

  Product expression{};
  {
    ExprPtr g_ptr = std::make_shared<Tensor>(g);
    ExprPtr t2_ptr = std::make_shared<Tensor>(T2);
    ExprPtr t3_ptr = std::make_shared<Tensor>(T3);
    ExprPtr t4_ptr = std::make_shared<Tensor>(T4);

    expression.append(1, g_ptr);
    expression.append(1, t2_ptr);
    expression.append(1, t3_ptr);
    expression.append(1,t4_ptr);
    // std::wcout << "expression: " << expression.to_latex() << std::endl;
  }
  ExprPtr expr_ptr = std::make_shared<Product>(expression);
  Sum expr_sum{};
  expr_sum.append(expr_ptr);
  ExprPtr temp_expr = std::make_shared<Sum>(expr_sum);
  auto function_check = spintrace(temp_expr);

  // std::wcout << "CCSD R2 =" << to_latex_align(ccsd[2], 20) << std::endl;
  // function_check = spintrace(ccsd[2]);

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

  return 0;
}
