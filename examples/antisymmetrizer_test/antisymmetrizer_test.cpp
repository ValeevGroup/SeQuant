
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/rdm.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>

#include <clocale>
#include <iostream>

using namespace sequant;

void try_main() {
  using namespace sequant::mbpt;
  std::wcout << "START ANTISYMM_TEST: " << std::endl;
  const auto cumulant = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                                   WstrList{L"a_1"}, WstrList{L"i_1"});
  // const auto a =ex<Tensor>(L"a",WstrList{L"i_2", L"i_3"},WstrList{L"a_2",
  // L"a_3"});
  const auto a = ex<FNOperator>(
      std::initializer_list<Index>({Index(L"i_2"), Index(L"i_3")}),
      std::initializer_list<Index>({Index(L"a_2"), Index(L"a_3")}));
  auto a_cumulant = cumulant * a;
  std::wcout << "a_cumulant " << to_latex_align(a_cumulant) << std::endl;
  antisymmetrize _a_cumulant(a_cumulant);
  std::wcout << to_latex_align(_a_cumulant.result) << std::endl;

  auto cumulant2 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                              WstrList{L"a_2"}, WstrList{L"i_2"});
  auto cumulant3 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                              WstrList{L"a_3"}, WstrList{L"i_3"});
  auto cumulant_3x = cumulant * cumulant2 * cumulant3;
  std::wcout << "cumulant_3x " << to_latex_align(cumulant_3x) << std::endl;
  antisymmetrize _cumulant_3x(cumulant_3x);
  std::wcout << to_latex_align(_cumulant_3x.result) << std::endl;

  auto a1 = ex<FNOperator>(std::initializer_list<Index>({Index(L"i_1")}),
                           std::initializer_list<Index>({Index(L"a_1")}));
  auto a1_cumu1_cumu2 = a1 * cumulant2 * cumulant3;
  std::wcout << "a1 y1 y2 " << to_latex_align(a1_cumu1_cumu2) << std::endl;
  antisymmetrize _a1_cumu1_cumu2(a1_cumu1_cumu2);
  std::wcout << to_latex_align(_a1_cumu1_cumu2.result) << std::endl;

  auto two_body_cumu =
      ex<Tensor>(optype2label.at(OpType::RDMCumulant), WstrList{L"a_2", L"a_3"},
                 WstrList{L"i_2", L"i_3"});
  auto a1_cumu2 = a1 * two_body_cumu;
  std::wcout << " a1 y2 " << to_latex_align(a1_cumu2) << std::endl;
  antisymmetrize _a1_cumu2(a1_cumu2);
  std::wcout << to_latex_align(_a1_cumu2.result) << std::endl;

  auto cumu1_cumu2 = cumulant * two_body_cumu;
  std::wcout << " y1 y2 " << to_latex_align(cumu1_cumu2) << std::endl;
  antisymmetrize _cumu1_cumu2(cumu1_cumu2);
  std::wcout << to_latex_align(_cumu1_cumu2.result) << std::endl;

  auto cumu3 = ex<Tensor>(optype2label.at(OpType::RDMCumulant),
                          WstrList{L"a_1", L"a_2", L"a_3"},
                          WstrList{L"i_1", L"i_2", L"i_3"});
  std::wcout << " y3 " << to_latex_align(cumu3) << std::endl;
  antisymmetrize _cumu3(cumu3);
  std::wcout << to_latex_align(_cumu3.result) << std::endl;

  auto antisymm_test =
      ex<Constant>(-1) * _cumu1_cumu2.result + _a1_cumu2.result +
      ex<Constant>(-2) * _a1_cumu1_cumu2.result +
      ex<Constant>(4) * _cumulant_3x.result +
      _a_cumulant
          .result;  // expression from Cite as: J. Chem. Phys. 127, 104107
                    // (2007); https://doi.org/10.1063/1.2761870 eqn 26.
  std::wcout << "END ANTISYMM TEST: " << std::endl << std::endl << std::endl;

  auto a3 = ex<FNOperator>(std::initializer_list<Index>(
                               {Index(L"i_1"), Index(L"i_2"), Index(L"i_3")}),
                           std::initializer_list<Index>(
                               {Index(L"a_1"), Index(L"a_2"), Index(L"a_3")}));
  auto new_a3 = decompositions::three_body_substitution(a3, 2);
  std::wcout << "Now TESTING SPIN SUMMATION";
  // auto result = new_a3 - antisymm_test;
  // simplify(result);
  std::wcout << to_latex_align(new_a3, 20, 7) << std::endl;
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
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(Context(Vacuum::Physical, IndexSpaceMetric::Unit,
                                       BraKetSymmetry::conjugate,
                                       SPBasis::spinfree));
  mbpt::set_default_convention();

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  // WARNING some code is not thread safe ...
  // set_num_threads(1);

  try {
    try_main();
  } catch (std::exception& ex) {
    std::cerr << "caught a std::exception: " << ex.what();
  } catch (...) {
    std::cerr << "caught an unknown exception, ouch";
  }

  return 0;
}
