
#include <SeQuant/core/density.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/rdm.hpp>

#include <clocale>
#include <iostream>

using namespace sequant;

void try_main() {
  using namespace sequant::mbpt;
  std::cout << "START ANTISYMM_TEST: " << std::endl;
  const auto cumulant = density::make_cumulant(bra{"a_1"}, ket{"i_1"});
  // const auto a =ex<Tensor>("a",bra{"i_2", "i_3"},ket{"a_2",
  // "a_3"});
  const auto a = ex<FNOperator>(cre{Index("i_2"), Index("i_3")},
                                ann{Index("a_2"), Index("a_3")});
  auto a_cumulant = cumulant * a;
  std::cout << "a_cumulant " << to_latex_align(a_cumulant) << std::endl;
  antisymmetrize _a_cumulant(a_cumulant);
  std::cout << to_latex_align(_a_cumulant.result) << std::endl;

  auto cumulant2 = density::make_cumulant(bra{"a_2"}, ket{"i_2"});
  auto cumulant3 = density::make_cumulant(bra{"a_3"}, ket{"i_3"});
  auto cumulant_3x = cumulant * cumulant2 * cumulant3;
  std::cout << "cumulant_3x " << to_latex_align(cumulant_3x) << std::endl;
  antisymmetrize _cumulant_3x(cumulant_3x);
  std::cout << to_latex_align(_cumulant_3x.result) << std::endl;

  auto a1 = ex<FNOperator>(cre{Index("i_1")}, ann{Index("a_1")});
  auto a1_cumu1_cumu2 = a1 * cumulant2 * cumulant3;
  std::cout << "a1 y1 y2 " << to_latex_align(a1_cumu1_cumu2) << std::endl;
  antisymmetrize _a1_cumu1_cumu2(a1_cumu1_cumu2);
  std::cout << to_latex_align(_a1_cumu1_cumu2.result) << std::endl;

  auto two_body_cumu =
      density::make_cumulant(bra{"a_2", "a_3"}, ket{"i_2", "i_3"});
  auto a1_cumu2 = a1 * two_body_cumu;
  std::cout << " a1 y2 " << to_latex_align(a1_cumu2) << std::endl;
  antisymmetrize _a1_cumu2(a1_cumu2);
  std::cout << to_latex_align(_a1_cumu2.result) << std::endl;

  auto cumu1_cumu2 = cumulant * two_body_cumu;
  std::cout << " y1 y2 " << to_latex_align(cumu1_cumu2) << std::endl;
  antisymmetrize _cumu1_cumu2(cumu1_cumu2);
  std::cout << to_latex_align(_cumu1_cumu2.result) << std::endl;

  auto cumu3 = density::make_cumulant(bra{"a_1", "a_2", "a_3"},
                                      ket{"i_1", "i_2", "i_3"});
  std::cout << " y3 " << to_latex_align(cumu3) << std::endl;
  antisymmetrize _cumu3(cumu3);
  std::cout << to_latex_align(_cumu3.result) << std::endl;

  auto antisymm_test =
      ex<Constant>(-1) * _cumu1_cumu2.result + _a1_cumu2.result +
      ex<Constant>(-2) * _a1_cumu1_cumu2.result +
      ex<Constant>(4) * _cumulant_3x.result +
      _a_cumulant
          .result;  // expression from Cite as: J. Chem. Phys. 127, 104107
                    // (2007); https://doi.org/10.1063/1.2761870 eqn 26.
  std::cout << "END ANTISYMM TEST: " << std::endl << std::endl << std::endl;

  auto a3 = ex<FNOperator>(cre{Index("i_1"), Index("i_2"), Index("i_3")},
                           ann{Index("a_1"), Index("a_2"), Index("a_3")});
  auto new_a3 = decompositions::three_body_substitution(a3, 2);
  std::cout << "Now TESTING SPIN SUMMATION";
  // auto result = new_a3 - antisymm_test;
  // simplify(result);
  std::cout << to_latex_align(new_a3, 20, 7) << std::endl;
}

int main() {
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  std::cerr.precision(std::numeric_limits<double>::max_digits10);
  sequant::set_locale();
  sequant::set_default_context(
      {.index_basis_registry_shared_ptr = sequant::mbpt::make_min_sr_spaces(),
       .spbasis = SPBasis::Spinfree});

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
