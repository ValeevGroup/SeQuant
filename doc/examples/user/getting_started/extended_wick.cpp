#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <iostream>

int main() {
  using namespace sequant;

  // start-snippet-1
  // a multireference vocabulary of index spaces: core (o, i), active (u),
  // virtual (a, g); the reference is a general state in the active space
  set_default_context(
      Context({.index_basis_registry_shared_ptr = mbpt::make_mr_spaces(),
               .vacuum = Vacuum::MultiProduct}));

  // a product of two generalized-normal-ordered one-body operators
  auto expr = ex<FNOperator>(cre({L"u_1"}), ann({L"u_2"})) *
              ex<FNOperator>(cre({L"u_3"}), ann({L"u_4"}));

  // its reference expectation value: a γη pair plus a 2-body cumulant κ
  auto vev = FWickTheorem{expr}.compute();
  std::wcout << to_latex(vev) << std::endl;

  // its generalized-normal-ordered form, cumulants truncated at rank 2
  auto gno = FWickTheorem{expr}
                 .full_contractions(false)
                 .max_cumulant_rank(2)
                 .compute();
  std::wcout << to_latex(gno) << std::endl;
  // end-snippet-1

  SEQUANT_ASSERT(vev->size() == 2);

  return 0;
}
