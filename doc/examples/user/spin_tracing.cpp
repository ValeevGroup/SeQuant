#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;

  // start-snippet-1
  // register the standard occupied/virtual spaces used in single-reference
  // quantum chemistry (see the User Guide's "Context and configuration" page)
  load(Convention::SR);

  Index i1(L"i_1"), a1(L"a_1");
  TensorSymmetries particle_symmetric{.column = ColumnSymmetry::Symm};

  // a spin-orbital expression: indices carry no explicit spin label, so each
  // one implicitly ranges over both alpha and beta spin-orbitals
  auto expr = ex<Tensor>(L"t", bra{i1}, ket{a1}, particle_symmetric) *
              ex<Tensor>(L"F", bra{a1}, ket{i1}, particle_symmetric);

  // spintrace() sums a spin-orbital expression over spin, producing its
  // spin-free (spatial-orbital) form; the summation over the 2 spin cases
  // that agree with each other produces the factor of 2 below
  auto expr_st = spintrace(expr);
  simplify(expr_st);

  std::wcout << to_latex(expr_st) << std::endl;
  // end-snippet-1

  SEQUANT_ASSERT(expr_st.is<Product>());
  SEQUANT_ASSERT(expr_st.as<Product>().scalar() == 2);

  // start-snippet-2
  // Spin quantum numbers can also be inspected/attached to individual
  // indices directly, e.g. when interpreting an already spin-orbital
  // expression by hand
  Index i1_up = make_spinalpha(i1);
  SEQUANT_ASSERT(to_spin(i1_up.space().qns()) == Spin::alpha);
  // end-snippet-2

  return 0;
}
