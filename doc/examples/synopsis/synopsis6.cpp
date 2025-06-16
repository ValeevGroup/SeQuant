#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#if __cplusplus >= 202002L
inline auto commutator(auto op1, auto op2) { return op1 * op2 - op2 * op1; }
#else
inline auto commutator(sequant::ExprPtr op1, sequant::ExprPtr op2) {
  return op1 * op2 - op2 * op1;
}
#endif

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context(Context(make_min_sr_spaces(), Vacuum::SingleProduct));

  auto hbar =
      H(2) + commutator(H(2), T_(2)) +
      ex<Constant>(rational(1, 2)) * commutator(commutator(H(2), T_(2)), T_(2));
  auto ccd_eq = op::vac_av(P(nₚ(2)) * hbar);
  std::wcout << "<" << to_latex(P(nₚ(2)) * hbar) << "> = " << to_latex(ccd_eq)
             << std::endl;

  return 0;
}
