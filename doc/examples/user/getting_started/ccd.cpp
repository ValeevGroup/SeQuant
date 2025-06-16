#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

inline auto commutator(auto op1, auto op2) { return op1 * op2 - op2 * op1; }

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context(Context(make_min_sr_spaces(), Vacuum::SingleProduct));

  auto hbar =
      H(2) + commutator(H(2), T_(2)) +
      ex<Constant>(rational(1, 2)) * commutator(commutator(H(2), T_(2)), T_(2));
  auto ccd_eq = vac_av(P(2) * hbar);
  std::wcout << "<" << to_latex(P(2) * hbar) << "> = " << to_latex(ccd_eq)
             << std::endl;

  return 0;
}
