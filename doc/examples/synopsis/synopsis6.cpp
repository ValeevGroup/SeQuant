#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

inline auto commutator(auto op1, auto op2) { return op1 * op2 - op2 * op1; }

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_space_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  set_default_mbpt_context({.op_registry_ptr = make_legacy_registry()});

  auto hbar =
      H(2) + commutator(H(2), T_(2)) +
      ex<Constant>(rational(1, 2)) * commutator(commutator(H(2), T_(2)), T_(2));
  auto ccd_eq = op::vac_av(P(nₚ(2)) * hbar);
  std::wcout << "<" << to_latex(P(nₚ(2)) * hbar) << "> = " << to_latex(ccd_eq)
             << std::endl;

  return 0;
}
