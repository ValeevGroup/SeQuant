#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <iostream>

// each consumed backend: its exported target must provide the backend's
// headers, its own prerequisite's headers, and the SEQUANT_HAS_<BACKEND> macro
#ifdef SEQUANT_CONSUMER_TA
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
static_assert(SEQUANT_HAS_TILEDARRAY == 1);
#endif
#ifdef SEQUANT_CONSUMER_BTAS
#include <SeQuant/core/eval/backends/btas/result.hpp>
static_assert(SEQUANT_HAS_BTAS == 1);
#endif
#ifdef SEQUANT_CONSUMER_TAPP
#include <SeQuant/core/eval/backends/tapp/result.hpp>
static_assert(SEQUANT_HAS_TAPP == 1);
#endif

int main() {
  using namespace sequant;
  using namespace sequant::mbpt;
  set_default_context({.index_basis_registry_shared_ptr = make_min_sr_spaces(),
                       .vacuum = Vacuum::SingleProduct});
  set_default_mbpt_context({.op_registry_ptr = make_legacy_registry()});

  const auto result = ref_av(H() * T(2));
  std::wcout << result.to_latex() << '\n';
  return 0;
}
