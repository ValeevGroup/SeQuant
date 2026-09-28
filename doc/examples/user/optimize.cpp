#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  // start-snippet-1
  // registers the occupied ("i") / virtual ("a") spaces with concrete
  // extents, which optimize() needs to estimate contraction cost
  mbpt::load(mbpt::Convention::SR);

  // a chain of 3 tensors written as one flat Product; SeQuant does not
  // commit to a pairwise evaluation order until one is asked for
  auto expr = ex<Tensor>(L"A", bra{L"a_1"}, ket{L"i_1"}) *
              ex<Tensor>(L"B", bra{L"i_1"}, ket{L"i_2"}) *
              ex<Tensor>(L"C", bra{L"i_2"}, ket{L"a_2"});

  // optimize() picks a pairwise contraction order (and, for a Sum, reorders
  // summands to share intermediates) that minimizes a cost metric -- by
  // default the total floating-point operation count, using
  // IndexSpace::approximate_size() for index extents
  auto optimized = optimize(expr);

  std::wcout << to_latex(optimized) << std::endl;
  // end-snippet-1

  SEQUANT_ASSERT(optimized.is<Product>());

  return 0;
}
