//
// Kramers tracer (spinor.{hpp,cpp}) unit tests.
//

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/spinor.hpp>

#include <cstddef>
#include <memory>

TEST_CASE("kramers_trace", "[spinor]") {
  using namespace sequant;
  using namespace sequant::mbpt;

  auto ctx = get_default_context();
  ctx.set(CanonicalizeOptions{.method = CanonicalizationMethod::Complete});
  auto _ = set_scoped_default_context(ctx);
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // counts RealPart nodes anywhere in the expression tree
  auto count_realparts = [](const ExprPtr& expr) {
    std::size_t n = 0;
    expr->visit(
        [&n](const ExprPtr& current) {
          if (current->is<RealPart>()) ++n;
        },
        /* atoms_only = */ false);
    return n;
  };

  SECTION("MP2 energy") {
    // E = 1/4 g-bar^{a1 a2}_{i1 i2} t-bar^{i1 i2}_{a1 a2}  (Kramers-free)
    const auto E = ex<Constant>(rational{1, 4}) *
                   ex<Tensor>(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                              Symmetry::Antisymm) *
                   ex<Tensor>(L"t", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"},
                              Symmetry::Antisymm);

    ExprPtr result;
    REQUIRE_NOTHROW(result = closed_shell_kramers_trace(E));
    REQUIRE(result);

    INFO("closed_shell_kramers_trace(E_MP2) =\n" << toUtf8(to_latex(result)));

    REQUIRE(result->is<Sum>());

    // Each summand is `c * Re[1/2 g . t-bar]`, one per TRS-canonical Kramers
    // representative. With the global-T (whole-config) fold and sigma handled
    // by canonicalize — but the internal-T-reach (self-complementary block
    // fold) deferred — the MP2 energy yields 7 representatives: the 6-term form
    // plus the unfolded self-complementary pair (g^{a^ a_} and -g^{a_ a^}).
    REQUIRE(result->size() == 7);
    REQUIRE(count_realparts(result) == 7);

    // Coefficient bookkeeping: every summand is `Constant * RealPart`, and the
    // outer coefficients sum to 2^n = 16 (all Kramers configurations accounted
    // for: T-fold contributes x2, sigma multiplicity the rest).
    Constant::scalar_type coeff_sum = 0;
    for (const auto& term : *result) {
      REQUIRE(term->is<Product>());
      coeff_sum += term->as<Product>().scalar();
    }
    REQUIRE(coeff_sum == Constant::scalar_type{16});
  }
}
