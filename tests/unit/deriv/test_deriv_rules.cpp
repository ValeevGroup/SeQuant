#include <SeQuant/core/deriv/deriv.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/rational.hpp>

#include <catch2/catch_test_macros.hpp>
#include "catch2_sequant.hpp"

#include <string_view>

using namespace sequant;

TEST_CASE("deriv/rules/scalars", "[deriv][rules]") {
  const auto diff = [](const ExprPtr& E, const ExprPtr& x) {
    return sequant::deriv::diff(E, {x});
  };
  auto T = [](std::wstring_view s) { return deserialize(s); };
  const auto x = ex<Variable>(L"x");
  const auto x_star = adjoint(x);

  REQUIRE(diff(ex<Constant>(3), x)->is_zero());
  REQUIRE(diff(T(L"y z"), x)->is_zero());
  REQUIRE_THAT(diff(x, x), EquivalentTo(L"1"));
  REQUIRE_THAT(diff(T(L"3 x^(2) y + x"), x), EquivalentTo(L"6 x y + 1"));
  REQUIRE_THAT(diff(T(L"x y x"), x), EquivalentTo(L"2 x y"));
  REQUIRE_THAT(diff(T(L"x^(1/2)"), x), EquivalentTo(L"1/2 x^(-1/2)"));
  REQUIRE_THAT(diff(T(L"x^(-1)"), x), EquivalentTo(L"-1 x^(-2)"));
  // x and x^* are independent
  REQUIRE_THAT(diff(x * x_star, x), EquivalentTo(x_star));
  REQUIRE_THAT(diff(x * x_star, x_star), EquivalentTo(x));
  // a conjugated power is a power of x^*
  auto p = ex<Power>(x, rational(3));
  p->adjoint();
  auto q = ex<Power>(x, rational(2));
  q->adjoint();
  REQUIRE(diff(p, x)->is_zero());
  REQUIRE_THAT(diff(p, x_star), EquivalentTo(ex<Constant>(3) * q));
  REQUIRE_THAT(diff(q, x_star), EquivalentTo(ex<Constant>(2) * x_star));
  // a power of x^* is one too
  const auto r = ex<Power>(x_star, rational(3));
  REQUIRE(diff(r, x)->is_zero());
  REQUIRE_THAT(diff(r, x_star),
               EquivalentTo(ex<Constant>(3) * ex<Power>(x_star, rational(2))));
  // a nested product with its own scalar
  REQUIRE_THAT(
      diff(ex<Product>(2,
                       ExprPtrList{x, ex<Product>(3, ExprPtrList{x, T(L"y")})}),
           x),
      EquivalentTo(L"12 x y"));
  // a deeper nested sum
  REQUIRE_THAT(diff(x * (T(L"y") + x * (x + T(L"z"))), x),
               EquivalentTo(L"y + 3 x^(2) + 2 x z"));
  // tensors are constants
  REQUIRE_THAT(diff(T(L"x^(3) f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"), x),
               EquivalentTo(L"3 x^(2) f{i_1;a_1}:N-N-N t{a_1;i_1}:N-N-N"));
  // an input that is not expanded
  REQUIRE_THAT(diff((x + T(L"y")) * (x + T(L"z")), x),
               EquivalentTo(L"2 x + y + z"));
  // any order, left to right; no variable is the function itself
  const auto y = ex<Variable>(L"y");
  REQUIRE_THAT(sequant::deriv::diff(T(L"x^(3) y^(2)"), {x, x, y}),
               EquivalentTo(L"12 x y"));
  REQUIRE_THAT(sequant::deriv::diff(T(L"x^(3) y^(2)"), {}),
               EquivalentTo(L"x^(3) y^(2)"));
}
