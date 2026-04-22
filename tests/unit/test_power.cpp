//
// Created by Ajay Melekamburath on 4/18/26.
//

#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

TEST_CASE("power", "[elements]") {
  using namespace sequant;

  SECTION("constructors") {
    const auto c2 = ex<Constant>(rational{1, 2});
    const auto vx = ex<Variable>(L"x");

    REQUIRE_NOTHROW(Power(c2, rational{1, 2}));
    REQUIRE_NOTHROW(Power(vx, rational{3, 1}));

    // convenience ctors:
    REQUIRE(Power(L"x", 2) == Power(vx, rational{2}));
    REQUIRE(Power(L"x", rational{1, 2}) == Power(vx, rational{1, 2}));
    REQUIRE(Power(2, 3) == Power(ex<Constant>(2), rational{3}));
    REQUIRE(Power(rational{2, 3}, 2) ==
            Power(ex<Constant>(rational{2, 3}), rational{2}));
    // power-of-power flattens: Power(Power(b, e1), e2) -> Power(b, e1*e2)
    const auto inner = ex<Power>(c2, rational{1, 2});
    Power outer(inner, rational{2, 3});
    REQUIRE(outer.exponent() == rational{1, 3});  // 1/2 * 2/3 = 1/3
    REQUIRE(outer.base() == c2);

    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      // base must be a Constant, Variable, or Power
      auto bad_base = ex<Product>(Product{});
      REQUIRE_THROWS(Power(bad_base, rational{2}));

      // 0^n is defined only for n >= 0
      REQUIRE_THROWS(Power(ex<Constant>(0), rational{-1}));
      REQUIRE_THROWS(Power(ex<Constant>(0), rational{-1, 2}));
    }
  }

  SECTION("accessors") {
    const auto c2 = ex<Constant>(rational{1, 2});
    Power p(c2, rational{1, 2});
    REQUIRE(p.base() == ex<Constant>(rational{1, 2}));
    REQUIRE(p.exponent() == rational{1, 2});

    // is_zero: base == 0, exponent > 0
    Power pz(ex<Constant>(0), rational{2});
    REQUIRE(pz.is_zero());
    // 0^0 is not zero by our convention
    Power pz2(ex<Constant>(0), rational{0});
    REQUIRE(!pz2.is_zero());
    REQUIRE(!p.is_zero());
  }

  SECTION("comparison") {
    const auto c2 = ex<Constant>(rational{1, 2});
    const auto vx = ex<Variable>(L"x");

    Power p1(c2, rational{1, 2});
    Power p2(c2, rational{1, 2});
    REQUIRE(p1 == p2);

    Power p3(c2, rational{1, 3});
    REQUIRE(!(p1 == p3));

    Power p4(vx, rational{1, 2});
    REQUIRE(!(p1 == p4));

    // static_less_than: same base, compare by exponent
    Power plt_a(vx, rational{1, 2});
    Power plt_b(vx, rational{3, 4});
    REQUIRE(plt_a < plt_b);
    Power plt_c(vx, rational{1, 2});
    REQUIRE(!(plt_a < plt_c));
  }

  SECTION("hashing") {
    const auto c2 = ex<Constant>(rational{1, 2});
    const auto vx = ex<Variable>(L"x");

    Power p1(c2, rational{1, 2});
    Power p2(c2, rational{1, 2});
    Power p3(c2, rational{1, 3});  // different exponent
    Power p4(vx, rational{1, 2});  // different base
    REQUIRE(sequant::hash::value(p1) == sequant::hash::value(p2));
    REQUIRE(sequant::hash::value(p1) != sequant::hash::value(p3));
    REQUIRE(sequant::hash::value(p1) != sequant::hash::value(p4));
  }

  SECTION("adjoint") {
    // adjoint conjugates the base; exponent unchanged
    Power pv(ex<Variable>(L"z"), rational{1, 2});
    pv.adjoint();
    REQUIRE(pv.base()->as<Variable>().conjugated());
    REQUIRE(pv.exponent() == rational{1, 2});

    // adjoint must not mutate a base shared with another expression
    auto shared_base = ex<Variable>(L"w");
    Power ps1(shared_base, rational{1, 2});
    Power ps2(shared_base, rational{1, 3});
    ps1.adjoint();
    REQUIRE(ps1.base()->as<Variable>().conjugated());
    REQUIRE(!ps2.base()->as<Variable>().conjugated());
    REQUIRE(!shared_base->as<Variable>().conjugated());
  }

  SECTION("latex") {
    const auto c2 = ex<Constant>(rational{1, 2});

    Power p1(c2, rational{1, 2});
    REQUIRE(p1.to_latex() == L"{{{\\frac{1}{2}}}}^{\\frac{1}{2}}");

    Power p_exp1(c2, rational{1});
    REQUIRE(p_exp1.to_latex() == L"{{{\\frac{1}{2}}}}");

    Power pv(ex<Variable>(L"x"), rational{2, 1});
    REQUIRE(pv.to_latex() == L"{x}^{2}");
  }

  SECTION("flatten") {
    // Power::flatten: Constant base + integer exponent folds in place
    auto pf1 = ex<Power>(2, 3);
    Power::flatten(pf1);
    REQUIRE(pf1 == ex<Constant>(rational{8}));

    auto pf2 = ex<Power>(2, -3);
    Power::flatten(pf2);
    REQUIRE(pf2 == ex<Constant>(rational{1, 8}));

    auto pf3 = ex<Power>(rational{2, 3}, 0);
    Power::flatten(pf3);
    REQUIRE(pf3 == ex<Constant>(rational{1}));

    // non-integer exponent on a Constant base is a no-op
    auto pf4 = ex<Power>(2, rational{1, 2});
    Power::flatten(pf4);
    REQUIRE(pf4->is<Power>());

    // Variable base is a no-op
    auto pf5 = ex<Power>(L"x", 2);
    Power::flatten(pf5);
    REQUIRE(pf5->is<Power>());

    // large exponents
    auto pf6 = ex<Power>(2, 1000000);
    Power::flatten(pf6);
    REQUIRE(pf6->is<Constant>());

    // 2^20 = 1048576
    auto pf7 = ex<Power>(2, 20);
    Power::flatten(pf7);
    REQUIRE(pf7 == ex<Constant>(rational{1048576}));

    // 2^(-20) = 1/1048576
    auto pf8 = ex<Power>(2, -20);
    Power::flatten(pf8);
    REQUIRE(pf8 == ex<Constant>(rational{1, 1048576}));

    // simplify folds a foldable Power-in-Product into the Product scalar
    auto pwf = ex<Power>(2, 2) * ex<Variable>(L"x");
    simplify(pwf);
    REQUIRE(pwf->is<Product>());
    REQUIRE(pwf->as<Product>().scalar() == 4);

    // non-foldable Power stays as a factor; scalar is 1
    auto pwnf = ex<Power>(2, rational{1, 2}) * ex<Variable>(L"x");
    simplify(pwnf);
    REQUIRE(pwnf->is<Product>());
    REQUIRE(pwnf->as<Product>().scalar() == 1);
  }

  SECTION("operator*=") {
    const auto vx = ex<Variable>(L"x");

    // b^e1 *= b^e2 -> b^(e1+e2)
    Power pa(vx, rational{1, 2});
    Power pb(vx, rational{1, 3});
    pa *= pb;
    REQUIRE(pa.exponent() == rational{5, 6});  // 1/2 + 1/3

    // bare base: b^e *= b -> b^(e+1)
    Power pc(vx, rational{1, 2});
    pc *= static_cast<const Expr &>(*vx);
    REQUIRE(pc.exponent() == rational{3, 2});

    // 2^{1/2} * 2^{1/2} = 2
    Power pe(ex<Constant>(2), rational{1, 2});
    Power pf(ex<Constant>(2), rational{1, 2});
    pe *= pf;
    REQUIRE(pe.exponent() == rational{1});
    REQUIRE(pe.to_latex() == Constant(2).to_latex());
  }

  SECTION("in Product") {
    const auto vx = ex<Variable>(L"x");

    // Power should NOT be absorbed into Product::scalar_
    auto pw = ex<Power>(vx, rational{1, 2});
    auto prod = ex<Product>(Product{});
    prod->as<Product>().append(1, pw, Product::Flatten::Yes);
    REQUIRE(prod->as<Product>().factors().size() == 1);
    REQUIRE(prod->as<Product>().scalar() == 1);
  }

  SECTION("in Sum") {
    const auto vx = ex<Variable>(L"x");

    // x^{1/2} + x^{1/2} = 2 * x^{1/2}
    auto pw1 = ex<Power>(vx, rational{1, 2});
    auto pw2 = ex<Power>(vx, rational{1, 2});
    auto sum_expr = pw1 + pw2;
    simplify(sum_expr);
    REQUIRE(sum_expr->is<Product>());
    REQUIRE(sum_expr->as<Product>().scalar() == 2);

    // 2 x^{1/2} + 3 x^{1/2} = 5 x^{1/2}
    auto s1 = ex<Constant>(2) * ex<Power>(vx, rational{1, 2});
    auto s2 = ex<Constant>(3) * ex<Power>(vx, rational{1, 2});
    auto sum2 = s1 + s2;
    simplify(sum2);
    REQUIRE(sum2->is<Product>());
    REQUIRE(sum2->as<Product>().scalar() == 5);
  }
}
