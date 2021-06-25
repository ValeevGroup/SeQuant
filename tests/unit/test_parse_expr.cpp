#include "catch.hpp"

#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>

// parse a non-symmetric tensor
sequant::ExprPtr  parse(std::wstring_view raw) {
  return sequant::parse_expr(raw, sequant::Symmetry::nonsymm);
}

TEST_CASE("TEST_PARSE_EXPR", "[parse_expr]"){
  using namespace sequant;
  SECTION("Tensor"){
    auto expr = parse(L"t{i1;a1}");
    REQUIRE(expr->is<Tensor>());
    REQUIRE(expr->as<Tensor>().label() == L"t");
    REQUIRE(expr->as<Tensor>().bra().size() == 1);
    REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
    REQUIRE(expr->as<Tensor>().ket().size() == 1);
    REQUIRE(expr->as<Tensor>().ket().at(0) == L"a_1");

    REQUIRE(*expr == *parse(L"t_{i1}^{a1}"));
    REQUIRE(*expr == *parse(L"t^{a1}_{i1}"));
    REQUIRE(*expr == *parse(L"t{i_1; a_1}"));
    REQUIRE(*expr == *parse(L"t_{i_1}^{a_1}"));

     expr = parse(L"t{i1,i2;a1,a2}");
    REQUIRE(expr->as<Tensor>().bra().size() == 2);
    REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
    REQUIRE(expr->as<Tensor>().bra().at(1).label() == L"i_2");
    REQUIRE(expr->as<Tensor>().ket().size() == 2);
    REQUIRE(expr->as<Tensor>().ket().at(0).label() == L"a_1");
    REQUIRE(expr->as<Tensor>().ket().at(1).label() == L"a_2");

    REQUIRE(*expr == *parse(L"+t{i1, i2; a1, a2}"));
    REQUIRE(parse(L"-t{i1;a1}")->is<Product>());
    REQUIRE(*expr == *parse(L"t{\ti1, \ti2; \na1,\t a2 \t}"));
  }

  SECTION("Constant") {
      REQUIRE(parse(L"1/2")->is<Constant>());
      REQUIRE(parse(L"0.1")->is<Constant>());
      REQUIRE(parse(L"0.1/0.2")->is<Constant>());
      REQUIRE(parse(L".1/.2")->is<Constant>());
  }

  SECTION("Product") {
    auto expr = parse(L"-1/2 g{i2,i3; i1,a2} t{a1,a2; i2,i3}");
    REQUIRE(expr->is<Product>());

    auto const& prod = expr->as<Product>();
    REQUIRE(prod.scalar() == -0.5);
    REQUIRE(*prod.factor(0) == *parse(L"g_{i_2, i_3}^{i_1, a_2}"));
    REQUIRE(*prod.factor(1) == *parse(L"t^{i2, i3}_{a1, a2}"));
  }

  SECTION("Sum") {
    auto expr = parse(L"f{a1;i1}"
                "- 1/2*g{i2,a1; a2,a3}t{a2,a3; i1,i2}");
    REQUIRE(expr->is<Sum>());

    auto const& sum = expr->as<Sum>();
    REQUIRE(*sum.summand(0) == *parse(L"f{a1;i1}"));
    REQUIRE(*sum.summand(1) == *parse(L"- 1/2 * g{i2,a1; a2,a3} * t{a2,a3; i1,i2}"));
  }

  SECTION("Parentheses"){
    auto expr1 = parse(L"-1/2 g{i2,i3; a2,a3} * ( t{a1,a3; i2,i3} * t{a2;i1} )");
    REQUIRE(expr1->is<Product>());

    auto const& prod1 = expr1->as<Product>();
    REQUIRE(prod1.size() == 2);
    REQUIRE(prod1.scalar() == -0.5);
    REQUIRE(prod1.factor(0)->is<Tensor>());
    REQUIRE(prod1.factor(1)->is<Product>());
    REQUIRE(prod1.factor(1)->size() == 2);

    auto expr2 = parse(L"(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (t{a2;i1})");
    REQUIRE(expr2->is<Product>());

    auto const& prod2 = expr2->as<Product>();
    REQUIRE(prod2.size() == 2);
    REQUIRE(prod2.scalar() == -0.5);
    REQUIRE(prod2.factor(0)->is<Product>());
    REQUIRE(*prod2.factor(0)->at(0) == *parse(L"g{i2,i3; a2,a3}"));
    REQUIRE(*prod2.factor(0)->at(1) == *parse(L"t{a1,a3; i2,i3}"));
    REQUIRE(*prod2.factor(1) == *parse(L"t{a2;i1}"));

    auto expr3 = parse(L"(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (1/2) * ((t{a2;i1}))");
    REQUIRE(expr3->is<Product>());

    auto const& prod3 = expr3->as<Product>();
    REQUIRE(prod3.size() == 2);
    REQUIRE(prod3.scalar() == -0.25);
    REQUIRE(prod3.factor(0)->is<Product>());
    REQUIRE(*prod3.factor(0)->at(0) == *parse(L"g{i2,i3; a2,a3}"));
    REQUIRE(*prod3.factor(0)->at(1) == *parse(L"t{a1,a3; i2,i3}"));
    REQUIRE(*prod3.factor(1) == *parse(L"t{a2;i1}"));
  }

  SECTION("Mixed") {
    auto expr = parse(L"1/4 g{a1,a2; i1,i2}"
        "+ 1/4 g{i3,i4; a3,a4} (t{a3;i1} * t{a4;i2}) * (t{a1;i3} * t{a2;i4})");
    REQUIRE(expr->is<Sum>());
    auto const& sum = expr->as<Sum>();

    REQUIRE(sum.size() == 2);

    REQUIRE(sum.summand(0)->is<Product>());
    REQUIRE(sum.summand(0)->as<Product>().scalar() == 0.25);
    REQUIRE(sum.summand(0)->size() == 1);
    REQUIRE(*sum.summand(0)->at(0) == *parse(L"g{a1,a2; i1,i2}"));

    REQUIRE(sum.summand(1)->is<Product>());
    auto const& prod = sum.summand(1)->as<Product>();

    REQUIRE(prod.scalar() == 0.25);
    REQUIRE(prod.size() == 3);
    REQUIRE(*prod.factor(0) == *parse(L"g{i3,i4; a3,a4}"));

    REQUIRE(prod.factor(1)->is<Product>());
    REQUIRE(*prod.factor(1)->at(0) == *parse(L"t{a3;i1}"));
    REQUIRE(*prod.factor(1)->at(1) == *parse(L"t{a4;i2}"));

    REQUIRE(prod.factor(2)->is<Product>());
    REQUIRE(*prod.factor(2)->at(0) == *parse(L"t{a1;i3}"));
    REQUIRE(*prod.factor(2)->at(1) == *parse(L"t{a2;i4}"));
  }

}