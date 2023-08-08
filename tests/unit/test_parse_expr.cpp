#include "catch.hpp"

#include <SeQuant/core/parse/regex_sequant.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <boost/regex.hpp>
#include <locale>

// parse a non-symmetric tensor
sequant::ExprPtr parse(std::wstring_view raw) {
  return sequant::parse_expr(raw, sequant::Symmetry::nonsymm);
}

TEST_CASE("TEST_REGEX", "[parse_expr]") {
  using namespace sequant;

  // set up regex objects
  // ==============
  std::locale prev_locale{};
  std::locale::global(std::locale::classic());
  auto const rgx_label = boost::wregex{parse::regex_patterns::label().data()};
  auto const rgx_pure_index =
      boost::wregex{parse::regex_patterns::pure_index()};
  auto const rgx_index = boost::wregex{parse::regex_patterns::index_capture(),
                                       boost::regex::nosubs};
  auto const rgx_tensor_exp = boost::wregex{
      parse::regex_patterns::tensor_expanded().data(), boost::regex::nosubs};
  auto const rgx_tensor_terse = boost::wregex{
      parse::regex_patterns::tensor_terse().data(), boost::regex::nosubs};
  auto const rgx_ratio = boost::wregex{
      parse::regex_patterns::abs_real_frac().data(), boost::regex::nosubs};
  std::locale::global(prev_locale);
  // done setting up regex objects.
  // --------------

  SECTION("Labels") {
    for (auto const& l : std::initializer_list<std::wstring>{
             L"i", L"i1", L"a", L"a12", L"t⁔1", L"t¹", L"i↑", L"a↓", L"x⁺",
             L"x⁻", L"xₐ", L"xₑ", L"x₊", L"x₋", L"x₌"})
      REQUIRE(boost::regex_match(l, rgx_label));
  }

  SECTION("Index") {
    // pure index
    for (auto const& i : std::initializer_list<std::wstring>{
             L"i1", L"i_1", L"a1", L"a_1", L"i⁺_12", L"i⁻_34", L"i↑_0",
             L"a↓_24"}) {
      REQUIRE(boost::regex_match(i, rgx_pure_index));
      REQUIRE(boost::regex_match(i, rgx_index));
    }

    // index with proto indices
    for (auto const& i : std::initializer_list<std::wstring>{
             L"i1<a1 , a2>", L"i_1< a_1, a_2 >", L"i1 < a1,a2 >"})
      REQUIRE(boost::regex_match(i, rgx_index));

    // invalid index spec
    for (auto const& i :
         std::initializer_list<std::wstring>{L"i", L"a", L"i1<>", L"i_"})
      REQUIRE(!boost::regex_match(i, rgx_index));
  }

  SECTION("Tensor") {
    REQUIRE(boost::regex_match(L"t_{i1}^{a1}", rgx_tensor_exp));
    REQUIRE(boost::regex_match(L"t_{i1,i2}^{a1,a2}", rgx_tensor_exp));
    REQUIRE(boost::regex_match(L"t^{i1,i2}_{a1,a2}", rgx_tensor_exp));
    REQUIRE(boost::regex_match(L"t{i1;a1}", rgx_tensor_terse));
    REQUIRE(boost::regex_match(L"t{i1,i2;a1,a2}", rgx_tensor_terse));
    REQUIRE(boost::regex_match(L"t{ i1<a1, a2>, i2<a1, a2>; a1, a2 }",
                               rgx_tensor_terse));
  }

  SECTION("Ratio") {
    for (auto const& r : std::initializer_list<std::wstring>{
             L"1", L"123", L"1.", L"01.00", L"0 / 10", L"0.5/.25", L".4",
             L".1/.2"})
      REQUIRE(boost::regex_match(r, rgx_ratio));
  }

  SECTION("sequant::Variable") {
    // sequant variable is just a label
    for (auto const& v : std::initializer_list<std::wstring>{
             L"a", L"α", L"b", L"β", L"γ", L"λ", L"δ"})
      REQUIRE(boost::regex_match(v, rgx_label));
  }
}

TEST_CASE("TEST_PARSE_EXPR", "[parse_expr]") {
  using namespace sequant;
  SECTION("Tensor") {
    auto expr = parse_expr(L"t{i1;a1}");
    REQUIRE(expr->is<Tensor>());
    REQUIRE(expr->as<Tensor>().label() == L"t");
    REQUIRE(expr->as<Tensor>().bra().size() == 1);
    REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
    REQUIRE(expr->as<Tensor>().ket().size() == 1);
    REQUIRE(expr->as<Tensor>().ket().at(0) == L"a_1");

    REQUIRE(*expr == *parse_expr(L"t_{i1}^{a1}"));
    REQUIRE(*expr == *parse_expr(L"t^{a1}_{i1}"));
    REQUIRE(*expr == *parse_expr(L"t{i_1; a_1}"));
    REQUIRE(*expr == *parse_expr(L"t_{i_1}^{a_1}"));

    expr = parse_expr(L"t{i1,i2;a1,a2}");
    REQUIRE(expr->as<Tensor>().bra().size() == 2);
    REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
    REQUIRE(expr->as<Tensor>().bra().at(1).label() == L"i_2");
    REQUIRE(expr->as<Tensor>().ket().size() == 2);
    REQUIRE(expr->as<Tensor>().ket().at(0).label() == L"a_1");
    REQUIRE(expr->as<Tensor>().ket().at(1).label() == L"a_2");

    REQUIRE(*expr == *parse_expr(L"+t{i1, i2; a1, a2}"));
    REQUIRE(parse_expr(L"-t{i1;a1}")->is<Product>());
    REQUIRE(*expr == *parse_expr(L"t{\ti1, \ti2; \na1,\t a2 \t}"));
    REQUIRE_NOTHROW(parse_expr(L"t⁔1{a1;i1}"));
    REQUIRE_NOTHROW(parse_expr(L"t¹{a1;i1}"));
  }

  SECTION("Tensor with symmetry annotation") {
    auto expr1 = parse_expr(L"t{a1;i1}:A");
    auto expr2 = parse_expr(L"t{a1;i1}:S");
    auto expr3 = parse_expr(L"t{a1;i1}:N");
    REQUIRE(expr1->as<Tensor>().symmetry() == sequant::Symmetry::antisymm);
    REQUIRE(expr2->as<Tensor>().symmetry() == sequant::Symmetry::symm);
    REQUIRE(expr3->as<Tensor>().symmetry() == sequant::Symmetry::nonsymm);
  }

  SECTION("Constant") {
    REQUIRE(parse_expr(L"1/2")->is<Constant>());
    REQUIRE(parse_expr(L"0/2")->is<Constant>());
    REQUIRE(!parse_expr(L"-1/2")->is<Constant>());
    REQUIRE(!parse_expr(L"-0/2")->is<Constant>());
  }

  SECTION("Product") {
    auto expr = parse_expr(L"-1/2 g{i2,i3; i1,a2} t{a1,a2; i2,i3}");
    REQUIRE(expr->is<Product>());

    auto const& prod = expr->as<Product>();
    REQUIRE(prod.scalar() == rational{-1, 2});
    REQUIRE(*prod.factor(0) == *parse_expr(L"g_{i_2, i_3}^{i_1, a_2}"));
    REQUIRE(*prod.factor(1) == *parse_expr(L"t^{i2, i3}_{a1, a2}"));
    REQUIRE(parse_expr(L"-1/2 * δ * t{i1;a1}") ==
            parse_expr(L"-1/2  δ  t{i1;a1}"));
    auto const prod2 = parse_expr(L"-1/2 * δ * γ * t{i1;a1}")->as<Product>();
    REQUIRE(prod2.scalar() == rational{-1,2});
    REQUIRE(prod2.factor(0) == ex<Variable>(L"δ"));
    REQUIRE(prod2.factor(1) == ex<Variable>(L"γ"));
    REQUIRE(prod2.factor(2)->is<Tensor>());
  }

  SECTION("Sum") {
    auto expr = parse_expr(
        L"f{a1;i1}"
        "- 1/2*g{i2,a1; a2,a3}t{a2,a3; i1,i2}");
    REQUIRE(expr->is<Sum>());

    auto const& sum = expr->as<Sum>();
    REQUIRE(*sum.summand(0) == *parse_expr(L"f{a1;i1}"));
    REQUIRE(*sum.summand(1) ==
            *parse_expr(L"- 1/2 * g{i2,a1; a2,a3} * t{a2,a3; i1,i2}"));
  }

  SECTION("Parentheses") {
    auto expr1 =
        parse_expr(L"-1/2 g{i2,i3; a2,a3} * ( t{a1,a3; i2,i3} * t{a2;i1} )");
    REQUIRE(expr1->is<Product>());

    auto const& prod1 = expr1->as<Product>();
    REQUIRE(prod1.size() == 2);
    REQUIRE(prod1.scalar() == rational{-1, 2});
    REQUIRE(prod1.factor(0)->is<Tensor>());
    REQUIRE(prod1.factor(1)->is<Product>());
    REQUIRE(prod1.factor(1)->size() == 2);

    auto expr2 = parse_expr(
        L"(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (t{a2;i1})");
    REQUIRE(expr2->is<Product>());

    auto const& prod2 = expr2->as<Product>();
    REQUIRE(prod2.size() == 2);
    REQUIRE(prod2.scalar() == rational{-1, 2});
    REQUIRE(prod2.factor(0)->is<Product>());
    REQUIRE(*prod2.factor(0)->at(0) == *parse_expr(L"g{i2,i3; a2,a3}"));
    REQUIRE(*prod2.factor(0)->at(1) == *parse_expr(L"t{a1,a3; i2,i3}"));
    REQUIRE(*prod2.factor(1) == *parse_expr(L"t{a2;i1}"));

    auto expr3 = parse_expr(
        L"(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (1/2) * ((t{a2;i1}))");
    REQUIRE(expr3->is<Product>());

    auto const& prod3 = expr3->as<Product>();
    REQUIRE(prod3.size() == 2);
    REQUIRE(prod3.scalar() == rational{-1, 4});
    REQUIRE(prod3.factor(0)->is<Product>());
    REQUIRE(*prod3.factor(0)->at(0) == *parse_expr(L"g{i2,i3; a2,a3}"));
    REQUIRE(*prod3.factor(0)->at(1) == *parse_expr(L"t{a1,a3; i2,i3}"));
    REQUIRE(*prod3.factor(1) == *parse_expr(L"t{a2;i1}"));
  }

  SECTION("Mixed") {
    auto expr = parse_expr(
        L"1/4 g{a1,a2; i1,i2}"
        "+ 1/4 g{i3,i4; a3,a4} (t{a3;i1} * t{a4;i2}) * (t{a1;i3} * t{a2;i4})");
    REQUIRE(expr->is<Sum>());
    auto const& sum = expr->as<Sum>();

    REQUIRE(sum.size() == 2);

    REQUIRE(sum.summand(0)->is<Product>());
    REQUIRE(sum.summand(0)->as<Product>().scalar() == rational{1, 4});
    REQUIRE(sum.summand(0)->size() == 1);
    REQUIRE(*sum.summand(0)->at(0) == *parse_expr(L"g{a1,a2; i1,i2}"));

    REQUIRE(sum.summand(1)->is<Product>());
    auto const& prod = sum.summand(1)->as<Product>();

    REQUIRE(prod.scalar() == rational{1, 4});
    REQUIRE(prod.size() == 3);
    REQUIRE(*prod.factor(0) == *parse_expr(L"g{i3,i4; a3,a4}"));

    REQUIRE(prod.factor(1)->is<Product>());
    REQUIRE(*prod.factor(1)->at(0) == *parse_expr(L"t{a3;i1}"));
    REQUIRE(*prod.factor(1)->at(1) == *parse_expr(L"t{a4;i2}"));

    REQUIRE(prod.factor(2)->is<Product>());
    REQUIRE(*prod.factor(2)->at(0) == *parse_expr(L"t{a1;i3}"));
    REQUIRE(*prod.factor(2)->at(1) == *parse_expr(L"t{a2;i4}"));
  }
}
