#include "catch.hpp"

#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <algorithm>
#include <locale>
#include <sstream>
#include <string>

namespace Catch {
template <>
struct StringMaker<sequant::ParseError> {
  static std::string convert(const sequant::ParseError& error) {
    return "ParseError{offset: " + std::to_string(error.offset) +
           ", length: " + std::to_string(error.length) + ", what(): '" +
           error.what() + "'}";
  }
};
}  // namespace Catch

struct ParseErrorMatcher : Catch::MatcherBase<sequant::ParseError> {
  std::size_t offset;
  std::size_t length;
  std::string messageFragment;

  ParseErrorMatcher(std::size_t offset, std::size_t length,
                    std::string messageFragment = "")
      : offset(offset),
        length(length),
        messageFragment(std::move(messageFragment)) {}

  bool match(const sequant::ParseError& exception) const override {
    if (exception.offset != offset) {
      return false;
    }
    if (exception.length != length) {
      return false;
    }
    if (!messageFragment.empty()) {
      std::string message(exception.what());
      auto iter = message.find(messageFragment);

      return iter != std::string::npos;
    }
    return true;
  }

  std::string describe() const override {
    return "-- expected {offset: " + std::to_string(offset) +
           ", length: " + std::to_string(length) +
           (messageFragment.empty()
                ? std::string("}")
                : ", what() references '" + messageFragment + "'}");
  }
};

ParseErrorMatcher parseErrorMatches(std::size_t offset, std::size_t length,
                                    std::string messageFragment = "") {
  return ParseErrorMatcher{offset, length, std::move(messageFragment)};
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
    REQUIRE_NOTHROW(parse_expr(L"α{a1;i1}"));
    REQUIRE_NOTHROW(parse_expr(L"γ_1{a1;i1}"));
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
    REQUIRE(parse_expr(L"-1/2")->is<Constant>());
    REQUIRE(parse_expr(L"-0/2")->is<Constant>());
    REQUIRE(parse_expr(L"1")->is<Constant>());
    REQUIRE(parse_expr(L"123")->is<Constant>());
    REQUIRE(parse_expr(L"1.")->is<Constant>());
    REQUIRE(parse_expr(L"01.00")->is<Constant>());
    REQUIRE(parse_expr(L"0 / 10")->is<Constant>());
    REQUIRE(parse_expr(L"0.5/0.25")->is<Constant>());
    REQUIRE(parse_expr(L".4")->is<Constant>());
  }

  SECTION("Variable") {
    // sequant variable is just a label followed by an optional ^*
    // to denote if the variable is conjugated
    REQUIRE(parse_expr(L"a")->is<Variable>());
    REQUIRE(parse_expr(L"α")->is<Variable>());
    REQUIRE(parse_expr(L"β")->is<Variable>());
    REQUIRE(parse_expr(L"γ")->is<Variable>());
    REQUIRE(parse_expr(L"λ")->is<Variable>());
    REQUIRE(parse_expr(L"δ")->is<Variable>());
    REQUIRE(parse_expr(L"a^*")->is<Variable>());
    REQUIRE(parse_expr(L"α^*")->is<Variable>());
    REQUIRE(parse_expr(L"β^*")->is<Variable>());
    REQUIRE(parse_expr(L"b^*")->is<Variable>());
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
    REQUIRE(prod2.scalar() == rational{-1, 2});
    REQUIRE(prod2.factor(0) == ex<Variable>(L"δ"));
    REQUIRE(prod2.factor(1) == ex<Variable>(L"γ"));
    REQUIRE(prod2.factor(2)->is<Tensor>());
  }

  SECTION("Sum") {
    auto expr1 = parse_expr(
        L"f{a1;i1}"
        "- 1/2*g{i2,a1; a2,a3}t{a2,a3; i1,i2}");
    REQUIRE(expr1->is<Sum>());

    auto const& sum1 = expr1->as<Sum>();
    REQUIRE(*sum1.summand(0) == *parse_expr(L"f{a1;i1}"));
    REQUIRE(*sum1.summand(1) ==
            *parse_expr(L"- 1/2 * g{i2,a1; a2,a3} * t{a2,a3; i1,i2}"));

    auto expr2 = parse_expr(L"a - 4");
    REQUIRE(expr2->is<Sum>());

    auto const& sum2 = expr2->as<Sum>();
    REQUIRE(*sum2.summand(0) == *parse_expr(L"a"));
    REQUIRE(*sum2.summand(1) == *parse_expr(L"-4"));
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

    auto expr4 = parse_expr(L"1/2 (a + b) * c");
    REQUIRE(expr4->is<Product>());

    const auto& prod4 = expr4->as<Product>();
    REQUIRE(prod4.size() == 2);
    REQUIRE(prod4.scalar() == rational{1, 2});
    REQUIRE(prod4.factor(1)->is<Variable>());
    REQUIRE(prod4.factor(1)->as<Variable>().label() == L"c");
    REQUIRE(prod4.factor(0)->is<Sum>());
    const auto& nestedSum = prod4.factor(0)->as<Sum>();
    REQUIRE(nestedSum.size() == 2);
    REQUIRE(nestedSum.summand(0)->is<Variable>());
    REQUIRE(nestedSum.summand(0)->as<Variable>().label() == L"a");
    REQUIRE(nestedSum.summand(1)->is<Variable>());
    REQUIRE(nestedSum.summand(1)->as<Variable>().label() == L"b");
  }

  SECTION("Mixed") {
    auto expr = parse_expr(
        L"0.25 g{a1,a2; i1,i2}"
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

  SECTION("Empty input") { REQUIRE(parse_expr(L"") == nullptr); }

  SECTION("Error handling") {
    SECTION("Exception type") {
      std::vector<std::wstring> inputs = {L"t^",
                                          L"a + + b"
                                          L"1/t"
                                          L"T{}"
                                          L"T^{i1}{a1}"};

      for (const std::wstring& current : inputs) {
        REQUIRE_THROWS_AS(parse_expr(current), ParseError);
      }
    }

    SECTION("Invalid index") {
      REQUIRE_THROWS_MATCHES(parse_expr(L"t{i1<az1>;}"), ParseError,
                             parseErrorMatches(5, 3, "proto"));
      REQUIRE_THROWS_MATCHES(parse_expr(L"t{i1;az3}"), ParseError,
                             parseErrorMatches(5, 3, "Unknown index space"));
    }

    SECTION("Invalid symmetry") {
      REQUIRE_THROWS_MATCHES(
          parse_expr(L"t{i1;a3}:P"), ParseError,
          parseErrorMatches(9, 1, "Invalid symmetry specifier"));
    }
  }
}

TEST_CASE("TEST_DEPARSE_EXPR", "[parse_expr]") {
  using namespace sequant;

  std::vector<std::wstring> expressions = {
      L"t{a_1,a_2;a_3,a_4}:N",
      L"42",
      L"1/2",
      L"-1/4 t{a_1,i_1<a_1>;a_2,i_2}:S",
      L"a + b - 4 specialVariable",
      L"variable + A{a_1;i_1}:N * B{i_1;a_1}:A",
      L"1/2 (a + b) * c"};

  for (const std::wstring& current : expressions) {
    ExprPtr expression = parse_expr(current);

    REQUIRE(deparse_expr(expression, true) == current);
  }
}
