#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/tensor.hpp>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <algorithm>
#include <cstddef>
#include <locale>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

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

struct ParseErrorMatcher : Catch::Matchers::MatcherBase<sequant::ParseError> {
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

TEST_CASE("parsing", "[parse]") {
  SECTION("parse_expr") {
    using namespace sequant;

    auto ctx_resetter = set_scoped_default_context(
        Context(mbpt::make_sr_spaces(), Vacuum::SingleProduct));

    SECTION("Scalar tensor") {
      auto expr = parse_expr(L"t{}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().empty());
      REQUIRE(expr->as<Tensor>().ket().empty());
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == parse_expr(L"t{;}"));
      REQUIRE(expr == parse_expr(L"t{;;}"));
      REQUIRE(expr == parse_expr(L"t^{}_{}"));
      REQUIRE(expr == parse_expr(L"t_{}^{}"));
    }
    SECTION("Tensor") {
      auto expr = parse_expr(L"t{i1;a1}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().label() == L"t");
      REQUIRE(expr->as<Tensor>().bra().size() == 1);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
      REQUIRE(expr->as<Tensor>().ket().size() == 1);
      REQUIRE(expr->as<Tensor>().ket().at(0) == L"a_1");
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == parse_expr(L"t_{i1}^{a1}"));
      REQUIRE(expr == parse_expr(L"t^{a1}_{i1}"));
      REQUIRE(expr == parse_expr(L"t{i_1; a_1}"));
      REQUIRE(expr == parse_expr(L"t{i_1; a_1;}"));
      REQUIRE(expr == parse_expr(L"t_{i_1}^{a_1}"));

      expr = parse_expr(L"t{i1,i2;a1,a2}");
      REQUIRE(expr->as<Tensor>().bra().size() == 2);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
      REQUIRE(expr->as<Tensor>().bra().at(1).label() == L"i_2");
      REQUIRE(expr->as<Tensor>().ket().size() == 2);
      REQUIRE(expr->as<Tensor>().ket().at(0).label() == L"a_1");
      REQUIRE(expr->as<Tensor>().ket().at(1).label() == L"a_2");
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == parse_expr(L"+t{i1, i2; a1, a2}"));
      REQUIRE(parse_expr(L"-t{i1;a1}")->is<Product>());
      REQUIRE(expr == parse_expr(L"t{\ti1, \ti2; \na1,\t a2 \t}"));

      // Tensor labels including underscores
      REQUIRE(parse_expr(L"T_1{i_1;a_1}")->as<Tensor>().label() == L"T_1");

      // "Non-standard" tensor labels
      REQUIRE(parse_expr(L"α{a1;i1}")->as<Tensor>().label() == L"α");
      REQUIRE(parse_expr(L"γ_1{a1;i1}")->as<Tensor>().label() == L"γ_1");
      REQUIRE(parse_expr(L"t⁔1{a1;i1}")->as<Tensor>().label() == L"t⁔1");
      REQUIRE(parse_expr(L"t¹{a1;i1}")->as<Tensor>().label() == L"t¹");
      REQUIRE(parse_expr(L"t⁸{a1;i1}")->as<Tensor>().label() == L"t⁸");
      REQUIRE(parse_expr(L"t⁻{a1;i1}")->as<Tensor>().label() == L"t⁻");
      REQUIRE(parse_expr(L"tₐ{a1;i1}")->as<Tensor>().label() == L"tₐ");
      REQUIRE(parse_expr(L"t₋{a1;i1}")->as<Tensor>().label() == L"t₋");
      REQUIRE(parse_expr(L"t₌{a1;i1}")->as<Tensor>().label() == L"t₌");
      REQUIRE(parse_expr(L"t↓{a1;i1}")->as<Tensor>().label() == L"t↓");
      REQUIRE(parse_expr(L"t↑{a1;i1}")->as<Tensor>().label() == L"t↑");

      // "Non-standard" index names
      auto expr1 = parse_expr(L"t{a↓1;i↑1}");
      REQUIRE(expr1->as<Tensor>().bra().at(0).label() == L"a↓_1");
      REQUIRE(expr1->as<Tensor>().ket().at(0).label() == L"i↑_1");

      // Auxiliary indices
      expr = parse_expr(L"t{;;i1}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().empty());
      REQUIRE(expr->as<Tensor>().ket().empty());
      REQUIRE(expr->as<Tensor>().aux().size() == 1);
      REQUIRE(expr->as<Tensor>().aux()[0].label() == L"i_1");

      // All index groups at once
      expr = parse_expr(L"t{i1,i2;a1;x1,x2}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().size() == 2);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == L"i_1");
      REQUIRE(expr->as<Tensor>().bra().at(1).label() == L"i_2");
      REQUIRE(expr->as<Tensor>().ket().size() == 1);
      REQUIRE(expr->as<Tensor>().ket().at(0).label() == L"a_1");
      REQUIRE(expr->as<Tensor>().aux().size() == 2);
      REQUIRE(expr->as<Tensor>().aux().at(0).label() == L"x_1");
      REQUIRE(expr->as<Tensor>().aux().at(1).label() == L"x_2");
    }

    SECTION("Tensor with symmetry annotation") {
      auto expr1 = parse_expr(L"t{a1;i1}:A");
      auto expr2 = parse_expr(L"t{a1;i1}:S-C");
      auto expr3 = parse_expr(L"t{a1;i1}:N-S-N");

      const Tensor& t1 = expr1->as<Tensor>();
      const Tensor& t2 = expr2->as<Tensor>();
      const Tensor& t3 = expr3->as<Tensor>();

      REQUIRE(t1.symmetry() == Symmetry::antisymm);

      REQUIRE(t2.symmetry() == Symmetry::symm);
      REQUIRE(t2.braket_symmetry() == BraKetSymmetry::conjugate);

      REQUIRE(t3.symmetry() == Symmetry::nonsymm);
      REQUIRE(t3.braket_symmetry() == BraKetSymmetry::symm);
      REQUIRE(t3.particle_symmetry() == ParticleSymmetry::nonsymm);
    }

    SECTION("NormalOperator") {
      {
        using NOp = FNOperator;
        auto expr = parse_expr(L"a{i1;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[0]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == L"a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == L"i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::Physical);
      }
      {
        using NOp = FNOperator;
        auto expr = parse_expr(L"ã{i1;}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[1]);
        REQUIRE(expr->as<NOp>().creators().size() == 0);
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == L"i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::SingleProduct);
      }

      {
        using NOp = BNOperator;
        auto expr = parse_expr(L"b{i1;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[0]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == L"a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == L"i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::Physical);
      }
      {
        using NOp = BNOperator;
        auto expr = parse_expr(L"b̃{;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[1]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == L"a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 0);
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::SingleProduct);
      }
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
      // SeQuant variable is just a label followed by an optional ^*
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
      REQUIRE(parse_expr(L"b^*")->as<Variable>().conjugated());
      REQUIRE(parse_expr(L"b^*")->as<Variable>().label() == L"b");
    }

    SECTION("Product") {
      auto expr = parse_expr(L"-1/2 g{i2,i3; i1,a2} t{a1,a2; i2,i3}");
      REQUIRE(expr->is<Product>());

      auto const& prod = expr->as<Product>();
      REQUIRE(prod.scalar() == rational{-1, 2});
      REQUIRE(prod.factor(0) == parse_expr(L"g_{i_2, i_3}^{i_1, a_2}"));
      REQUIRE(prod.factor(1) == parse_expr(L"t^{i2, i3}_{a1, a2}"));
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
      REQUIRE(sum1.summand(0) == parse_expr(L"f{a1;i1}"));
      REQUIRE(sum1.summand(1) ==
              parse_expr(L"- 1/2 * g{i2,a1; a2,a3} * t{a2,a3; i1,i2}"));

      auto expr2 = parse_expr(L"a - 4");
      REQUIRE(expr2->is<Sum>());

      auto const& sum2 = expr2->as<Sum>();
      REQUIRE(sum2.summand(0) == parse_expr(L"a"));
      REQUIRE(sum2.summand(1) == parse_expr(L"-4"));
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
      REQUIRE(prod2.factor(0)->at(0) == parse_expr(L"g{i2,i3; a2,a3}"));
      REQUIRE(prod2.factor(0)->at(1) == parse_expr(L"t{a1,a3; i2,i3}"));
      REQUIRE(prod2.factor(1) == parse_expr(L"t{a2;i1}"));

      auto expr3 = parse_expr(
          L"(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (1/2) * "
          L"((t{a2;i1}))");
      REQUIRE(expr3->is<Product>());

      auto const& prod3 = expr3->as<Product>();
      REQUIRE(prod3.size() == 2);
      REQUIRE(prod3.scalar() == rational{-1, 4});
      REQUIRE(prod3.factor(0)->is<Product>());
      REQUIRE(prod3.factor(0)->at(0) == parse_expr(L"g{i2,i3; a2,a3}"));
      REQUIRE(prod3.factor(0)->at(1) == parse_expr(L"t{a1,a3; i2,i3}"));
      REQUIRE(prod3.factor(1) == parse_expr(L"t{a2;i1}"));

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
          "+ 1/4 g{i3,i4; a3,a4} (t{a3;i1} * t{a4;i2}) * (t{a1;i3} * "
          "t{a2;i4})");
      REQUIRE(expr->is<Sum>());
      auto const& sum = expr->as<Sum>();

      REQUIRE(sum.size() == 2);

      REQUIRE(sum.summand(0)->is<Product>());
      REQUIRE(sum.summand(0)->as<Product>().scalar() == rational{1, 4});
      REQUIRE(sum.summand(0)->size() == 1);
      REQUIRE(sum.summand(0)->at(0) == parse_expr(L"g{a1,a2; i1,i2}"));

      REQUIRE(sum.summand(1)->is<Product>());
      auto const& prod = sum.summand(1)->as<Product>();

      REQUIRE(prod.scalar() == rational{1, 4});
      REQUIRE(prod.size() == 3);
      REQUIRE(prod.factor(0) == parse_expr(L"g{i3,i4; a3,a4}"));

      REQUIRE(prod.factor(1)->is<Product>());
      REQUIRE(prod.factor(1)->at(0) == parse_expr(L"t{a3;i1}"));
      REQUIRE(prod.factor(1)->at(1) == parse_expr(L"t{a4;i2}"));

      REQUIRE(prod.factor(2)->is<Product>());
      REQUIRE(prod.factor(2)->at(0) == parse_expr(L"t{a1;i3}"));
      REQUIRE(prod.factor(2)->at(1) == parse_expr(L"t{a2;i4}"));
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

  SECTION("parse_result_expr") {
    using namespace sequant;

    SECTION("constant") {
      ResultExpr result = parse_result_expr(L"A = 3");

      REQUIRE(result.has_label());
      REQUIRE(result.label() == L"A");
      REQUIRE(result.bra().empty());
      REQUIRE(result.ket().empty());
      REQUIRE(result.symmetry() == Symmetry::nonsymm);
      REQUIRE(result.braket_symmetry() == BraKetSymmetry::nonsymm);
      REQUIRE(result.particle_symmetry() == ParticleSymmetry::nonsymm);

      REQUIRE(result.expression().is<Constant>());
      REQUIRE(result.expression().as<Constant>().value<int>() == 3);
    }
    SECTION("contraction") {
      ResultExpr result =
          parse_result_expr(L"R{i1,i2;e1,e2}:A = f{e2;e3} t{e1,e3;i1,i2}");

      REQUIRE(result.has_label());
      REQUIRE(result.label() == L"R");
      REQUIRE(result.bra().size() == 2);
      REQUIRE(result.bra()[0].full_label() == L"i_1");
      REQUIRE(result.bra()[1].full_label() == L"i_2");
      REQUIRE(result.ket().size() == 2);
      REQUIRE(result.ket()[0].full_label() == L"e_1");
      REQUIRE(result.ket()[1].full_label() == L"e_2");
      REQUIRE(result.symmetry() == Symmetry::antisymm);
      REQUIRE(result.braket_symmetry() ==
              get_default_context().braket_symmetry());
      REQUIRE(result.particle_symmetry() == ParticleSymmetry::symm);

      REQUIRE(result.expression().is<Product>());
      const Product& prod = result.expression().as<Product>();
      REQUIRE(prod.size() == 2);
      REQUIRE(prod.factor(0).is<Tensor>());
      REQUIRE(prod.factor(0).as<Tensor>().label() == L"f");
      REQUIRE(prod.factor(1).is<Tensor>());
      REQUIRE(prod.factor(1).as<Tensor>().label() == L"t");
    }
  }

  SECTION("deparse") {
    using namespace sequant;

    std::vector<std::wstring> expressions = {
        L"t{a_1,a_2;a_3,a_4}:N-C-S",
        L"42",
        L"1/2",
        L"-1/4 t{a_1,i_1<a_1>;a_2,i_2}:S-N-S",
        L"a + b - 4 specialVariable",
        L"variable + A{a_1;i_1}:N-N-S * B{i_1;a_1}:A-C-S",
        L"1/2 (a + b) * c",
        L"T1{}:N-N-N + T2{;;x_1}:N-N-N * T3{;;x_1}:N-N-N + T4{a_1;;x_2}:S-C-S "
        L"* "
        L"T5{;a_1;x_2}:S-S-S",
        L"q1 * q2^* * q3",
        L"1/2 ã{i_1;i_2} * b̃{i_3;i_4}"};

    for (const std::wstring& current : expressions) {
      ExprPtr expression = parse_expr(current);

      REQUIRE(deparse(expression, true) == current);
    }

    SECTION("result_expressions") {
      std::vector<std::wstring> expressions = {
          L"A = 5",
          L"A = g{i_1,i_2;e_1,e_2}:S-N-S * t{e_1,e_2;i_1,i_2}:N-N-S",
          L"R{i_1,i_2;e_1,e_2}:A-N-S = f{e_2;e_3}:A-N-S * "
          L"t{e_1,e_3;i_1,i_2}:A-N-S + "
          L"g{i_1,i_2;e_1,e_2}:A-N-S",
      };

      for (const std::wstring& current : expressions) {
        ResultExpr result = parse_result_expr(current);

        REQUIRE(deparse(result, true) == current);
      }
    }
  }
}
