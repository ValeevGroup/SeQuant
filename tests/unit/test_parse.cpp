#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/rational.hpp>

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

namespace Catch {
template <>
struct StringMaker<sequant::io::serialization::SerializationError> {
  static std::string convert(
      const sequant::io::serialization::SerializationError& error) {
    return "io::serialization::SerializationError{offset: " +
           std::to_string(error.offset) +
           ", length: " + std::to_string(error.length) + ", what(): '" +
           error.what() + "'}";
  }
};
}  // namespace Catch

struct SerializationErrorMatcher
    : Catch::Matchers::MatcherBase<
          sequant::io::serialization::SerializationError> {
  std::size_t offset;
  std::size_t length;
  std::string messageFragment;

  SerializationErrorMatcher(std::size_t offset, std::size_t length,
                            std::string messageFragment = "")
      : offset(offset),
        length(length),
        messageFragment(std::move(messageFragment)) {}

  bool match(const sequant::io::serialization::SerializationError& exception)
      const override {
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

SerializationErrorMatcher serializationErrorMatches(
    std::size_t offset, std::size_t length, std::string messageFragment = "") {
  return SerializationErrorMatcher{offset, length, std::move(messageFragment)};
}

TEST_CASE("serialization", "[serialization]") {
  SECTION("deserialize<ExprPtr>") {
    using namespace sequant;

    auto ctx = get_default_context();
    ctx.set(mbpt::make_sr_spaces());
    auto ctx_resetter = set_scoped_default_context(ctx);

    SECTION("default symmetries follow the Context") {
      // the deserializer sources unspecified tensor symmetries from the active
      // Context (unlike the programmatic Tensor ctor, which uses fixed
      // defaults -- see test_tensor "programmatic ctor defaults are
      // Context-independent"). Flip the Context column-symmetry default and
      // confirm parsed tensors follow it.
      auto ctx_symm = get_default_context();
      ctx_symm.set(ColumnSymmetry::Symm);
      {
        auto resetter = set_scoped_default_context(ctx_symm);
        REQUIRE(deserialize<ExprPtr>("t{i1,i2;a1,a2}")
                    ->as<Tensor>()
                    .column_symmetry() == ColumnSymmetry::Symm);
      }
      auto ctx_nonsymm = get_default_context();
      ctx_nonsymm.set(ColumnSymmetry::Nonsymm);
      {
        auto resetter = set_scoped_default_context(ctx_nonsymm);
        REQUIRE(deserialize<ExprPtr>("t{i1,i2;a1,a2}")
                    ->as<Tensor>()
                    .column_symmetry() == ColumnSymmetry::Nonsymm);
      }
    }

    SECTION("Scalar tensor") {
      auto expr = deserialize<ExprPtr>("t{}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().empty());
      REQUIRE(expr->as<Tensor>().ket().empty());
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == deserialize<ExprPtr>("t{;}"));
      REQUIRE(expr == deserialize<ExprPtr>("t{;;}"));
      REQUIRE(expr == deserialize<ExprPtr>("t^{}_{}"));
      REQUIRE(expr == deserialize<ExprPtr>("t_{}^{}"));
    }
    SECTION("Tensor") {
      auto expr = deserialize<ExprPtr>("t{i1;a1}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().label() == "t");
      REQUIRE(expr->as<Tensor>().bra().size() == 1);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == "i_1");
      REQUIRE(expr->as<Tensor>().ket().size() == 1);
      REQUIRE(expr->as<Tensor>().ket().at(0) == "a_1");
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == deserialize<ExprPtr>("t_{i1}^{a1}"));
      REQUIRE(expr == deserialize<ExprPtr>("t^{a1}_{i1}"));
      REQUIRE(expr == deserialize<ExprPtr>("t{i_1; a_1}"));
      REQUIRE(expr == deserialize<ExprPtr>("t{i_1; a_1;}"));
      REQUIRE(expr == deserialize<ExprPtr>("t_{i_1}^{a_1}"));

      expr = deserialize<ExprPtr>("t{i1,i2;a1,a2}");
      REQUIRE(expr->as<Tensor>().bra().size() == 2);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == "i_1");
      REQUIRE(expr->as<Tensor>().bra().at(1).label() == "i_2");
      REQUIRE(expr->as<Tensor>().ket().size() == 2);
      REQUIRE(expr->as<Tensor>().ket().at(0).label() == "a_1");
      REQUIRE(expr->as<Tensor>().ket().at(1).label() == "a_2");
      REQUIRE(expr->as<Tensor>().aux().empty());

      REQUIRE(expr == deserialize<ExprPtr>("+t{i1, i2; a1, a2}"));
      REQUIRE(deserialize<ExprPtr>("-t{i1;a1}")->is<Product>());
      REQUIRE(expr == deserialize<ExprPtr>("t{\ti1, \ti2; \na1,\t a2 \t}"));

      // Tensor labels including underscores
      REQUIRE(deserialize<ExprPtr>("T_1{i_1;a_1}")->as<Tensor>().label() ==
              "T_1");

      // "Non-standard" tensor labels
      REQUIRE(deserialize<ExprPtr>("α{a1;i1}")->as<Tensor>().label() == "α");
      REQUIRE(deserialize<ExprPtr>("γ_1{a1;i1}")->as<Tensor>().label() ==
              "γ_1");
      REQUIRE(deserialize<ExprPtr>("t⁔1{a1;i1}")->as<Tensor>().label() ==
              "t⁔1");
      REQUIRE(deserialize<ExprPtr>("t¹{a1;i1}")->as<Tensor>().label() == "t¹");
      REQUIRE(deserialize<ExprPtr>("t⁸{a1;i1}")->as<Tensor>().label() == "t⁸");
      REQUIRE(deserialize<ExprPtr>("t⁻{a1;i1}")->as<Tensor>().label() == "t⁻");
      REQUIRE(deserialize<ExprPtr>("tₐ{a1;i1}")->as<Tensor>().label() == "tₐ");
      REQUIRE(deserialize<ExprPtr>("t₋{a1;i1}")->as<Tensor>().label() == "t₋");
      REQUIRE(deserialize<ExprPtr>("t₌{a1;i1}")->as<Tensor>().label() == "t₌");
      REQUIRE(deserialize<ExprPtr>("t↓{a1;i1}")->as<Tensor>().label() == "t↓");
      REQUIRE(deserialize<ExprPtr>("t↑{a1;i1}")->as<Tensor>().label() == "t↑");

      // "Non-standard" index names
      auto expr1 = deserialize<ExprPtr>("t{a↓1;i↑1}");
      REQUIRE(expr1->as<Tensor>().bra().at(0).label() == "a↓_1");
      REQUIRE(expr1->as<Tensor>().ket().at(0).label() == "i↑_1");

      // Auxiliary indices
      expr = deserialize<ExprPtr>("t{;;i1}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().empty());
      REQUIRE(expr->as<Tensor>().ket().empty());
      REQUIRE(expr->as<Tensor>().aux().size() == 1);
      REQUIRE(expr->as<Tensor>().aux()[0].label() == "i_1");

      // All index groups at once
      expr = deserialize<ExprPtr>("t{i1,i2;a1;x1,x2}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(expr->as<Tensor>().bra().size() == 2);
      REQUIRE(expr->as<Tensor>().bra().at(0).label() == "i_1");
      REQUIRE(expr->as<Tensor>().bra().at(1).label() == "i_2");
      REQUIRE(expr->as<Tensor>().ket().size() == 1);
      REQUIRE(expr->as<Tensor>().ket().at(0).label() == "a_1");
      REQUIRE(expr->as<Tensor>().aux().size() == 2);
      REQUIRE(expr->as<Tensor>().aux().at(0).label() == "x_1");
      REQUIRE(expr->as<Tensor>().aux().at(1).label() == "x_2");
    }

    SECTION("Tensor with symmetry annotation") {
      auto expr1 = deserialize<ExprPtr>("t{a1;i1}:A");
      auto expr2 = deserialize<ExprPtr>("t{a1;i1}:S-C");
      auto expr3 = deserialize<ExprPtr>("t{a1;i1}:N-S-N");
      auto expr4 = deserialize<ExprPtr>("t{a1;i1}:N-H-N");
      auto expr5 = deserialize<ExprPtr>("t{a1;i1}:N-A-N");
      auto expr6 = deserialize<ExprPtr>("t{a1;i1}:N-N-N");

      const Tensor& t1 = expr1->as<Tensor>();
      const Tensor& t2 = expr2->as<Tensor>();
      const Tensor& t3 = expr3->as<Tensor>();
      const Tensor& t4 = expr4->as<Tensor>();
      const Tensor& t5 = expr5->as<Tensor>();
      const Tensor& t6 = expr6->as<Tensor>();

      REQUIRE(t1.symmetry() == Symmetry::Antisymm);

      REQUIRE(t2.symmetry() == Symmetry::Symm);
      REQUIRE(t2.braket_symmetry() == BraKetSymmetry::Conjugate);

      REQUIRE(t3.symmetry() == Symmetry::Nonsymm);
      REQUIRE(t3.braket_symmetry() == BraKetSymmetry::Symm);
      REQUIRE(t3.column_symmetry() == ColumnSymmetry::Nonsymm);

      REQUIRE(t4.hermiticity() == Hermiticity::Hermitian);
      REQUIRE(t5.hermiticity() == Hermiticity::AntiHermitian);
      REQUIRE(t6.hermiticity() == Hermiticity::NonHermitian);
    }

    SECTION("NormalOperator") {
      {
        using NOp = FNOperator;
        auto expr = deserialize<ExprPtr>("a{i1;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[0]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == "a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == "i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::Physical);
      }
      {
        using NOp = FNOperator;
        auto expr = deserialize<ExprPtr>("ã{i1;}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[1]);
        REQUIRE(expr->as<NOp>().creators().size() == 0);
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == "i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::SingleProduct);
      }

      {
        using NOp = BNOperator;
        auto expr = deserialize<ExprPtr>("b{i1;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[0]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == "a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 1);
        REQUIRE(expr->as<NOp>().annihilators().at(0).index() == "i_1");
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::Physical);
      }
      {
        using NOp = BNOperator;
        auto expr = deserialize<ExprPtr>("b̃{;a1}");
        REQUIRE(expr->is<NOp>());
        REQUIRE(expr->as<NOp>().label() == NOp::labels()[1]);
        REQUIRE(expr->as<NOp>().creators().size() == 1);
        REQUIRE(expr->as<NOp>().creators().at(0).index().label() == "a_1");
        REQUIRE(expr->as<NOp>().annihilators().size() == 0);
        REQUIRE(expr->as<NOp>().vacuum() == Vacuum::SingleProduct);
      }
      {
        // a tilde operator takes the Context's (non-Physical) vacuum
        auto ctx_mp = get_default_context();
        ctx_mp.set(Vacuum::MultiProduct);
        auto resetter = set_scoped_default_context(ctx_mp);
        REQUIRE(deserialize<ExprPtr>("ã{i1;i2}")->as<FNOperator>().vacuum() ==
                Vacuum::MultiProduct);
        REQUIRE(deserialize<ExprPtr>("b̃{i1;i2}")->as<BNOperator>().vacuum() ==
                Vacuum::MultiProduct);
      }
      {
        auto ctx_phys = get_default_context();
        ctx_phys.set(Vacuum::Physical);
        auto resetter = set_scoped_default_context(ctx_phys);
        REQUIRE(deserialize<ExprPtr>("ã{i1;i2}")->as<FNOperator>().vacuum() ==
                Vacuum::SingleProduct);
        REQUIRE(deserialize<ExprPtr>("b̃{i1;i2}")->as<BNOperator>().vacuum() ==
                Vacuum::SingleProduct);
      }
    }

    SECTION("Constant") {
      REQUIRE(deserialize<ExprPtr>("1/2")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("0/2")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("-1/2")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("-0/2")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("1")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("123")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("1.")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("01.00")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("0 / 10")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>("0.5/0.25")->is<Constant>());
      REQUIRE(deserialize<ExprPtr>(".4")->is<Constant>());
    }

    SECTION("Variable") {
      // SeQuant variable is just a label followed by an optional ^*
      // to denote if the variable is conjugated
      REQUIRE(deserialize<ExprPtr>("a")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("α")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("β")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("γ")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("λ")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("δ")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("a^*")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("α^*")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("β^*")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("b^*")->is<Variable>());
      REQUIRE(deserialize<ExprPtr>("b^*")->as<Variable>().conjugated());
      REQUIRE(deserialize<ExprPtr>("b^*")->as<Variable>().label() == "b");
    }

    SECTION("Power") {
      auto half = deserialize<ExprPtr>("2^(1/2)");
      REQUIRE(half->is<Power>());
      REQUIRE(half->as<Power>().base()->is<Constant>());
      REQUIRE(half->as<Power>().base()->as<Constant>().value() == rational{2});
      REQUIRE(half->as<Power>().exponent() == rational{1, 2});

      auto x_square = deserialize<ExprPtr>("x^(2)");
      REQUIRE(x_square->is<Power>());
      REQUIRE(x_square->as<Power>().base()->is<Variable>());
      REQUIRE(x_square->as<Power>().base()->as<Variable>().label() == "x");
      REQUIRE(x_square->as<Power>().exponent() == rational{2});

      // whitespace
      REQUIRE(deserialize<ExprPtr>("x ^ ( 1/2 )") ==
              deserialize<ExprPtr>("x^(1/2)"));

      // Power composes inside Products
      auto prod = deserialize<ExprPtr>("3 * x^(1/2) * t{i1;a1}");
      REQUIRE(prod->is<Product>());
      REQUIRE(prod->as<Product>().scalar() == rational{3});
      REQUIRE(prod->as<Product>().factor(0)->is<Power>());
      REQUIRE(prod->as<Product>().factor(1)->is<Tensor>());

      auto neg_half = ex<Power>(ex<Constant>(rational{-2}), rational{1, 2});
      const auto neg_half_str = "(-2)^(1/2)";
      REQUIRE(serialize(neg_half) == neg_half_str);
      REQUIRE(deserialize<ExprPtr>("(-2)^(1/2)") == neg_half);

      auto frac_half = ex<Power>(ex<Constant>(rational{3, 2}), rational{1, 2});
      const auto frac_half_str = "(3/2)^(1/2)";
      REQUIRE(serialize(frac_half) == frac_half_str);
      REQUIRE(deserialize<ExprPtr>("(3/2)^(1/2)") == frac_half);

      // conjugated Variable base must be parenthesized
      auto conj_var = ex<Variable>("x");
      conj_var->as<Variable>().conjugate();
      auto pw_conj_var = ex<Power>(conj_var, rational{2});
      const auto pw_conj_var_str = "(x^*)^(2)";
      REQUIRE(serialize(pw_conj_var) == pw_conj_var_str);
      REQUIRE(deserialize<ExprPtr>(pw_conj_var_str) == pw_conj_var);

      // conjugated power with conjugated Variable base: ((x^*)^(2))^*
      auto pw_conj_var_conj = pw_conj_var->clone();
      pw_conj_var_conj->as<Power>().conjugate();
      const auto pw_conj_var_conj_str = "((x^*)^(2))^*";
      REQUIRE(serialize(pw_conj_var_conj) == pw_conj_var_conj_str);
      REQUIRE(deserialize<ExprPtr>(pw_conj_var_conj_str) == pw_conj_var_conj);
    }

    SECTION("Product") {
      auto expr = deserialize<ExprPtr>("-1/2 g{i2,i3; i1,a2} t{a1,a2; i2,i3}");
      REQUIRE(expr->is<Product>());

      auto const& prod = expr->as<Product>();
      REQUIRE(prod.scalar() == rational{-1, 2});
      REQUIRE(prod.factor(0) ==
              deserialize<ExprPtr>("g_{i_2, i_3}^{i_1, a_2}"));
      REQUIRE(prod.factor(1) == deserialize<ExprPtr>("t^{i2, i3}_{a1, a2}"));
      REQUIRE(deserialize<ExprPtr>("-1/2 * δ * t{i1;a1}") ==
              deserialize<ExprPtr>("-1/2  δ  t{i1;a1}"));
      auto const prod2 =
          deserialize<ExprPtr>("-1/2 * δ * γ * t{i1;a1}")->as<Product>();
      REQUIRE(prod2.scalar() == rational{-1, 2});
      REQUIRE(prod2.factor(0) == ex<Variable>("δ"));
      REQUIRE(prod2.factor(1) == ex<Variable>("γ"));
      REQUIRE(prod2.factor(2)->is<Tensor>());
    }

    SECTION("Sum") {
      auto expr1 = deserialize<ExprPtr>(
          "f{a1;i1}"
          "- 1/2*g{i2,a1; a2,a3}t{a2,a3; i1,i2}");
      REQUIRE(expr1->is<Sum>());

      auto const& sum1 = expr1->as<Sum>();
      REQUIRE(sum1.summand(0) == deserialize<ExprPtr>("f{a1;i1}"));
      REQUIRE(
          sum1.summand(1) ==
          deserialize<ExprPtr>("- 1/2 * g{i2,a1; a2,a3} * t{a2,a3; i1,i2}"));

      auto expr2 = deserialize<ExprPtr>("a - 4");
      REQUIRE(expr2->is<Sum>());

      auto const& sum2 = expr2->as<Sum>();
      REQUIRE(sum2.summand(0) == deserialize<ExprPtr>("a"));
      REQUIRE(sum2.summand(1) == deserialize<ExprPtr>("-4"));
    }

    SECTION("Parentheses") {
      auto expr1 = deserialize<ExprPtr>(
          "-1/2 g{i2,i3; a2,a3} * ( t{a1,a3; i2,i3} * t{a2;i1} )");
      REQUIRE(expr1->is<Product>());

      auto const& prod1 = expr1->as<Product>();
      REQUIRE(prod1.size() == 2);
      REQUIRE(prod1.scalar() == rational{-1, 2});
      REQUIRE(prod1.factor(0)->is<Tensor>());
      REQUIRE(prod1.factor(1)->is<Product>());
      REQUIRE(prod1.factor(1)->size() == 2);

      auto expr2 = deserialize<ExprPtr>(
          "(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (t{a2;i1})");
      REQUIRE(expr2->is<Product>());

      auto const& prod2 = expr2->as<Product>();
      REQUIRE(prod2.size() == 2);
      REQUIRE(prod2.scalar() == rational{-1, 2});
      REQUIRE(prod2.factor(0)->is<Product>());
      REQUIRE(prod2.factor(0)->at(0) ==
              deserialize<ExprPtr>("g{i2,i3; a2,a3}"));
      REQUIRE(prod2.factor(0)->at(1) ==
              deserialize<ExprPtr>("t{a1,a3; i2,i3}"));
      REQUIRE(prod2.factor(1) == deserialize<ExprPtr>("t{a2;i1}"));

      auto expr3 = deserialize<ExprPtr>(
          "(-1/2) ( g{i2,i3; a2,a3} * t{a1,a3; i2,i3} ) * (1/2) * "
          "((t{a2;i1}))");
      REQUIRE(expr3->is<Product>());

      auto const& prod3 = expr3->as<Product>();
      REQUIRE(prod3.size() == 2);
      REQUIRE(prod3.scalar() == rational{-1, 4});
      REQUIRE(prod3.factor(0)->is<Product>());
      REQUIRE(prod3.factor(0)->at(0) ==
              deserialize<ExprPtr>("g{i2,i3; a2,a3}"));
      REQUIRE(prod3.factor(0)->at(1) ==
              deserialize<ExprPtr>("t{a1,a3; i2,i3}"));
      REQUIRE(prod3.factor(1) == deserialize<ExprPtr>("t{a2;i1}"));

      auto expr4 = deserialize<ExprPtr>("1/2 (a + b) * c");
      REQUIRE(expr4->is<Product>());

      const auto& prod4 = expr4->as<Product>();
      REQUIRE(prod4.size() == 2);
      REQUIRE(prod4.scalar() == rational{1, 2});
      REQUIRE(prod4.factor(1)->is<Variable>());
      REQUIRE(prod4.factor(1)->as<Variable>().label() == "c");
      REQUIRE(prod4.factor(0)->is<Sum>());
      const auto& nestedSum = prod4.factor(0)->as<Sum>();
      REQUIRE(nestedSum.size() == 2);
      REQUIRE(nestedSum.summand(0)->is<Variable>());
      REQUIRE(nestedSum.summand(0)->as<Variable>().label() == "a");
      REQUIRE(nestedSum.summand(1)->is<Variable>());
      REQUIRE(nestedSum.summand(1)->as<Variable>().label() == "b");
    }

    SECTION("Mixed") {
      auto expr = deserialize<ExprPtr>(
          "0.25 g{a1,a2; i1,i2}"
          "+ 1/4 g{i3,i4; a3,a4} (t{a3;i1} * t{a4;i2}) * (t{a1;i3} * "
          "t{a2;i4})");
      REQUIRE(expr->is<Sum>());
      auto const& sum = expr->as<Sum>();

      REQUIRE(sum.size() == 2);

      REQUIRE(sum.summand(0)->is<Product>());
      REQUIRE(sum.summand(0)->as<Product>().scalar() == rational{1, 4});
      REQUIRE(sum.summand(0)->size() == 1);
      REQUIRE(sum.summand(0)->at(0) == deserialize<ExprPtr>("g{a1,a2; i1,i2}"));

      REQUIRE(sum.summand(1)->is<Product>());
      auto const& prod = sum.summand(1)->as<Product>();

      REQUIRE(prod.scalar() == rational{1, 4});
      REQUIRE(prod.size() == 3);
      REQUIRE(prod.factor(0) == deserialize<ExprPtr>("g{i3,i4; a3,a4}"));

      REQUIRE(prod.factor(1)->is<Product>());
      REQUIRE(prod.factor(1)->at(0) == deserialize<ExprPtr>("t{a3;i1}"));
      REQUIRE(prod.factor(1)->at(1) == deserialize<ExprPtr>("t{a4;i2}"));

      REQUIRE(prod.factor(2)->is<Product>());
      REQUIRE(prod.factor(2)->at(0) == deserialize<ExprPtr>("t{a1;i3}"));
      REQUIRE(prod.factor(2)->at(1) == deserialize<ExprPtr>("t{a2;i4}"));
    }

    SECTION("Empty input") { REQUIRE(deserialize<ExprPtr>("") == nullptr); }

    // every registered space label is an index space name of the grammar
    SECTION("Index in a primed space") {
      auto f12 = get_default_context();
      f12.set(mbpt::make_F12_sr_spaces());
      auto f12_resetter = set_scoped_default_context(f12);
      const auto expr = deserialize<ExprPtr>("t{α'_1;i_1}");
      REQUIRE(expr->is<Tensor>());
      REQUIRE(
          expr->as<Tensor>().bra().at(0).space() ==
          get_default_context().index_basis_registry()->retrieve_space("α'"));
      REQUIRE(*deserialize<ExprPtr>(serialize(expr)) == *expr);
    }

    SECTION("Unicode labels") {
      // labels consisting of multi-byte UTF-8 sequences, incl. characters
      // outside of the BMP (𝛼 = U+1D6FC) and combining characters (f̃)
      for (const std::string label : {"λ", "𝛼", "𝛼β", "f̃", "Ж₂"}) {
        CAPTURE(label);
        const ExprPtr var = deserialize<ExprPtr>(label);
        REQUIRE(var->is<Variable>());
        REQUIRE(var->as<Variable>().label() == label);
        REQUIRE(serialize(var) == label);

        const ExprPtr tensor = deserialize<ExprPtr>(label + "{i1;a1}");
        REQUIRE(tensor->is<Tensor>());
        REQUIRE(tensor->as<Tensor>().label() == label);
        REQUIRE(deserialize<ExprPtr>(serialize(tensor)) == tensor);
      }
    }

    SECTION("Error handling") {
      SECTION("Exception type") {
        std::vector<std::string> inputs = {"t^",
                                           "a + + b"
                                           "1/t"
                                           "T{}"
                                           "T^{i1}{a1}"};

        for (const std::string& current : inputs) {
          REQUIRE_THROWS_AS(deserialize<ExprPtr>(current),
                            io::serialization::SerializationError);
        }
      }

      SECTION("Invalid index") {
        REQUIRE_THROWS_MATCHES(deserialize<ExprPtr>("t{i1<az1>;}"),
                               io::serialization::SerializationError,
                               serializationErrorMatches(5, 3, "proto"));
        REQUIRE_THROWS_MATCHES(
            deserialize<ExprPtr>("t{i1;az3}"),
            io::serialization::SerializationError,
            serializationErrorMatches(5, 3, "Unknown index space"));
      }

      SECTION("Offsets are given in bytes") {
        // "α" occupies 2 bytes
        REQUIRE_THROWS_MATCHES(
            deserialize<ExprPtr>("α{i1;az3}"),
            io::serialization::SerializationError,
            serializationErrorMatches(6, 3, "Unknown index space"));
        // "𝛼" occupies 4 bytes
        REQUIRE_THROWS_MATCHES(
            deserialize<ExprPtr>("𝛼{i1;a3}:P"),
            io::serialization::SerializationError,
            serializationErrorMatches(12, 1, "Invalid symmetry specifier"));
      }

      SECTION("Invalid UTF-8") {
        using namespace std::string_literals;
        REQUIRE_THROWS_MATCHES(deserialize<ExprPtr>("t{i1;a1} + a\xFF"s),
                               io::serialization::SerializationError,
                               serializationErrorMatches(12, 1, "UTF-8"));
        // truncated multi-byte sequence
        REQUIRE_THROWS_MATCHES(
            deserialize<ExprPtr>(std::string_view("αβ").substr(0, 3)),
            io::serialization::SerializationError,
            serializationErrorMatches(2, 1, "UTF-8"));
        REQUIRE_THROWS_MATCHES(deserialize<ResultExpr>("R = a\xFF"s),
                               io::serialization::SerializationError,
                               serializationErrorMatches(5, 1, "UTF-8"));
      }

      SECTION("Invalid symmetry") {
        REQUIRE_THROWS_MATCHES(
            deserialize<ExprPtr>("t{i1;a3}:P"),
            io::serialization::SerializationError,
            serializationErrorMatches(9, 1, "Invalid symmetry specifier"));
      }

      SECTION("an overlap or Kronecker delta is between two indices") {
        // rejected by the Tensor ctor
        for (const char* bad : {"s{i1,i2;a1,a2}", "δ{i1,i2;a1,a2}", "s{;i1}",
                                "δ{i1;}", "s{i1;a1;x}"})
          REQUIRE_THROWS_AS(deserialize<ExprPtr>(bad), Exception);
        REQUIRE_NOTHROW(deserialize<ExprPtr>("s{i1;a1} * δ{i2;a2}"));
      }

      SECTION("(anti)symmetrization operators reject braket symmetry") {
        // a reserved (anti)symmetrizer must be braket-Nonsymm; an explicit
        // non-Nonsymm braket spec is rejected by the Tensor ctor
        REQUIRE_THROWS(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}:N-C-S"));
        REQUIRE_THROWS(deserialize<ExprPtr>("Â{i1,i2;a1,a2}:A-C-S"));
        // the default (braket-Nonsymm) form is accepted
        REQUIRE_NOTHROW(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}"));
        REQUIRE_NOTHROW(deserialize<ExprPtr>("Â{i1,i2;a1,a2}"));
        // ... and stays accepted even when the Context defaults to Hermitian:
        // reserved (anti)symmetrizers must not inherit the Context braket
        // default (which would otherwise derive a non-Nonsymm braket and throw)
        {
          auto ctx = get_default_context();
          ctx.set(Hermiticity::Hermitian);
          auto resetter = set_scoped_default_context(ctx);
          REQUIRE_NOTHROW(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}"));
          REQUIRE_NOTHROW(deserialize<ExprPtr>("Â{i1,i2;a1,a2}"));
        }
      }

      SECTION("(anti)symmetrization operators reject contradicting perm") {
        // an explicitly spelled perm symmetry that contradicts the reserved
        // label's defining one is rejected rather than silently overwritten,
        // as for the braket and column specs
        REQUIRE_THROWS(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}:A-N-S"));
        REQUIRE_THROWS(deserialize<ExprPtr>("Â{i1,i2;a1,a2}:N-N-S"));
        // the matching spellings stay accepted
        REQUIRE_NOTHROW(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}:N-N-S"));
        REQUIRE_NOTHROW(deserialize<ExprPtr>("Â{i1,i2;a1,a2}:A-N-S"));
      }

      SECTION("(anti)symmetrization operators have fixed perm symmetry") {
        // the defining bra/ket permutational symmetry of a reserved
        // (anti)symmetrizer is not a free parameter either, so it must not be
        // inherited from the Context; a parsed operator must match the one
        // built by make_(anti)symmetrizer()
        auto ctx = get_default_context();
        ctx.set(Symmetry::Antisymm);
        auto resetter = set_scoped_default_context(ctx);
        const auto S = deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}");
        REQUIRE(S->as<Tensor>().symmetry() == Symmetry::Nonsymm);
        REQUIRE(S == make_symmetrizer(bra{Index("i_1"), Index("i_2")},
                                      ket{Index("a_1"), Index("a_2")}));
        const auto A = deserialize<ExprPtr>("Â{i1,i2;a1,a2}");
        REQUIRE(A->as<Tensor>().symmetry() == Symmetry::Antisymm);
        REQUIRE(A == make_antisymmetrizer(bra{Index("i_1"), Index("i_2")},
                                          ket{Index("a_1"), Index("a_2")}));
      }

      SECTION("(anti)symmetrization operators are always column symmetric") {
        // ... including under a Context whose column default is the library
        // default Nonsymm: the deserializer must supply Symm rather than pass
        // the (contradicting) Context default to the Tensor ctor
        auto ctx = get_default_context();
        ctx.set(ColumnSymmetry::Nonsymm);
        auto resetter = set_scoped_default_context(ctx);
        for (const auto& input :
             {"Ŝ{i1,i2;a1,a2}", "Â{i1,i2;a1,a2}", "Ŝ{i1,i2;a1,a2}:N-N-S",
              "Â{i1,i2;a1,a2}:A-N-S"}) {
          ExprPtr expr;
          REQUIRE_NOTHROW(expr = deserialize<ExprPtr>(input));
          REQUIRE(expr->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);
        }
        // a *spelled-out* column symmetry that contradicts the defining one is
        // rejected, exactly as the equivalent programmatic ctor call is --
        // only a defaulted one is silently replaced above
        REQUIRE_THROWS(deserialize<ExprPtr>("Ŝ{i1,i2;a1,a2}:N-N-N"));
        REQUIRE_THROWS(deserialize<ExprPtr>("Â{i1,i2;a1,a2}:A-N-N"));
      }

      SECTION("(anti)symmetric bra/ket implies column symmetry") {
        // an unspecified column symmetry follows the implication rather than
        // the Context default, ...
        auto ctx = get_default_context();
        ctx.set(ColumnSymmetry::Nonsymm);
        auto resetter = set_scoped_default_context(ctx);
        for (const auto& input : {"t{i1,i2;a1,a2}:A", "t{i1,i2;a1,a2}:S"}) {
          ExprPtr expr;
          REQUIRE_NOTHROW(expr = deserialize<ExprPtr>(input));
          REQUIRE(expr->as<Tensor>().column_symmetry() == ColumnSymmetry::Symm);
        }
        // ... while a spelled-out one that contradicts it is rejected
        REQUIRE_THROWS(deserialize<ExprPtr>("t{i1,i2;a1,a2}:A-N-N"));
        REQUIRE_THROWS(deserialize<ExprPtr>("t{i1,i2;a1,a2}:S-N-N"));
        // a Nonsymm bra/ket implies nothing, so the Context default stands
        REQUIRE(deserialize<ExprPtr>("t{i1,i2;a1,a2}:N")
                    ->as<Tensor>()
                    .column_symmetry() == ColumnSymmetry::Nonsymm);
      }
    }
  }

  SECTION("deserialize<ResultExpr>") {
    using namespace sequant;

    SECTION("constant") {
      ResultExpr result = deserialize<ResultExpr>("A = 3");

      REQUIRE(result.has_label());
      REQUIRE(result.label() == "A");
      REQUIRE(result.bra().empty());
      REQUIRE(result.ket().empty());
      REQUIRE(result.symmetry() == Symmetry::Nonsymm);
      REQUIRE(result.braket_symmetry() == BraKetSymmetry::Nonsymm);
      REQUIRE(result.column_symmetry() == ColumnSymmetry::Nonsymm);

      REQUIRE(result.expression().is<Constant>());
      REQUIRE(result.expression().as<Constant>().value<int>() == 3);
    }
    SECTION("contraction") {
      ResultExpr result =
          deserialize<ResultExpr>("R{i1,i2;e1,e2}:A = f{e2;e3} t{e1,e3;i1,i2}");

      REQUIRE(result.has_label());
      REQUIRE(result.label() == "R");
      REQUIRE(result.bra().size() == 2);
      REQUIRE(result.bra()[0].full_label() == "i_1");
      REQUIRE(result.bra()[1].full_label() == "i_2");
      REQUIRE(result.ket().size() == 2);
      REQUIRE(result.ket()[0].full_label() == "e_1");
      REQUIRE(result.ket()[1].full_label() == "e_2");
      REQUIRE(result.symmetry() == Symmetry::Antisymm);
      REQUIRE(result.braket_symmetry() == BraKetSymmetry::Nonsymm);
      // This column symmetry is implicit with Symmetry::Antisymm
      REQUIRE(result.column_symmetry() == ColumnSymmetry::Symm);

      REQUIRE(result.expression().is<Product>());
      const Product& prod = result.expression().as<Product>();
      REQUIRE(prod.size() == 2);
      REQUIRE(prod.factor(0).is<Tensor>());
      REQUIRE(prod.factor(0).as<Tensor>().label() == "f");
      REQUIRE(prod.factor(1).is<Tensor>());
      REQUIRE(prod.factor(1).as<Tensor>().label() == "t");
    }
  }

  SECTION("serialize") {
    using namespace sequant;

    std::vector<std::string> expressions = {
        "t{a_1,a_2;a_3,a_4}:N-C-S",
        "42",
        "1/2",
        "-1/4 t{a_1,i_1<a_1>;a_2,i_2}:S-N-S",
        "a + b - 4 specialVariable",
        "variable + A{a_1;i_1}:N-N-S * B{i_1;a_1}:A-C-S",
        "1/2 (a + b) * c",
        "T1{}:N-N-N + T2{;;x_1}:N-N-N * T3{;;x_1}:N-N-N "
        "+ T4{a_1;;x_2}:S-C-S * T5{;a_1;x_2}:S-S-S",
        "q1 * q2^* * q3",
        "1/2 ã{i_1;i_2} * b̃{i_3;i_4}"};

    for (const std::string& current : expressions) {
      ExprPtr expression = deserialize<ExprPtr>(current);

      REQUIRE(serialize(expression, {.annot_symm = true}) == current);
    }

    SECTION("nested Sum") {
      Sum::summands_type summands{
          ex<Variable>("a"),
          ex<Sum>(ExprPtrList{ex<Variable>("b"), ex<Variable>("c")})};
      auto nested = ex<Sum>(std::move(summands), Sum::move_only_tag{});
      REQUIRE(serialize(nested) == "a + (b + c)");
    }

    SECTION("result_expressions") {
      std::vector<std::string> expressions = {
          "A = 5",
          "A = g{i_1,i_2;e_1,e_2}:S-N-S * t{e_1,e_2;i_1,i_2}:N-N-S",
          "R{i_1,i_2;e_1,e_2}:A-N-S = f{e_2;e_3}:A-N-S * "
          "t{e_1,e_3;i_1,i_2}:A-N-S + "
          "g{i_1,i_2;e_1,e_2}:A-N-S",
      };

      for (const std::string& current : expressions) {
        ResultExpr result = deserialize<ResultExpr>(current);

        REQUIRE(serialize(result, {.annot_symm = true}) == current);
      }
    }
  }
}
