#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/polynomial_division.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string>
#include <tuple>
#include <vector>

namespace {

bool is_zero(const sequant::ExprPtr &expr) {
  return expr.is<sequant::Constant>() && expr.as<sequant::Constant>().is_zero();
}

}  // namespace

TEST_CASE("polynomial_division", "[utilities][polynomial_division]") {
  using namespace sequant;

  SECTION("scalars") {
    SECTION("exact") {
      const ExprPtr dividend = deserialize(L"u w + u z + v w + v z");
      const ExprPtr divisor = deserialize(L"u + v");

      for (DivisionMode mode : {DivisionMode::Strict, DivisionMode::Split,
                                DivisionMode::Remainder}) {
        DivisionResult result =
            divide(dividend, divisor, mode, IndexMatching::Exact);

        REQUIRE(result.exact());
        REQUIRE(result.groups.size() == 1);
        REQUIRE(result.groups[0].relabeling.empty());
        REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"w + z"));
        REQUIRE_THAT(result.groups[0].divisor, EquivalentTo(L"u + v"));
        REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
      }
    }

    SECTION("coefficients") {
      const ExprPtr dividend =
          deserialize(L"2 u w + 1/2 u z + 6 v w + 3/2 v z");
      const ExprPtr divisor = deserialize(L"1/2 u + 3/2 v");

      DivisionResult result =
          divide(dividend, divisor, DivisionMode::Strict, IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"4 w + z"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("powers") {
      const ExprPtr dividend = deserialize(L"x^(3) y + x");
      const ExprPtr divisor = deserialize(L"x");

      DivisionResult result =
          divide(dividend, divisor, DivisionMode::Strict, IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"x^(2) y + 1"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));

      // repeated factors
      result = divide(deserialize(L"x x y"), divisor, DivisionMode::Strict,
                      IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"x y"));

      // insufficient power
      result = divide(deserialize(L"x y"), deserialize(L"x^(2)"),
                      DivisionMode::Remainder, IndexMatching::Exact);
      REQUIRE(result.groups.empty());
      REQUIRE_THAT(result.remainder, EquivalentTo(L"x y"));
    }

    SECTION("constant divisor") {
      const ExprPtr dividend = deserialize(L"u + 2 v");

      DivisionResult result =
          divide(dividend, deserialize(L"2"), DivisionMode::Strict,
                 IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"1/2 u + v"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("divisor with constant term") {
      const ExprPtr dividend = deserialize(L"u w + w + u + 1");
      const ExprPtr divisor = deserialize(L"u + 1");

      DivisionResult result =
          divide(dividend, divisor, DivisionMode::Strict, IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"w + 1"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("zero dividend") {
      DivisionResult result =
          divide(ex<Constant>(0), deserialize(L"u + v"), DivisionMode::Strict,
                 IndexMatching::Exact);

      REQUIRE(result.exact());
      REQUIRE(result.groups.empty());
      REQUIRE(is_zero(result.to_expr()));
    }

    SECTION("zero divisor") {
      if (assert_behavior() != AssertBehavior::Abort) {
        REQUIRE_THROWS(divide(deserialize(L"u"), ex<Constant>(0),
                              DivisionMode::Remainder, IndexMatching::Exact));
      }
    }
  }

  SECTION("modes") {
    // (u + v) w + z = u w + v w + z
    const ExprPtr divisor = deserialize(L"u + v");

    SECTION("Strict") {
      REQUIRE_THROWS_AS(divide(deserialize(L"u w + z"), divisor,
                               DivisionMode::Strict, IndexMatching::Exact),
                        Exception);
    }

    SECTION("Split") {
      DivisionResult result = divide(deserialize(L"u w + z"), divisor,
                                     DivisionMode::Split, IndexMatching::Exact);

      REQUIRE(result.groups.empty());
      REQUIRE_THAT(result.remainder, EquivalentTo(L"u w + z"));

      // only the cleanly factorizable group is factored
      const ExprPtr dividend = deserialize(L"u w + v w + u y + z");
      result =
          divide(dividend, divisor, DivisionMode::Split, IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"w"));
      REQUIRE_THAT(result.remainder, EquivalentTo(L"u y + z"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));

      // coefficients of the remaining terms may change but no new terms are
      // introduced
      const ExprPtr dividend2 = deserialize(L"u w + 3 v w + z");
      result =
          divide(dividend2, divisor, DivisionMode::Split, IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      // depending on the pivot, this is either 2 v w + z or z - 2 u w
      REQUIRE(result.remainder.is<Sum>());
      REQUIRE(result.remainder.as<Sum>().summands().size() == 2);
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend2));
    }

    SECTION("Remainder") {
      const ExprPtr dividend = deserialize(L"u w + z");
      DivisionResult result = divide(dividend, divisor, DivisionMode::Remainder,
                                     IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"w"));
      REQUIRE_THAT(result.remainder, EquivalentTo(L"z - v w"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }
  }

  SECTION("tensors") {
    const ExprPtr divisor = deserialize(L"f{a1;a2} t{a2;i1} + g{a1;i1}");

    SECTION("same externals") {
      const ExprPtr dividend =
          deserialize(L"f{a1;a3} t{a3;i1} X{i2;a2} + g{a1;i1} X{i2;a2}");

      for (IndexMatching matching :
           {IndexMatching::Exact, IndexMatching::Relabel}) {
        DivisionResult result =
            divide(dividend, divisor, DivisionMode::Strict, matching);

        REQUIRE(result.groups.size() == 1);
        REQUIRE(result.groups[0].relabeling.empty());
        REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"X{i2;a2}"));
        REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
      }
    }

    SECTION("bound to dummies") {
      const ExprPtr dividend =
          deserialize(L"f{a3;a4} t{a4;i3} Y{i3;a3} + g{a4;i4} Y{i4;a4}");

      for (IndexMatching matching :
           {IndexMatching::Exact, IndexMatching::Relabel}) {
        DivisionResult result =
            divide(dividend, divisor, DivisionMode::Strict, matching);

        REQUIRE(result.groups.size() == 1);
        REQUIRE(result.groups[0].relabeling.empty());
        REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"Y{i1;a1}"));
        REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
      }
    }

    SECTION("dummy clashing with divisor externals") {
      // a1 is used as a dummy in the dividend, while a3 is bound to the
      // divisor's external index a1
      const ExprPtr dividend =
          deserialize(L"f{a3;a1} t{a1;i3} Y{i3;a3} + g{a3;i3} Y{i3;a3}");

      DivisionResult result =
          divide(dividend, divisor, DivisionMode::Strict, IndexMatching::Exact);

      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"Y{i1;a1}"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("different externals") {
      const ExprPtr dividend =
          deserialize(L"f{a2;a3} t{a3;i2} X{i1;a1} + g{a2;i2} X{i1;a1}");

      // Exact matching doesn't allow relabeling the divisor's externals
      DivisionResult result = divide(dividend, divisor, DivisionMode::Remainder,
                                     IndexMatching::Exact);
      REQUIRE(result.groups.empty());
      REQUIRE_THAT(result.remainder, EquivalentTo(dividend));

      REQUIRE_THROWS_AS(
          divide(dividend, divisor, DivisionMode::Strict, IndexMatching::Exact),
          Exception);

      // Relabel matching does
      result = divide(dividend, divisor, DivisionMode::Strict,
                      IndexMatching::Relabel);
      REQUIRE(result.groups.size() == 1);
      REQUIRE(result.groups[0].relabeling ==
              container::map<Index, Index>{{Index{L"a_1"}, Index{L"a_2"}},
                                           {Index{L"i_1"}, Index{L"i_2"}}});
      REQUIRE_THAT(result.groups[0].divisor,
                   EquivalentTo(L"f{a2;a3} t{a3;i2} + g{a2;i2}"));
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"X{i1;a1}"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("multiple groups") {
      const ExprPtr dividend = deserialize(
          L"(f{a1;a3} t{a3;i1} + g{a1;i1}) X{i2;a2} + "
          L"(f{a2;a3} t{a3;i2} + g{a2;i2}) X{i1;a1}");

      DivisionResult result = divide(dividend, divisor, DivisionMode::Strict,
                                     IndexMatching::Relabel);
      REQUIRE(result.groups.size() == 2);
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));

      // With exact matching, only the identically indexed group is found
      result =
          divide(dividend, divisor, DivisionMode::Split, IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"X{i2;a2}"));
      REQUIRE_THAT(result.remainder,
                   EquivalentTo(L"(f{a2;a3} t{a3;i2} + g{a2;i2}) X{i1;a1}"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("partial match") {
      const ExprPtr dividend =
          deserialize(L"f{a1;a3} t{a3;i1} X{i2;a2} + h{a1;i1} X{i2;a2}");

      DivisionResult result =
          divide(dividend, divisor, DivisionMode::Split, IndexMatching::Exact);
      REQUIRE(result.groups.empty());
      REQUIRE_THAT(result.remainder, EquivalentTo(dividend));

      result = divide(dividend, divisor, DivisionMode::Remainder,
                      IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"X{i2;a2}"));
      REQUIRE_THAT(result.remainder,
                   EquivalentTo(L"h{a1;i1} X{i2;a2} - g{a1;i1} X{i2;a2}"));
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("repeated tensors") {
      const ExprPtr dividend = deserialize(L"t{a1;i1} t{a2;i2}");

      DivisionResult result =
          divide(dividend, deserialize(L"t{a1;i1}"), DivisionMode::Strict,
                 IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"t{a2;i2}"));

      result = divide(dividend, deserialize(L"t{a3;i3}"), DivisionMode::Strict,
                      IndexMatching::Relabel);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
    }

    SECTION("antisymmetric tensors") {
      const ExprPtr divisor2 =
          deserialize(L"g{a1,a2;i1,i2}:A + f{a1;a3} t{a3,a2;i1,i2}:A");
      const ExprPtr dividend = deserialize(
          L"(g{a4,a3;i3,i4}:A + f{a4;a5} t{a5,a3;i3,i4}:A) Y{i3,i4;a3,a4}");

      for (IndexMatching matching :
           {IndexMatching::Exact, IndexMatching::Relabel}) {
        DivisionResult result =
            divide(dividend, divisor2, DivisionMode::Strict, matching);
        REQUIRE(result.groups.size() == 1);
        REQUIRE_THAT(result.to_expr(), EquivalentTo(dividend));
      }
    }

    SECTION("scalars and tensors") {
      const ExprPtr dividend =
          deserialize(L"x f{a1;a2} t{a2;i1} y + x g{a1;i1} y + z");
      const ExprPtr divisor2 = deserialize(L"x f{a1;a2} t{a2;i1} + x g{a1;i1}");

      DivisionResult result = divide(
          dividend, divisor2, DivisionMode::Remainder, IndexMatching::Exact);
      REQUIRE(result.groups.size() == 1);
      REQUIRE_THAT(result.groups[0].quotient, EquivalentTo(L"y"));
      REQUIRE_THAT(result.remainder, EquivalentTo(L"z"));
    }

    SECTION("operators") {
      const ExprPtr dividend = deserialize(L"t{a1;i1}") *
                               ex<FNOperator>(cre({L"a_1"}), ann({L"i_1"}));

      DivisionResult result =
          divide(dividend, deserialize(L"t{a2;i2}"), DivisionMode::Strict,
                 IndexMatching::Relabel);
      REQUIRE(result.groups.size() == 1);
      REQUIRE(result.groups[0].quotient.is<FNOperator>());

      if (assert_behavior() != AssertBehavior::Abort) {
        REQUIRE_THROWS(divide(dividend,
                              ex<FNOperator>(cre({L"a_1"}), ann({L"i_1"})),
                              DivisionMode::Remainder, IndexMatching::Exact));
      }
    }

    SECTION("inconsistent divisor") {
      if (assert_behavior() != AssertBehavior::Abort) {
        REQUIRE_THROWS(divide(deserialize(L"u"),
                              deserialize(L"t{a1;i1} + t{a2;i2}"),
                              DivisionMode::Remainder, IndexMatching::Exact));
      }
    }
  }
}

TEST_CASE("factor_out", "[utilities][polynomial_division]") {
  using namespace sequant;

  SECTION("scalars") {
    const ExprPtr expr = deserialize(L"u w + u z + v w + v z + y");
    const ExprPtr factor = deserialize(L"u + v");

    ExprPtr keep = factor_out(expr, factor, DivisionMode::Split,
                              IndexMatching::Exact, FactorStyle::Keep);
    REQUIRE(keep.is<Sum>());
    REQUIRE(keep.as<Sum>().summands().size() == 2);
    REQUIRE_THAT(keep, EquivalentTo(expr));

    ExprPtr drop = factor_out(expr, factor, DivisionMode::Split,
                              IndexMatching::Exact, FactorStyle::Drop);
    REQUIRE_THAT(drop, EquivalentTo(L"w + z"));

    ExprPtr substituted = factor_out(expr, factor, *deserialize(L"s"),
                                     DivisionMode::Split, IndexMatching::Exact);
    REQUIRE_THAT(substituted, EquivalentTo(L"s w + s z + y"));
  }

  SECTION("intermediates") {
    const ResultExpr intermediate =
        deserialize<ResultExpr>(L"I{a1;i1} = f{a1;a2} t{a2;i1} + g{a1;i1}");

    SECTION("single group") {
      const ExprPtr expr = deserialize(
          L"f{a1;a3} t{a3;i1} X{i2;a2} + g{a1;i1} X{i2;a2} + h{a1,a2;i1,i2}");

      REQUIRE_THAT(factor_out(expr, intermediate, DivisionMode::Split,
                              IndexMatching::Exact),
                   EquivalentTo(L"I{a1;i1} X{i2;a2} + h{a1,a2;i1,i2}"));

      const ResultExpr result_expr = deserialize<ResultExpr>(
          L"R{a1,a2;i1,i2} = f{a1;a3} t{a3;i1} X{i2;a2} + g{a1;i1} X{i2;a2}");
      REQUIRE_THAT(factor_out(result_expr, intermediate, DivisionMode::Strict,
                              IndexMatching::Exact),
                   EquivalentTo(L"R{a1,a2;i1,i2} = I{a1;i1} X{i2;a2}"));
    }

    SECTION("multiple groups") {
      const ExprPtr expr = deserialize(
          L"(f{a1;a3} t{a3;i1} + g{a1;i1}) X{i2;a2} + "
          L"(f{a2;a3} t{a3;i2} + g{a2;i2}) X{i1;a1}");

      REQUIRE_THAT(factor_out(expr, intermediate, DivisionMode::Strict,
                              IndexMatching::Relabel),
                   EquivalentTo(L"I{a1;i1} X{i2;a2} + I{a2;i2} X{i1;a1}"));

      if (assert_behavior() != AssertBehavior::Abort) {
        REQUIRE_THROWS_AS(
            factor_out(expr, intermediate.expression(), DivisionMode::Strict,
                       IndexMatching::Relabel, FactorStyle::Drop),
            Exception);
      }
    }

    SECTION("contracted intermediate") {
      const ExprPtr expr =
          deserialize(L"f{a3;a4} t{a4;i3} Y{i3;a3} + g{a4;i4} Y{i4;a4}");

      REQUIRE_THAT(factor_out(expr, intermediate, DivisionMode::Strict,
                              IndexMatching::Exact),
                   EquivalentTo(L"I{a1;i1} Y{i1;a1}"));
    }
  }
}
