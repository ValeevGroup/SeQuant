//
// Created by Bimal Gaudel on 6/12/21.
//

#include "catch.hpp"

#include <SeQuant/core/parse/token_sequant.hpp>

TEST_CASE("TEST_TOKEN", "[parse_expr]") {
  using namespace sequant::parse;
  using sequant::Constant;
  using sequant::ex;
  using sequant::ExprPtr;
  using sequant::IndexList;
  using sequant::Tensor;

  SECTION("Ctor") {
    REQUIRE_NOTHROW(LeftParenthesis{});
    REQUIRE_NOTHROW(RightParenthesis{});
    REQUIRE_NOTHROW(OperatorTimes{});
    REQUIRE_NOTHROW(OperatorPlus{});
    REQUIRE_NOTHROW(OperatorMinus{});
    REQUIRE_NOTHROW(OperandTensor{L"t", IndexList{L"a_1"}, IndexList{L"i_1"}});
    REQUIRE_NOTHROW(OperandConstant{1});
  }

  SECTION("Id") {
    auto const lparen = LeftParenthesis{};
    REQUIRE(lparen.is<LeftParenthesis>());
    REQUIRE_FALSE(lparen.is<Operand>());
    REQUIRE_FALSE(lparen.is<Operator>());

    auto const rparen = RightParenthesis{};
    REQUIRE(rparen.is<RightParenthesis>());
    REQUIRE_FALSE(rparen.is<Operand>());
    REQUIRE_FALSE(rparen.is<Operator>());

    auto const times = OperatorTimes{};
    REQUIRE(times.is<OperatorTimes>());
    REQUIRE(times.is<Operator>());
    REQUIRE_FALSE(times.is<Operand>());

    auto const plus = OperatorPlus{};
    REQUIRE(plus.is<OperatorPlus>());
    REQUIRE(plus.is<Operator>());
    REQUIRE_FALSE(plus.is<Operand>());

    auto const minus = OperatorMinus{};
    REQUIRE(minus.is<OperatorMinus>());
    REQUIRE(minus.is<Operator>());
    REQUIRE_FALSE(minus.is<Operand>());

    auto const uplus = OperatorPlusUnary{};
    REQUIRE(uplus.is<OperatorPlusUnary>());
    REQUIRE(uplus.is<Operator>());
    REQUIRE_FALSE(uplus.is<Operand>());

    auto const uminus = OperatorMinusUnary{};
    REQUIRE(uminus.is<OperatorMinusUnary>());
    REQUIRE(uminus.is<Operator>());
    REQUIRE_FALSE(uminus.is<Operand>());

    auto const fraction = OperandConstant{42};
    REQUIRE(fraction.is<OperandConstant>());
    REQUIRE(fraction.is<Operand>());
    REQUIRE_FALSE(fraction.is<Operator>());

    auto const tensor =
        OperandTensor{L"t", IndexList{L"a_1"}, IndexList{L"i_1"}};
    REQUIRE(tensor.is<OperandTensor>());
    REQUIRE(tensor.is<Operand>());
    REQUIRE_FALSE(tensor.is<Operator>());
  }

  SECTION("Unique pointer and casting") {
    auto const pl = token<LeftParenthesis>();
    REQUIRE(pl->is<LeftParenthesis>());

    auto const pr = token<RightParenthesis>();
    REQUIRE(pr->is<RightParenthesis>());

    auto const ot = token<OperatorTimes>();
    REQUIRE(ot->is<Operator>());
    REQUIRE(ot->is<OperatorTimes>());
    REQUIRE_NOTHROW(ot->as<Operator>());
    REQUIRE_NOTHROW(ot->as<OperatorTimes>());

    auto const op = token<OperatorPlus>();
    REQUIRE(op->is<Operator>());
    REQUIRE(op->is<OperatorPlus>());
    REQUIRE_NOTHROW(op->as<Operator>());
    REQUIRE_NOTHROW(op->as<OperatorPlus>());

    auto const om = token<OperatorMinus>();
    REQUIRE(om->is<Operator>());
    REQUIRE(om->is<OperatorMinus>());
    REQUIRE_NOTHROW(om->as<Operator>());
    REQUIRE_NOTHROW(om->as<OperatorMinus>());

    auto const uop = token<OperatorPlusUnary>();
    REQUIRE(uop->is<Operator>());
    REQUIRE(uop->is<OperatorPlusUnary>());
    REQUIRE_NOTHROW(uop->as<Operator>());
    REQUIRE_NOTHROW(uop->as<OperatorPlusUnary>());

    auto const uom = token<OperatorMinusUnary>();
    REQUIRE(uom->is<Operator>());
    REQUIRE(uom->is<OperatorMinusUnary>());
    REQUIRE_NOTHROW(uom->as<Operator>());
    REQUIRE_NOTHROW(uom->as<OperatorMinusUnary>());

    auto const f1 = token<OperandConstant>(2);
    REQUIRE(f1->is<Operand>());
    REQUIRE(f1->is<OperandConstant>());
    REQUIRE_NOTHROW(f1->as<Operand>());
    REQUIRE_NOTHROW(f1->as<OperandConstant>());
    REQUIRE(f1->as<OperandConstant>().clone() == ex<Constant>(2));

    auto const t1 =
        token<OperandTensor>(L"t", IndexList{L"i_1"}, IndexList{L"a_1"});
    REQUIRE(t1->is<Operand>());
    REQUIRE(t1->is<OperandTensor>());
    REQUIRE_NOTHROW(t1->as<Operand>());
    REQUIRE_NOTHROW(t1->as<OperandTensor>());
    REQUIRE(t1->as<OperandTensor>().clone() ==
            ex<Tensor>(L"t", IndexList{L"i_1"}, IndexList{L"a_1"}));
  }

  SECTION("Operator precedence") {
    REQUIRE(OperatorTimes{}.precedence() > OperatorPlus{}.precedence());
    REQUIRE(OperatorTimes{}.precedence() > OperatorMinus{}.precedence());
    REQUIRE(OperatorPlus{}.precedence() < OperatorTimes{}.precedence());
    REQUIRE(OperatorMinus{}.precedence() < OperatorTimes{}.precedence());
    REQUIRE(OperatorMinus{}.precedence() == OperatorPlus{}.precedence());
  }
}
