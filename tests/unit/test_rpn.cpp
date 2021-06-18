//
// Created by Bimal Gaudel on 6/13/21.
//

#include <SeQuant/core/parse/rpn.hpp>
#include <SeQuant/core/parse/token_sequant.hpp>
#include "catch.hpp"

TEST_CASE("TEST_RPN", "[parse_expr]") {
  using namespace sequant::parse;
  using Number = OperandConstant;
  using Times = OperatorTimes;
  using Plus = OperatorPlus;
  using Minus = OperatorMinus;
  using sequant::Constant;
  using sequant::Tensor;
  using sequant::ex;

  auto rpn = ReversePolishNotation{};

  SECTION("close method") {
    REQUIRE(rpn.close());
    rpn.add_left_parenthesis();
    REQUIRE_FALSE(rpn.close());
  }

  SECTION("tokens method") {
    rpn.reset();
    REQUIRE(rpn.tokens().empty());
  }

  SECTION("add tokens: simple") {
    rpn.reset();
    rpn.add_operand(token<Number>(Constant{0.5}));
    rpn.add_operator(token<Times>());
    rpn.add_operand(token<Number>(Constant{0.3}));
    // expression: '1/2 * 1/3'
    // post-fix notation (aka Reverse Polish Notation):
    //                   '1/2 1/3 *'
    REQUIRE(rpn.close());
    REQUIRE(rpn.tokens().size() == 3);
    REQUIRE(rpn.tokens().at(0)->is<Operand>());
    REQUIRE(rpn.tokens().at(0)->as<Number>() ==
            sequant::ex<Constant>(0.5));
    REQUIRE(rpn.tokens().at(1)->is<Operand>());
    REQUIRE(rpn.tokens().at(1)->as<OperandSequant>() ==
                sequant::ex<Constant>(0.3));
    REQUIRE(rpn.tokens().at(2)->is<Times>());
  }

  SECTION("add tokens: equal precedence operators") {
    rpn.reset();
    rpn.add_operand(token<Number>(Constant{2}));
    rpn.add_operator(token<Plus>());
    rpn.add_operand(token<Number>(Constant{3}));
    rpn.add_operator(token<Minus>());
    rpn.add_operand(token<Number>(Constant{1}));
    // expression: '2 + 2/3 - 1'
    // pos-fix notation:
    //                 '2 2/3 + 1 -'
    REQUIRE(rpn.close());
    REQUIRE(rpn.tokens().at(0)->as<Number>() ==
            ex<Constant>(2));
    REQUIRE(rpn.tokens().at(1)->as<Number>() ==
            ex<Constant>(3));
    REQUIRE(rpn.tokens().at(2)->is<Plus>());
    REQUIRE(rpn.tokens().at(3)->as<Number>() ==
            ex<Constant>(1));
    REQUIRE(rpn.tokens().at(4)->is<Minus>());
  }

  SECTION("add tokens: varying precedence operators") {
    rpn.reset();
    rpn.add_operand(token<Number>(Constant{1}));
    rpn.add_operator(token<Times>());
    rpn.add_operand(token<Number>(Constant{2}));
    rpn.add_operator(token<Plus>());
    rpn.add_operand(token<Number>(Constant{3}));
    rpn.add_operator(token<Minus>());
    rpn.add_operand(token<Number>(Constant{3}));
    rpn.add_operator(token<Times>());
    rpn.add_operand(token<Number>(Constant{4}));
    // expression: '1 * 2 + 3 - 3 * 4'
    // post-fix notation:
    //             '1 2 * 3 + 3 4 * -'
    REQUIRE(rpn.close());
    REQUIRE(rpn.tokens().at(0)->as<Number>() == ex<Constant>(1));
    REQUIRE(rpn.tokens().at(1)->as<Number>() == ex<Constant>(2));
    REQUIRE(rpn.tokens().at(2)->is<Times>());
    REQUIRE(rpn.tokens().at(3)->as<Number>() == ex<Constant>(3));
    REQUIRE(rpn.tokens().at(4)->is<Plus>());
    REQUIRE(rpn.tokens().at(5)->as<Number>() == ex<Constant>(3));
    REQUIRE(rpn.tokens().at(6)->as<Number>() == ex<Constant>(4));
    REQUIRE(rpn.tokens().at(7)->is<Times>());
    REQUIRE(rpn.tokens().at(8)->is<Minus>());
  }

  SECTION("add tokens: parenthesis") {
    // expression: '(1 + 2) * 3'
    rpn.reset();
    rpn.add_left_parenthesis();
    rpn.add_operand(token<Number>(Constant{1}));
    rpn.add_operator(token<Plus>());
    rpn.add_operand(token<Number>(Constant{2}));
    REQUIRE(rpn.add_right_parenthesis());
    rpn.add_operator(token<Times>());
    rpn.add_operand(token<Number>(Constant{3}));
    // post-fix notation: '1 2 + 3 *'
    REQUIRE(rpn.close());
    REQUIRE(rpn.tokens().at(0)->as<Number>() == ex<Constant>(1));
    REQUIRE(rpn.tokens().at(1)->as<Number>() == ex<Constant>(2));
    REQUIRE(rpn.tokens().at(2)->is<Plus>());
    REQUIRE(rpn.tokens().at(3)->as<Number>() == ex<Constant>(3));
    REQUIRE(rpn.tokens().at(4)->is<Times>());

    // expression: '(3 * (4 + 9 - 10) + 2) - 12'
    rpn.reset();
    rpn.add_left_parenthesis();
    rpn.add_operand(token<Number>(Constant{3}));
    rpn.add_operator(token<Times>());
    rpn.add_left_parenthesis();
    rpn.add_operand(token<Number>(Constant{4}));
    rpn.add_operator(token<Plus>());
    rpn.add_operand(token<Number>(Constant{9}));
    rpn.add_operator(token<Minus>());
    rpn.add_operand(token<Number>(Constant{10}));
    REQUIRE(rpn.add_right_parenthesis());
    rpn.add_operator(token<Plus>());
    rpn.add_operand(token<Number>(Constant{2}));
    REQUIRE(rpn.add_right_parenthesis());
    rpn.add_operator(token<Minus>());
    rpn.add_operand(token<Number>(Constant{12}));
    // post-fix notation:
    //             '3 4 9 + 10 - * 2 + 12 -'
    REQUIRE(rpn.close());
    REQUIRE(rpn.tokens().at(0)->as<Number>() == ex<Constant>(3));
    REQUIRE(rpn.tokens().at(1)->as<Number>() == ex<Constant>(4));
    REQUIRE(rpn.tokens().at(2)->as<Number>() == ex<Constant>(9));
    REQUIRE(rpn.tokens().at(4)->as<Number>() == ex<Constant>(10));
    REQUIRE(rpn.tokens().at(7)->as<Number>() == ex<Constant>(2));
    REQUIRE(rpn.tokens().at(9)->as<Number>() == ex<Constant>(12));
    REQUIRE(rpn.tokens().at(3)->is<Plus>());
    REQUIRE(rpn.tokens().at(5)->is<Minus>());
    REQUIRE(rpn.tokens().at(6)->is<Times>());
    REQUIRE(rpn.tokens().at(8)->is<Plus>());
    REQUIRE(rpn.tokens().at(10)->is<Minus>());
  }
}
