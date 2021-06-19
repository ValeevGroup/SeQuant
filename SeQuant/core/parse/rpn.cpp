//
// Created by Bimal Gaudel on 6/13/21.
//

#include "rpn.hpp"

#include <cassert>

namespace sequant::parse {

bool ReversePolishNotation::add_operator(token_ptr oprtr) {
  assert(oprtr->is<Operator>());
  while (!operator_stack_.empty()
         && !operator_stack_.top()->is<LeftParenthesis>()
         && operator_stack_.top()->as<Operator>().precedence()
         >= oprtr->as<Operator>().precedence()) {
    token_vec_.emplace_back(std::move(operator_stack_.top()));
    operator_stack_.pop();
  }
  operator_stack_.emplace(std::move(oprtr));
  return true;
}

bool ReversePolishNotation::add_operand(token_ptr oprnd) {
  assert(oprnd->is<Operand>());
  token_vec_.emplace_back(std::move(oprnd));
  return true;
}

bool ReversePolishNotation::add_left_parenthesis() {
  operator_stack_.emplace(token<LeftParenthesis>());
  return true;
}
bool ReversePolishNotation::add_right_parenthesis() {
  while (!operator_stack_.empty()
         && !operator_stack_.top()->is<LeftParenthesis>()){
    token_vec_.emplace_back(std::move(operator_stack_.top()));
    operator_stack_.pop();
  }
  if (operator_stack_.empty()) return false;
  operator_stack_.pop();
  return true;
}

bool ReversePolishNotation::close() {
  while (!operator_stack_.empty()){
    if (operator_stack_.top()->is<LeftParenthesis>()) return false;
    token_vec_.emplace_back(std::move(operator_stack_.top()));
    operator_stack_.pop();
  }
  return true;
}
const std::vector<ReversePolishNotation::token_ptr>&
    ReversePolishNotation::tokens() {
  return token_vec_;
}

void ReversePolishNotation::reset() {
  token_vec_.clear();
  operator_stack_ = decltype(operator_stack_){};
}

} // namespace