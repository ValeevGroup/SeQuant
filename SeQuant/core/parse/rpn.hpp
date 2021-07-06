//
// Created by Bimal Gaudel on 6/13/21.
//

#ifndef CPPLAB_RPN_HPP
#define CPPLAB_RPN_HPP

#include "token.hpp"

#include <stack>
#include <vector>

namespace sequant::parse {

class ReversePolishNotation {
 public:
  using token_ptr = std::unique_ptr<Token>;

  ReversePolishNotation() = default;

  [[maybe_unused]] bool add_left_parenthesis();

  [[nodiscard]] bool add_right_parenthesis();

  [[maybe_unused]] bool add_operator(token_ptr oprtr);

  [[maybe_unused]] bool add_operand(token_ptr oprnd);

  [[nodiscard]] bool close();

  [[nodiscard]] std::vector<token_ptr> const& tokens();

  void reset();

 private:
  std::stack<token_ptr> operator_stack_;

  std::vector<token_ptr> token_vec_;
};

} // namespace

#endif  // CPPLAB_RPN_HPP
