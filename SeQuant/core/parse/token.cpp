//
// Created by Bimal Gaudel on 6/12/21.
//
#include "token.hpp"
#include <atomic>

namespace sequant::parse {

Token::type_id_type Token::get_next_type_id() {
  static std::atomic<type_id_type> grand_type_id = 0;
  return ++grand_type_id;
}

Token::type_id_type LeftParenthesis::type_id() const {
  return Token::get_type_id<LeftParenthesis>();
}

Token::type_id_type RightParenthesis::type_id() const {
  return Token::get_type_id<RightParenthesis>();
}

}  // namespace sequant::parse
