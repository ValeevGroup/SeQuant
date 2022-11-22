//
// Created by Bimal Gaudel on 6/12/21.
//

#ifndef SEQUANT_PARSE_EXPR_TOKEN_HPP
#define SEQUANT_PARSE_EXPR_TOKEN_HPP

#include <cassert>
#include <memory>
#include <string>

namespace sequant::parse {

class Operator;

class Operand;

class Token {
 public:
  using type_id_type = int;

  virtual ~Token() = default;

 protected:
  Token() = default;

  [[nodiscard]] virtual type_id_type type_id() const = 0;

  static type_id_type get_next_type_id();

  /// sets (unique) type id of class T
  /// @param id the value of type id of class T
  template <typename T>
  static type_id_type &type_id_accessor() {
    static type_id_type type_id = get_next_type_id();
    return type_id;
  };

  /// @return (unique) type id of class T
  template <typename T>
  static type_id_type get_type_id() {
    return type_id_accessor<T>();
  }

 public:
  /// @tparam T a Token type
  /// @return true if this object is of type @c T
  template <typename T>
  [[nodiscard]] bool is() const {
    if constexpr (std::is_same_v<std::decay_t<T>, Token>)
      return true;
    else
      return this->type_id() == get_type_id<std::decay_t<T>>() ||
             dynamic_cast<std::decay_t<T> const *>(this);
  }

  /// @tparam T a Token type
  /// @return this object cast to type @c T
  template <typename T>
  const T &as() const {
    assert(this->is<std::decay_t<T>>());  // so that as<const T>() works fine
    return static_cast<const T &>(*this);
  }
};

class Operator : public Token {
 public:
  ~Operator() override = default;

  [[nodiscard]] virtual int precedence() const = 0;

 protected:
  Operator() = default;
};

class Operand : public Token {
 public:
  ~Operand() override = default;

 protected:
  Operand() = default;
};

struct LeftParenthesis : public Token {
  [[nodiscard]] Token::type_id_type type_id() const override;
};

struct RightParenthesis : public Token {
  [[nodiscard]] Token::type_id_type type_id() const override;
};

struct TokenDummy: public Token {
  [[nodiscard]] Token::type_id_type type_id() const override;
};

template <typename T, typename... Args>
std::unique_ptr<Token> token(Args &&...args) {
  return std::make_unique<T>(std::forward<Args>(args)...);
}

}  // namespace sequant::parse

#endif  // SEQUANT_PARSE_EXPR_TOKEN_HPP
