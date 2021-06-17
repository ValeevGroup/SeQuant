//
// Created by Bimal Gaudel on 6/15/21.
//

#ifndef SEQUANT_PARSE_TOKEN_SEQUANT_HPP
#define SEQUANT_PARSE_TOKEN_SEQUANT_HPP

#include "token.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant::parse {

struct OperatorTimes : public Operator {
  [[nodiscard]] int precedence() const override;

  [[nodiscard]] Token::type_id_type type_id() const override;
};

struct OperatorPlus : public Operator {
  [[nodiscard]] int precedence() const override;

  [[nodiscard]] Token::type_id_type type_id() const override;
};

struct OperatorMinus : public Operator {
  [[nodiscard]] int precedence() const override;

  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperandSequant: public Operand {
 public:
  ~OperandSequant() override = default;

  [[nodiscard]] virtual operator ExprPtr() const = 0;

 protected:
  OperandSequant() = default;
};


class OperandConstant final: public OperandSequant {
 public:
  explicit OperandConstant(Constant const& expr);

  [[nodiscard]] operator ExprPtr() const override;

 private:
  ExprPtr expr_;

  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperandTensor final : public OperandSequant {
 public:
  explicit OperandTensor(Tensor const & expr);

  [[nodiscard]] operator ExprPtr() const override;

 private:
  ExprPtr expr_;

  [[nodiscard]] Token::type_id_type type_id() const override;
};

} // namespace

#endif  // SEQUANT_PARSE_TOKEN_SEQUANT_HPP
