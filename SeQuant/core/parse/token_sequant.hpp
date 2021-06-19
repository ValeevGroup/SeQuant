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

class OperandConstant final: public Operand, public Constant {
 public:
  using Constant::Constant;
  using Operand::is;
  using Operand::as;

 private:
  Token::type_id_type type_id() const override;
};

class OperandTensor final: public Operand, public Tensor {
 public:
  using Tensor::Tensor;
  using Operand::is;
  using Operand::as;

 private:
  Token::type_id_type type_id() const override;
};

} // namespace

#endif  // SEQUANT_PARSE_TOKEN_SEQUANT_HPP
