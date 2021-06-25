//
// Created by Bimal Gaudel on 6/15/21.
//

#ifndef SEQUANT_PARSE_TOKEN_SEQUANT_HPP
#define SEQUANT_PARSE_TOKEN_SEQUANT_HPP

#include "token.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant::parse {

class OperatorTimes : public Operator {
 public:
  [[nodiscard]] int precedence() const override;

 private:
  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperatorPlus : public Operator {
 public:
  [[nodiscard]] int precedence() const override;

 private:
  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperatorMinus : public Operator {
 public:
  [[nodiscard]] int precedence() const override;

 private:
  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperatorPlusUnary: public Operator {
 public:
  [[nodiscard]] int precedence() const override;

 private:
  [[nodiscard]] Token::type_id_type type_id() const override;
};

class OperatorMinusUnary: public Operator {
 public:
  [[nodiscard]] int precedence() const override;

 private:
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
