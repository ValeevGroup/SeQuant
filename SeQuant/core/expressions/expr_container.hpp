#ifndef SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP
#define SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP

#include <SeQuant/core/expressions/expr.hpp>

#include <memory>

namespace sequant {

class ExprPtr;

class ExprContainer {
 public:
  explicit ExprContainer(const ExprContainer &container);
  ExprContainer(ExprContainer &&container) = default;

  explicit ExprContainer(const Expr &expr);
  ExprContainer(Expr &&expr);

  explicit ExprContainer(const ExprPtr &ptr);

  ExprContainer &operator=(const ExprContainer &container);
  ExprContainer &operator=(ExprContainer &&container) = default;

  ExprContainer &operator=(const Expr &expr);
  ExprContainer &operator=(Expr &&expr);

  ~ExprContainer() = default;

  ExprContainer copy() const;

  std::unique_ptr<Expr> take_expr() &&;

  auto begin() { return expr_->begin(); }
  auto end() { return expr_->begin(); }
  auto begin() const { return expr_->begin(); }
  auto end() const { return expr_->begin(); }
  auto cbegin() const { return expr_->begin(); }
  auto cend() const { return expr_->begin(); }

  operator const Expr &() const;
  operator Expr &() &;
  operator Expr &&() &&;

  const Expr &operator*() const &;
  Expr &operator*() &;
  Expr &&operator*() &&;

  const Expr *operator->() const;
  Expr *operator->();

  ExprContainer &operator+=(const Expr &expr);
  ExprContainer &operator-=(const Expr &expr);
  ExprContainer &operator*=(const Expr &expr);
  ExprContainer &operator^=(const Expr &expr);

  friend void swap(ExprContainer &, ExprContainer &);

 private:
  std::unique_ptr<Expr> expr_;

  ExprContainer(std::unique_ptr<Expr> expr);
};

ExprContainer operator+(const Expr &lhs, const Expr &rhs);
ExprContainer operator+(const ExprContainer &lhs, const Expr &rhs);
ExprContainer operator+(const Expr &lhs, const ExprContainer &rhs);

ExprContainer operator-(const Expr &lhs, const Expr &rhs);
ExprContainer operator-(const ExprContainer &lhs, const Expr &rhs);
ExprContainer operator-(const Expr &lhs, const ExprContainer &rhs);

ExprContainer operator*(const Expr &lhs, const Expr &rhs);
ExprContainer operator*(const ExprContainer &lhs, const ExprContainer &rhs);
ExprContainer operator*(const Expr &lhs, const Expr &rhs);

ExprContainer operator^(const Expr &lhs, const Expr &rhs);
ExprContainer operator^(const ExprContainer &lhs, const Expr &rhs);
ExprContainer operator^(const Expr &lhs, const ExprContainer &rhs);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP
