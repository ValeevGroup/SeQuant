#ifndef SEQUANT_EXPRESSIONS_EXPR_HOLDER_HPP
#define SEQUANT_EXPRESSIONS_EXPR_HOLDER_HPP

#include <SeQuant/core/expressions/expr_iterator.hpp>

#include <initializer_list>
#include <memory>

namespace sequant {

class Expr;
class ExprPtr;

/// Holds a single expression with value semantics: copying an ExprHolder
/// copies the held expression (deep copy) and no two ExprHolder objects ever
/// share an expression. ExprHolder is a regular, std::copyable type and can
/// therefore be used in standard containers and ranges.
class ExprHolder {
 public:
  ExprHolder(const ExprHolder &holder);
  ExprHolder(ExprHolder &&holder) = default;

  explicit ExprHolder(const Expr &expr);
  ExprHolder(Expr &&expr);

  explicit ExprHolder(const ExprPtr &ptr);

  ExprHolder &operator=(const ExprHolder &holder);
  ExprHolder &operator=(ExprHolder &&holder) = default;

  ExprHolder &operator=(const Expr &expr);
  ExprHolder &operator=(Expr &&expr);

  ~ExprHolder() = default;

  ExprHolder copy() const;

  std::unique_ptr<Expr> take_expr() &&;

  ExprIterator begin();
  ExprIterator end();
  ConstExprIterator begin() const;
  ConstExprIterator end() const;
  ConstExprIterator cbegin() const;
  ConstExprIterator cend() const;

  operator const Expr &() const;
  operator Expr &() &;
  operator Expr &&() &&;

  const Expr &operator*() const &;
  Expr &operator*() &;
  Expr &&operator*() &&;

  const Expr *operator->() const;
  Expr *operator->();

  ExprHolder &operator+=(const Expr &expr);
  ExprHolder &operator-=(const Expr &expr);
  ExprHolder &operator*=(const Expr &expr);
  ExprHolder &operator^=(const Expr &expr);

  friend void swap(ExprHolder &, ExprHolder &);

 private:
  std::unique_ptr<Expr> expr_;

  ExprHolder(std::unique_ptr<Expr> expr);
};

using ExprHolderList = std::initializer_list<ExprHolder>;

ExprHolder operator+(const Expr &lhs, const Expr &rhs);
ExprHolder operator+(const ExprHolder &lhs, const Expr &rhs);
ExprHolder operator+(const Expr &lhs, const ExprHolder &rhs);

ExprHolder operator-(const Expr &lhs, const Expr &rhs);
ExprHolder operator-(const ExprHolder &lhs, const Expr &rhs);
ExprHolder operator-(const Expr &lhs, const ExprHolder &rhs);

ExprHolder operator*(const Expr &lhs, const Expr &rhs);
ExprHolder operator*(const ExprHolder &lhs, const ExprHolder &rhs);
ExprHolder operator*(const Expr &lhs, const Expr &rhs);

ExprHolder operator^(const Expr &lhs, const Expr &rhs);
ExprHolder operator^(const ExprHolder &lhs, const Expr &rhs);
ExprHolder operator^(const Expr &lhs, const ExprHolder &rhs);

bool operator==(const ExprHolder &lhs, const ExprPtr &rhs);
bool operator==(const ExprPtr &lhs, const ExprHolder &rhs);

ExprHolder adjoint(const ExprHolder &expr);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_HOLDER_HPP
