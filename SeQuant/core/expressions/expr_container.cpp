#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_container.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <memory>

namespace sequant {

ExprPtr to_expr_ptr(ExprContainer &&container) {
  return std::shared_ptr<Expr>(std::move(container.expr_));
}

ExprContainer::ExprContainer(const ExprContainer &container)
    : ExprContainer(container->unique_copy()) {}

ExprContainer::ExprContainer(std::unique_ptr<Expr> expr)
    : expr_(std::move(expr)) {}

ExprContainer::ExprContainer(const Expr &expr)
    : ExprContainer(expr.unique_copy()) {}

ExprContainer::ExprContainer(Expr &&expr)
    : ExprContainer(std::move(expr).unique_copy()) {}

ExprContainer &ExprContainer::operator=(const ExprContainer &container) {
  // Copy-and-swap
  ExprContainer copy(container);

  swap(*this, copy);

  return *this;
}

ExprContainer &ExprContainer::operator=(const Expr &expr) {
  // Copy-and-swap
  ExprContainer copy(expr);

  swap(*this, copy);

  return *this;
}

ExprContainer &ExprContainer::operator=(Expr &&expr) {
  // Copy-and-swap
  ExprContainer copy(std::move(expr));

  swap(*this, copy);

  return *this;
}

ExprContainer ExprContainer::copy() const { return expr_->unique_copy(); }

ExprContainer::operator const Expr &() const { return *expr_; }

ExprContainer::operator Expr &() & { return *expr_; }
ExprContainer::operator Expr &&() && { return std::move(*expr_); }

const Expr &ExprContainer::operator*() const { return *expr_; }

Expr &ExprContainer::operator*() { return *expr_; }

const Expr *ExprContainer::operator->() const { return expr_.get(); }

Expr *ExprContainer::operator->() { return expr_.get(); }

ExprContainer &ExprContainer::operator+=(const Expr &expr) {
  if (expr_->is<Sum>()) {
    expr_->as<Sum>() += expr;
  } else if (expr_->is<Constant>() && expr.is<Constant>()) {
    *this = expr_->as<Constant>() + expr.as<Constant>();
  } else {
    *this = Sum(ExprPtrList{to_expr_ptr(std::move(*this)), expr.clone()});
  }

  return *this;
}

ExprContainer &ExprContainer::operator-=(const Expr &expr) {
  if (expr_->is<Sum>()) {
    expr_->as<Sum>() -= expr;
  } else if (expr_->is<Constant>() && expr.is<Constant>()) {
    *this = expr_->as<Constant>() - expr.as<Constant>();
  } else {
    *this = Sum(ExprPtrList{to_expr_ptr(std::move(*this)),
                            ex<Product>(-1, ExprPtrList{expr.clone()})});
  }

  return *this;
}

ExprContainer &ExprContainer::operator*=(const Expr &expr) {
  if (expr_->is<Product>()) {
    expr_->as<Product>() *= expr;
  } else if (expr_->is<Constant>() && expr.is<Constant>()) {
    *this = expr_->as<Constant>() * expr.as<Constant>();
  } else {
    *this = Product(ExprPtrList{to_expr_ptr(std::move(*this)), expr.clone()});
  }
  return *this;
}

void swap(ExprContainer &lhs, ExprContainer &rhs) {
  std::swap(lhs.expr_, rhs.expr_);
}

ExprContainer operator+(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont += rhs;

  return cont;
}

ExprContainer operator-(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont -= rhs;

  return cont;
}

ExprContainer operator*(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont *= rhs;

  return cont;
}

}  // namespace sequant
