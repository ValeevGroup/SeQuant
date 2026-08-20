#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_container.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <memory>

namespace sequant {

ExprPtr to_expr_ptr(ExprContainer &&container) {
  return std::move(container).take_expr();
}

ExprContainer::ExprContainer(const ExprContainer &container)
    : ExprContainer(container->unique_copy()) {}

ExprContainer::ExprContainer(std::unique_ptr<Expr> expr)
    : expr_(std::move(expr)) {
  SEQUANT_ASSERT(expr_ != nullptr);
}

ExprContainer::ExprContainer(const Expr &expr)
    : ExprContainer(expr.unique_copy()) {}

ExprContainer::ExprContainer(Expr &&expr)
    : ExprContainer(std::move(expr).unique_copy()) {}

ExprContainer::ExprContainer(const ExprPtr &ptr)
    : ExprContainer(ptr->unique_copy()) {}

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

std::unique_ptr<Expr> ExprContainer::take_expr() && { return std::move(expr_); }

ExprIterator ExprContainer::begin() {
  SEQUANT_ASSERT(expr_);
  return expr_->begin();
}

ExprIterator ExprContainer::end() {
  SEQUANT_ASSERT(expr_);
  return expr_->end();
}

ConstExprIterator ExprContainer::begin() const {
  SEQUANT_ASSERT(expr_);
  return std::as_const(*expr_).begin();
}

ConstExprIterator ExprContainer::end() const {
  SEQUANT_ASSERT(expr_);
  return std::as_const(*expr_).end();
}

ConstExprIterator ExprContainer::cbegin() const { return begin(); }

ConstExprIterator ExprContainer::cend() const { return end(); }

ExprContainer::operator const Expr &() const { return *expr_; }

ExprContainer::operator Expr &() & { return *expr_; }
ExprContainer::operator Expr &&() && { return std::move(*expr_); }

const Expr &ExprContainer::operator*() const & { return *expr_; }

Expr &ExprContainer::operator*() & { return *expr_; }

Expr &&ExprContainer::operator*() && { return std::move(*expr_); }

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

ExprContainer &ExprContainer::operator^=(const Expr &other) {
  auto this_is_product = expr_->is<Product>();
  auto other_is_product = other.is<Product>();
  if (!this_is_product && !other_is_product) {
    *this =
        NCProduct(ExprPtrList{to_expr_ptr(std::move(*this)), other.clone()});
  } else if (this_is_product) {
    *this = NCProduct(std::move(expr_->as<Product>()));
    expr_->as<NCProduct>().append(1, other.clone());
  } else {  // other_is_product
    NCProduct result(other.clone().as<Product>());
    result.prepend(1, to_expr_ptr(std::move(*this)));
    *this = std::move(result);
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

ExprContainer operator+(const ExprContainer &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) + rhs;
}

ExprContainer operator+(const Expr &lhs, const ExprContainer &rhs) {
  return lhs + static_cast<const Expr &>(rhs);
}

ExprContainer operator-(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont -= rhs;

  return cont;
}

ExprContainer operator-(const ExprContainer &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) - rhs;
}

ExprContainer operator-(const Expr &lhs, const ExprContainer &rhs) {
  return lhs - static_cast<const Expr &>(rhs);
}

ExprContainer operator*(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont *= rhs;

  return cont;
}

ExprContainer operator*(const ExprContainer &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) * rhs;
}

ExprContainer operator*(const Expr &lhs, const ExprContainer &rhs) {
  return lhs * static_cast<const Expr &>(rhs);
}

ExprContainer operator^(const Expr &lhs, const Expr &rhs) {
  ExprContainer cont(lhs);
  cont ^= rhs;

  return cont;
}

ExprContainer operator^(const ExprContainer &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) ^ rhs;
}

ExprContainer operator^(const Expr &lhs, const ExprContainer &rhs) {
  return lhs ^ static_cast<const Expr &>(rhs);
}

bool operator==(const ExprContainer &lhs, const ExprPtr &rhs) {
  return *lhs == *rhs;
}

bool operator==(const ExprPtr &lhs, const ExprContainer &rhs) {
  return *lhs == *rhs;
}

ExprContainer adjoint(const ExprContainer &cont) {
  ExprContainer copy = cont.copy();
  copy->adjoint();
  return copy;
}

}  // namespace sequant
