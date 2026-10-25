#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_holder.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <memory>

namespace sequant {

ExprPtr to_expr_ptr(ExprHolder &&holder) {
  return std::move(holder).take_expr();
}

ExprHolder::ExprHolder(const ExprHolder &holder)
    : ExprHolder(holder->unique_copy()) {}

ExprHolder::ExprHolder(std::unique_ptr<Expr> expr) : expr_(std::move(expr)) {
  SEQUANT_ASSERT(expr_ != nullptr);
}

ExprHolder::ExprHolder(const Expr &expr) : ExprHolder(expr.unique_copy()) {}

ExprHolder::ExprHolder(Expr &&expr)
    : ExprHolder(std::move(expr).unique_move()) {}

ExprHolder::ExprHolder(const ExprPtr &ptr) : ExprHolder(ptr->unique_copy()) {}

ExprHolder &ExprHolder::operator=(const ExprHolder &holder) {
  // Copy-and-swap
  ExprHolder copy(holder);

  swap(*this, copy);

  return *this;
}

ExprHolder &ExprHolder::operator=(const Expr &expr) {
  // Copy-and-swap
  ExprHolder copy(expr);

  swap(*this, copy);

  return *this;
}

ExprHolder &ExprHolder::operator=(Expr &&expr) {
  // Copy-and-swap
  ExprHolder copy(std::move(expr));

  swap(*this, copy);

  return *this;
}

ExprHolder ExprHolder::copy() const { return expr_->unique_copy(); }

std::unique_ptr<Expr> ExprHolder::take_expr() && { return std::move(expr_); }

ExprIterator ExprHolder::begin() {
  SEQUANT_ASSERT(expr_);
  return expr_->begin();
}

ExprIterator ExprHolder::end() {
  SEQUANT_ASSERT(expr_);
  return expr_->end();
}

ConstExprIterator ExprHolder::begin() const {
  SEQUANT_ASSERT(expr_);
  return std::as_const(*expr_).begin();
}

ConstExprIterator ExprHolder::end() const {
  SEQUANT_ASSERT(expr_);
  return std::as_const(*expr_).end();
}

ConstExprIterator ExprHolder::cbegin() const { return begin(); }

ConstExprIterator ExprHolder::cend() const { return end(); }

ExprHolder::operator const Expr &() const { return *expr_; }

ExprHolder::operator Expr &() & { return *expr_; }
ExprHolder::operator Expr &&() && { return std::move(*expr_); }

const Expr &ExprHolder::operator*() const & { return *expr_; }

Expr &ExprHolder::operator*() & { return *expr_; }

Expr &&ExprHolder::operator*() && { return std::move(*expr_); }

const Expr *ExprHolder::operator->() const { return expr_.get(); }

Expr *ExprHolder::operator->() { return expr_.get(); }

ExprHolder &ExprHolder::operator+=(const Expr &expr) {
  if (expr_->is<Sum>()) {
    expr_->as<Sum>() += expr;
  } else if (expr_->is<Constant>() && expr.is<Constant>()) {
    *this = expr_->as<Constant>() + expr.as<Constant>();
  } else {
    *this = Sum(ExprPtrList{to_expr_ptr(std::move(*this)), expr.clone()});
  }

  return *this;
}

ExprHolder &ExprHolder::operator-=(const Expr &expr) {
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

ExprHolder &ExprHolder::operator*=(const Expr &expr) {
  if (expr_->is<Product>()) {
    expr_->as<Product>() *= expr;
  } else if (expr_->is<Constant>() && expr.is<Constant>()) {
    *this = expr_->as<Constant>() * expr.as<Constant>();
  } else {
    *this = Product(ExprPtrList{to_expr_ptr(std::move(*this)), expr.clone()});
  }
  return *this;
}

ExprHolder &ExprHolder::operator^=(const Expr &other) {
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

void swap(ExprHolder &lhs, ExprHolder &rhs) { std::swap(lhs.expr_, rhs.expr_); }

ExprHolder operator+(const Expr &lhs, const Expr &rhs) {
  ExprHolder holder(lhs);
  holder += rhs;

  return holder;
}

ExprHolder operator+(const ExprHolder &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) + rhs;
}

ExprHolder operator+(const Expr &lhs, const ExprHolder &rhs) {
  return lhs + static_cast<const Expr &>(rhs);
}

ExprHolder operator-(const Expr &lhs, const Expr &rhs) {
  ExprHolder holder(lhs);
  holder -= rhs;

  return holder;
}

ExprHolder operator-(const ExprHolder &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) - rhs;
}

ExprHolder operator-(const Expr &lhs, const ExprHolder &rhs) {
  return lhs - static_cast<const Expr &>(rhs);
}

ExprHolder operator*(const Expr &lhs, const Expr &rhs) {
  ExprHolder holder(lhs);
  holder *= rhs;

  return holder;
}

ExprHolder operator*(const ExprHolder &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) * rhs;
}

ExprHolder operator*(const Expr &lhs, const ExprHolder &rhs) {
  return lhs * static_cast<const Expr &>(rhs);
}

ExprHolder operator^(const Expr &lhs, const Expr &rhs) {
  ExprHolder holder(lhs);
  holder ^= rhs;

  return holder;
}

ExprHolder operator^(const ExprHolder &lhs, const Expr &rhs) {
  return static_cast<const Expr &>(lhs) ^ rhs;
}

ExprHolder operator^(const Expr &lhs, const ExprHolder &rhs) {
  return lhs ^ static_cast<const Expr &>(rhs);
}

bool operator==(const ExprHolder &lhs, const ExprPtr &rhs) {
  return *lhs == *rhs;
}

bool operator==(const ExprPtr &lhs, const ExprHolder &rhs) {
  return *lhs == *rhs;
}

ExprHolder adjoint(const ExprHolder &holder) {
  ExprHolder copy = holder.copy();
  copy->adjoint();
  return copy;
}

}  // namespace sequant
