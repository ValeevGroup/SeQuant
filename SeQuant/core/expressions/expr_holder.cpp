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
#include <utility>

namespace sequant {

namespace {

const Expr &deref(const ExprPtr &ptr) {
  SEQUANT_ASSERT(ptr != nullptr);
  return *ptr;
}

}  // namespace

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

ExprHolder::ExprHolder(const ExprPtr &ptr)
    : ExprHolder(deref(ptr).unique_copy()) {}

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

ExprHolder ExprHolder::copy() const { return get().unique_copy(); }

std::unique_ptr<Expr> ExprHolder::take_expr() && {
  SEQUANT_ASSERT(expr_ != nullptr);
  return std::move(expr_);
}

Expr &ExprHolder::get() {
  SEQUANT_ASSERT(expr_ != nullptr &&
                 "Attempted to access the expression of a moved-from "
                 "ExprHolder");
  return *expr_;
}

const Expr &ExprHolder::get() const {
  SEQUANT_ASSERT(expr_ != nullptr &&
                 "Attempted to access the expression of a moved-from "
                 "ExprHolder");
  return *expr_;
}

ExprIterator ExprHolder::begin() { return get().begin(); }

ExprIterator ExprHolder::end() { return get().end(); }

ConstExprIterator ExprHolder::begin() const {
  return std::as_const(get()).begin();
}

ConstExprIterator ExprHolder::end() const { return std::as_const(get()).end(); }

ConstExprIterator ExprHolder::cbegin() const { return begin(); }

ConstExprIterator ExprHolder::cend() const { return end(); }

ExprHolder::operator const Expr &() const { return get(); }

ExprHolder::operator Expr &() & { return get(); }
ExprHolder::operator Expr &&() && { return std::move(get()); }

const Expr &ExprHolder::operator*() const & { return get(); }

Expr &ExprHolder::operator*() & { return get(); }

Expr &&ExprHolder::operator*() && { return std::move(get()); }

const Expr *ExprHolder::operator->() const { return &get(); }

Expr *ExprHolder::operator->() { return &get(); }

// N.B. in the operators below, the held expression is kept alive (in held)
// until the result has been formed, as the other operand may refer to (a part
// of) it

ExprHolder &ExprHolder::operator+=(const Expr &expr) {
  if (get().is<Sum>()) {
    get().as<Sum>() += expr;
  } else if (get().is<Constant>() && expr.is<Constant>()) {
    *this = get().as<Constant>() + expr.as<Constant>();
  } else {
    const ExprPtr held = to_expr_ptr(std::move(*this));
    Sum sum(Sum::summands_type{held}, Sum::move_only_tag{});
    sum += expr;
    *this = std::move(sum);
  }

  return *this;
}

ExprHolder &ExprHolder::operator-=(const Expr &expr) {
  if (get().is<Sum>()) {
    get().as<Sum>() -= expr;
  } else if (get().is<Constant>() && expr.is<Constant>()) {
    *this = get().as<Constant>() - expr.as<Constant>();
  } else {
    const ExprPtr held = to_expr_ptr(std::move(*this));
    Sum sum(Sum::summands_type{held}, Sum::move_only_tag{});
    sum -= expr;
    *this = std::move(sum);
  }

  return *this;
}

ExprHolder &ExprHolder::operator*=(const Expr &expr) {
  if (get().is<Product>()) {
    get().as<Product>() *= expr;
  } else if (get().is<Constant>() && expr.is<Constant>()) {
    *this = get().as<Constant>() * expr.as<Constant>();
  } else {
    const ExprPtr held = to_expr_ptr(std::move(*this));
    Product product(ExprPtrList{held});
    product *= expr;
    *this = std::move(product);
  }
  return *this;
}

ExprHolder &ExprHolder::operator^=(const Expr &other) {
  const ExprPtr held = to_expr_ptr(std::move(*this));
  if (held->is<Product>()) {
    // N.B. the (shallow) copy shares the factors with held
    NCProduct result(std::as_const(held->as<Product>()));
    result.append(1, other);
    *this = std::move(result);
  } else if (other.is<Product>()) {
    NCProduct result(other.as<Product>().deep_copy());
    result.prepend(1, *held);
    *this = std::move(result);
  } else {
    NCProduct result(ExprPtrList{held});
    result.append(1, other);
    *this = std::move(result);
  }

  return *this;
}

void swap(ExprHolder &lhs, ExprHolder &rhs) { std::swap(lhs.expr_, rhs.expr_); }

ExprHolder operator+(const Expr &lhs, const Expr &rhs) {
  return ExprHolder(lhs) + rhs;
}

ExprHolder operator+(ExprHolder lhs, const Expr &rhs) {
  lhs += rhs;

  return lhs;
}

ExprHolder operator+(const Expr &lhs, const ExprHolder &rhs) {
  return lhs + *rhs;
}

ExprHolder operator+(ExprHolder lhs, const ExprHolder &rhs) {
  return std::move(lhs) + *rhs;
}

ExprHolder operator-(const Expr &lhs, const Expr &rhs) {
  return ExprHolder(lhs) - rhs;
}

ExprHolder operator-(ExprHolder lhs, const Expr &rhs) {
  lhs -= rhs;

  return lhs;
}

ExprHolder operator-(const Expr &lhs, const ExprHolder &rhs) {
  return lhs - *rhs;
}

ExprHolder operator-(ExprHolder lhs, const ExprHolder &rhs) {
  return std::move(lhs) - *rhs;
}

ExprHolder operator*(const Expr &lhs, const Expr &rhs) {
  return ExprHolder(lhs) * rhs;
}

ExprHolder operator*(ExprHolder lhs, const Expr &rhs) {
  lhs *= rhs;

  return lhs;
}

ExprHolder operator*(const Expr &lhs, const ExprHolder &rhs) {
  return lhs * *rhs;
}

ExprHolder operator*(ExprHolder lhs, const ExprHolder &rhs) {
  return std::move(lhs) * *rhs;
}

ExprHolder operator^(const Expr &lhs, const Expr &rhs) {
  return ExprHolder(lhs) ^ rhs;
}

ExprHolder operator^(ExprHolder lhs, const Expr &rhs) {
  lhs ^= rhs;

  return lhs;
}

ExprHolder operator^(const Expr &lhs, const ExprHolder &rhs) {
  return lhs ^ *rhs;
}

ExprHolder operator^(ExprHolder lhs, const ExprHolder &rhs) {
  return std::move(lhs) ^ *rhs;
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
