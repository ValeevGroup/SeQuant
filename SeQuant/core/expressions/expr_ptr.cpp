#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_container.hpp>
#include <SeQuant/core/expressions/expr_operators.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>

namespace sequant {

ExprPtr::ExprPtr(const ExprContainer &container) : ExprPtr(container.copy()) {}

ExprPtr::ExprPtr(ExprContainer &&container)
    : ExprPtr(std::move(container).take_expr()) {}

ExprPtr ExprPtr::clone() const & {
  if (!*this) return {};
  return ExprPtr(as_shared_ptr()->clone());
}

ExprPtr ExprPtr::clone() && noexcept { return std::move(*this); }

ExprPtr::base_type &ExprPtr::as_shared_ptr() & {
  return static_cast<base_type &>(*this);
}
const ExprPtr::base_type &ExprPtr::as_shared_ptr() const & {
  return static_cast<const base_type &>(*this);
}
ExprPtr::base_type &&ExprPtr::as_shared_ptr() && {
  return static_cast<base_type &&>(*this);
}

Expr &ExprPtr::operator*() & {
  SEQUANT_ASSERT(this->operator bool());
  return *(this->get());
}

const Expr &ExprPtr::operator*() const & {
  SEQUANT_ASSERT(this->operator bool());
  return *(this->get());
}

Expr &&ExprPtr::operator*() && {
  SEQUANT_ASSERT(this->operator bool());
  return std::move(*(this->get()));
}

ExprPtr &ExprPtr::operator+=(const ExprPtr &other) {
  if (!other) return *this;

  if (!*this) {
    *this = other.clone();
  } else if (as_shared_ptr()->is<Sum>()) {
    as<Sum>() += *other;
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() +
                         other->as<Constant>().value());
  } else {
    *this = ex<Sum>(ExprPtrList{*this, other});
  }
  return *this;
}

ExprPtr &ExprPtr::operator-=(const ExprPtr &other) {
  if (!other) return *this;

  if (!*this) {
    *this = ex<Constant>(-1) * other.clone();
  } else if (as_shared_ptr()->is<Sum>()) {
    as<Sum>() -= *other;
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() -
                         other->as<Constant>().value());
  } else {
    *this = ex<Sum>(ExprPtrList{*this, ex<Product>(-1, ExprPtrList{other})});
  }
  return *this;
}

ExprPtr &ExprPtr::operator*=(const ExprPtr &other) {
  if (!other) return *this;

  if (!*this) {
    *this = other.clone();
  } else if (as_shared_ptr()->is<Product>()) {
    as<Product>() *= *other;
  } else if (as_shared_ptr()->is<Constant>() && other->is<Constant>()) {
    *this = ex<Constant>(this->as<Constant>().value() *
                         other->as<Constant>().value());
  } else {
    *this = ex<Product>(ExprPtrList{*this, other});
  }
  return *this;
}

std::size_t ExprPtr::size() const { return this->get()->size(); }

std::wstring ExprPtr::to_latex() const { return as_shared_ptr()->to_latex(); }

ExprPtr adjoint(const ExprPtr &expr) {
  auto result = expr->clone();
  result->adjoint();
  return result;
}

bool operator==(const ExprPtr &left, const ExprPtr &right) {
  return *left == *right;
}

ExprPtr operator*(const ExprPtr &left, const ExprPtr &right) {
  if (left.is<Constant>() && right.is<Constant>()) {
    auto c_ = left->clone();
    auto &c = c_.as<Constant>();
    c *= right.as<Constant>();
    return c_;
  }

  auto left_is_product = left->is<Product>();
  auto right_is_product = right->is<Product>();
  if (!left_is_product && !right_is_product) {
    return ex<Product>(ExprPtrList{left, right});
  } else if (left_is_product) {
    auto result = std::static_pointer_cast<Product>(left->clone());
    result->append(1, right);
    return result;
  } else {  // right_is_product
    auto result = std::static_pointer_cast<Product>(right->clone());
    result->prepend(1, left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

/// Unlike @code operator*(const ExprPtr&, const ExprPtr&) @endcode this
/// produces a non-commutative product (i.e. NCProduct)
ExprPtr operator^(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_product = left->is<Product>();
  auto right_is_product = right->is<Product>();
  if (!left_is_product && !right_is_product) {
    return ex<NCProduct>(ExprPtrList{left, right});
  } else if (left_is_product) {
    auto result = std::make_shared<NCProduct>(left->clone().as<Product>());
    result->append(1, right);
    return result;
  } else {  // right_is_product
    auto result = std::make_shared<NCProduct>(right->clone().as<Product>());
    result->prepend(1, left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

ExprPtr operator+(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  auto right_is_sum = right->is<Sum>();
  if (!left_is_sum && !right_is_sum) {
    return ex<Sum>(ExprPtrList{left, right});
  } else if (left_is_sum) {
    auto result = std::static_pointer_cast<Sum>(left->clone());
    result->append(right);
    return result;
  } else {  // right_is_sum
    auto result = std::static_pointer_cast<Sum>(right->clone());
    result->prepend(left);
    return result;
  }

  SEQUANT_UNREACHABLE;
}

ExprPtr operator-(const ExprPtr &left, const ExprPtr &right) {
  auto left_is_sum = left->is<Sum>();
  if (!left_is_sum) {
    return ex<Sum>(ExprPtrList{
        left,
        (right->is<Constant>() ? ex<Constant>(-right->as<Constant>().value())
                               : ex<Product>(-1, ExprPtrList{right}))});
  } else if (left_is_sum) {
    auto result = std::static_pointer_cast<Sum>(left->clone());
    if (right->is<Constant>())
      result->append(ex<Constant>(-right->as<Constant>().value()));
    else
      result->append(ex<Product>(-1, ExprPtrList{right}));
    return result;
  }

  SEQUANT_UNREACHABLE;
}

}  // namespace sequant
