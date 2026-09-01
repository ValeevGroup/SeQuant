#include <SeQuant/core/expressions/complex.hpp>

namespace sequant {

ExprPtr RealPart::clone() const { return ex<RealPart>(inner_->clone()); }

std::wstring RealPart::to_latex() const {
  return L"\\Re\\left[" + inner_->to_latex() + L"\\right]";
}

Expr::hash_type RealPart::memoizing_hash() const {
  auto compute = [this]() {
    auto v = hash::value(*inner_);
    hash::combine(v, std::size_t{0xC0FFEE01ull});
    return v;
  };
  if (!hash_value_) hash_value_ = compute();
  return *hash_value_;
}

bool RealPart::static_equal(const Expr& that) const {
  return *inner_ == *static_cast<const RealPart&>(that).inner_;
}

bool RealPart::static_less_than(const Expr& that) const {
  return *inner_ < *static_cast<const RealPart&>(that).inner_;
}

ExprPtr ImagPart::clone() const { return ex<ImagPart>(inner_->clone()); }

std::wstring ImagPart::to_latex() const {
  return L"\\Im\\left[" + inner_->to_latex() + L"\\right]";
}

Expr::hash_type ImagPart::memoizing_hash() const {
  auto compute = [this]() {
    auto v = hash::value(*inner_);
    hash::combine(v, std::size_t{0xC0FFEE02ull});
    return v;
  };
  if (!hash_value_) hash_value_ = compute();
  return *hash_value_;
}

bool ImagPart::static_equal(const Expr& that) const {
  return *inner_ == *static_cast<const ImagPart&>(that).inner_;
}

bool ImagPart::static_less_than(const Expr& that) const {
  return *inner_ < *static_cast<const ImagPart&>(that).inner_;
}

namespace detail {
ExprPtr strip_scalar(const Product& prod) {
  auto rest = std::make_shared<Product>();
  for (const auto& f : prod.factors()) rest->append(1, f, Product::Flatten::No);
  return rest;
}
}  // namespace detail

ExprPtr real_part(ExprPtr expr) {
  SEQUANT_ASSERT(expr);
  if (expr->is<Constant>())
    return ex<Constant>(expr->as<Constant>().value().real());
  if (expr->is<RealPart>() || expr->is<ImagPart>()) return expr;
  if (expr->is<Product>()) {
    const auto& prod = expr->as<Product>();
    const auto c = prod.scalar();
    if (c.imag() == 0 && c.real() != 1)
      return ex<Constant>(c) * real_part(detail::strip_scalar(prod));
    if (c.real() == 0 && c.imag() != 0)
      return ex<Constant>(-c.imag()) *
             imaginary_part(detail::strip_scalar(prod));
  }
  return ex<RealPart>(std::move(expr));
}

ExprPtr imaginary_part(ExprPtr expr) {
  SEQUANT_ASSERT(expr);
  if (expr->is<Constant>())
    return ex<Constant>(expr->as<Constant>().value().imag());
  if (expr->is<RealPart>() || expr->is<ImagPart>()) return ex<Constant>(0);
  if (expr->is<Product>()) {
    const auto& prod = expr->as<Product>();
    const auto c = prod.scalar();
    if (c.imag() == 0 && c.real() != 1)
      return ex<Constant>(c) * imaginary_part(detail::strip_scalar(prod));
    if (c.real() == 0 && c.imag() != 0)
      return ex<Constant>(c.imag()) * real_part(detail::strip_scalar(prod));
  }
  return ex<ImagPart>(std::move(expr));
}

}  // namespace sequant
