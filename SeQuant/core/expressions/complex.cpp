#include <SeQuant/core/expressions/complex.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_operators.hpp>

namespace sequant {

namespace {
/// Takes a real scalar out of @p inner in place, as `Re(c X) = c Re(X)` and
/// `Im(c X) = c Im(X)` allow for a real @c c. A complex scalar stays put:
/// neither wrapper is linear over it.
/// @return the scalar taken out, or nullptr if there was none to take
ExprPtr hoist_real_scalar(ExprPtr& inner) {
  if (!inner->is<Product>()) return {};
  const auto& prod = inner->as<Product>();
  const auto scalar = prod.scalar();
  if (scalar.imag() != 0 || scalar.real() == 1) return {};
  inner = detail::strip_scalar(prod);
  return ex<Constant>(scalar);
}
}  // namespace

ExprPtr RealPart::clone() const { return ex<RealPart>(inner_->clone()); }

std::wstring RealPart::to_latex() const {
  return L"\\Re\\left[" + inner_->to_latex() + L"\\right]";
}

ExprPtr RealPart::canonicalize(CanonicalizeOptions opts) {
  // a Constant byproduct of the inner canonicalization multiplies the wrapped
  // expression back, and the eager hoist then takes a real scalar out again:
  // `Re(c X)` is `c Re(X)` for a real c, so that scalar leaves as this node's
  // own byproduct and no canonical wrapper holds a scaled Product. A complex
  // scalar stays inside, the wrapper not being linear over it
  if (auto byproduct = inner_->canonicalize(opts);
      byproduct && byproduct->is<Constant>())
    inner_ = byproduct * inner_;
  ExprPtr hoisted = hoist_real_scalar(inner_);
  reset_hash_value();
  return hoisted;
}

ExprPtr RealPart::rapid_canonicalize(CanonicalizeOptions opts) {
  SEQUANT_ASSERT(opts.method == CanonicalizationMethod::Rapid);
  if (auto byproduct = inner_->rapid_canonicalize(opts);
      byproduct && byproduct->is<Constant>())
    inner_ = byproduct * inner_;
  ExprPtr hoisted = hoist_real_scalar(inner_);
  reset_hash_value();
  return hoisted;
}

ExprIterator RealPart::begin_subexpr() {
  // N.B. a mutable iterator into inner_ invalidates the memoized hash
  reset_hash_value();
  return ExprIterator{&inner_};
}

ExprIterator RealPart::end_subexpr() {
  reset_hash_value();
  return ExprIterator{&inner_ + 1};
}

ConstExprIterator RealPart::begin_subexpr() const {
  return ConstExprIterator{&inner_};
}

ConstExprIterator RealPart::end_subexpr() const {
  return ConstExprIterator{&inner_ + 1};
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

ExprPtr ImagPart::canonicalize(CanonicalizeOptions opts) {
  // a Constant byproduct of the inner canonicalization multiplies the wrapped
  // expression back, and the eager hoist then takes a real scalar out again:
  // `Im(c X)` is `c Im(X)` for a real c, so that scalar leaves as this node's
  // own byproduct and no canonical wrapper holds a scaled Product. A complex
  // scalar stays inside, the wrapper not being linear over it
  if (auto byproduct = inner_->canonicalize(opts);
      byproduct && byproduct->is<Constant>())
    inner_ = byproduct * inner_;
  ExprPtr hoisted = hoist_real_scalar(inner_);
  reset_hash_value();
  return hoisted;
}

ExprPtr ImagPart::rapid_canonicalize(CanonicalizeOptions opts) {
  SEQUANT_ASSERT(opts.method == CanonicalizationMethod::Rapid);
  if (auto byproduct = inner_->rapid_canonicalize(opts);
      byproduct && byproduct->is<Constant>())
    inner_ = byproduct * inner_;
  ExprPtr hoisted = hoist_real_scalar(inner_);
  reset_hash_value();
  return hoisted;
}

ExprIterator ImagPart::begin_subexpr() {
  // N.B. a mutable iterator into inner_ invalidates the memoized hash
  reset_hash_value();
  return ExprIterator{&inner_};
}

ExprIterator ImagPart::end_subexpr() {
  reset_hash_value();
  return ExprIterator{&inner_ + 1};
}

ConstExprIterator ImagPart::begin_subexpr() const {
  return ConstExprIterator{&inner_};
}

ConstExprIterator ImagPart::end_subexpr() const {
  return ConstExprIterator{&inner_ + 1};
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
  // a lone factor stands on its own: wrapping it in a one-factor Product
  // would make `Re[c X]` and `c Re[X]` different expressions spelled alike
  if (prod.factors().size() == 1) return prod.factors().front();
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
