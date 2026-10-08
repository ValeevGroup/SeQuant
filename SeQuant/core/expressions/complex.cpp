#include <SeQuant/core/expressions/complex.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_operators.hpp>

namespace sequant {

namespace detail {
ExprPtr hoist_real_scalar(ExprPtr& inner) {
  if (!inner->is<Product>()) return {};
  const auto& prod = inner->as<Product>();
  const auto scalar = prod.scalar();
  if (scalar.imag() != 0 || scalar.real() == 1) return {};
  inner = strip_scalar(prod);
  return ex<Constant>(scalar);
}

ExprPtr canonicalize_projected(ExprPtr& inner, CanonicalizeOptions opts,
                               bool rapid) {
  // a Constant byproduct of the inner canonicalization multiplies the wrapped
  // expression back, and the eager hoist then takes a real scalar out again:
  // `Re(c X)` is `c Re(X)` for a real c, so that scalar leaves as the node's
  // own byproduct and no canonical wrapper holds a scaled Product. A complex
  // scalar stays inside, the wrapper not being linear over it
  if (auto byproduct =
          rapid ? inner->rapid_canonicalize(opts) : inner->canonicalize(opts);
      byproduct && byproduct->is<Constant>())
    inner = byproduct * inner;
  return hoist_real_scalar(inner);
}
}  // namespace detail

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
