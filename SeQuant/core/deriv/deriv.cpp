#include <SeQuant/core/deriv/deriv.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cstddef>
#include <tuple>
#include <utility>

namespace sequant::deriv {

namespace {

/// @return the derivative of the factor @p f by the scalar variable @p y, or
/// nullptr if @p f is a constant with respect to @p y
ExprPtr scalar_rule(const ExprPtr& f, const ExprPtr& y) {
  if (f->is<Variable>()) return *f == *y ? ex<Constant>(1) : nullptr;
  if (!f->is<Power>()) return nullptr;
  const auto& power = f->as<Power>();
  auto base = power.base()->clone();
  if (power.conjugated()) base->adjoint();
  if (*base != *y) return nullptr;
  ExprPtr lowered = ex<Power>(power.base(), power.exponent() - 1);
  if (power.conjugated()) lowered->as<Power>().conjugate();
  Power::flatten(lowered);
  return ex<Product>(power.exponent(), ExprPtrList{lowered});
}

/// @return the scalar and the factors of @p term , a term of an expanded
/// expression
std::pair<Constant::scalar_type, container::svector<ExprPtr>> split(
    const ExprPtr& term) {
  Constant::scalar_type scalar = 1;
  container::svector<ExprPtr> factors;
  const auto walk = [&](const auto& self, const ExprPtr& e) -> void {
    if (e->is<Product>()) {
      scalar *= e->as<Product>().scalar();
      for (const auto& f : e->as<Product>().factors()) self(self, f);
    } else if (e->is<Constant>()) {
      scalar *= e->as<Constant>().value();
    } else {
      SEQUANT_ASSERT(!e->is<Sum>());
      factors.push_back(e);
    }
  };
  walk(walk, term);
  return {scalar, std::move(factors)};
}

/// @return the terms of the expanded @p E
container::svector<ExprPtr> terms_of(const ExprPtr& E) {
  if (!E->is<Sum>()) return {E};
  const auto& summands = E->as<Sum>().summands();
  return {summands.begin(), summands.end()};
}

/// a term of an expanded expression: its scalar and its factors
using Term = std::tuple<Constant::scalar_type, container::svector<ExprPtr>>;

/// the terms of the expanded @p E
container::svector<Term> split_terms(const ExprPtr& E) {
  container::svector<Term> result;
  for (const auto& term : terms_of(E)) {
    auto [scalar, factors] = split(term);
    result.emplace_back(scalar, std::move(factors));
  }
  return result;
}

/// D_y(E): the product rule over the factors of each of the @p terms of E,
/// with the scalar rule for a scalar @p y
ExprPtr one_way(const container::svector<Term>& terms, const ExprPtr& y) {
  container::vector<ExprPtr> result;
  for (const auto& [scalar, factors] : terms) {
    const auto product = [&](std::size_t skip) {
      auto p = ex<Product>(scalar, ExprPtrList{});
      for (std::size_t j = 0; j != factors.size(); ++j)
        if (j != skip) p->as<Product>().append(1, factors[j]);
      return p;
    };
    for (std::size_t k = 0; k != factors.size(); ++k) {
      const auto& f = factors[k];
      auto d = scalar_rule(f, y);
      if (!d) continue;
      auto p = product(k);
      p->as<Product>().append(1, d);
      result.push_back(p);
    }
  }
  return result.empty() ? ex<Constant>(0)
                        : ex<Sum>(result.begin(), result.end());
}

/// the derivative of @p E by the one variable @p x
ExprPtr by_one(const ExprPtr& E, const ExprPtr& x) {
  auto e = E->clone();
  expand(e);
  auto terms = split_terms(e);
  auto result = one_way(terms, x);
  rapid_simplify(result);
  return result;
}

}  // namespace

ExprPtr diff(const ExprPtr& E, const container::svector<ExprPtr>& variables) {
  if (variables.empty()) return E->clone();
  // by_one() works on a clone of its input
  auto result = E;
  for (const auto& x : variables) result = by_one(result, x);
  return result;
}

}  // namespace sequant::deriv
