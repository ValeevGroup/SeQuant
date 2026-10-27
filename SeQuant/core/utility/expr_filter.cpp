#include <SeQuant/core/utility/expr_filter.hpp>

#include <SeQuant/core/expr.hpp>

#include <algorithm>
#include <memory>
#include <regex>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {

bool ExpressionFilter::matches(const Expr &expr) const {
  for (const std::unique_ptr<Rule> &current : rules_) {
    if (current->matches(expr) == !current->negate) {
      if (!require_all_) {
        return true;
      }
    } else if (require_all_) {
      return false;
    }
  }

  // In "any" mode, getting here means that no rule matched
  return require_all_;
}

void ExpressionFilter::set_require_all(bool require) { require_all_ = require; }

void ExpressionFilter::add_rule(std::unique_ptr<Rule> rule) {
  rules_.emplace_back(std::move(rule));
}

namespace {

/// Checks @p expr and, recursively, every descendant subexpression against
/// @p predicate, stopping at the first match.
template <typename Predicate>
bool matches_recursively(const Expr &expr, Predicate &&predicate) {
  if (predicate(expr)) {
    return true;
  }

  for (const ExprPtr &sub : expr) {
    if (matches_recursively(*sub, predicate)) {
      return true;
    }
  }

  return false;
}

}  // namespace

ContainsRule::ContainsRule(ExprMatcher matcher)
    : matcher_(std::move(matcher)) {}

bool ContainsRule::matches(const Expr &expr) const {
  return matches_recursively(expr,
                             [&](const Expr &sub) { return sub == matcher_; });
}

LabelRule::LabelRule(std::vector<std::wregex> patterns)
    : patterns_(std::move(patterns)) {}

bool LabelRule::matches(const Expr &expr) const {
  return matches_recursively(expr, [&](const Expr &sub) {
    if (!sub.is<Labeled>()) {
      return false;
    }

    std::wstring_view label = sub.as<Labeled>().label();

    return std::ranges::any_of(patterns_, [&](const std::wregex &pattern) {
      return std::regex_match(label.begin(), label.end(), pattern);
    });
  });
}

/// Recursively rebuilds @p expr, keeping only the parts relevant to
/// @p filter: for a Sum, keeps (and further prunes) only those summands
/// that (recursively) contain a match, dropping the rest; for a Product,
/// keeps every factor (dropping one would change what the product
/// computes), pruning each factor in turn; leaves are returned unchanged.
ExprPtr prune_to_matches(const Expr &expr, const ExpressionFilter &filter) {
  if (expr.is<Sum>()) {
    Sum::summands_type kept;

    for (const ExprPtr &summand : expr.as<Sum>().summands()) {
      if (filter.matches(*summand)) {
        kept.push_back(prune_to_matches(*summand, filter));
      }
    }

    return ex<Sum>(std::move(kept));
  }

  if (expr.is<Product>()) {
    const Product &prod = expr.as<Product>();
    Product::factors_type pruned_factors;
    pruned_factors.reserve(prod.factors().size());

    for (const ExprPtr &factor : prod.factors()) {
      pruned_factors.push_back(prune_to_matches(*factor, filter));
    }

    return ex<Product>(prod.scalar(), std::move(pruned_factors),
                       Product::Flatten::No);
  }

  return expr.clone();
}

}  // namespace sequant
