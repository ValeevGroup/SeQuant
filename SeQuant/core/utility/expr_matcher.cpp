#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/tensor.hpp>

#include <compare>

namespace sequant {

ExprMatcher::ExprMatcher(const Expr &expr, ExprMatcherOptions opts)
    : expr_(expr.clone()), opts_(std::move(opts)) {
  if (!expr.empty()) {
    throw Exception("ExprMatcher does not (yet) support composite expressions");
  }
}

bool ExprMatcher::is_equal(const ExprPtr &expr) const {
  if (!expr) {
    return false;
  }

  return is_equal(*expr);
}

bool ExprMatcher::is_equal(const Expr &other) const {
  if (other.type_id() != expr_->type_id()) {
    return false;
  }

  if (opts_.tensor_cmp == TensorComparison::Block && other.is<Tensor>()) {
    TensorBlockEqualComparator cmp;

    return cmp(expr_->as<Tensor>(), other.as<Tensor>());
  } else {
    return *expr_ == other;
  }
}

std::partial_ordering ExprMatcher::compare(const ExprPtr &other) const {
  if (!other) {
    return std::partial_ordering::unordered;
  }

  return compare(*other);
}

std::partial_ordering ExprMatcher::compare(const Expr &other) const {
  if (!opts_.cross_comparisons && other.type_id() != expr_->type_id()) {
    return std::partial_ordering::unordered;
  }

  bool is_less = false;
  if (opts_.tensor_cmp == TensorComparison::Block && expr_->is<Tensor>() &&
      other.is<Tensor>()) {
    TensorBlockLessThanComparator cmp;

    is_less = cmp(expr_->as<Tensor>(), other.as<Tensor>());
  } else {
    is_less = *expr_ < other;
  }

  if (is_less) {
    return std::partial_ordering::less;
  } else if (is_equal(other)) {
    return std::partial_ordering::equivalent;
  } else {
    return std::partial_ordering::greater;
  }
}

bool ExprMatcher::operator==(const ExprMatcher &other) const {
  return opts_ == other.opts_ && expr_ == other;
}

std::partial_ordering ExprMatcher::operator<=>(const ExprMatcher &other) const {
  if (opts_ != other.opts_) {
    return std::partial_ordering::unordered;
  }

  return expr_ <=> other;
}

bool operator==(const ExprMatcher &matcher, const Expr &expr) {
  return matcher.is_equal(expr);
}
bool operator==(const ExprMatcher &matcher, const ExprPtr &expr) {
  return matcher.is_equal(expr);
}
bool operator==(const Expr &expr, const ExprMatcher &matcher) {
  return matcher.is_equal(expr);
}
bool operator==(const ExprPtr &expr, const ExprMatcher &matcher) {
  return matcher.is_equal(expr);
}

std::partial_ordering operator<=>(const ExprMatcher &matcher,
                                  const Expr &expr) {
  return matcher.compare(expr);
}
std::partial_ordering operator<=>(const ExprMatcher &matcher,
                                  const ExprPtr &expr) {
  return matcher.compare(expr);
}
std::partial_ordering operator<=>(const Expr &expr,
                                  const ExprMatcher &matcher) {
  std::partial_ordering result = matcher.compare(expr);

  // Take care of the fact that the matcher is actually on the RHS
  if (result == std::partial_ordering::less) {
    result = std::partial_ordering::greater;
  } else if (result == std::partial_ordering::greater) {
    result = std::partial_ordering::less;
  }

  return result;
}
std::partial_ordering operator<=>(const ExprPtr &expr,
                                  const ExprMatcher &matcher) {
  return *expr <=> matcher;
}

}  // namespace sequant
