#ifndef SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP
#define SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

#include <compare>

namespace sequant {

enum class TensorComparison {
  /// Tensors are compared including their specific index names
  Identity,
  /// Tensor indices are only compared by their index space
  Shape,
  /// @sa Shape
  Block = Shape,
};

/// A wrapper around an expression that allows customization on how
/// exactly two expressions are supposed to be compared.
/// It is directly comparable to expression objects and can therefore
/// be used for expression lookup in containers (assuming the container
/// allows for transparent comparators).
struct ExprMatcherOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;

  /// The way tensor-valued expressions shall be compared
  TensorComparison tensor_cmp = TensorComparison::Identity;
  /// Whether to compare expressions across different expression  types
  bool cross_comparisons = true;

  bool operator==(const ExprMatcherOptions &) const = default;
};

class ExprMatcher {
 public:
  explicit ExprMatcher(const Expr &expr, ExprMatcherOptions opts = {});

  bool is_equal(const ExprPtr &other) const;
  bool is_equal(const Expr &other) const;

  std::partial_ordering compare(const ExprPtr &other) const;
  std::partial_ordering compare(const Expr &other) const;

  bool operator==(const ExprMatcher &other) const;
  std::partial_ordering operator<=>(const ExprMatcher &other) const;

 private:
  ExprPtr expr_;
  ExprMatcherOptions opts_;
};

bool operator==(const ExprMatcher &matcher, const Expr &expr);
bool operator==(const ExprMatcher &matcher, const ExprPtr &expr);
bool operator==(const Expr &expr, const ExprMatcher &matcher);
bool operator==(const ExprPtr &exp, const ExprMatcher &matcherr);

std::partial_ordering operator<=>(const ExprMatcher &matcher, const Expr &expr);
std::partial_ordering operator<=>(const ExprMatcher &matcher,
                                  const ExprPtr &expr);
std::partial_ordering operator<=>(const Expr &expr, const ExprMatcher &matcher);
std::partial_ordering operator<=>(const ExprPtr &expr,
                                  const ExprMatcher &matcher);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP
