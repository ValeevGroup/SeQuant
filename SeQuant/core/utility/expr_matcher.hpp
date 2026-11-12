#ifndef SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP
#define SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP

#include <SeQuant/api.hpp>
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

class SEQUANT_SYMB_INTERNAL_API ExprMatcher {
 public:
  explicit ExprMatcher(const Expr &expr, ExprMatcherOptions opts = {});

  bool is_equal(const ExprPtr &other) const;
  bool is_equal(const Expr &other) const;

  std::partial_ordering compare(const ExprPtr &other) const;
  std::partial_ordering compare(const Expr &other) const;

  const Expr &expr() const;

  bool operator==(const ExprMatcher &other) const;
  std::partial_ordering operator<=>(const ExprMatcher &other) const;

 private:
  ExprPtr expr_;
  ExprMatcherOptions opts_;
};

SEQUANT_SYMB_INTERNAL_API bool operator==(const ExprMatcher &matcher,
                                          const Expr &expr);
SEQUANT_SYMB_INTERNAL_API bool operator==(const ExprMatcher &matcher,
                                          const ExprPtr &expr);
SEQUANT_SYMB_INTERNAL_API bool operator==(const Expr &expr,
                                          const ExprMatcher &matcher);
SEQUANT_SYMB_INTERNAL_API bool operator==(const ExprPtr &exp,
                                          const ExprMatcher &matcherr);

SEQUANT_SYMB_INTERNAL_API std::partial_ordering operator<=>(
    const ExprMatcher &matcher, const Expr &expr);
SEQUANT_SYMB_INTERNAL_API std::partial_ordering operator<=>(
    const ExprMatcher &matcher, const ExprPtr &expr);
SEQUANT_SYMB_INTERNAL_API std::partial_ordering operator<=>(
    const Expr &expr, const ExprMatcher &matcher);
SEQUANT_SYMB_INTERNAL_API std::partial_ordering operator<=>(
    const ExprPtr &expr, const ExprMatcher &matcher);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_EXPRMATCHER_HPP
