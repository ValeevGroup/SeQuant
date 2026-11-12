#ifndef SEQUANT_CORE_UTILITY_EXPRFILTER_HPP
#define SEQUANT_CORE_UTILITY_EXPRFILTER_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>

#include <cstddef>
#include <limits>
#include <memory>
#include <regex>
#include <vector>

namespace sequant {

/// A predicate on expressions, composed of rules that are combined either
/// conjunctively (default) or disjunctively
class ExpressionFilter {
 public:
  struct Rule {
    bool negate = false;

    virtual ~Rule() = default;

    virtual bool matches(const Expr &expr) const = 0;
  };

  ExpressionFilter() = default;

  /// @returns Whether the given expression satisfies this filter. Without any
  /// rules, this is the case if and only if all rules are required.
  bool matches(const Expr &) const;

  void set_require_all(bool require);
  void add_rule(std::unique_ptr<Rule> rule);

 private:
  bool require_all_ = true;
  std::vector<std::unique_ptr<Rule>> rules_;
};

/// Rule matching expressions that contain (at any depth) a subexpression
/// matching the given matcher
struct ContainsRule : ExpressionFilter::Rule {
  explicit ContainsRule(ExprMatcher matcher);

  bool matches(const Expr &expr) const override;

  ExprMatcher matcher_;
};

/// Rule matching expressions that contain (at any depth) a labeled
/// subexpression whose label fully matches one of the given patterns
struct LabelRule : ExpressionFilter::Rule {
  explicit LabelRule(std::vector<std::wregex> patterns);

  bool matches(const Expr &expr) const override;

  std::vector<std::wregex> patterns_;
};

/// Bounds on the number of a tensor's indices that live in a given space
struct IndexSpaceCount {
  IndexSpace space;
  std::size_t min = 0;
  std::size_t max = std::numeric_limits<std::size_t>::max();
};

/// Rule matching expressions that contain (at any depth) a tensor whose
/// number of (non-null) indices in each of the given spaces lies within the
/// respective bounds. Indices in other spaces are not constrained. If
/// patterns are given, the tensor's label additionally has to fully match one
/// of them.
struct IndexSpaceCountRule : ExpressionFilter::Rule {
  explicit IndexSpaceCountRule(std::vector<IndexSpaceCount> counts,
                               std::vector<std::wregex> label_patterns = {});

  bool matches(const Expr &expr) const override;

  std::vector<IndexSpaceCount> counts_;
  std::vector<std::wregex> label_patterns_;
};

/// Rule matching expressions that contain (at any depth) a tensor with exactly
/// the given total number of (non-null) indices, i.e. modes
struct IndexCountRule : ExpressionFilter::Rule {
  explicit IndexCountRule(std::size_t count);

  bool matches(const Expr &expr) const override;

  std::size_t count_;
};

/// Recursively rebuilds @p expr, keeping only the parts relevant to
/// @p filter: for a Sum, keeps (and further prunes) only those summands
/// that (recursively) contain a match, dropping the rest; for a Product,
/// keeps every factor (dropping one would change what the product
/// computes), pruning each factor in turn; leaves are returned unchanged.
ExprPtr prune_to_matches(const Expr &expr, const ExpressionFilter &filter);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_EXPRFILTER_HPP
