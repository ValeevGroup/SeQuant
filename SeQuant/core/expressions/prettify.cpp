#include <SeQuant/core/expressions/prettify.hpp>

#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/lexicographic_order.hpp>

#include <memory>

namespace sequant {

ExprPtr prettify(const ExprPtr& expr, const NamedIndexSpec& named_indices,
                 const PrettifyOptions& opts) {
  if (!expr) return expr;
  ExprPtr result = expr->clone();

  if (opts.recanonicalize) {
    canonicalize(result, named_indices,
                 CanonicalizeOptions{
                     .method = CanonicalizationMethod::Complete,
                     .presentation =
                         std::make_shared<const CanonicalizationPresentation>(
                             opts.presentation)});
  }

  // the canonical marks need no reset: canonicalization with presentation
  // settings does not mark, and a mark copied by the clone is not valid
  // anymore once the expression has changed
  sort_expr(*result, LexicographicExprCompare(LexicographicCompareOptions{
                         .cardinal_tensor_labels =
                             opts.presentation.cardinal_tensor_labels.value_or(
                                 default_cardinal_tensor_labels())}));
  return result;
}

ResultExpr prettify(const ResultExpr& expr, const PrettifyOptions& opts) {
  // shares the expression with expr until it is replaced by the prettified
  // copy
  ResultExpr result = expr;
  result = prettify(expr.expression(), named_indices(expr), opts);
  return result;
}

}  // namespace sequant
