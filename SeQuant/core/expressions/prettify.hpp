#ifndef SEQUANT_EXPRESSIONS_PRETTIFY_HPP
#define SEQUANT_EXPRESSIONS_PRETTIFY_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/options.hpp>
#include <SeQuant/core/utility/macros.hpp>

namespace sequant {

/// options of prettify()
struct PrettifyOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;

  /// whether to canonicalize (Complete) with the #presentation settings before
  /// sorting; if false, the summands and the commuting factors of the
  /// expression are only sorted
  bool recanonicalize = true;
  /// the presentation settings: the canonicalizers and comparers used to
  /// recanonicalize, and the cardinal tensor labels used to recanonicalize
  /// and to sort (unset members default to those of the library's canonical
  /// form)
  CanonicalizationPresentation presentation = {};
};

/// @brief Brings an expression into a form meant for reading
///
/// The library canonicalizes to one fixed internal form, which is meant to
/// identify equivalent expressions, not to be read. This returns an equivalent
/// expression in a form chosen by @p opts : (optionally) canonicalized with
/// the presentation settings, then ordered lexicographically (see sort_expr()
/// and LexicographicExprCompare), which, unlike canonicalization, does not
/// depend on the hashing algorithm. The named indices are never renamed.
/// @note the result is meant for presentation only, i.e. it is not
/// canonical in the sense of Expr::is_canonical(); pass it to
/// canonicalize()/simplify() before relying on its form
/// @param expr the expression to prettify; it is not modified
/// @param named_indices the indices whose labels are meaningful (ordinarily
/// the external indices of @p expr ); see NamedIndexSpec
/// @param opts the options
/// @return the prettified copy of @p expr
[[nodiscard]] ExprPtr prettify(const ExprPtr& expr,
                               const NamedIndexSpec& named_indices,
                               const PrettifyOptions& opts = {});

/// @brief Brings a ResultExpr into a form meant for reading, with the indices
/// of its result as the named indices
/// @sa prettify(const ExprPtr&, const NamedIndexSpec&, const PrettifyOptions&)
/// @param expr the expression to prettify; it is not modified
/// @param opts the options
/// @return the prettified copy of @p expr
[[nodiscard]] ResultExpr prettify(const ResultExpr& expr,
                                  const PrettifyOptions& opts = {});

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_PRETTIFY_HPP
