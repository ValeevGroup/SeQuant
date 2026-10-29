#ifndef SEQUANT_CORE_UTILITY_POLYNOMIAL_DIVISION_HPP
#define SEQUANT_CORE_UTILITY_POLYNOMIAL_DIVISION_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

namespace sequant {

/// Determines how terms of the dividend that are not (cleanly) divisible by
/// the divisor are handled by divide()
enum class DivisionMode {
  /// The dividend must be exactly divisible by the divisor, i.e. the
  /// remainder must vanish. Otherwise, an Exception is thrown.
  Strict,
  /// Only those quotient terms q are accepted for which every term of
  /// divisor * q is already present in the dividend. Hence, the remainder
  /// only consists of monomials that are part of the dividend (though their
  /// coefficients may be different) and no cancellation terms are introduced.
  Split,
  /// Standard multivariate polynomial division: every term divisible by the
  /// pivot term of the divisor contributes to the quotient. The remainder
  /// contains all non-divisible terms of the dividend as well as any
  /// cancellation terms that are required to keep the overall expression
  /// unchanged.
  Remainder,
};

/// Determines how the external indices of a tensorial divisor are matched
/// against the terms of the dividend
enum class IndexMatching {
  /// The divisor's external indices are taken literally: only dummy indices
  /// of the dividend may be renamed to match the divisor. The result is of the
  /// form B = A * Q + R.
  Exact,
  /// The divisor acts as a template whose external indices may bind to any
  /// (compatible) index of the dividend. This yields one group per distinct
  /// binding σ of the divisor's external indices to the dividend's external
  /// indices, i.e. B = Σ_σ A_σ * Q_σ + R.
  Relabel,
};

/// One group of the result of a polynomial division
struct DivisionGroup {
  /// The divisor with its external indices relabeled according to relabeling
  ExprPtr divisor;
  /// The quotient associated with divisor
  ExprPtr quotient;
  /// Maps external indices of the original divisor to the indices they are
  /// bound to in this group. Indices that are mapped onto themselves are not
  /// listed (thus, this is always empty for IndexMatching::Exact).
  container::map<Index, Index> relabeling;
};

/// The result of a polynomial division of B by A, representing
/// B = Σ_σ A_σ * Q_σ + R
struct DivisionResult {
  /// The individual (divisor, quotient) groups. Empty, if no term of the
  /// dividend was divisible.
  container::svector<DivisionGroup> groups;
  /// The remainder R (Constant 0 if the division is exact)
  ExprPtr remainder;

  /// @returns Whether the remainder vanishes
  [[nodiscard]] bool exact() const;

  /// @returns The expression Σ_σ A_σ * (Q_σ) + R in which the divisors are
  /// kept as (parenthesized) factors. Dummy indices of the divisors are renamed
  /// as necessary to not clash with the indices of the respective quotient.
  [[nodiscard]] ExprPtr to_expr() const;
};

/// Performs a polynomial division of @p dividend by @p divisor, i.e. it
/// determines quotient(s) and remainder such that
/// dividend = Σ_σ divisor_σ * quotient_σ + remainder
///
/// Both expressions are brought into a simplified (expanded and canonical)
/// form before performing the division. Terms of the dividend may be
/// arbitrary products (including non-commuting operators, which simply become
/// part of the quotient), whereas the terms of the divisor may only consist
/// of tensors and scalars (constants, variables and powers thereof). All terms
/// of the divisor must have the same external indices.
///
/// Since the result of multivariate polynomial division depends on the choice
/// of the pivot ("leading") term of the divisor, every term of the divisor is
/// tried as the pivot and the result with the fewest remainder terms is
/// chosen.
///
/// @param dividend The expression to divide
/// @param divisor The expression to divide by
/// @param mode The DivisionMode to use
/// @param matching The IndexMatching to use
/// @returns The result of the division
/// @throws Exception if mode is DivisionMode::Strict and @p dividend is not
/// divisible by @p divisor
[[nodiscard]] DivisionResult divide(const ExprPtr &dividend,
                                    const ExprPtr &divisor, DivisionMode mode,
                                    IndexMatching matching);

/// Determines in which form factor_out() returns the factored expression
enum class FactorStyle {
  /// Keep the factor as a (parenthesized) factor, i.e. A (Q) + R
  Keep,
  /// Drop the factor and the remainder, i.e. only return Q. This requires that
  /// at most one group exists in the division result.
  Drop,
};

/// Factors @p factor out of @p expr
///
/// @param expr The expression to factor
/// @param factor The factor to factor out
/// @param mode The DivisionMode to use
/// @param matching The IndexMatching to use
/// @param style The FactorStyle in which to produce the result
/// @returns The factored expression
/// @throws Exception if @p style is FactorStyle::Drop and the division
/// produced more than one group
/// @sa divide
[[nodiscard]] ExprPtr factor_out(const ExprPtr &expr, const ExprPtr &factor,
                                 DivisionMode mode, IndexMatching matching,
                                 FactorStyle style);

/// Factors @p factor out of @p expr and replaces every occurrence of
/// (a relabeled version of) @p factor with (a correspondingly relabeled
/// version of) @p substitute, i.e. it produces Σ_σ S_σ Q_σ + R
///
/// Indices of @p substitute that coincide with external indices of @p factor
/// are relabeled in the same way as the external indices of @p factor in
/// every group. All other indices of @p substitute are left as-is.
///
/// @param expr The expression to factor
/// @param factor The factor to factor out
/// @param substitute The expression to replace (relabeled) factors with
/// @param mode The DivisionMode to use
/// @param matching The IndexMatching to use
/// @returns The factored expression
/// @sa divide
[[nodiscard]] ExprPtr factor_out(const ExprPtr &expr, const ExprPtr &factor,
                                 const Expr &substitute, DivisionMode mode,
                                 IndexMatching matching);

/// Factors the expression of @p factor out of @p expr and replaces it with the
/// result (tensor or variable) of @p factor, e.g. factoring I{a1;i1} = f{a1;a2}
/// t{a2;i1} + g{a1;i1} out of an expression produces terms of the form I{a1;i1}
/// Q + R
///
/// @param expr The expression to factor
/// @param factor The factor (and its result) to factor out
/// @param mode The DivisionMode to use
/// @param matching The IndexMatching to use
/// @returns The factored expression
/// @sa divide
[[nodiscard]] ExprPtr factor_out(const ExprPtr &expr, const ResultExpr &factor,
                                 DivisionMode mode, IndexMatching matching);

/// Like factor_out(const ExprPtr&, const ExprPtr&, DivisionMode,
/// IndexMatching, FactorStyle) with FactorStyle::Keep but for a ResultExpr,
/// whose result is retained
[[nodiscard]] ResultExpr factor_out(const ResultExpr &expr,
                                    const ExprPtr &factor, DivisionMode mode,
                                    IndexMatching matching);

/// Like factor_out(const ExprPtr&, const ResultExpr&, DivisionMode,
/// IndexMatching) but for a ResultExpr, whose result is retained
[[nodiscard]] ResultExpr factor_out(const ResultExpr &expr,
                                    const ResultExpr &factor, DivisionMode mode,
                                    IndexMatching matching);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_POLYNOMIAL_DIVISION_HPP
