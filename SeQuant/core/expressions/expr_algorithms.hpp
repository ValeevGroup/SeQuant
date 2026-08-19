//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_EXPRESSIONS_ALGORITHMS_HPP
#define SEQUANT_EXPRESSIONS_ALGORITHMS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/access.hpp>
#include <range/v3/view/transform.hpp>

#include <functional>
#include <string>

namespace sequant {

/// splits long outer sum into a multiline align
/// @param exprptr the expression to be converted to a string
/// @param max_lines_per_align the maximum number of lines in the align before
/// starting new align block (if zero, will produce single align block)
/// @param max_terms_per_line the maximum number of terms per line
std::wstring to_latex_align(const ExprPtr& exprptr,
                            size_t max_lines_per_align = 0,
                            size_t max_terms_per_line = 1);

template <typename Sequence>
std::decay_t<Sequence> clone(Sequence&& exprseq) {
  auto cloned_seq = exprseq | ranges::views::transform([](const ExprPtr& ptr) {
                      return ptr ? ptr->clone() : nullptr;
                    });
  return std::decay_t<Sequence>(ranges::begin(cloned_seq),
                                ranges::end(cloned_seq));
}

/// @param[in] expr an expression
/// @return number of subexpressions in @p expr, i.e., 0 for atoms (Constant,
/// Variable, Tensor, etc.), >0 for nontrivial Product or Sum
std::size_t size(const Expr& expr);

/// @param[in] exprptr (a pointer to) an expression
/// @return number of subexpressions in @p exprptr , i.e., 0 if @p exprptr is
/// null or an atom (Constant, Variable, Tensor, etc.), >0 for nontrivial
/// Product or Sum
std::size_t size(const ExprPtr& exprptr);

/// @param[in] exprptr (a pointer to) an expression
/// @return begin iterator to the expression range
inline decltype(auto) begin(const ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::begin(*exprptr);
}

/// @param[in] exprptr (a pointer to) an expression
/// @return begin iterator to the expression range
inline decltype(auto) begin(ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::begin(*exprptr);
}

/// @param[in] exprptr (a pointer to) an expression
/// @return begin iterator to the expression range
inline decltype(auto) cbegin(const ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::cbegin(*exprptr);
}

/// @param[in] exprptr (a pointer to) an expression
/// @return end iterator to the expression range
inline decltype(auto) end(const ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::end(*exprptr);
}

/// @param[in] exprptr (a pointer to) an expression
/// @return end iterator to the expression range
inline decltype(auto) end(ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::end(*exprptr);
}

/// @param[in] exprptr (a pointer to) an expression
/// @return end iterator to the expression range
inline decltype(auto) cend(const ExprPtr& exprptr) {
  SEQUANT_ASSERT(exprptr);
  return ranges::cend(*exprptr);
}

template <typename T>
bool ExprPtr::is() const {
  return as_shared_ptr()->is<T>();
}

template <typename T>
const T& ExprPtr::as() const {
  return as_shared_ptr()->as<T>();
}

template <typename T>
T& ExprPtr::as() {
  return as_shared_ptr()->as<T>();
}

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return \p expr to facilitate chaining
ExprPtr& canonicalize(
    ExprPtr& expr,
    CanonicalizeOptions opts = CanonicalizeOptions::default_options());

/// Recursively canonicalizes an Expr; like mutating canonicalize() but works
/// for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be canonicalized
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return canonicalized form of \p expr_rv
ExprPtr canonicalize(
    ExprPtr&& expr_rv,
    CanonicalizeOptions opts = CanonicalizeOptions::default_options());

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return \p expr to facilitate chaining
ResultExpr& canonicalize(
    ResultExpr& expr,
    CanonicalizeOptions opts = CanonicalizeOptions::default_options());

/// Recursively canonicalizes an Expr; like mutating canonicalize() but works
/// for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be canonicalized
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return canonicalized form of \p expr_rv
[[nodiscard]] ResultExpr& canonicalize(
    ResultExpr&& expr,
    CanonicalizeOptions opts = CanonicalizeOptions::default_options());

/// Folds complex-conjugate-related summand pairs of a sum whose VALUE the
/// caller asserts to be real.
///
/// For a real-valued sum, Re(s + s*) = Re(2 s). The complex conjugate of a
/// scalar (fully contracted) summand is its adjoint (conjugated scalar,
/// reversed adjoint factors; a BraKetSymmetry::Conjugate tensor's adjoint is
/// its bra<->ket-swapped orientation), so every summand pair {s, adjoint(s)}
/// collapses to 2*s without changing the sum's (real) value. This is the
/// symbolic-layer exploitation of conjugate braket symmetry: the eval-layer
/// Conjugate fold folds conjugate-related LEAVES onto one cache
/// slot, while this folds conjugate-related TERMS out of the sum entirely.
/// Summands whose adjoint is not present among the other summands --
/// including self-adjoint (manifestly real) summands -- are left untouched.
///
/// @warning The caller asserts the sum's VALUE is real (e.g. an expectation
/// value consumed through its real part); the folded expression's imaginary
/// part differs from the input's (both are discarded by that assertion).
///
/// @param[in] expr the sum to fold; returned unchanged if not a Sum
/// @param[in] opts canonicalization options used to identify pairs (named
///            index labels are always treated as meaningful, as in
///            Sum::canonicalize_impl)
/// @param[in] conjugate_op optional map from a summand to an expression the
///            caller asserts to EQUAL the summand's complex conjugate in
///            value. Defaults to the algebraic adjoint. Supply a custom map
///            when a domain identity relates the conjugate to a different
///            symbolic form than the adjoint (e.g. a symmetry of the leaf
///            tensors expressed as an index relabeling), so conjugate pairs
///            written in that form can be recognized.
/// @return the folded expression
ExprPtr fold_conjugate_pairs_of_real_sum(
    ExprPtr const& expr,
    CanonicalizeOptions opts = CanonicalizeOptions::default_options(),
    std::function<ExprPtr(ExprPtr const&)> conjugate_op = {});

/// Recursively expands products of sums
/// @param[in,out] expr expression to be expanded
/// @return \p expr to facilitate chaining
ExprPtr& expand(ExprPtr& expr);

/// Recursively expands products of sums
/// @param[in,out] expr expression to be expanded
/// @return \p expr to facilitate chaining
ExprPtr expand(ExprPtr&& expr);

/// Recursively expands products of sums
/// @param[in,out] expr expression to be expanded
/// @return \p expr to facilitate chaining
ResultExpr& expand(ResultExpr& expr);

/// Recursively expands products of sums
/// @param[in,out] expr expression to be expanded
/// @return The expanded expression
[[nodiscard]] ResultExpr& expand(ResultExpr&& expr);

/// Recursively flattens Sum of Sum's and Product of Product's
/// @param[in,out] expr expression to be flattened
/// @return \p expr to facilitate chaining
ExprPtr& flatten(ExprPtr& expr);

/// Recursively flattens Sum of Sum's and Product of Product's
/// @param[in,out] expr expression to be flattened
/// @return \p expr to facilitate chaining
ExprPtr flatten(ExprPtr&& expr);

/// Recursively flattens Sum of Sum's and Product of Product's
/// @param[in,out] expr expression to be flattened
/// @return \p expr to facilitate chaining
ResultExpr& flatten(ResultExpr& expr);

/// Recursively flattens Sum of Sum's and Product of Product's
/// @param[in,out] expr expression to be flattened
/// @return The expanded expression
[[nodiscard]] ResultExpr& flatten(ResultExpr&& expr);

/// Simplifies an Expr by applying cheap transformations (e.g. eliminating
/// trivial math, flattening sums and products, etc.)
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @sa simplify()
/// @return \p expr to facilitate chaining
ExprPtr& rapid_simplify(
    ExprPtr& expr, SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by applying cheap transformations (e.g. eliminating
/// trivial math, flattening sums and products, etc.)
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @sa simplify()
/// @return \p expr to facilitate chaining
ResultExpr& rapid_simplify(
    ResultExpr& expr,
    SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by applying cheap transformations (e.g. eliminating
/// trivial math, flattening sums and products, etc.)
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @sa simplify()
/// @return \p expr to facilitate chaining
[[nodiscard]] ResultExpr& rapid_simplify(
    ResultExpr&& expr,
    SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @sa rapid_simplify()
/// @return \p expr to facilitate chaining
ExprPtr& simplify(ExprPtr& expr,
                  SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify; like mutating simplify() but works for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be simplified
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return simplified form of \p expr_rv
ExprPtr simplify(ExprPtr&& expr_rv,
                 SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @sa rapid_simplify()
/// @return \p expr to facilitate chaining
ResultExpr& simplify(ResultExpr& expr,
                     SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify; like mutating simplify() but works for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be simplified
/// @param[in] opts canonicalization options (if not given, uses
///            CanonicalizeOptions::default_options() to obtain the default)
/// @return simplified form of \p expr_rv
[[nodiscard]] ResultExpr& simplify(
    ResultExpr&& expr,
    SimplifyOptions opts = SimplifyOptions::default_options());

/// Simplifies an Expr by a combination of expansion and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @sa simplify()
/// @return \p expr to facilitate chaining
ExprPtr& non_canon_simplify(ExprPtr& expr);

/// Simplifies an Expr by a combination of expansion and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @sa simplify()
/// @return \p expr to facilitate chaining
ResultExpr& non_canon_simplify(ResultExpr& expr);

/// Simplifies an Expr by a combination of expansion and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @sa simplify()
/// @return Simplified expression
[[nodiscard]] ResultExpr non_canon_simplify(ResultExpr&& expr);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_ALGORITHMS_HPP
