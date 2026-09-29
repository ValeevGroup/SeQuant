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
/// _replaced_ (i.e. `&expr` may be mutated by call); a Sum left with at most
/// one summand, or a Product with one factor and a unit scalar, is replaced
/// by that summand or factor (or 0), canonicalized on its own
/// @note the canonicalization options are those of the default context
/// @return \p expr to facilitate chaining
ExprPtr& canonicalize(ExprPtr& expr);

/// Recursively canonicalizes an Expr; like mutating canonicalize() but works
/// for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be canonicalized
/// @note the canonicalization options are those of the default context
/// @return canonicalized form of \p expr_rv
ExprPtr canonicalize(ExprPtr&& expr_rv);

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @note the canonicalization options are those of the default context
/// @return \p expr to facilitate chaining
ResultExpr& canonicalize(ResultExpr& expr);

/// Recursively canonicalizes an Expr; like mutating canonicalize() but works
/// for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be canonicalized
/// @note the canonicalization options are those of the default context
/// @return canonicalized form of \p expr_rv
[[nodiscard]] ResultExpr& canonicalize(ResultExpr&& expr);

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
/// @sa simplify()
/// @return \p expr to facilitate chaining
ExprPtr& rapid_simplify(ExprPtr& expr);

/// Simplifies an Expr by applying cheap transformations (e.g. eliminating
/// trivial math, flattening sums and products, etc.)
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @sa simplify()
/// @return \p expr to facilitate chaining
ResultExpr& rapid_simplify(ResultExpr& expr);

/// Simplifies an Expr by applying cheap transformations (e.g. eliminating
/// trivial math, flattening sums and products, etc.)
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @sa simplify()
/// @return \p expr to facilitate chaining
[[nodiscard]] ResultExpr& rapid_simplify(ResultExpr&& expr);

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @note the canonicalization options are those of the default context
/// @sa rapid_simplify()
/// @return \p expr to facilitate chaining
ExprPtr& simplify(ExprPtr& expr);

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify; like mutating simplify() but works for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be simplified
/// @note the canonicalization options are those of the default context
/// @return simplified form of \p expr_rv
ExprPtr simplify(ExprPtr&& expr_rv);

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @note the canonicalization options are those of the default context
/// @sa rapid_simplify()
/// @return \p expr to facilitate chaining
ResultExpr& simplify(ResultExpr& expr);

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify; like mutating simplify() but works for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be simplified
/// @note the canonicalization options are those of the default context
/// @return simplified form of \p expr_rv
[[nodiscard]] ResultExpr& simplify(ResultExpr&& expr);

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

/// @brief the complex conjugate of the value of a c-number expression
///
/// `conj <p|O|q> = <q|O⁺|p>`, so on c-number content this is
/// sequant::adjoint(const ExprPtr&), one implementation: every `Tensor` is
/// adjointed (over a real basis the coset rule spells that adjoint as the
/// `꙳` state with the slots in place), `Constant`, `Variable` and `Power`
/// conjugate, `RealPart`/`ImagPart` are invariant, a `Sum` maps its
/// summands, and a `Product` conjugates its scalar and reverses its factors
/// -- a reversal that is not observable, the factors of a c-number product
/// commuting. For an expression with open indices the head's bra and ket are
/// exchanged, so a `ResultExpr` caller exchanges the head as well.
///
/// The result is assembled with the Sum and Product constructors' default
/// flattening, as kconjugate()'s is: a nested product is spliced into the
/// rebuilt product and a factor's sign byproduct folds into its scalar, and
/// a nested sum is spliced into the rebuilt sum.
///
/// Conjugation is an involution up to that flattening:
/// `conjugate(conjugate(e))` equals `e` for an @p expr whose Products and
/// Sums are already flat, and equals its flattened form otherwise.
///
/// @param expr a c-number expression
/// @return a new expression denoting `conj(expr)`
/// @throw sequant::Exception for operator-valued content, which has no value
///        to conjugate; sequant::kconjugate() is the conjugation of an
///        operator
/// @sa kconjugate(const ExprPtr&)
[[nodiscard]] ExprPtr conjugate(const ExprPtr& expr);

/// @brief the K-conjugate of an expression, `K E K⁻¹`
///
/// Dispatches to Expr::kconjugate() on a clone: the complex conjugate of a
/// scalar, the K-conjugated state of a Tensor, the identity on an operator
/// string. A sign byproduct that the conjugated object cannot hold (a Tensor
/// holds no scalar) becomes a scalar factor here.
///
/// The result is assembled with the Sum and Product constructors' default
/// flattening, as conjugate()'s is:
/// - a Product is rebuilt with that flattening, so a nested product is
///   spliced into it and a factor's sign byproduct folds into the product's
///   scalar;
/// - a Sum is rebuilt with that flattening, so a nested sum is spliced in,
///   but a summand whose K-conjugate carries a sign stays a
///   `Product{-1, summand}` -- a Sum has no scalar to fold into, so the
///   sign is carried by that summand's own scalar.
///
/// @param expr an expression
/// @return a new expression denoting `K expr K⁻¹`
/// @throw sequant::Exception if @p expr is operator-valued and a normal
///        operator in it acts on an index space that `K` does not close
///        (only a real-field space, with a real-field proto-index closure,
///        is `K`-closed today). The check covers `NormalOperator` and
///        `NormalOperatorSequence`; two operator kinds are accepted without
///        it: a non-normal-ordered `Operator<S>` (`FOperator`/`BOperator`),
///        which the check does not reach, and an `mbpt::Operator`, whose
///        K-conjugate is the operator itself -- complex conjugation acts
///        through the coefficients of its tensor form.
/// @sa conjugate(const ExprPtr&)
[[nodiscard]] ExprPtr kconjugate(const ExprPtr& expr);

/// Folds complex-conjugate-related summand pairs of a sum, exactly.
///
/// A summand pair {s, s*} contributes s + s* = 2 Re(s), and a pair
/// {s, -s*} contributes s - s* = 2i Im(s): both identities are
/// unconditional (no reality assumption on the sum), the imaginary/real
/// parts being carried symbolically by RealPart/ImagPart nodes. Pair
/// detection goes through canonical forms (robust to dummy renaming and
/// factor reordering), matched via a hash-bucketed lookup. Summands whose
/// conjugate is not present -- including self-conjugate (manifestly real)
/// summands -- are left untouched, as are operator-valued summands (the
/// fold applies to c-number content only).
///
/// A folded pair's Re/Im wrapper carries the pair's canonical
/// representative, which may be the other member of the pair and has its
/// dummy indices relabeled; that is what makes the fold independent of the
/// order the pair was written in, up to a canonicalization of the result.
///
/// @param[in] expr the sum to fold; returned unchanged if not a Sum
/// @param[in] conjugate_op optional map from a summand to an expression the
///            caller asserts to _equal_ the summand's complex conjugate in
///            value. Defaults to the algebraic adjoint (for a fully
///            contracted c-number summand the adjoint IS its complex
///            conjugate). Supply a custom map when a domain identity
///            relates the conjugate to a different symbolic form than the
///            adjoint (e.g. a symmetry of the leaf tensors expressed as an
///            index relabeling), so conjugate pairs written in that form
///            can be recognized.
/// @return the folded expression
/// @note the pairs are identified by canonical form under the context's
///       canonicalization options, with named index labels always treated as
///       meaningful, as in Sum::canonicalize_impl
[[nodiscard]] ExprPtr fold_conjugate_pairs(
    ExprPtr const& expr,
    std::function<ExprPtr(ExprPtr const&)> conjugate_op = {});

/// Back-compat variant of fold_conjugate_pairs() for a sum whose _value_ the
/// caller asserts to be real: a pair {s, s*} folds to 2*s (the imaginary
/// parts of the folded and input expressions differ; both are discarded by
/// the caller's reality assertion). Difference pairs {s, -s*} are _not_
/// folded. Prefer fold_conjugate_pairs(), which needs no assertion.
[[deprecated(
    "fragile: asserts (unverifiably) that the sum's value is real, and "
    "silently leaves {s, -s*} difference pairs unfolded; use "
    "fold_conjugate_pairs()")]] [[nodiscard]] ExprPtr
fold_conjugate_pairs_of_real_sum(
    ExprPtr const& expr,
    std::function<ExprPtr(ExprPtr const&)> conjugate_op = {});

/// @return whether the c-number expression @p expr denotes a Hermitian
///         network: its _value_ equals that of its adjoint (conjugate
///         transpose), decided by comparing canonical forms. For a closed
///         (fully contracted, scalar-valued) network this is reality
///         recognition: N == conj(N). This derived recognition is what the
///         time-reversal-symmetry folding builds on; subexpressions carry no
///         first-class hermiticity tag today (a cached tag on subnetworks is
///         a possible later extension). For an open network the comparison
///         answers the strict expression-level question (adjoint exchanges
///         the named bra/ket slots), not block hermiticity under a slot
///         pairing -- that refinement also belongs to the time-reversal
///         work.
/// @pre `expr->is_cnumber()`
[[nodiscard]] bool is_hermitian_network(ExprPtr const& expr);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_ALGORITHMS_HPP
