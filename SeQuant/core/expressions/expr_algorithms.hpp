//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_EXPRESSIONS_ALGORITHMS_HPP
#define SEQUANT_EXPRESSIONS_ALGORITHMS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/print_options.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/access.hpp>
#include <range/v3/view/transform.hpp>

#include <string>
#include <utility>
#include <vector>

namespace sequant {

/// splits long outer sum into a multiline align
/// @param exprptr the expression to be converted to a string
/// @param max_lines_per_align the maximum number of lines in the align before
/// starting new align block (if zero, will produce single align block)
/// @param max_terms_per_line the maximum number of terms per line
/// @note uses the PrintOptions of the default context
std::wstring to_latex_align(const ExprPtr& exprptr,
                            size_t max_lines_per_align = 0,
                            size_t max_terms_per_line = 1);

/// splits long outer sum into a multiline align
/// @param exprptr the expression to be converted to a string
/// @param opts the PrintOptions
/// @param max_lines_per_align the maximum number of lines in the align before
/// starting new align block (if zero, will produce single align block)
/// @param max_terms_per_line the maximum number of terms per line
std::wstring to_latex_align(const ExprPtr& exprptr, const PrintOptions& opts,
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

/// Recursively canonicalizes an Expr and replaces it as needed
/// @param[in,out] expr expression to be canonicalized; may be
/// _replaced_ (i.e. `&expr` may be mutated by call); a Sum left with at most
/// one summand, or a Product with one factor and a unit scalar, is replaced
/// by that summand or factor (or 0), canonicalized on its own
/// @param[in] named_indices the indices whose labels are meaningful, i.e. that
/// canonicalization must not rename (ordinarily the external indices of
/// @p expr); see NamedIndexSpec
/// @param[in] opts the canonicalization options
/// @return \p expr to facilitate chaining
ExprPtr& canonicalize(ExprPtr& expr, const NamedIndexSpec& named_indices,
                      const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Recursively canonicalizes an Expr; like mutating canonicalize() but works
/// for temporary expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be canonicalized
/// @param[in] named_indices the indices that canonicalization must not rename
/// @param[in] opts the canonicalization options
/// @return canonicalized form of \p expr_rv
ExprPtr canonicalize(ExprPtr&& expr_rv, const NamedIndexSpec& named_indices,
                     const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// @return the named indices of @p expr : the indices of its result
NamedIndexSpec named_indices(const ResultExpr& expr);

namespace detail {

/// @return the permutations of the indices of @p result that its symmetries
/// relate, as index replacement maps, each with its phase (the identity, with
/// an empty map, first); only those that map every index to one of the same
/// space, i.e. that relate elements of the block of the result
std::vector<std::pair<container::map<Index, Index>, int>> result_symmetry_group(
    const ResultExpr& result);

}  // namespace detail

/// Recursively canonicalizes the expression of a ResultExpr, with the indices
/// of the result as the named indices
/// @param[in,out] expr expression to be canonicalized
/// @param[in] opts the canonicalization options
/// @return \p expr to facilitate chaining
ResultExpr& canonicalize(
    ResultExpr& expr, const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Recursively canonicalizes the expression of a ResultExpr; like mutating
/// canonicalize() but works for temporary expressions
/// @param[in] expr expression to be canonicalized
/// @param[in] opts the canonicalization options
/// @return \p expr
[[nodiscard]] ResultExpr& canonicalize(
    ResultExpr&& expr, const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Simplifies an Expr by a combination of expansion, canonicalization, and
/// rapid_simplify
/// @param[in,out] expr expression to be simplified; may be
/// _replaced_ (i.e. `&expr` may be mutated by call)
/// @param[in] named_indices the indices that canonicalization must not rename
/// (ordinarily the external indices of @p expr); see NamedIndexSpec
/// @param[in] opts the canonicalization options; must produce the library's
/// canonical form (`opts.is_internal_form()`)
/// @sa rapid_simplify()
/// @return \p expr to facilitate chaining
ExprPtr& simplify(ExprPtr& expr, const NamedIndexSpec& named_indices,
                  const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Simplifies an Expr like simplify(ExprPtr&, ...), but works for temporary
/// expressions
/// @param[in] expr_rv rvalue-ref-to-expression to be simplified
/// @param[in] named_indices the indices that canonicalization must not rename
/// @param[in] opts the canonicalization options
/// @return simplified form of \p expr_rv
ExprPtr simplify(ExprPtr&& expr_rv, const NamedIndexSpec& named_indices,
                 const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Simplifies the expression of a ResultExpr, with the indices of the result
/// as the named indices
/// @param[in,out] expr expression to be simplified
/// @param[in] opts the canonicalization options
/// @return \p expr to facilitate chaining
ResultExpr& simplify(ResultExpr& expr,
                     const CanonicalizeOptions& opts = CanonicalizeOptions{});

/// Simplifies the expression of a ResultExpr; like mutating simplify() but
/// works for temporary expressions
/// @param[in] expr expression to be simplified
/// @param[in] opts the canonicalization options
/// @return \p expr
[[nodiscard]] ResultExpr& simplify(
    ResultExpr&& expr, const CanonicalizeOptions& opts = CanonicalizeOptions{});

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
