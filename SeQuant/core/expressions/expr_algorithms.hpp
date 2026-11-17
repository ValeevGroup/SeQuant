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

#include <functional>
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

/// @return whether @p expr is a tensor (any AbstractTensor), or a (nested)
/// product whose factors are scalars (see Expr::is_scalar()) and tensor trees
bool is_tensor_tree(const ExprPtr& expr);

/// permutations of indices, as index replacement maps, each with its phase
using IndexPermutationGroup =
    std::vector<std::pair<container::map<Index, Index>, int>>;

/// @param symmetrizer an antisymmetrizer (labeled reserved::antisymm_label())
/// or a symmetrizer (labeled reserved::symm_label())
/// @return the permutations of the indices of @p symmetrizer that it absorbs,
/// i.e. the P for which S * X = sign(P) S * P(X) (the identity, with an empty
/// map, first): for Â the permutations of its bra and of its ket indices, the
/// phase being the product of their parities (see mbpt::A_maps()); for Ŝ the
/// permutations of its columns (the bra and the ket permuted together), with
/// phase 1 (see mbpt::S_replacement_maps()); only those that map every index
/// to one of the same space, i.e. that relate elements of the same block
IndexPermutationGroup symmetrizer_group(const Tensor& symmetrizer);

/// @return the product of the groups of @p symmetrizers (see
/// symmetrizer_group()), which act on disjoint sets of indices
IndexPermutationGroup symmetrizers_group(
    const container::svector<Tensor>& symmetrizers);

/// Applies @p f to the parts of @p term that may be replaced by their images
/// under the group of the (anti)symmetrizers of @p result that multiply them
/// (see symmetrizers_group()), as S * X = sign(P) S * P(X); @p f receives
/// the part and the (anti)symmetrizers. The group is that
/// of the symmetrizers that are present, not that of the declared symmetry of
/// @p result : that the result has a symmetry does not make its terms have
/// it, so nothing else may be replaced.
/// - S * (X_1 + X_2 + ...), S being one or more (anti)symmetrizers of the
///   result: @p f is applied to each X_k ;
/// - any other product S * X with (anti)symmetrizers S of the result as
///   factors: @p f is applied to X, the product of the other factors (with
///   the scalar), and @p term is reassembled.
///
/// The (anti)symmetrizers themselves are never passed to @p f : the image of
/// X is sign(P) P(X), and S is not permuted. Several (anti)symmetrizers are
/// only handled if their index sets are disjoint; their group is then the
/// product of their groups.
///
/// An (anti)symmetrizer of the result is a tensor labeled
/// reserved::symm_label() or reserved::antisymm_label() whose indices are
/// indices of @p result .
/// @return true if @p f was applied
bool for_each_symmetrized_term(
    ExprPtr& term, const ResultExpr& result,
    const std::function<void(
        ExprPtr& part, const container::svector<Tensor>& symmetrizers)>& f);

/// @return whether @p tensor is a projector, i.e. is labeled
/// reserved::antisymm_label() or reserved::symm_label()
bool is_projector(const Tensor& tensor);

/// Renames the indices of a term P' * @p rest , P' being the projector
/// @p projector , so that P' becomes @p target : the slots of @p projector are
/// matched to those of @p target as the symmetry of the projector allows
/// (see symmetrizer_group()), each index of @p projector is renamed to the
/// index of the matched slot of @p target , and the indices of @p rest that
/// have the labels of the indices of @p target , but are not renamed so, are
/// renamed to fresh labels. This is only valid if the indices of
/// @p projector are dummies of the term, i.e. occur nowhere else than in
/// P' * @p rest .
/// @return the renamed @p rest , with the phase of the permutation of the
/// slots of @p projector , such that P' * @p rest = @p target * result, or
/// nullopt if @p target does not have the label, the symmetries and the slot
/// spaces of @p projector (up to its symmetry) or the slots are not of the
/// supported form (aux, null or protoindexed indices)
std::optional<ExprPtr> rename_to_projector(const ExprPtr& rest,
                                           const Tensor& projector,
                                           const Tensor& target);

/// @return the distinct projectors (see is_projector()) of @p expr
container::svector<Tensor> projectors(const ExprPtr& expr);

/// Renames the indices of each term of @p expr (each summand, or @p expr
/// itself) that has a projector as a factor so that the projector is one of
/// @p targets (see rename_to_projector()). The canonicalization of a
/// projected expression (whose projector indices are not named, so that terms
/// that differ only by a relabeling of them merge, see
/// without_projector_indices()) names the projector indices as it names the
/// other dummies, so that the projector of each term carries labels of its
/// own; this gives them back the labels of the projectors of the input.
/// @param expr an expression whose projector indices are dummies
/// @param targets the projectors to give the terms, e.g. projectors() of the
/// expression before canonicalization
/// @throw Exception if a projector of a term matches no target, or more than
/// one that differ (and is none of them)
void restore_projectors(ExprPtr& expr,
                        const container::svector<Tensor>& targets);

/// @return @p projector with its indices numbered from 1 within each space
/// in the order of its slots (bra, then ket), the protoindices of an index
/// renamed as the indices of the projector they are (e.g. Â{i_7,i_3;a_9,a_4}
/// becomes Â{i_1,i_2;a_1,a_2})
Tensor conventional_projector(const Tensor& projector);

/// Renames the indices of each term of @p expr that has a projector as a
/// factor so that its projector is conventional_projector() of itself (see
/// restore_projectors()); the convention by which the projector of an
/// expression derived from operators, e.g. a projected expectation value,
/// whose projectors have temporary indices of their own, is named
void name_projectors(ExprPtr& expr);

/// restore_projectors() with @p targets if they have only ordinary (not
/// temporary, see Index::min_tmp_index()) indices, else name_projectors(): the
/// labels of projectors made with temporary indices (e.g. by operators) are
/// not meant to be kept
void restore_or_name_projectors(ExprPtr& expr,
                                const container::svector<Tensor>& targets);

/// simplify() of a projected expression whose projector indices are not
/// named, followed by restore_or_name_projectors() with the projectors of
/// @p expr as they were
/// @param expr an expression
/// @param named_indices the named indices; should not include the indices of
/// the projectors (see without_projector_indices())
ExprPtr& simplify_projected(ExprPtr& expr, const NamedIndexSpec& named_indices);

/// canonicalize() counterpart of simplify_projected()
ExprPtr& canonicalize_projected(ExprPtr& expr,
                                const NamedIndexSpec& named_indices);

}  // namespace detail

/// Recursively canonicalizes the expression of a ResultExpr, with the indices
/// of the result as the named indices; terms that an explicit
/// (anti)symmetrizer of the result multiplies and that differ by a permutation
/// of its indices that it absorbs are brought into the same form (see
/// detail::for_each_symmetrized_term())
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
