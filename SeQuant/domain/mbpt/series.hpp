#ifndef SEQUANT_DOMAIN_MBPT_SERIES_HPP
#define SEQUANT_DOMAIN_MBPT_SERIES_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>

#include <cstddef>
#include <functional>
#include <string_view>

namespace sequant {
class Index;
class IndexBasis;
}  // namespace sequant

namespace sequant::mbpt {

/// @return @p expr where every tensor labelled @p label , or its adjoint
/// label, is replaced by its Taylor series in @p strength through @p order :
/// p → Σ_k (1/k!) ε^k p⁽ᵏ⁾, the coefficient p⁽ᵏ⁾ being p labelled by
/// detail::decorate_with_pert_order(label, k), with its slots and symmetries.
/// For a strength tensor ε{;;z} the k-th term is
/// (1/k!) p⁽ᵏ⁾{…;z_1…z_k} ε{;;z_1}…ε{;;z_k}, the z_j temporary indices
/// (Index::make_tmp_index) named unlike any index of @p expr , and appended to
/// the leaf's auxiliary indices. The coefficient of p⁺ is the adjoint of that
/// of p, p being bra-ket Nonsymm as adjoint() labels no other tensor with ⁺.
/// Nothing is expanded or truncated. A strength that already occurs in
/// @p expr is accepted: it is the same variable.
/// @param label a nonempty label that is not an adjoint label
/// @param strength an unconjugated Variable, or a tensor with one auxiliary
/// index, without protoindices, and no other slot; the index space of that
/// index is the caller's choice. It is labelled neither @p label nor like a
/// coefficient
/// @throw Exception for a null argument, another @p label or @p strength ,
/// @p order above 9, an mbpt::Operator or a Variable labelled @p label in
/// @p expr , a leaf of @p label under a TensorPower, or a tensor in @p expr
/// labelled like a coefficient p⁽ᵏ⁾ (or its adjoint), 1 ≤ k ≤ @p order
[[nodiscard]] ExprPtr series(const ExprPtr& expr, std::wstring_view label,
                             const ExprPtr& strength, std::size_t order);

namespace detail {

/// calls @p f with each atom of @p e , @p e itself if it is one; the base of
/// a TensorPower is not an atom of its own
template <typename E, typename F>
void for_each_atom(E& e, const F& f) {
  if (e->is_atom())
    f(e);
  else
    e->visit(f, /*atoms_only=*/true);
}

/// @return the base of @p atom if it is a TensorPower or a Power, else @p atom
[[nodiscard]] const ExprPtr& base_of(const ExprPtr& atom);

/// @return whether @p e is an mbpt::Operator
[[nodiscard]] bool is_operator(const Expr& e);

/// @return a tensor labelled @p label with the bra @p b , the ket @p k , the
/// auxiliary indices @p a and the symmetries of @p model
[[nodiscard]] ExprPtr tensor_like(const Tensor& model, std::wstring_view label,
                                  container::svector<Index> b,
                                  container::svector<Index> k,
                                  container::svector<Index> a);

/// @return a function that returns, at each call, a new temporary index of
/// the basis it is given (Index::make_tmp_index), named unlike every index of
/// @p expr
[[nodiscard]] std::function<Index(const IndexBasis&)> fresh_indices(
    const ExprPtr& expr);

/// @throw Exception what series(expr, label, strength, order) throws, except,
/// if @p expanded is false, for a Variable labelled @p label or a leaf of
/// @p label under a TensorPower: leaves of @p label are not expanded
void check_series(const ExprPtr& expr, std::wstring_view label,
                  const ExprPtr& strength, std::size_t order,
                  bool expanded = true);

/// @return the term of order @p k of the series of the leaf @p t in
/// @p strength , as series() writes it, @p t itself for k = 0; @p label is
/// the base label of @p t , the arguments are those check_series() accepts,
/// and @p fresh makes the perturbation indices
[[nodiscard]] ExprPtr series_term(
    const Tensor& t, std::wstring_view label, const ExprPtr& strength,
    std::size_t k, const std::function<Index(const IndexBasis&)>& fresh);

/// @return a copy of @p expr where every tensor leaf for which @p make
/// returns non-null is replaced by what it returns; @p expr is not modified.
/// @p make is called for every Tensor leaf, symmetrizers and Kronecker deltas
/// included, but not for the base of a TensorPower; its argument is valid only
/// during the call. It returns null to keep the leaf, else a new expression
/// that no other call returns; indices it adds must be fresh in all of @p expr
[[nodiscard]] ExprPtr substituted(
    const ExprPtr& expr, const std::function<ExprPtr(const Tensor&)>& make);

}  // namespace detail

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SERIES_HPP
