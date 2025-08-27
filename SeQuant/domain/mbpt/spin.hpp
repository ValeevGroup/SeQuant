//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_DOMAIN_MBPT_SPIN_HPP
#define SEQUANT_DOMAIN_MBPT_SPIN_HPP

#include <SeQuant/domain/mbpt/fwd.hpp>

#include <SeQuant/domain/mbpt/space_qns.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace sequant {

namespace mbpt {

// Spin is a scoped enum, hence not implicitly convertible to
// QuantumNumbersAttr::bitset_t
static_assert(!std::is_convertible_v<sequant::mbpt::Spin, bitset_t>);
// QuantumNumbersAttr::bitset_t cannot be constructed from Spin
static_assert(!std::is_constructible_v<bitset_t, sequant::mbpt::Spin>);
// but Spin can be cast to QuantumNumbersAttr::bitset_t
static_assert(meta::is_statically_castable_v<sequant::mbpt::Spin, bitset_t>);
// Spin cannot be cast to nonsense ...
static_assert(
    !meta::is_statically_castable_v<sequant::mbpt::Spin, std::string>);

inline Spin operator~(Spin s) {
  return static_cast<Spin>(~(static_cast<bitset_t>(s)));
}
inline Spin operator|(Spin s1, Spin s2) {
  return static_cast<Spin>(static_cast<bitset_t>(s1) |
                           static_cast<bitset_t>(s2));
}
inline Spin operator&(Spin s1, Spin s2) {
  return static_cast<Spin>(static_cast<bitset_t>(s1) &
                           static_cast<bitset_t>(s2));
}

/// converts QuantumNumbersAttr to Spin
/// @note this filters out all bits not used in Spin
inline Spin to_spin(const QuantumNumbersAttr& t) {
  assert((t.to_int32() & static_cast<int>(Spin::mask)) != 0);
  return static_cast<Spin>(static_cast<Spin>(t.to_int32()) & Spin::mask);
}

/// removes spin annotation in QuantumNumbersAttr by unsetting the bits used by
/// Spin
inline QuantumNumbersAttr spinannotation_remove(const QuantumNumbersAttr& t) {
  return t.intersection(QuantumNumbersAttr(~Spin::mask));
}

/// removes spin annotation, if any
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring spinannotation_remove(WS&& label) {
  auto view = to_basic_string_view(label);
  assert(!ranges::contains(view, L'_'));
  const auto has_annotation = view.back() == L'↑' || view.back() == L'↓';
  return std::wstring{view.data(),
                      view.data() + view.size() - (has_annotation ? 1 : 0)};
}

/// adds spin annotation to IndexSpace base label or Index (full) label
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring spinannotation_add(WS&& label, Spin s) {
  [[maybe_unused]] auto view = to_basic_string_view(label);
  assert(!ranges::contains(view, L'_'));
  assert(view.back() != L'↑' && view.back() != L'↓');
  switch (s) {
    case Spin::any:
      return to_wstring(std::forward<WS>(label));
    case Spin::alpha:
      return to_wstring(std::forward<WS>(label)) + L'↑';
    case Spin::beta:
      return to_wstring(std::forward<WS>(label)) + L'↓';
    case Spin::null:
      SEQUANT_ABORT("Invalid spin quantum number");
  }

  SEQUANT_UNREACHABLE;
}

/// replaces spin annotation to
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring spinannotation_replacе(WS&& label, Spin s) {
  auto label_sf = spinannotation_remove(std::forward<WS>(label));
  return spinannotation_add(label_sf, s);
}

}  // namespace mbpt

// make alpha-spin idx
[[nodiscard]] Index make_spinalpha(const Index& idx);

// make beta-spin idx
[[nodiscard]] Index make_spinbeta(const Index& idx);

// make null-spin idx
[[nodiscard]] Index make_spinfree(const Index& idx);

/// @brief Preserving particle symmetry, swaps bra and ket labels on all tensors
/// in an expression
/// @param expr ExprPtr to transform
/// @return transformed expression
ExprPtr swap_bra_ket(const ExprPtr& expr);

/// @brief Adds spins to indices in an expression with a replacement map
/// @param expr an expression pointer
/// @param index_replacements a map of pairs containing the index and the
/// corresponding replacement
/// @return expr the ExprPtr with substituted indices
ExprPtr append_spin(const ExprPtr& expr,
                    const container::map<Index, Index>& index_replacements);

/// @brief Removes spin label from all indices in an expression
/// @param expr an ExprPtr with spin indices
/// @return expr an ExprPtr with spin labels removed
ExprPtr remove_spin(const ExprPtr& expr);

/// @brief Checks the spin symmetry of a pair of indices corresponding to a
/// particle in tensor notation
/// @param tensor a tensor with indices containing spin labels
/// @return true if spin symmetry matches for all pairs of indices
bool spin_symm_tensor(const Tensor& tensor);

/// @brief Checks if spin labels on a tensor are same
/// @param tensor a tensor with indices containing spin labels
/// @return true if all spin labels on a tensor are the same
bool same_spin_tensor(const Tensor& tensor);

/// @brief Check if the number of alpha spins in bra and ket are equal;
/// beta spins will match if total number of indices is the same
/// @param tensor with spin indices
/// @return true if number of alpha spins match in bra and ket
bool can_expand(const Tensor& tensor);

/// @brief expand an antisymmetric tensor
///
/// @details For spin-indices, the tensor is NOT expanded if all spin-labels
/// are either alpha or beta
/// @param[in] tensor a tensor
/// @param[in] skip_spinsymm is true, will not expand tensors whose indices
///            all have the same spin [default=false]
/// @return an ExprPtr containing the sum of expanded terms, if antisymmetric
ExprPtr expand_antisymm(const Tensor& tensor, bool skip_spinsymm = false);

/// @brief expands all antisymmetric tensors in an expression
/// @param expr an expression to expand
/// @param skip_spinsymm is true, will not expand tensors whose indices all have
/// same spin [default=false]
/// @return an expression pointer with expanded tensors as a sum
ExprPtr expand_antisymm(const ExprPtr& expr, bool skip_spinsymm = false);

/// @brief Check if a tensor with a certain label is present in an expression
/// @param expr input expression
/// @param label tensor label to find in the expression
/// @return true if tensor with given label is found
bool has_tensor(const ExprPtr& expr, std::wstring label);

/// @brief Generates a vector of replacement maps for antisymmetrization (A)
/// tensor
/// @param A An antisymmetrizer tensor (A) (with > 2 particle indices)
/// @return Vector of replacement maps
container::svector<container::map<Index, Index>> A_maps(const Tensor& A);

/// @brief Removes tensor with a certain label from product
/// @param product A product expression
/// @param label Label of the tensor to remove
/// @return ExprPtr with the tensor removed
ExprPtr remove_tensor(const Product& product, std::wstring label);

/// @brief Removes tensor with a certain label from an expression
/// @param expr An expression pointer
/// @param label Label of the tensor to remove
/// @return ExprPtr with the tensor removed
ExprPtr remove_tensor(const ExprPtr& expr, std::wstring label);

/// @brief Expand a product containing the antisymmetrization (A) tensor
/// @param product a Product that may or may not include the antisymmetrizer
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expand_A_op(const Product& product);

/// @brief Write expression in terms of Symmetrizer (S operator)
/// @param product
/// @return expression pointer with Symmstrizer operator
ExprPtr symmetrize_expr(const Product& product);

/// @brief Expand an expression containing the antisymmetrization (A) tensor
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr symmetrize_expr(const ExprPtr& expr);

/// @brief Expand an expression containing the antisymmetrization (A) tensor
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expand_A_op(const ExprPtr& expr);

/// @brief Generates a vector of replacement maps for particle permutation
/// operator
/// @param P a particle permutation operator (with > 2 particle indices)
/// @return Vector of replacement maps
container::svector<container::map<Index, Index>> P_maps(const Tensor& P);

/// @brief Expand a product containing the particle permutation (P) tensor
/// @param product a Product that may or may not contain P tensor
/// @return an ExprPtr containing sum of expanded terms if P is present
ExprPtr expand_P_op(const Product& product);

/// @brief Expand an expression containing the particle permutation (P) tensor
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if P is present
ExprPtr expand_P_op(const ExprPtr& expr);

container::svector<container::map<Index, Index>> S_replacement_maps(
    const Tensor& S);

/// @brief Expand S operator
ExprPtr S_maps(const ExprPtr& expr);

/// @brief Filters and compacts expression based on a hash-based heuristic.
/// @details This function processes a sum expression, grouping product terms by
/// hash of their canonicalized tensor network forms. For each group, it
/// retains only the terms with the largest absolute scalar coefficient. This is
/// intended to compact large expressions by keeping only the dominant terms
/// after a hash-based grouping.
/// This is particularly useful to eliminate the redudancy caused by
/// biorthogonal trnasformation in coupled-cluster equations.
/// @param expr The input expression, expected to be a `Sum` of `Product` terms.
/// @param ext_idxs A vector of external index groups. The function will not
/// apply the filtering logic if `ext_idxs.size()` is 2 or less.
/// @return A new `ExprPtr` representing the filtered and compacted expression.
ExprPtr hash_filter_compact_set(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs);

/// @brief Transforms an expression from spin orbital to spatial orbitals
/// @details This functions is designed for integrating spin out of expression
/// with Coupled Cluster equations in mind.
/// @attention This function may fail on arbitrarily created expressions that
/// lacks proper index attributes.
/// @param expr ExprPtr with spin orbital indices
/// @param ext_index_groups groups of external indices
/// @param full_expansion if true, we first fully expand the
/// anti-symmetrizer, which makes spintracing expensive because of the large
/// number of terms. If false, we exapnd it in terms of the symmetrizer,
/// which results in a partial expansion.
/// @return an expression with spin integrated/adapted
ExprPtr closed_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups = {},
    bool full_expansion = false);

container::svector<ResultExpr> closed_shell_spintrace(
    const ResultExpr& expr, bool full_expansion = false);

/// @brief Transforms Coupled cluster from spin orbital to spatial orbitals
/// @details The external indices are deduced from Antisymmetrization operator
/// @param expr ExprPtr to Sum type with spin orbital indices
/// @return an expression with spin integrated/adapted
ExprPtr closed_shell_CC_spintrace(ExprPtr const& expr);

/// \brief Same as \c closed_shell_CC_spintrace except it treats expression
/// differently. It generates a compact set expression with spin integration.
/// It applies a specific filtering and compaction procedure. It internally
/// expands the symmetrizer to access all individual terms. The function
/// then uses a heuristic hash filter method to identify and retain only the
/// most significant terms, thereby reducing the total number of terms that
/// need to be evaluated.
/// After the evaluation of the final expression, a `mbpt-cleanup` function
/// is applied to restore the effects of deleted terms.
/// To work correctly with the cleanup function, the intermediate expression
/// is internally multiplied by a rescaling factor of $n!/(n!-1)$.
/// If you need to print the compact-set equations, you must manually compensate
/// for this factor and multiply the expression by $(n!-1)/n!$.
/// @param expr ExprPtr to Sum type with spin orbital indices
/// @return An expression pointer representing the most compact set of
/// spin-integrated terms.
ExprPtr closed_shell_CC_spintrace_compact_set(ExprPtr const& expr);

/// \brief Same as \c closed_shell_CC_spintrace except internally uses
///        \c sequant::spintrace instead of sequant::closed_shell_spintrace.
ExprPtr closed_shell_CC_spintrace_rigorous(ExprPtr const& expr);

/// Collect all indices from an expression
container::set<Index, Index::LabelCompare> index_list(const ExprPtr& expr);

/// @brief Swap spin labels in a tensor
Tensor swap_spin(const Tensor& t);

/// @brief Swap spin labels in an expression
ExprPtr swap_spin(const ExprPtr& expr);

/// @brief Merge operators into a single operator (designed for P operator)
/// @warning Repetition of indices is allowed in a bra or a ket
ExprPtr merge_tensors(const Tensor& O1, const Tensor& O2);

/// @brief Vector of Anti-symmetrizers for spin-traced open-shell expr
std::vector<ExprPtr> open_shell_A_op(const Tensor& A);

/// @brief Generate a vector of permutation operators for partial expansion of
/// antisymmstrizer
/// @details The antisymmetrizer need not be fully expanded to spin-trace for
/// open-shell case. By expanding only the unlike spin terms, the
/// antisymmetrizer is pereserved for same-spin particle indices.
/// @param A Antisymmetrizer tensor produced in Coupled Cluster
/// @return a vector of expression pointers containing permutation operators as
/// a sum
/// @warning This function assumes the antisymmetrizer (A) has a canonical form
std::vector<ExprPtr> open_shell_P_op_vector(const Tensor& A);

/// @brief Generates spin expressions to be used for open-shell coupled cluster
/// @details Every spin combination of external indices will have all spin
/// combinations of internal indices.
/// @param expr ExprPtr with spin orbital indices
/// @param ext_index_groups groups of external indices
/// @param single_spin_case Calculate open-shell expression for a specific spin
/// case
/// @return a vector of expr ptrs with spin expressions
std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    const int single_spin_case = 0);

/// @brief Transforms Coupled cluster from spin orbital to spatial orbitals
/// @param expr ExprPtr to Sum type with spin orbital indices
/// @return a vector of spin expressions for open-shell reference
std::vector<ExprPtr> open_shell_CC_spintrace(const ExprPtr& expr);

/// @brief Transforms an expression from spin orbital to spin-free (spatial)
/// orbital form
/// @details Given an expression, this function extracts all indices and adds a
/// spin attribute to all the indices in the expression. A map is generated with
/// all possible spin permutations and substituted in the expression. Based on
/// spin symmetry of particle indices: the non-zero terms are expanded, the spin
/// labels removed and a sum of all non-zero expressions is returned.
/// @param expr ExprPtr with spin orbital indices
/// @param ext_index_groups groups of external indices
/// @param spinfree_index_spaces if true, will assume that all index spaces are
/// spin-free and will remove spin from all indices of the result before
/// returning
/// @return an expression with spin integrated/adapted
/// @warning The result of this function is not simplified since this is a
/// building block for more specialized spin-tracing functions
ExprPtr spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups = {},
    bool spinfree_index_spaces = true);

container::svector<ResultExpr> spintrace(const ResultExpr& expr,
                                         bool spinfree_index_spaces = true);

/// @brief Factorize S out of terms
/// @details Given an expression, permute indices and check if a given product
/// @param expression Expression pointer
/// @param fast_method use hash maps (memory intensive) for faster evaluation
/// @param ext_index_groups External index groups to generate the S operator
/// @return ExprPtr with terms with S operator as a factor
ExprPtr factorize_S(const ExprPtr& expression,
                    std::initializer_list<IndexList> ext_index_groups,
                    bool fast_method = true);

}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_SPIN_HPP
