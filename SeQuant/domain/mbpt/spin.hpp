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

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace sequant::mbpt {

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
  SEQUANT_ASSERT((t.to_int32() & mask_v<Spin>) != 0);
  return static_cast<Spin>(t.to_int32() & mask_v<Spin>);
}

/// removes spin annotation in QuantumNumbersAttr by unsetting the bits used by
/// Spin
inline QuantumNumbersAttr spinannotation_remove(const QuantumNumbersAttr& t) {
  static_assert((~(~mask_v<Spin> & ~bitset::reserved) & ~bitset::reserved) ==
                    mask_v<Spin>,
                "Spin bitmask uses reserved bits");
  return t.intersection(QuantumNumbersAttr(~mask_v<Spin> & ~bitset::reserved));
}

/// removes spin annotation, if any
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring spinannotation_remove(WS&& label) {
  auto view = to_basic_string_view(label);
  SEQUANT_ASSERT(!ranges::contains(view, L'_'));
  const auto has_annotation = view.back() == L'↑' || view.back() == L'↓';
  return std::wstring{view.data(),
                      view.data() + view.size() - (has_annotation ? 1 : 0)};
}

/// adds spin annotation to IndexSpace base label or Index (full) label
template <typename WS, typename = std::enable_if_t<
                           meta::is_wstring_or_view_v<std::decay_t<WS>>>>
std::wstring spinannotation_add(WS&& label, Spin s) {
  [[maybe_unused]] auto view = to_basic_string_view(label);
  SEQUANT_ASSERT(!ranges::contains(view, L'_'));
  SEQUANT_ASSERT(view.back() != L'↑' && view.back() != L'↓');
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

/// @brief Checks that columns conserve Ms (azimuthal spin qn); only
/// filled columns (with 2 non-null indices) are considered
/// @param tensor a tensor
/// @return true if  columns conserve Ms
bool ms_conserving_columns(const AbstractTensor& tensor);

/// @brief Checks if all indices have same ms
/// @param tensor a tensor
/// @return true if all indices have same ms
bool ms_uniform_tensor(const AbstractTensor& tensor);

/// @brief Check if the number of alpha spins in bra and ket are equal;
/// beta spins will match if total number of indices is the same
/// @param tensor with spin indices
/// @return true if number of alpha spins match in bra and ket
bool can_expand(const AbstractTensor& tensor);

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

/// @brief Generates a vector of replacement maps for antisymmetrization (A)
/// tensor
/// @param A An antisymmetrizer tensor (A) (with > 2 particle indices)
/// @return Vector of replacement maps
container::svector<container::map<Index, Index>> A_maps(const Tensor& A);

/// @brief Expand a product containing the antisymmetrization (A) tensor
/// @param product a Product that may or may not include the antisymmetrizer
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expand_A_op(const ProductPtr& product);

/// @brief Write expression in terms of Symmetrizer (S operator)
/// @param product
/// @return expression pointer with Symmstrizer operator
ExprPtr symmetrize_expr(const ProductPtr& product);

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
ExprPtr expand_P_op(const ProductPtr& product);

/// @brief Expand an expression containing the particle permutation (P) tensor
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if P is present
ExprPtr expand_P_op(const ExprPtr& expr);

container::svector<container::map<Index, Index>> S_replacement_maps(
    const Tensor& S);

/// @brief Expand S operator
ExprPtr S_maps(const ExprPtr& expr);

/// @brief filters out the nonunique terms in Wang-Knizia biorthogonalization

/// WK biorthogonalization rewrites biorthogonal expressions as a projector
/// onto non-null-space (NNS)
/// applied to the biorothogonal expressions where out of each
/// group of terms related by permutation of external indices
/// those with the largest coefficients are selected.
/// This function performs the selection by forming groups of terms that
/// are equivalent modulo external index permutation (all terms in a group
/// have identical graph hashes).
/// @details This function processes a sum expression, grouping product terms by
/// hash of their canonicalized tensor network forms. For each group, it
/// retains only the terms with the largest absolute scalar coefficient.
/// @param expr The input expression, expected to be a `Sum` of `Product` terms.
/// @param ext_idxs A vector of external index groups. The function will not
/// apply the filtering logic if `ext_idxs.size()` is 2 or less.
/// @return A new `ExprPtr` representing the filtered and compacted expression.
ExprPtr WK_biorthogonalization_filter(
    ExprPtr expr,
    const container::svector<container::svector<Index>>& ext_idxs);

// clang-format off
/// @brief Traces out spin degrees of freedom from fermionic operator moments
/// @details This function is designed for integrating spin out of
/// expressions obtained as (k,k)-moment of a fermionic operator:
/// A_{p_1 .. p_k}^{q_1 .. q_k} <a^{p_1 .. p_k}_{q_1 .. q_k} op>
/// with A antisymmetric tensor, op an arbitrary fermionic operator, and
/// the expectation value is with respect to a closed-shell Fermi vacuum.
/// @param expr an input expression
/// @param ext_index_groups groups of external indices
/// @param full_expansion if true, we first fully expand the
/// antisymmetrizer, which makes spintracing expensive because of the large
/// number of terms. If false, we expand it in terms of the symmetrizer,
/// which results in a partial expansion.
/// @return the spin-free form of expr
/// @warning the "antisymmetrizer" tensor A is assumed to be at the front of each tensor
/// network, hence must use "complete" canonicalization to produce the input expression.
// clang-format on
ExprPtr closed_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups = {},
    bool full_expansion = false);

container::svector<ResultExpr> closed_shell_spintrace(
    const ResultExpr& expr, bool full_expansion = false);

// clang-format off
/// biorthogonalization variants
enum class BiorthogonalizationMethod {
  /// standard biorthogonalization method
  V1,
  /// improved Wang-Knizia biorthogonalization (DOI 10.48550/arXiv.1805.00565), with factored out the Non-Null-Space (NNS) projector (expressed as a linear combination of permutations)
  V2
};
// clang-format on

/// controls behavior of biorthogonal closed-shell spin-tracing
struct ClosedShellCCSpintraceOptions {
  BiorthogonalizationMethod method = BiorthogonalizationMethod::V2;
  /// set to true to use sequant::spintrace which does not assume closed-shell
  /// (spin-free) basis and thus has an exponential cost;
  /// the default is to use closed_shell_spintrace, which is more efficient
  bool naive_spintrace = false;
};

// clang-format off
/// @brief like closed_shell_spintrace but transforms spin-free moments to biorthogonal form
/// The algorithm (most steps are the same in V1 and V2):
/// - factor out symmetrizer from the antisymmetrizer,
///   "set it aside" (remove from the expression), and
///   expand the rest of the antisymmetrizer.
/// - spin-trace the resulting expression assuming spin-restricted basis and
///   closed-shell reference
/// - biorthogonalize
/// - V2 only: collect terms that differ only by external index permutation
/// (i.e. have same colored graph hash), for each group select the terms with
/// the largest coefficients)
/// - apply the particle symmetrizer and simplify
///
/// The V2 method produces an expression that becomes equivalent that of V1 by
/// application of if the biorthogonal NNS projector (particular linear combination of
/// permutation operators) is applied. The biorthogonal NNS projection should be
/// performed numerically.
/// @warning the antisymmetrizer is assumed to be at the front of each tensor
/// network, hence must use "complete" canonicalization to produce the input expression.
/// @param expr input expression
/// @param options optional parameter controlling the algorithm selection and
///                other traits; the default is to use
///                the V2 method
// clang-format on
ExprPtr closed_shell_CC_spintrace(ExprPtr const& expr,
                                  ClosedShellCCSpintraceOptions options = {});

/// @sa closed_shell_CC_spintrace
ExprPtr closed_shell_CC_spintrace_v1(
    ExprPtr const& expr,
    ClosedShellCCSpintraceOptions options = {
        .method = BiorthogonalizationMethod::V1, .naive_spintrace = false});

/// @sa closed_shell_CC_spintrace
ExprPtr closed_shell_CC_spintrace_v2(
    ExprPtr const& expr,
    ClosedShellCCSpintraceOptions options = {
        .method = BiorthogonalizationMethod::V2, .naive_spintrace = false});

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

// clang-format off
/// @brief Traces out spin degrees of freedom from fermionic operator moments
/// @details This function is designed for integrating spin out of
/// expressions obtained as (k,k)-moment of a fermionic operator:
/// A_{p_1 .. p_k}^{q_1 .. q_k} <a^{p_1 .. p_k}_{q_1 .. q_k} op>
/// with A antisymmetric tensor, op an arbitrary fermionic operator, and
/// the expectation value is with respect to an open-shell Fermi vacuum.
/// @param expr an input expression
/// @param ext_index_groups groups of external indices
/// @param target_spin_case if non-null specifies the target spin case of the external indices, else produces equations for all spin cases
/// @return vector of spin-traced expressions for each spincase
/// @warning the "antisymmetrizer" tensor A is assumed to be at the front of each tensor
/// network, hence must use "complete" canonicalization to produce the input expression.
/// @note this performs full expansion of the antisymmetrizer
/// @sa open_shell_CC_spintrace
// clang-format on
std::vector<ExprPtr> open_shell_spintrace(
    const ExprPtr& expr,
    const container::svector<container::svector<Index>>& ext_index_groups,
    std::optional<int> target_spin_case = std::nullopt);

// clang-format off
/// @brief Like open_shell_spintrace but uses minimal expansion of the antisymmetrizer
/// @details This function is designed for integrating spin out of
/// expressions obtained as (k,k)-moment of a fermionic operator:
/// A_{p_1 .. p_k}^{q_1 .. q_k} <a^{p_1 .. p_k}_{q_1 .. q_k} op>
/// with A antisymmetric tensor, op an arbitrary fermionic operator, and
/// the expectation value is with respect to an open-shell Fermi vacuum.
/// Antisymmetrizer is expanded partially to produce antisymmetrizer for
/// spin-up and spin-down columns.
/// @param expr the input expression
/// @return vector of spin-traced expressions for each spincase
/// @warning the "antisymmetrizer" tensor A is assumed to be at the front of each tensor
/// network, hence must use "complete" canonicalization to produce the input expression.
// clang-format on
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

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPIN_HPP
