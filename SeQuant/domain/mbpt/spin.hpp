//
// Created by Eduard Valeyev on 2019-02-27.
//

#ifndef SEQUANT_SPIN_HPP
#define SEQUANT_SPIN_HPP

#include <SeQuant/core/tensor_network.hpp>
#include <unordered_map>
#include "SeQuant/core/tensor.hpp"

namespace sequant {

/// @brief Applies index replacement rules to an ExprPtr
/// @param expr ExprPtr to transform
/// @param index_replacements index replacement map
/// @param scaling_factor to scale the result
/// @return a substituted and scaled expression pointer
ExprPtr transform_expression(const ExprPtr& expr,
                             const std::map<Index, Index>& index_replacements,
                             double scaling_factor = 1.0);

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
ExprPtr append_spin(ExprPtr& expr,
                    const std::map<Index, Index>& index_replacements);

/// @brief Removes spin label from all indices in an expression
/// @param expr an ExprPtr with spin indices
/// @return expr an ExprPtr with spin labels removed
ExprPtr remove_spin(ExprPtr& expr);

/// @brief Checks the spin symmetry of a pair of indices corresponding to a
/// particle in tensor notation
/// @param tensor a tensor with indices containing spin labels
/// @return true if spin symmetry matches for all pairs of indices
bool is_tensor_spin_symm(const Tensor& tensor);

/// @brief Check if the number of alpha spins in bra and ket are equal;
/// beta spins will match if total number of indices is the same
/// @param tensor with spin indices
/// @return true if number of alpha spins match in bra and ket
bool can_expand(const Tensor& tensor);

/// @brief expand an antisymmetric tensor
/// @param tensor a tensor from a product
/// @return an ExprPtr containing the sum of expanded terms if antisymmetric OR
/// @return an ExprPtr containing the tensor otherwise
ExprPtr expand_antisymm(const Tensor& tensor);

// TODO: Correct this function
/// @brief expands all antisymmetric tensors in a product
/// @param expr an expression pointer to expand
/// @return an expression pointer with expanded tensors as a sum
ExprPtr expand_antisymm(const ExprPtr& expr);

/// @brief Check if label A is present in an expression pointer
/// @detail This function assumes canonical ordering of tensors
/// in product for performance
/// @input expr an expression pointer
/// @return true if A is present in input
bool has_A_label(const ExprPtr& expr);

/// @brief Check if a tensor with a certain label is present in an expression
/// @param expr input expression
/// @param label tensor label to find in the expression
/// @return true if tensor with given label is found
bool has_tensor_label(const ExprPtr& expr, std::wstring label);

/// @brief Check if an operator with a certain label is present in an expression
/// @detailed Specifically designed for
/// @param expr input expression
/// @param label tensor label to find in the expression
/// @return true if tensor with given label is found
bool has_operator_label(const ExprPtr& expr, std::wstring label);

/// @brief Generates a vector of replacement maps for Antisymmetrizer operator
/// @param A An antisymmetrizer tensor (A) (with > 2 particle indices)
/// @return Vector of replacement maps
std::vector<std::map<Index, Index>> A_replacement_map(const Tensor& A);

/// @brief Removes tensor with a certain label from product
/// @param product A product expression
/// @param label Label of the tensor to remove
/// @return ExprPtr with the tensor removed
ExprPtr remove_tensor_from_product(const Product& product, std::wstring label);

/// @brief Expand a product containing the Antisymmetrization (A) operator
/// @param A product term with/without A operator
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expand_A_operator(const Product& product);

/// @brief Write expression in terms of Symmetrizer (S operator)
/// @param product
/// @return expression pointer with Symmstrizer operator
ExprPtr expr_symmetrize(const Product& product);

/// @brief Expand an expression containing the Antisymmetrization (A) operator
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expr_symmetrize(const ExprPtr& expr);

/// @brief Expand an expression containing the Antisymmetrization (A) operator
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if A is present
ExprPtr expand_A_operator(const ExprPtr& expr);

/// @brief Generates a vector of replacement maps for particle permutation
/// operator
/// @param P A particle permutation operator (with > 2 particle indices)
/// @return Vector of replacement maps
std::vector<std::map<Index, Index>> P_replacement_map(const Tensor& P);

/// @brief Expand a product containing the particle permutation (P) operator
/// @param A product term with/without P operator
/// @return an ExprPtr containing sum of expanded terms if P is present
ExprPtr expand_P_operator(const Product& product);

/// @brief Expand an expression containing the particle permutation (P) operator
/// @param expr any ExprPtr
/// @return an ExprPtr containing sum of expanded terms if P is present
ExprPtr expand_P_operator(const ExprPtr& expr);

std::vector<std::map<Index, Index>> S_replacement_maps(const Tensor& S);

/// @brief Expand S operator
ExprPtr expand_S_operator(const ExprPtr& expr);

/// @brief Returns the number of cycles
/// @detailed
/// @param
/// @return
int count_cycles(const container::svector<int, 6>& vec1,
                 const container::svector<int, 6>& vec2);

/// @brief Transforms an expression from spin orbital to spatial orbitals
/// @detailed This functions is designed for integrating spin out of expression
/// with Coupled Cluster equations in mind.
/// @attention This function may fail on arbitrarily created expressions that
/// lacks proper index attributes.
/// @param expr ExprPtr with spin orbital indices
/// @param ext_index_groups groups of external indices
/// @return an expression with spin integrated/adapted
ExprPtr closed_shell_spintrace(
    const ExprPtr& expression,
    const container::vector<container::vector<Index>> ext_index_groups = {{}});

/// @brief Transforms Coupled cluster from spin orbital to spatial orbitals
/// @detailed The external indices are deduced from Antisymmetrization operator
/// @param expr ExprPtr with spin orbital indices
/// @return an expression with spin integrated/adapted
ExprPtr closed_shell_cc_spintrace(const ExprPtr& expr);

/// Collect all indices from an expression
auto index_list(const ExprPtr& expr);

std::vector<ExprPtr> open_shell_spintrace(const ExprPtr& expr,
    const std::vector<std::vector<Index>> ext_index_groups = {{}});

/// @brief Transforms an expression from spin orbital to spatial orbitals
/// @detailed Given an expression, this function extracts all indices and adds a
/// spin attribute to all the indices in the expression. A map is generated with
/// all possible spin permutations and substituted in the expression. Based on
/// spin symmetry of particle indices: the non-zero terms are expanded, the spin
/// labels removed and a sum of all non-zero expressions is returned.
/// @param expr ExprPtr with spin orbital indices
/// @param ext_index_groups groups of external indices
/// @return an expression with spin integrated/adapted
ExprPtr spintrace(
    ExprPtr expression,
    container::vector<container::vector<Index>> ext_index_groups = {{}});

/// @brief Factorize S out of terms
/// @detailed Given an expression, permute indices and check if a given product
/// @param expression Expression pointer
/// @param fast_method use hash maps (memory intensive) for faster evaluation
/// @param ext_index_groups External index groups to geenrate S operator
/// @return ExprPtr with terms with S operator as a factor
ExprPtr factorize_S_operator(
    const ExprPtr& expression,
    const std::initializer_list<IndexList> ext_index_groups = {{}},
    const bool fast_method = true);

}  // namespace sequant

#endif  // SEQUANT_SPIN_HPP
