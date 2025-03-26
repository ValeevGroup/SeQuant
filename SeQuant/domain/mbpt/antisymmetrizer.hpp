//
// created by Conner Masteran 06/1/2021
//

#ifndef SEQUANT_DOMAIN_MBPT_ANTISYMMETRIZER_HPP
#define SEQUANT_DOMAIN_MBPT_ANTISYMMETRIZER_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <algorithm>
#include <iostream>
#include <vector>

namespace sequant {

/// @brief generates all unique permutations of a product where products
/// differing only by internal tensor antisymmetry are non-unique. i.e.
/// \f$ a^{i_1 i_2}_{a_1 a_2} = - a^{i_2 i_1}_{a_1 a_2} \f$.
/// RHS is non-unique in this context.
class antisymm_element {
  using IndexGroup = std::pair<size_t, size_t>;

 private:
  std::vector<IndexGroup>
      index_group;  // where each tensor begins and ends. assumes particle
                    // conserving ops. needed to keep track of canonical
                    // ordering
  std::vector<std::pair<int, std::vector<Index>>>
      unique_bras_list;  // list of unique bra orderings with associated integer
                         // for the sign
  std::vector<std::pair<int, std::vector<Index>>>
      unique_kets_list;  // list of unique ket orderings with associated integer
                         // for the sign
  ExprPtr current_product;  // used to keep track of the original expression
                            // recieved by the constructor

  /// generates all possible permutations while observing the canonical ordering
  /// of each tensor/operator function kept general to work with other data
  /// types, algorithm does not require sequant::Index objects
  /// @param ordered_indices a list of T objects/indices in the canonical
  /// ordering
  /// @return a list of all permutations each with an associated sign
  template <typename T>
  std::vector<std::pair<int, std::vector<T>>> gen_antisymm_unique(
      std::vector<T> ordered_indices);

 public:
  // takes a sequant::ExprPtr and generates all antisymmetric unique
  // permutations of that expression. requires that ex_ is a product expression
  // at this point
  // @param ex_ as product
  // populates a result ExprPtr that the user can grab. result is in general a
  // Sum.
  antisymm_element(ExprPtr ex_);

  std::vector<Index> sorted_bra_indices;  // The original order of the upper
                                          // indices on a given term
  std::vector<Index> sorted_ket_indices;  // the original order of the lower
                                          // indices on a given term
  ExprPtr result;
};

/// @brief simple class to call antisymm_element on only products
class antisymmetrize {
 public:
  ExprPtr result = ex<Constant>(0);

  antisymmetrize(ExprPtr s);
};

namespace antisymm {

/// function which counts number of closed loops from the starting order of
/// upper and lower indices and the contracted indices. since the ordering of
/// the new contracted indices is arbitrary, the algorithm searched for the
/// upper index which would connect to the lower index checks if contracted
/// lower index closes the loop, if not, continue searching until
/// the corresponding upper index is not present, or the loop closes.
/// keep track of which indices are
/// used so that loops are not double counted
/// @param init_upper initial order of upper indices before contraction,
/// @param init_lower initial
/// order of lower indices before contraction,
/// @param new_upper set of contracted upper
/// indices,
/// @param new_lower set of lower contracted indices.
/// @out the number of loops.
// TODO Test this function extensively and add more asserts
int num_closed_loops(std::vector<Index> init_upper,
                     std::vector<Index> init_lower,
                     std::vector<Index> new_upper,
                     std::vector<Index> new_lower);

/// for the mnemonic rules of spin summing to work, each individual
/// tensor/FNOperator needs to maximally resemble the original tensor,
/// so indices may need to be swapped (and the sign changed).
/// this does not matter for spin-orbital expressions, but the rules become
/// incredibly simple if things stay maximally the same.
/// @param original_upper original upper indices
/// @param original_lower original lower indices
/// @param expression starting expression
ExprPtr max_similarity(const std::vector<Index>& original_upper,
                       const std::vector<Index>& original_lower,
                       ExprPtr expression);

/// not a general spin-summing procedure, implementation for a known singlet
/// state for the prefactor rules to apply.
/// @param original_upper an antisymm_element object strictly so
/// the original ordering of the indices is known
/// @param original_lower bool singlet_state? the looping rules and
/// contraction prefactors are a direct result of the singlet state
/// approximation to densities.
// TODO: use a generalized spin summation for non-singlet states
ExprPtr spin_sum(std::vector<Index> original_upper,
                 std::vector<Index> original_lower, ExprPtr expression,
                 bool singlet_state = true);

}  // namespace antisymm
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_ANTISYMMETRIZER_HPP
