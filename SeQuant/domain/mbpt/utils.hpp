//
// Created by Ajay Melekamburath on 9/26/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_UTILS_HPP
#define SEQUANT_DOMAIN_MBPT_UTILS_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>

#include <SeQuant/domain/mbpt/op.hpp>

#include <range/v3/view.hpp>

namespace sequant::mbpt {

namespace detail {

/// @brief computes the commutator [a,b] = a b - b a
/// @param a first expression
/// @param b second expression
/// @return the non-canonicalized commutator expression
inline auto commutator(const ExprPtr& a, const ExprPtr& b) {
  auto result = a * b - b * a;
  return non_canon_simplify(result);
}

/// @brief computes the nested commutator up to given rank:
/// e.g., for rank=3, computes A + [A,B] + (1/2)[[A,B],B] + (1/3!)[[[A,B],B],B]
/// @param a first expression
/// @param b second expression
/// @param rank the rank of nested commutator
/// @return the non-canonicalized nested commutator expression
inline auto nested_commutator(const ExprPtr& a, const ExprPtr& b,
                              const int rank) {
  ExprPtr result = a;
  ExprPtr tmp = a;
  for (int i = 1; i <= rank; ++i) {
    auto comm = commutator(tmp, b);
    comm *= ex<Constant>(rational{1, i});
    result += comm;
    tmp = comm;
  }
  return non_canon_simplify(result);
}
}  // namespace detail

// clang-format off
/// @brief Computes the Lie similarity transformation, e^(-B) A e^B, using its Campbell expansion (DOI 10.1112/plms/s1-28.1.381) as a series of nested commutators:
/// `e^(-B) A e^B = A + [A,B] + (1/2!)[[A,B],B] + (1/3!)[[[A,B],B],B] + ...`
///
/// Notes:
/// - If \p unitary is true, the ansatz uses B - B^+ instead of B.
///
/// @param A Expression A (e.g., the Hamiltonian)
/// @param B Expression B (e.g., the cluster operator T)
/// @param commutator_rank The maximum order of nested commutators to retain (e.g. 4 for traditional coupled-cluster)
/// @param unitary If true, uses unitary generator with B-B^+
/// @param skip_clone if true, will not clone the input expression
/// @pre This function expects \p A and \p B to be composed of mbpt::Operators
// clang-format on
ExprPtr lst(ExprPtr A, ExprPtr B, size_t commutator_rank, bool unitary = false,
            bool skip_clone = false);

/// @brief Screens out terms in the expression \p expr that cannot contribute to
/// expectation value
/// @param expr input expression
/// @param skip_clone if true, will not clone the input expression
/// @return return screened expression
/// @code
/// // example usage:
/// auto expr1 = screen_vac_av(expr); // screens for <0| expr |0>
/// auto expr2 = screen_vac_av(P(2) * expr); // screens for <P(2)| expr |0>
/// @endcode
/// @pre This function expects \p input to be composed of mbpt::Operators
ExprPtr screen_vac_av(ExprPtr expr, bool skip_clone = false);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_UTILS_HPP
