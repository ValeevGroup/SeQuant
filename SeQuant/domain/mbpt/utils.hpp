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

/// @brief Computes the net quantum number change produced by a product of
/// mbpt::Operators
/// @param pdt a product of mbpt::op_t
/// @return net quantum number change produced by \p pdt
/// @pre This function expects \p pdt to be a Product of mbpt::op_t
inline qns_t compute_qnc(const Product& pdt) {
  qns_t qns;
  for (auto& factor : ranges::views::reverse(pdt.factors())) {
    if (factor->template is<Variable>()) continue;  // skip Variables
    assert(factor->template is<op_t>());  // everything else must be op_t
    const auto& op = factor->template as<op_t>();
    qns = op(qns);
  }
  return qns;
}

}  // namespace detail

// clang-format off
/// @brief Computes the similarity transformation e^(-B) A e^B using a commutator series truncated up to the specified order. Uses the Baker-Campbell-Hausdorff expansion.
/// The transformation is given by: e^(-B) A e^B = A + [A,B] + (1/2!)[[A,B],B] + (1/3!)[[[A,B],B],B] + ...
///
/// Notes:
/// - This implementation is specific to the exponential ansatz B = e^X, not the general similarity transform B^{-1} A B.
/// - If \p unitary is true, the ansatz uses B - B^+ instead of B, corresponding to a unitary coupled-cluster type transformation.
///
/// @param A Expression A (e.g., the Hamiltonian)
/// @param B Expression B (e.g., the cluster operator T)
/// @param commutator_rank The maximum order of nested commutators to retain (e.g. 4 for traditional coupled-cluster)
/// @param unitary If true, uses unitary ansatz with B-B^+
/// @param skip_clone if true, will not clone the input expression
/// @pre This function expects \p A and \p B to be composed of mbpt::Operators
// clang-format on
ExprPtr sim_tr(ExprPtr A, const ExprPtr& B, size_t commutator_rank,
               bool unitary = false, bool skip_clone = false);

/// @brief Screens out terms in the expression \p expr that cannot contribute to
/// expectation value
/// @param expr input expression
/// @param skip_clone if true, will not clone the input expression
/// @return return screened expression
/// @code
/// // example usage:
/// auto expr1 = screen_zero_terms(expr); // screens for <0| expr |0>
/// auto expr2 = screen_zero_terms(P(2) * expr); // screens for <P(2)| expr |0>
/// @endcode
/// @pre This function expects \p input to be composed of mbpt::Operators
ExprPtr screen_zero_terms(ExprPtr expr, bool skip_clone = false);

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_UTILS_HPP
