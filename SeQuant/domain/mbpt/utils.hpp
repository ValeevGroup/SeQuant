//
// Created by Ajay Melekamburath on 9/26/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_UTILS_HPP
#define SEQUANT_DOMAIN_MBPT_UTILS_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/aggregate.hpp>

#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant::mbpt {

/// Options for Lie Similarity Transformation. @see lst
struct LSTOptions {
  /// Present only to force designated initialization, e.g.
  /// `LSTOptions{.unitary = true}` rather than `LSTOptions{true, false}`.
  /// Renaming/reordering the fields below would otherwise silently reinterpret
  /// any positional-init caller instead of failing to compile.
  [[no_unique_address]] sequant::detail::designated_only designated_only_ = {};
  /// If true, uses unitary generator
  bool unitary = false;
  /// If true, uses connected products [A,B] = (AB)_c; otherwise uses explicit
  /// commutators [A,B] = AB - BA. The connected-product form is only
  /// equivalent if the caller supplies operator connectivity downstream, hence
  /// the default is the explicit form.
  bool use_connected_form = false;
  /// If true, will not clone the input expression
  bool skip_clone = false;
};

// clang-format off
/// @brief Computes the Lie similarity transformation, e^(-B) A e^B, using its Campbell expansion (DOI 10.1112/plms/s1-28.1.381) as a series of nested commutators:
/// `e^(-B) A e^B = A + [A,B] + (1/2!)[[A,B],B] + (1/3!)[[[A,B],B],B] + ...`
///
/// @param A Expression A (e.g., the Hamiltonian)
/// @param B Expression B (e.g., the cluster operator T)
/// @param commutator_rank The maximum order of nested commutators to retain (e.g. 4 for traditional coupled-cluster)
/// @param options Options controlling the transformation behavior, see LSTOptions.
/// @sa LSTOptions
/// Notes:
/// - If \p options.unitary is true, the ansatz uses B - B^+ instead of B.
/// - By default commutators are computed explicitly: [A,B] = AB - BA
/// - If \p options.use_connected_form is true, commutators are computed via connected products: [A,B] = (AB)_c ; this is only valid if the caller connects the operators (e.g. by passing OpConnections to vac_av/ref_av), so set it only if you do.
/// @pre This function expects \p A and \p B to be composed of mbpt::Operators
// clang-format on
ExprPtr lst(ExprPtr A, ExprPtr B, size_t commutator_rank,
            LSTOptions options = {});

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
