//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_DOMAIN_MBPT_MR_HPP
#define SEQUANT_DOMAIN_MBPT_MR_HPP

#include <SeQuant/domain/mbpt/op.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/interval.hpp>
#include <SeQuant/core/latex.hpp>
#include <SeQuant/core/space.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {
namespace mbpt {
namespace mr {

struct qns_tag;

// clang-format off
/// multireference operator algebra can be screened by tracking the number of creators and annihilators in the (active) occupied, active, and (active) unoccupied space
/// the order of of elements is {# of occupied creators, # of occupied annihilators, # of active creators, # of active annihilators, # of unoccupied creators, # of unoccupied annihilators}
/// \note use signed integer, although could use unsigned in this case, so that can represent quantum numbers and their changes by the same type
// clang-format on
using qns_t = mbpt::QuantumNumberChange<6, qns_tag, std::int64_t>;
/// changes in quantum number represented by quantum numbers themselves
using qnc_t = qns_t;
using op_t = mbpt::Operator<qnc_t>;

// clang-format off
/// @return the number of creators in \p qns acting on space \p s
/// @pre `(s.type()==IndexSpace::Type::active_occupied || s.type()==IndexSpace::Type::active_unoccupied)&&s.qns()==IndexSpace::null_qns`
// clang-format on
qninterval_t ncre(qns_t qns, const IndexSpace& s);

// clang-format off
/// @return the number of creators in \p qns acting on space \p s
/// @pre `s==IndexSpace::Type::active_occupied || s==IndexSpace::Type::active_unoccupied`
// clang-format on
qninterval_t ncre(qns_t qns, const IndexSpace::Type& s);

// clang-format off
/// @return the number of creators in \p qns acting on the occupied space
// clang-format on
qninterval_t ncre_occ(qns_t qns);

// clang-format off
/// @return the number of creators in \p qns acting on the active space
// clang-format on
qninterval_t ncre_act(qns_t qns);

// clang-format off
/// @return the number of creators in \p qns acting on the unoccupied space
// clang-format on
qninterval_t ncre_uocc(qns_t qns);

// clang-format off
/// @return the total number of creators in \p qns
// clang-format on
qninterval_t ncre(qns_t qns);

// clang-format off
/// @return the number of annihilators in \p qns acting on space \p s
/// @pre `(s.type()==IndexSpace::Type::active_occupied || s.type()==IndexSpace::Type::active_unoccupied)&&s.qns()==IndexSpace::null_qns`
// clang-format on
qninterval_t nann(qns_t qns, const IndexSpace& s);

// clang-format off
/// @return the number of annihilators in \p qns acting on space \p s
/// @pre `s==IndexSpace::Type::active_occupied || s==IndexSpace::Type::active_unoccupied`
// clang-format on
qninterval_t nann(qns_t qns, const IndexSpace::Type& s);

// clang-format off
/// @return the number of annihilators in \p qns acting on the occupied space
// clang-format on
qninterval_t nann_occ(qns_t qns);

// clang-format off
/// @return the number of annihilators in \p qns acting on the active space
// clang-format on
qninterval_t nann_act(qns_t qns);

// clang-format off
/// @return the number of annihilators in \p qns acting on the unoccupied space
// clang-format on
qninterval_t nann_uocc(qns_t qns);

// clang-format off
/// @return the total number of annihilators in \p qns
// clang-format on
qninterval_t nann(qns_t qns);

/// combines 2 sets of quantum numbers using Wick's theorem
qns_t combine(qns_t, qns_t);

}  // namespace mr
}  // namespace mbpt

/// @param qns the quantum numbers to adjoint
/// @return the adjoint of \p qns
mbpt::mr::qns_t adjoint(mbpt::mr::qns_t);

namespace mbpt {
namespace mr {

// clang-format off
/// @brief makes a tensor-level fermionic many-body operator for use in single-reference methods

/// A many-body operator has the following generic form:
/// \f$ \frac{1}{P} T_{b_1 b_2 \dots b_B}^{k_1 k_2 \dots k_K} A^{b_1 b_2 \dots b_B}_{k_1 k_2 \dots k_K} \f$
/// where \f$ \{B,K\} \f$ are number of bra/ket indices of \f$ T \f$ or, equivalently, the number of creators/annihilators
/// of normal-ordered (w.r.t. the default, not necessarily Fermi vacuum) operator \f$ A \f$.
/// Indices \f$ \{ b_i \} \f$ / \f$ \{ k_i \} \f$ are (active) unoccupied/occupied for (pure) excitation operators,
/// are occupied/unoccupied for deexcitation operators; for general operators complete basis indices are assumed by default,
/// unless overridden by user manually. \f$ P \f$ is the "normalization" factor and depends on the vacuum used to define \f$ A \f$,
/// and indices \f$ \{ b_i \} \f$ / \f$ \{ k_i \} \f$.
/// @note The choice of unoccupied indices/spaces can be controlled by the default Formalism:
/// - if `get_default_formalism().sum_over_uocc() == SumOverUocc::Complete` IndexSpace::complete_unoccupied will be used instead of IndexSpace::active_unoccupied
/// - if `get_default_formalism().csv() == CSVFormalism::CSV` will use cluster-specific (e.g., PNO) unoccupied indices
/// @warning Tensor \f$ T \f$ will be antisymmetrized if `get_default_context().spbasis() == SPBasis::spinorbital`, else it will be particle-symmetric; the latter is only valid if # of bra and ket indices coincide.
// clang-format on
class OpMaker : public mbpt::OpMaker<Statistics::FermiDirac> {
 public:
  using base_type = mbpt::OpMaker<Statistics::FermiDirac>;

  using base_type::base_type;

  // clang-format off
  /// @param[in] op the operator type:
  /// - if @p op is a (pure) excitation operator bra/ket indices
  ///   will be IndexSpace::active_unoccupied/IndexSpace::active_occupied,
  /// - for (pure) deexcitation @p op bra/ket will be IndexSpace::active_occupied/IndexSpace::active_unoccupied
  /// - for general @p op bra/ket will be IndexSpace::complete
  /// @param[in] nbra number of bra indices/creators
  /// @param[in] nket number of ket indices/annihilators; if not specified, will be set to @p nbra
  // clang-format on
  OpMaker(OpType op, std::size_t nbra,
          std::size_t nket = std::numeric_limits<std::size_t>::max());

  using base_type::operator();
};

#include <SeQuant/domain/mbpt/mr/op.impl.hpp>

/// @name tensor-level operators
/// @{

// clang-format off
/// @brief `k`-body contribution to the "generic" Hamiltonian (in normal order relative to the default vacuum)
/// @param[in] k the rank of the particle interactions; only `k<=2` is
/// supported
// clang-format on
ExprPtr H_(std::size_t k);

/// @brief total Hamiltonian including up to `k`-body interactions
/// @param[in] k the maximum rank of the particle interactions; only `k<=2` is
/// supported
ExprPtr H(std::size_t k = 2);

/// @brief Fock operator
/// @param use_f_tensor if true, will use Fock tensor, else will use tensors
/// used to define `H_(1)` and `H_(2)`
ExprPtr F(bool use_f_tensor = true);

/// @}

/// computes the vacuum expectation value (VEV)

/// @param[in] expr input expression
/// @param[in] nop_connections specifies the pairs of normal operators to be
/// connected
/// @param[in] use_top if true, topological equivalence will be utilized
/// @return the VEV
ExprPtr vac_av(ExprPtr expr,
               std::vector<std::pair<int, int>> nop_connections = {},
               bool use_top = true);

// these produce operator-level expressions
namespace op {

// clang-format off
/// @brief `k`-body contribution to the "generic" Hamiltonian (in normal order relative to the default vacuum)
/// @param[in] k the rank of the particle interactions; only `k<=2` is
/// supported
// clang-format on
ExprPtr H_(std::size_t k);

/// @brief total Hamiltonian including up to `k`-body interactions
/// @param[in] k the maximum rank of the particle interactions; only `k<=2` is
/// supported
ExprPtr H(std::size_t k = 2);

/// makes particle-conserving excitation operator of rank \p K
ExprPtr T_(std::size_t K);

/// makes sum of particle-conserving excitation operators of all ranks up to \p
/// K
ExprPtr T(std::size_t K);

/// makes particle-conserving deexcitation operator of rank \p K
ExprPtr Λ_(std::size_t K);

/// makes sum of particle-conserving deexcitation operators of all ranks up to
/// \p K
ExprPtr Λ(std::size_t K);

/// @}

/// @return true if \p op_or_op_product can change quantum numbers from \p source_qns to \p target_qns
bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t target_qns,
                    const qns_t source_qns = {});

#include <SeQuant/domain/mbpt/vac_av.hpp>

}  // namespace op

}  // namespace mr

extern template class Operator<mr::qns_t, Statistics::FermiDirac>;
extern template class Operator<mr::qns_t, Statistics::BoseEinstein>;

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_MR_HPP
