//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_DOMAIN_MBPT_MR_HPP
#define SEQUANT_DOMAIN_MBPT_MR_HPP

#include <limits>
#include <utility>

#include "SeQuant/core/expr_fwd.hpp"
#include "SeQuant/core/sequant.hpp"
#include "SeQuant/core/space.hpp"
#include "SeQuant/domain/mbpt/op.hpp"

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

ExprPtr H1();

ExprPtr H2();

ExprPtr F();

ExprPtr H();

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

ExprPtr H1();

ExprPtr H2();

}  // namespace op

}  // namespace mr

extern template class Operator<mr::qns_t, Statistics::FermiDirac>;
extern template class Operator<mr::qns_t, Statistics::BoseEinstein>;

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_MR_HPP
