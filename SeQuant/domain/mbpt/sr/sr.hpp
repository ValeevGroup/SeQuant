//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_SRCC_HPP
#define SEQUANT_SRCC_HPP

#include <limits>
#include <utility>

#include "SeQuant/core/expr_fwd.hpp"
#include "SeQuant/core/sequant.hpp"
#include "SeQuant/core/space.hpp"
#include "SeQuant/domain/mbpt/op.hpp"

namespace sequant {
namespace mbpt {
namespace sr {

/// single reference operator algebra can be screened by tracking total
/// _particle_ number and number of quasiparticles in occupied/unoccupied space
using qns_t = mbpt::ParticleNumberChange<2>;
using op_t = mbpt::Operator<qns_t>;

class make_op {
  std::size_t nbra_, nket_;
  OpType op_;

 public:
  make_op(std::size_t nbra, std::size_t nket, OpType op);

  ExprPtr operator()(IndexSpace::Type unocc, IndexSpace::Type occ) const;

  ExprPtr operator()() const;
};

make_op Op(OpType _Op, std::size_t Nbra,
           std::size_t Nket = std::numeric_limits<std::size_t>::max());

#include "sr_op.impl.hpp"

ExprPtr H1();

ExprPtr H2();

ExprPtr H0mp();
ExprPtr H1mp();

ExprPtr F();
ExprPtr W();

ExprPtr H();

/// computes the vacuum expectation value (VEV)

/// @param[in] expr input expression
/// @param[in] op_connections specifies the connectivity to be ensured
/// @param[in] use_top if true, topological equivalence will be utilized
/// @return the VEV
ExprPtr vac_av(ExprPtr expr,
               std::vector<std::pair<int, int>> op_connections = {},
               bool use_top = true);

// these produce operator-level expressions
namespace op {

ExprPtr H1();

ExprPtr H2();

ExprPtr H0mp();
ExprPtr H1mp();

ExprPtr F();
ExprPtr W();

ExprPtr H();

/// makes particle-conserving excitation operator of rank \p K
ExprPtr T_(std::size_t K);

/// makes sum of particle-conserving excitation operators of all ranks up to \p
/// K
ExprPtr T(std::size_t K);

/// makes particle-conserving deexcitation operator of rank \p K
ExprPtr Lambda_(std::size_t K);

/// makes sum of particle-conserving deexcitation operators of all ranks up to
/// \p
/// K
ExprPtr Lambda(std::size_t K);

/// makes deexcitation operator of rank \p K
ExprPtr A(std::size_t K);

/// @return true if \p op_or_op_product can produce determinant of excitation
/// rank \p k when applied to reference
bool contains_rank(const ExprPtr& op_or_op_product, const unsigned long k);

/// @return true if \p op_or_op_product can produce determinant of excitation
/// rank up to \p k when applied to reference
bool contains_up_to_rank(const ExprPtr& op_or_op_product,
                         const unsigned long k);

/// computes the vacuum expectation value (VEV)

/// @param[in] expr input expression
/// @param[in] op_connections list of pairs of operators to connect
/// @return the VEV
ExprPtr vac_av(
    ExprPtr expr,
    std::vector<std::pair<std::wstring, std::wstring>> op_connections = {
        {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"}});

}  // namespace op

}  // namespace sr
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_SRCC_HPP
