//
// Created by Eduard Valeyev on 3/9/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_SPACE_QNS_HPP
#define SEQUANT_DOMAIN_MBPT_SPACE_QNS_HPP

#include "SeQuant/core/fwd.hpp"

namespace sequant::mbpt {

// SeQuant/mbpt definitions of IndexSpace QNs

/// quantum numbers tags related to spin
/// \note spin quantum number takes 2 rightmost bits since there are 3 possible
/// states (any/no spin, spin-up, spin-down)
enum class Spin : bitset_t {
  alpha = 0b000001,  //!< spin-up (chemistry)
  beta = 0b000010,   //!< spin-down (chemistry)
  /// arbitrary spin-state represented by 2 bits
  /// so that overlap and union work as expected
  /// (`any & alpha = alpha`, `alpha | beta = any`, etc.)
  any = alpha | beta,
  // syntactic sugar
  free = any,   //!< spin-free means spin does not matter
  none = any,   //!< same as spin-free
  up = alpha,   //!< spin-up (physics)
  down = beta,  //!< spin-down (physics)
  mask = any,   //!< spin bit mask
  // only used by legacy code
  null = 0b000000
};

/// quantum numbers tags related to LCAO basis traits
/// \note sLCAO basis traits use 3rd and 4th rightmost bits
enum class LCAOQNS : bitset_t {
  ao = 0b000100,
  pao = 0b001000,  // projected AO space denotes unoccupied spaces made from AO
                   // basis and orthogonal to occupied orbitals
  mask = ao | pao
};

/// quantum numbers tags related to tensor factorization basis traits
/// \note TensorFactorization basis traits use 5th rightmost bits
enum class TensorFactorizationQNS : bitset_t { df = 0b010000, mask = df };

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_SPACE_QNS_HPP
