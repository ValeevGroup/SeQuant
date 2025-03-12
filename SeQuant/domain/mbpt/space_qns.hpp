//
// Created by Eduard Valeyev on 3/9/25.
//

#ifndef MPQC_DOMAIN_MBPT_SPACE_QNS_HPP
#define MPQC_DOMAIN_MBPT_SPACE_QNS_HPP

#include "SeQuant/core/fwd.hpp"

namespace sequant::mbpt {

// SeQuant/mbpt definitions of IndexSpace QNs

/// quantum numbers tags related to spin
/// \note spin quantum number takes 2 rightmost bits since there are 3 possible
/// states (any/no spin, spin-up, spin-down)
enum class Spin : bitset_t {
  alpha = 0b000001,
  beta = 0b000010,
  any = alpha | beta,  // using 2 bits so that overlap and union work expected
                       // (any & alpha = alpha, alpha | beta = any)
  // syntax sugar
  free = any,
  none = any,
  up = alpha,
  down = beta,
  mask = any
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

/// tags related to batching
/// \note BatchingQNS uses the 6th rightmost bit
enum class BatchingQNS : bitset_t {
  batch = 0b100000,  // for batching tensors
  mask = batch
};

}  // namespace sequant::mbpt

#endif  // MPQC_DOMAIN_MBPT_SPACE_QNS_HPP
