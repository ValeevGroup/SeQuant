//
// Created by Eduard Valeyev on 3/6/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_LCAO_HPP
#define SEQUANT_DOMAIN_MBPT_LCAO_HPP

#include <SeQuant/core/fwd.hpp>

namespace sequant::mbpt {

/// quantum numbers tags related to LCAO basis traits
/// \note sLCAO basis traits use 3rd and 4th rightmost bits
enum class LCAOQNS : bitset_t {
  ao = 0b000100,
  pao = 0b001000,  // projected AO space denotes unoccupied spaces made from AO
                   // basis and orthogonal to occupied orbitals
  mask = ao | pao
};

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_LCAO_HPP
