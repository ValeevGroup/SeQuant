//
// Created by Eduard Valeyev on 3/6/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_LCAO_HPP
#define SEQUANT_DOMAIN_MBPT_LCAO_HPP

#include <SeQuant/core/fwd.hpp>

namespace sequant::mbpt {

/// quantum numbers tags related to spin
/// \note spin quantum number takes 2 rightmost bits
enum class LCAOQNS : bitset_t {
  lcao = 0b000000,
  ao = 0b000100,
  mask = lcao | ao
};

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_LCAO_HPP
