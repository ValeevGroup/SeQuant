//
// Created by Eduard Valeyev on 7/18/23.
//

#ifndef SEQUANT_CORE_FWD_HPP
#define SEQUANT_CORE_FWD_HPP

/// @brief the main namespace of the SeQuant framework
namespace sequant {}

#include <SeQuant/core/expr_fwd.hpp>

namespace sequant {
namespace bitset {
using type = int32_t;
constexpr type reserved = 0x80000000;
constexpr type null = 0x00000000;
}  // namespace bitset

using bitset_t = bitset::type;
}  // namespace sequant

#endif  // SEQUANT_CORE_FWD_HPP
