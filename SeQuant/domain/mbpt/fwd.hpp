//
// Created by Eduard Valeyev on 7/18/23.
//

#ifndef SEQUANT_DOMAIN_MBPT_FWD_HPP
#define SEQUANT_DOMAIN_MBPT_FWD_HPP

#include <SeQuant/core/fwd.hpp>

#include <cstdint>

namespace sequant {

/// @brief the main namespace of the Many-Body Perturbation Theory (MBPT)
/// components of SeQuant
namespace mbpt {

/// @brief the namespace containing MBPT operators
inline namespace op {

/// @brief the namespace containing tensor form of MBPT operators
namespace tensor {}

}  // namespace op
}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_FWD_HPP
