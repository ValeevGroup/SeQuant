//
// Created by Eduard Valeyev on 7/18/23.
//

#ifndef SEQUANT_DOMAIN_MBPT_FWD_HPP
#define SEQUANT_DOMAIN_MBPT_FWD_HPP

#include "SeQuant/core/fwd.hpp"

#include <cstdint>

namespace sequant {

/// @brief the main namespace of the Many-Body Perturbation Theory (MBPT)
/// components of SeQuant
namespace mbpt {

/// @brief the namespace of the MBPT formalisms with respect to the single
/// determinant reference (SR) vacuum
namespace sr {}

/// @brief the namespace of the MBPT formalisms with respect to the multi
/// determinant reference (MR) wave function
namespace mr {}

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_DOMAIN_MBPT_FWD_HPP
