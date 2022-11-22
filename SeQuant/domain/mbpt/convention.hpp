//
// Created by Eduard Valeyev on 2019-04-01.
//

#ifndef SEQUANT_CONVENTION_HPP
#define SEQUANT_CONVENTION_HPP

namespace sequant {
namespace mbpt {

enum class Convention { QCiFS };

/// @brief Loads defaults for Convention @c conv

/// This registers IndexSpace objects standard for the chosen convention,
/// updates default context's IndexRegistry object, and
/// updates TensorCanonicalizer cardinal labels
/// @warning should be only called once
void set_default_convention(Convention conv = Convention::QCiFS);

}  // namespace mbpt
}  // namespace sequant

#endif  // SEQUANT_CONVENTION_HPP
