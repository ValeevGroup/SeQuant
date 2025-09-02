//
// Created by Eduard Valeyev on 3/9/25.
//

#include "SeQuant/core/index_space_registry.hpp"

namespace sequant {

void IndexSpaceRegistry::physical_particle_attribute_mask(bitset_t m) {
  physical_particle_attribute_mask_ = m;
}

bitset_t IndexSpaceRegistry::physical_particle_attribute_mask() const {
  return physical_particle_attribute_mask_;
}

IndexSpace::QuantumNumbers IndexSpaceRegistry::physical_particle_attributes(
    IndexSpace::QuantumNumbers qn) const {
  return static_cast<bitset_t>(qn) & physical_particle_attribute_mask_;
}

IndexSpace::QuantumNumbers IndexSpaceRegistry::other_attributes(
    IndexSpace::QuantumNumbers qn) const {
  return static_cast<bitset_t>(qn) & ~physical_particle_attribute_mask_;
}

IndexSpaceRegistry IndexSpaceRegistry::clone() const {
  IndexSpaceRegistry result(*this);
  result.spaces_ =
      std::make_shared<container::set<IndexSpace, IndexSpace::KeyCompare>>(
          *spaces_);
  return result;
}

}  // namespace sequant
