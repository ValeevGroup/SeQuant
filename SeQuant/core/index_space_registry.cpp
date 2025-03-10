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

}  // namespace sequant
