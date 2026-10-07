//
// Created by Eduard Valeyev on 3/9/25.
//

#include "SeQuant/core/index_basis_registry.hpp"

namespace sequant {

void IndexBasisRegistry::physical_particle_attribute_mask(bitset_t m) {
  physical_particle_attribute_mask_ = m;
}

bitset_t IndexBasisRegistry::physical_particle_attribute_mask() const {
  return physical_particle_attribute_mask_;
}

IndexSpace::QuantumNumbers IndexBasisRegistry::physical_particle_attributes(
    IndexSpace::QuantumNumbers qn) const {
  return to_bitset(qn) & physical_particle_attribute_mask_;
}

IndexSpace::QuantumNumbers IndexBasisRegistry::other_attributes(
    IndexSpace::QuantumNumbers qn) const {
  return to_bitset(qn) & ~physical_particle_attribute_mask_;
}

}  // namespace sequant
