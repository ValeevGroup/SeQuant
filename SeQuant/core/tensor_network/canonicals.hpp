//
// Created by Eduard Valeyev on 1/31/25.
//

#ifndef SEQUANT_TENSOR_NETWORK_CANONICALS_HPP
#define SEQUANT_TENSOR_NETWORK_CANONICALS_HPP

namespace sequant {

/// to produce the "canonical" layout of tensors/TNs this is used by
/// canonicalize_slots to sort BEFORE sorting by topological canonical
/// order produced by bliss::Graph:
/// - first order by # of protoindices
/// - then, order by index slot types: bra, ket, aux, spbundle (matches
/// the order of entries in IndexSlotType)
struct default_idxptr_slottype_lesscompare {
  template <typename IdxPtrPlusSlotType>
  bool operator()(const IdxPtrPlusSlotType& idxptr_slottype1,
                  const IdxPtrPlusSlotType& idxptr_slottype2) {
    const auto& [idxptr1, slottype1] = idxptr_slottype1;
    const auto& [idxptr2, slottype2] = idxptr_slottype2;
    if (idxptr1->proto_indices().size() != idxptr2->proto_indices().size())
      return idxptr1->proto_indices().size() < idxptr2->proto_indices().size();
    else {
      if (slottype1 != slottype2)
        return slottype1 < slottype2;
      else  // in same types of slots order by space
        return idxptr1->space() < idxptr2->space();
    }
  }
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_CANONICALS_HPP
