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
/// - then, order by index slot types: spbundle, bra, ket, aux(matches
/// the order of entries in IndexSlotType).
/// - then, by space
///
/// N.B. spbundle is at the front to make CSV tensor layout "natural", with all
/// occupieds (in particular the spbundle-only occupieds) first. E.g. we want
/// g^{a<ij> b<ij>}_{k l} to be laid out with ij first, since in
/// any tensor product involving such indices will treat ij as batch indices,
/// hence it's best to make them as "slow" as possible (in row-major layout)
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
