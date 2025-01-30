#ifndef SEQUANT_TENSOR_NETWORK_SLOT_H
#define SEQUANT_TENSOR_NETWORK_SLOT_H

namespace sequant {

/// types of slots that can host an index

/// @sa TensorNetwork
/// @sa TensorNetworkV2
enum class IndexSlotType {
  /// occupying tensor vector index slot(s)
  TensorBraKet,
  /// occupying tensor aux index slot(s)
  TensorAux,
  /// only part of proto index bundles
  SPBundle
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_SLOT_H
