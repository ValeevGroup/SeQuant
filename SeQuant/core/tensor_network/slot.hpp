#ifndef SEQUANT_TENSOR_NETWORK_SLOT_H
#define SEQUANT_TENSOR_NETWORK_SLOT_H

namespace sequant {

/// types of index slots that a tensor can have

/// @note the order orchestrated to produce intuitively "canonical" layout of
/// named indices of tensor networks
/// @sa TensorNetwork
/// @sa TensorNetworkV2
enum class TensorIndexSlotType {
  /// bra tensor vector index slot
  Bra = 0,
  /// ket tensor vector index slot
  Ket = 1,
  /// tensor aux index slot
  Aux = 2,

  Invalid
};

/// types of slots that can host an index

/// @note the order orchestrated to produce intuitively "canonical" layout of
/// named indices of tensor networks
/// @sa TensorNetwork
/// @sa TensorNetworkV2
enum class IndexSlotType {
  /// in bra tensor vector index slot
  TensorBra,
  /// in ket tensor vector index slot
  TensorKet,
  /// connecting two tensor vector index slots
  TensorBraKet,
  /// in tensor aux index slot
  TensorAux,
  /// connecting two tensor aux index slots
  TensorAuxAux,
  /// only part of proto index bundles
  SPBundle,

  Invalid
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_SLOT_H
