#ifndef SEQUANT_TENSOR_NETWORK_VERTEX_H
#define SEQUANT_TENSOR_NETWORK_VERTEX_H

namespace sequant {

/// types of vertices created on a colored graph representation of a tensor
/// network

/// @sa TensorNetwork
/// @sa TensorNetworkV2
/// @sa TensorNetworkV3
enum class VertexType {
  /// represent a Index object
  Index,
  /// represent a bundle of indices that serve as protoindices to Index objects
  IndexBundle,
  /// legacy alias for IndexBundle
  SPBundle = IndexBundle,
  /// represents a bra slot for an Index
  TensorBra,
  /// represents a ket slot for an Index
  TensorKet,
  /// represents an aux slot for an Index
  TensorAux,
  /// represents a group of bra slots (TensorBra)
  TensorBraBundle,
  /// represents a group of ket slots (TensorKet)
  TensorKetBundle,
  /// represents a group of aux slots (TensorAux)
  TensorAuxBundle,
  /// represents a slot for tensor label
  TensorCore,
  /// connects bra slots (or their bundle) to the matching ket slots (or their
  /// bundle):
  /// - for symmetric/antisymmetric case (i.e., tensor is
  /// symmetric/antisymmetric w.r.t. permutation of slots within bra or ket)
  /// bra and ket bundles connect to same braket vertex;
  /// - for particle-symmetric case (i.e., tensor is symmetric w.r.t.
  /// permutation of columns) there is one per column, all with the same color
  /// - for an asymmetric case there is also one per column, all with different
  /// colors
  TensorBraKet
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_VERTEX_H
