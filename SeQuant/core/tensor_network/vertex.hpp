#ifndef SEQUANT_TENSOR_NETWORK_VERTEX_H
#define SEQUANT_TENSOR_NETWORK_VERTEX_H

namespace sequant {

/// types of vertices created on a colored graph representation of a tensor
/// network

/// @sa TensorNetwork
/// @sa TensorNetworkV2
enum class VertexType {
  Index,
  SPBundle,
  TensorBra,
  TensorKet,
  TensorAux,
  TensorCore,
  Particle
};

}  // namespace sequant

#endif  // SEQUANT_TENSOR_NETWORK_VERTEX_H
