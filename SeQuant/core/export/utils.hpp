#ifndef SEQUANT_CORE_EXPORT_UTILS_HPP
#define SEQUANT_CORE_EXPORT_UTILS_HPP

#include <SeQuant/core/tensor.hpp>

namespace sequant {

/// Comparator that identifies Tensors only by their "block", which is defined
/// by its name, the amount of its indices as well as the space these indices
/// belong to. Note that it explicitly does not depend on the explicit index
/// labelling.
struct TensorBlockCompare {
  /// \returns Whether lhs compares less than rhs
  bool operator()(const Tensor &lhs, const Tensor &rhs) const;
};

}  // namespace sequant

#endif
