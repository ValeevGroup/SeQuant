#ifndef SEQUANT_CORE_EXPORT_COMPUTE_SELECTION_HPP
#define SEQUANT_CORE_EXPORT_COMPUTE_SELECTION_HPP

#include <SeQuant/api.hpp>

namespace sequant {

/// In a full binary tree, the ComputeSelection determines the left and/or right
/// subtree shall explicitly be considered to contribute to the computation of
/// the current node's result. If neither the left nor right subtree contribute,
/// that is ComputeSelection::None is used, the given node exists for tree
/// connectivity purposes only and shall be skipped during code generation.
enum class ComputeSelection {
  None = 0,
  Left = 0b01,
  Right = 0b10,
  Both = 0b11,
};

SEQUANT_EXPORT_API ComputeSelection operator|(ComputeSelection lhs,
                                              ComputeSelection rhs);
SEQUANT_EXPORT_API ComputeSelection operator&(ComputeSelection lhs,
                                              ComputeSelection rhs);
SEQUANT_EXPORT_API ComputeSelection& operator|=(ComputeSelection& lhs,
                                                ComputeSelection rhs);
SEQUANT_EXPORT_API ComputeSelection& operator&=(ComputeSelection& lhs,
                                                ComputeSelection rhs);

SEQUANT_EXPORT_API ComputeSelection operator~(ComputeSelection selection);

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_COMPUTE_SELECTION_HPP
