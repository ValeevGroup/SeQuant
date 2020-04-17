#ifndef SEQUANT_EVALUATE_EVAL_TENSOR_FWD_HPP
#define SEQUANT_EVALUATE_EVAL_TENSOR_FWD_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <memory>

namespace sequant::evaluate {

using HashType = size_t;

using ScalarType = double;

using IndexContainer = container::svector<Index, 4>;

using OpsCount = unsigned long long;

template <typename DataTensorType>
class EvalTree;

template <typename DataTensorType>
using EvalTreePtr = std::shared_ptr<EvalTree<DataTensorType>>;

enum class Operation {
  /// Represents the summation type binary evaluation
  SUM,
  /// Represents the contraction type binary evaluation
  PRODUCT,
  ///
  /// Represents the antisymmetrization type evaluation.
  ///
  /// An antisymmetrization type evaluation should be thought of as a binary
  /// operation between an antisymmetrization tensor (label "A" tensors in
  /// sequant) and the tensor to be antisymmetrized.
  ///
  ANTISYMMETRIZE,
  ///
  /// Represents the symmetrization type evaluation.
  ///
  /// An symmetrization type evaluation should be thought of as a binary
  /// operation between an symmetrization tensor (label "P" tensors in
  /// sequant) and the tensor to be symmetrized.
  SYMMETRIZE,
  /// Represents an invalid evaluation.
  INVALID
};

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_EVAL_TENSOR_FWD_HPP */
