#ifndef SEQUANT_FACTORIZER_EVAL_TENSOR_FWD_HPP
#define SEQUANT_FACTORIZER_EVAL_TENSOR_FWD_HPP

#include <SeQuant/core/container.hpp>
#include <cstddef>
#include <memory>

namespace sequant::factorize {

using HashType = size_t;

using IndexLabelContainer = container::svector<std::wstring_view, 4>;

using OpsCount = unsigned long long;

class EvalTensor;

using EvalTensorPtr = std::shared_ptr<EvalTensor>;

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
  /// Represents an invalid evaluation.
  INVALID
};

}  // namespace sequant::factorize

#endif /* ifndef SEQUANT_FACTORIZER_EVAL_TENSOR_FWD_HPP */
