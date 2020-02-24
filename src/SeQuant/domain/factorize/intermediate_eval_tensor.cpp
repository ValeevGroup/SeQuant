#include "intermediate_eval_tensor.hpp"
#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

#include <memory>

namespace sequant::factorize {

bool IntermediateEvalTensor::is_leaf() const { return false; }

void IntermediateEvalTensor::set_left_ptr(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  left_tensor_ptr_ = tensor_ptr;
}

const EvalTensorPtr& IntermediateEvalTensor::get_left_tensor() const {
  return left_tensor_ptr_;
}

void IntermediateEvalTensor::set_right_ptr(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  right_tensor_ptr_ = tensor_ptr;
}

const EvalTensorPtr& IntermediateEvalTensor::get_right_tensor() const {
  return right_tensor_ptr_;
}

void IntermediateEvalTensor::set_operation(Operation op) { op_type_ = op; }

const Operation IntermediateEvalTensor::get_operation() const {
  return op_type_;
}

}  // namespace sequant::factorize
