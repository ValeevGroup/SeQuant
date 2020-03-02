#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

#include <boost/any.hpp>
#include <memory>

namespace sequant::factorize {

// EvalTensor the base class
void EvalTensor::set_indices(const IndexLabelContainer& index_labels) {
  indices_ = index_labels;
}

const IndexLabelContainer& EvalTensor::get_indices() const { return indices_; }

void EvalTensor::set_hash_value(HashType hash_value) {
  hash_value_ = hash_value;
}

const HashType EvalTensor::get_hash_value() const { return hash_value_; }

void EvalTensor::set_ops_count(OpsCount count) { ops_count_ = count; }

const OpsCount EvalTensor::get_ops_count() const { return ops_count_; }

// EvalTensorIntermediate
bool EvalTensorIntermediate::is_leaf() const { return false; }

void EvalTensorIntermediate::set_left_ptr(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  left_tensor_ptr_ = tensor_ptr;
}

const EvalTensorPtr& EvalTensorIntermediate::get_left_tensor() const {
  return left_tensor_ptr_;
}

void EvalTensorIntermediate::set_right_ptr(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  right_tensor_ptr_ = tensor_ptr;
}

const EvalTensorPtr& EvalTensorIntermediate::get_right_tensor() const {
  return right_tensor_ptr_;
}

void EvalTensorIntermediate::set_operation(Operation op) { op_type_ = op; }

const Operation EvalTensorIntermediate::get_operation() const {
  return op_type_;
}

// EvalTensorLeaf
void EvalTensorLeaf::set_dtensor_ptr(
    const std::shared_ptr<boost::any>& dtensor_ptr) {
  dtensor_ptr_ = dtensor_ptr;
}

bool EvalTensorLeaf::is_leaf() const { return true; }

}  // namespace sequant::factorize
