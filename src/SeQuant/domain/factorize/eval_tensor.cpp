#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

namespace sequant::factorize {

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

}  // namespace sequant::factorize
