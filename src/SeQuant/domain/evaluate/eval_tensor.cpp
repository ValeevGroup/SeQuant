#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

#include <SeQuant/core/expr.hpp>

#include <boost/any.hpp>
#include <memory>

namespace sequant::evaluate {

// EvalTensor the base class
void EvalTensor::set_indices(const IndexLabelContainer& index_labels) {
  indices_ = index_labels;
}

const IndexLabelContainer& EvalTensor::get_indices() const { return indices_; }

void EvalTensor::set_hash_value(HashType hash_value) {
  hash_value_ = hash_value;
}

HashType EvalTensor::get_hash_value() const { return hash_value_; }

void EvalTensor::set_ops_count(OpsCount count) { ops_count_ = count; }

OpsCount EvalTensor::get_ops_count() const { return ops_count_; }

ScalarType EvalTensor::get_scalar() const { return scalar_; }

void EvalTensor::set_scalar(ScalarType scale) { scalar_ = scale; }

// EvalTensorIntermediate
const std::wstring EvalTensorIntermediate::SUM_COLOR = L"green";

const std::wstring EvalTensorIntermediate::PROD_COLOR = L"red";

const std::wstring EvalTensorIntermediate::ANTISYMM_COLOR = L"gray";

const std::wstring EvalTensorIntermediate::SYMM_COLOR = L"turquoise";

bool EvalTensorIntermediate::is_leaf() const { return false; }

void EvalTensorIntermediate::set_left_tensor(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  left_tensor_ = tensor_ptr;
}

const EvalTensorPtr& EvalTensorIntermediate::get_left_tensor() const {
  return left_tensor_;
}

void EvalTensorIntermediate::set_right_tensor(
    const std::shared_ptr<EvalTensor>& tensor_ptr) {
  right_tensor_ = tensor_ptr;
}

const EvalTensorPtr& EvalTensorIntermediate::get_right_tensor() const {
  return right_tensor_;
}

void EvalTensorIntermediate::set_operation(Operation op) { operation_ = op; }

Operation EvalTensorIntermediate::get_operation() const { return operation_; }

std::wstring EvalTensorIntermediate::to_digraph() const {
  std::wstring bra = L"";
  std::wstring ket = L"";
  for (auto i = 0; i < get_indices().size() / 2; ++i) {
    bra += L"{" + std::wstring(get_indices().at(i)) + L"}";
  }
  for (auto i = get_indices().size() / 2; i < get_indices().size(); ++i) {
    ket += L"{" + std::wstring(get_indices().at(i)) + L"}";
  }
  auto node_color = get_operation() == Operation::SUM
                        ? SUM_COLOR
                        : get_operation() == Operation::PRODUCT
                              ? PROD_COLOR
                              : get_operation() == Operation::ANTISYMMETRIZE
                                    ? ANTISYMM_COLOR
                                    : SYMM_COLOR;
  auto this_node = L"node" + std::to_wstring(get_hash_value()) +
                   L" [label=\"" + L"q^{" + ket + L"}_{" + bra + L"}\"" +
                   L", color=\"" + node_color + L"\", style = \"filled\"" + L", fontcolor=\"white\"];\n";

  auto left_node = get_left_tensor()->to_digraph();

  auto right_node = get_right_tensor()->to_digraph();
  // drawing the edges
  auto get_node_name = [](const std::wstring& str) {
    return str.substr(0, str.find(L" "));
  };
  auto parent_node_name = get_node_name(this_node);
  auto left_node_name = get_node_name(left_node);
  auto right_node_name = get_node_name(right_node);
  auto result = this_node + left_node + right_node;
  result += parent_node_name + L" -> {" + left_node_name + L"};\n";
  result += parent_node_name + L" -> {" + right_node_name + L"};\n";
  return result;
}

/// Visit the tree by pre-order traversal.
void EvalTensorIntermediate::visit(
    const std::function<void(const EvalTensor&)>& visitor) const {
  visitor(*this);
  get_left_tensor()->visit(visitor);
  get_right_tensor()->visit(visitor);
}

// EvalTensorLeaf

EvalTensorLeaf::EvalTensorLeaf(const ExprPtr& expr) : expr_{expr} {}

void EvalTensorLeaf::set_data_tensor(
    const std::shared_ptr<boost::any>& dtensor_ptr) {
  data_tensor_ = dtensor_ptr;
}

bool EvalTensorLeaf::is_leaf() const { return true; }

const ExprPtr& EvalTensorLeaf::get_expr() const { return expr_; }

std::wstring EvalTensorLeaf::to_digraph() const {
  return L"node" + std::to_wstring(get_hash_value()) + L" [label = \"" +
         get_expr()->to_latex() + L"\"];\n";
}

/// Visit the tree by pre-order traversal.
void EvalTensorLeaf::visit(
    const std::function<void(const EvalTensor&)>& visitor) const {
  visitor(*this);
}

}  // namespace sequant::evaluate
