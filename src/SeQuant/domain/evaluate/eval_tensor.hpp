#ifndef SEQUANT_EVALUATE_EVAL_TENSOR_HPP
#define SEQUANT_EVALUATE_EVAL_TENSOR_HPP

#include "eval_tensor_fwd.hpp"

#include <functional>
#include <memory>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant::evaluate {

namespace detail {
struct Perm {
  std::vector<size_t> sequence;
  int phase{1};
};

std::vector<Perm> perm_calc(std::vector<size_t> to_perm, size_t size,
                            size_t cswap = 0, size_t begin = 0) {
  if (begin + 1 == size) {
    // even permutations added. odds subtracted
    return std::vector<Perm>{Perm{to_perm, (cswap % 2 == 0) ? 1 : -1}};
  }

  auto result = std::vector<Perm>{};
  for (auto i = begin; i < size; ++i) {
    std::swap(to_perm[begin], to_perm[i]);
    auto more_result =
        perm_calc(to_perm, size, (begin == i) ? cswap : cswap + 1, begin + 1);
    for (auto& p : more_result) {
      result.push_back(p);
    }
    std::swap(to_perm[begin], to_perm[i]);
  }
  return result;
}
}  // namespace detail

///
/// @brief Representation of binary evaluation of sequant::Expr objects. A
/// context is needed for evaluation.
///
/// The atomic SeQuant Expr(s) can be thought of either as a representation of a
/// single data-tensor (eg. btas::Tensor<double>, TA::TArrayD), or binary
/// evaluations of such data-tensors. An evaluation could be summing two
/// data-tensors or taking a product of them. A result of a binary evaluation
/// can be an input for another evaluation, i.e. evaluations can be nested.
///
/// This class represents such evaluations. @note No evaluations on data-tensors
/// are actually performed. This is a tree data-structure only.
///
/// @author Bimal Gaudel
/// @date Feb 2020
///
template <typename DataTensorType>
class EvalTensor {
 private:
  /// The index labels of the tensor's bra and ket in that order.
  IndexLabelContainer indices_;

  /// A unique identifier of this evaluation.
  HashType hash_value_{0};

  /// The operations count that resulted during this evaluation.
  OpsCount ops_count_{0};

  /// The scalar to multiply this evaluation with.
  ScalarType scalar_{1};

 public:
  /// Setter method of the indices_ field.
  void set_indices(const IndexLabelContainer& index_labels) {
    indices_ = index_labels;
  }

  /// Getter method of the indices_ field.
  const IndexLabelContainer& get_indices() const { return indices_; }

  /// Setter method of the hash_value_ field.
  void set_hash_value(HashType hash_value) { hash_value_ = hash_value; }

  /// Getter method of the hash_value_ field.
  HashType get_hash_value() const { return hash_value_; }

  /// Setter method of the ops_count_ field.
  void set_ops_count(OpsCount count) { ops_count_ = count; }

  /// Getter method of the ops_count_ field.
  OpsCount get_ops_count() const { return ops_count_; }

  /// Getter for the scalar field.
  ScalarType get_scalar() const { return scalar_; }

  /// Setter for the scalar field.
  void set_scalar(ScalarType scale) { scalar_ = scale; }

  /// Check if this is the end of the evaluation.
  virtual bool is_leaf() const = 0;

  /// Get a directed graph definitions and paths in dot language.
  virtual std::wstring to_digraph() const = 0;

  /// Visit the tree by pre-order traversal.
  virtual void visit(
      const std::function<void(const EvalTensor&)>& visitor) const = 0;

  /// Evaluate the EvalTensor.
  virtual DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const = 0;
};

///
/// Non-leaf type EvalTensor node in a binary evaluation of the
/// sequant::Expr at some context.
///
template <typename DataTensorType>
class EvalTensorIntermediate : public EvalTensor<DataTensorType> {
 private:
  /// Evaluation node to the left.
  std::shared_ptr<EvalTensor<DataTensorType>> left_tensor_{nullptr};

  /// Evaluation node to the right
  std::shared_ptr<EvalTensor<DataTensorType>> right_tensor_{nullptr};

  /// The type of the binary evaluation.
  Operation operation_{Operation::INVALID};

 public:
  /// @return false
  bool is_leaf() const override { return false; }

  /// Setter method for the left evaluation node.
  void set_left_tensor(
      const std::shared_ptr<EvalTensor<DataTensorType>>& left_eval_tensor) {
    left_tensor_ = left_eval_tensor;
  }

  /// Getter method for the left evaluation node.
  const EvalTensorPtr<DataTensorType>& get_left_tensor() const {
    return left_tensor_;
  }

  /// Setter method for the right evaluation node.
  void set_right_tensor(
      const std::shared_ptr<EvalTensor<DataTensorType>>& right_eval_tensor) {
    right_tensor_ = right_eval_tensor;
  }

  /// Getter method for the left evaluation node.
  const EvalTensorPtr<DataTensorType>& get_right_tensor() const {
    return right_tensor_;
  }

  /// Setter method for the operation type.
  /// @param op An Operation.
  void set_operation(Operation op) { operation_ = op; }

  /// Getter method for the operation type.
  /// @return Operation for this intermediate tensor.
  Operation get_operation() const { return operation_; }

  /// Get a directed graph definitions and paths in dot language.
  std::wstring to_digraph() const override {
    std::wstring bra = L"";
    std::wstring ket = L"";
    for (auto i = 0; i < this->get_indices().size() / 2; ++i) {
      bra += L"{" + std::wstring(this->get_indices().at(i)) + L"}";
    }
    for (auto i = this->get_indices().size() / 2;
         i < this->get_indices().size(); ++i) {
      ket += L"{" + std::wstring(this->get_indices().at(i)) + L"}";
    }
    auto node_color =
        get_operation() == Operation::SUM
            ? L"green"
            : get_operation() == Operation::PRODUCT
                  ? L"red"
                  : get_operation() == Operation::ANTISYMMETRIZE
                        ? L"gray"
                        : L"turquoise";  // turquoise for symmetrization
    std::string scalar = std::to_string(this->get_scalar());
    scalar = scalar.substr(4);
    auto this_node =
        L"node" + std::to_wstring(this->get_hash_value()) + L" [label=\"" +
        ((this->get_scalar() != 1.)
             ? std::wstring(scalar.begin(), scalar.end()) + L"\\times "
             : L"");
    this_node += L"q^{" + ket + L"}_{" + bra + L"}\"" + L", color=\"" +
                 node_color + L"\", style = \"filled\"" +
                 L", fontcolor=\"white\"];\n";

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
  void visit(const std::function<void(const EvalTensor<DataTensorType>&)>&
                 visitor = {}) const override {
    /// Visit the tree by pre-order traversal.
    visitor(*this);
    get_left_tensor()->visit(visitor);
    get_right_tensor()->visit(visitor);
  }

  /// Evaluate the EvalTensor.
  DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const override {
    auto TA_annotation = [this](decltype(this->get_indices())& indices) {
      std::string annot = "";
      for (const auto& idx : indices)
        annot += std::string(idx.begin(), idx.end()) + ", ";

      annot.erase(annot.size() - 2);  // remove trailing L", "
      return annot;
    };

    auto left_annot = TA_annotation(get_left_tensor()->get_indices());
    auto right_annot = TA_annotation(get_right_tensor()->get_indices());
    auto this_annot = TA_annotation(this->get_indices());

    if (get_operation() == Operation::SUM) {
      DataTensorType result;
      result(this_annot) =
          get_left_tensor()->get_scalar() *
              get_left_tensor()->evaluate(context)(left_annot) +
          get_right_tensor()->get_scalar() *
              get_right_tensor()->evaluate(context)(right_annot);
      return result;
    } else if (get_operation() == Operation::PRODUCT) {
      DataTensorType result;
      result(this_annot) = get_left_tensor()->get_scalar() *
                           get_left_tensor()->evaluate(context)(left_annot) *
                           get_right_tensor()->get_scalar() *
                           get_right_tensor()->evaluate(context)(right_annot);
      return result;

    } else if (get_operation() == Operation::ANTISYMMETRIZE) {
      auto right_eval = get_right_tensor()->evaluate(context);

      size_t rank = right_eval.trange().rank();
      if (rank == 2)
        return right_eval;
      else if (rank % 2 != 0)
        throw std::domain_error("Can't handle odd-ordered tensors yet!");

      std::vector<size_t> to_perm;
      for (auto i = 0; i < (size_t)rank / 2; ++i) to_perm.push_back(i);
      auto vp = detail::perm_calc(to_perm, (size_t)rank / 2);

      auto ords_to_csv_str = [](const std::vector<size_t>& ords) {
        std::string str = "";
        for (const auto& o : ords) str += std::to_string(o) + ",";

        str.pop_back();  // remove the trailing comma(,)
        return str;
      };

      auto range_to_csv_str = [&ords_to_csv_str](const size_t& n) {
        std::vector<size_t> range_vec(n);
        for (auto i = 0; i < n; ++i) range_vec[i] = i;
        return ords_to_csv_str(range_vec);
      };

      auto inds = range_to_csv_str(rank);

      DataTensorType result(right_eval.world(), right_eval.trange());
      result.fill(0.);
      for (const auto& p : vp) {
        for (const auto& q : vp) {
          // permutation of the bra
          auto perm_vec = p.sequence;
          // permutation of the ket
          // ket indices = rank/2 + bra indices
          for (auto qq : q.sequence) perm_vec.push_back((size_t)rank / 2 + qq);

          // permute and add
          if (p.phase * q.phase == 1)
            result(inds) = result(inds) + right_eval(ords_to_csv_str(perm_vec));
          else
            result(inds) = result(inds) - right_eval(ords_to_csv_str(perm_vec));
        }  // for q: vp
      }    // for p: vp
      result(inds) = result(inds) * this->get_right_tensor()->get_scalar();
      return result;
    }

    throw std::domain_error("Functionality not yet implemented!");
  }
};

///
/// The leaf node of the EvalTensor tree.
///
/// Leaf type evaluations have an access to the appropriate data-tensor object
/// provided during evaluation. Consequently, evaluation should simply return
/// the data pointed by this EvalTensorLeaf.
///
template <typename DataTensorType>
class EvalTensorLeaf : public EvalTensor<DataTensorType> {
 private:
  /// The sequant Expression this leaf evaltensor corresponds to.
  ExprPtr expr_{nullptr};

 public:
  // Construct from a sequant tensor.
  EvalTensorLeaf(const ExprPtr& expr) : expr_{expr} {}

  /// @return True.
  bool is_leaf() const override { return true; }

  /// Getter of the sequant expr corresponding to this eval tensor.
  const ExprPtr& get_expr() const { return expr_; }

  /// Get a directed graph definitions and paths in dot language.
  std::wstring to_digraph() const override {
    return L"node" + std::to_wstring(this->get_hash_value()) + L" [label = \"" +
           get_expr()->to_latex() + L"\"];\n";
  }

  /// Visit the tree by pre-order traversal.
  void visit(const std::function<void(const EvalTensor<DataTensorType>&)>&
                 visitor = {}) const override {
    visitor(*this);
  }

  /// Evaluate the EvalTensor.
  DataTensorType evaluate(
      const container::map<HashType, std::shared_ptr<DataTensorType>>& context)
      const override {
    if (auto label = expr_->as<Tensor>().label();
        (label == L"A" || label == L"P"))
      throw std::logic_error(
          "(anti-)symmetrization tensors cannot be evaluated from here!");
    auto context_entry = context.find(this->get_hash_value());
    if (context_entry == context.end()) {
      throw std::runtime_error("Data tensor missing!");
    }
    return *context_entry->second;
  }
};

}  // namespace sequant::evaluate

#endif /* SEQUANT_EVALUATE_EVAL_TENSOR_HPP */
