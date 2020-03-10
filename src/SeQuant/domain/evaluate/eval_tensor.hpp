#ifndef SEQUANT_FACTORIZER_EVAL_TENSOR_HPP
#define SEQUANT_FACTORIZER_EVAL_TENSOR_HPP

#include <boost/any.hpp>
#include <memory>

#include "eval_tensor_fwd.hpp"

namespace sequant::factorize {
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

class EvalTensor {
 protected:
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
  void set_indices(const IndexLabelContainer&);

  /// Getter method of the indices_ field.
  const IndexLabelContainer& get_indices() const;

  /// Setter method of the hash_value_ field.
  void set_hash_value(HashType);

  /// Getter method of the hash_value_ field.
  HashType get_hash_value() const;

  /// Setter method of the ops_count_ field.
  void set_ops_count(OpsCount);

  /// Getter method of the ops_count_ field.
  OpsCount get_ops_count() const;

  /// Getter for the scalar field.
  ScalarType get_scalar() const;

  /// Setter for the scalar field.
  void set_scalar(ScalarType);

  /// Check if this is the end of the evaluation.
  virtual bool is_leaf() const = 0;
};

///
/// Non-leaf type EvalTensor node in a binary evaluation of the
/// sequant::Expr at some context.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
class EvalTensorIntermediate : public EvalTensor {
 private:
  /// Evaluation node to the left.
  std::shared_ptr<EvalTensor> left_tensor_{nullptr};

  /// Evaluation node to the right
  std::shared_ptr<EvalTensor> right_tensor_{nullptr};

  /// The type of the binary evaluation.
  Operation operation_;

 public:
  /// @return false
  bool is_leaf() const override;

  /// Setter method for the left evaluation node.
  void set_left_tensor(const std::shared_ptr<EvalTensor>&);

  /// Getter method for the left evaluation node.
  const EvalTensorPtr& get_left_tensor() const;

  /// Setter method for the right evaluation node.
  void set_right_tensor(const std::shared_ptr<EvalTensor>&);

  /// Getter method for the left evaluation node.
  const EvalTensorPtr& get_right_tensor() const;

  /// Setter method for the operation type.
  /// @param op An Operation.
  void set_operation(Operation op);

  /// Getter method for the operation type.
  /// @return Operation for this intermediate tensor.
  Operation get_operation() const;
};

///
/// The leaf node of the EvalTensor tree.
///
/// Leaf type evaluations have an access to the appropriate data-tensor object
/// provided during evaluation. Consequently, evaluation should simply return
/// the data pointed by this EvalTensorLeaf.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
class EvalTensorLeaf : public EvalTensor {
 private:
  /// The data-tensor pointer assigned during the evaluation time based on this
  /// EvalTensorLeaf's hash value. Eg. std::shared_ptr<btas::Tensor<double>>
  /// when using BTAS backend. See: https://github.com/BTAS/BTAS
  std::shared_ptr<boost::any> data_tensor_{nullptr};

 public:
  /// Setter method for the data_tensor_;
  /// @param dtensor_ptr A shared_ptr to the data-tensor associated with this
  /// leaf tensor.
  void set_data_tensor(const std::shared_ptr<boost::any>& dtensor_ptr);

  /// @return True.
  bool is_leaf() const override;
};

}  // namespace sequant::factorize

#endif /* SEQUANT_FACTORIZER_EVAL_TENSOR_HPP */
