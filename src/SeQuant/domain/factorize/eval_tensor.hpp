#ifndef SEQUANT_FACTORIZER_EVAL_TENSOR_HPP
#define SEQUANT_FACTORIZER_EVAL_TENSOR_HPP

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

 public:
  /// Default constructor.
  EvalTensor() = default;

  /// Default destructor.
  virtual ~EvalTensor() = default;

  /// Setter method of the indices_ field.
  void set_indices(const IndexLabelContainer &);

  /// Getter method of the indices_ field.
  const IndexLabelContainer &get_indices() const;

  /// Setter method of the hash_value_ field.
  void set_hash_value(HashType);

  /// Getter method of the hash_value_ field.
  const HashType get_hash_value() const;

  /// Setter method of the ops_count_ field.
  void set_ops_count(OpsCount);

  /// Getter method of the ops_count_ field.
  const OpsCount get_ops_count() const;

  /// Check if this is the end of the evaluation.
  virtual bool is_leaf() const = 0;
};

}  // namespace sequant::factorize

#endif /* SEQUANT_FACTORIZER_EVAL_TENSOR_HPP */
