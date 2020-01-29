//
// Created by Bimal Gaudel on 1/17/20.
//

#ifndef SEQUANT_EVAL_TENSOR_HPP
#define SEQUANT_EVAL_TENSOR_HPP

#include "eval_fwd.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>

#include <memory>
#include <string>

namespace sequant::evaluate {

///
/// @brief Binary evaluation of SeQuant Expr objects.
///
/// The atomic SeQuant Expr(s) can be thought of either as a representation of a
/// single data-tensor (eg. btas::Tensor<double>, TA::TArrayD), or binary
/// evaluations of such data-tensors. An evaluation could be summing two
/// data-tensors or taking a product of them. A result of a binary evaluation
/// can be an input for another evaluation, i.e. evaluations can be nested.
///
/// This class represents such evaluations. @note No evaluations on data-tensors
/// are actually performed. To do calculations on data-tensors look at eval_expr.hpp
/// @see eval_expr.hpp
/// @author Bimal Gaudel
/// @date Jan 2020
///
class EvalTensor {
 public:
  /// Represents a binary operation
  enum struct Operation {
    /// Represents summing of tensors
    Sum,
    /// Represents contraction of tensors
    Product,
    /// Represents the end of evaluation
    Eval
  };

  EvalTensor() = default;

  ///
  /// Constructs EvalTensor from a Tensor object.
  ///
  /// @note Tensor objects are leaf type sequant Exprs'.
  /// Thus the resulting EvalTensor will be constructed
  /// as leaf. That means, the indices and the hash value
  /// are filled during the construction of the object and
  /// are independent of the other EvalTensors as the object
  /// has no children.
  ///
  EvalTensor(const Tensor&);

  ///
  /// Construct from a general Expr.
  ///
  /// Calls BinaryOpTypeBuilder and EvalTensor(const Tensor&)
  /// recursively as needed.
  ///
  EvalTensor(const ExprPtr&);

  constant_type get_scalar() const;

  void set_scalar(constant_type);

  Operation get_op() const;

  void set_op(Operation);

  const label_container_type& indices() const;

  label_container_type& indices_mutable();

  hash_type get_hash_value() const;

  void set_hash_value(hash_type);

  const EvTensorPtr& left_tensor() const;

  void set_left_tensor(EvTensorPtr&);

  const EvTensorPtr& right_tensor() const;

  void set_right_tensor(EvTensorPtr&);

  /// @return true if the object is a leaf.
  ///
  /// i.e. left_tensor_ and right_tensor_ are null pointers.
  /// However, checking for nullity of either of them is
  /// sufficient as one of them canNOT be a null pointer
  /// while the other is non-null.
  bool is_leaf() const;

#ifdef SEQUANT_HAS_BTAS
  //< If indices empty, traverses the
  //< evaluation tree recursively and fills them up
  void fill_btas_indices();

  //< @return const reference to btas_indices_.
  const btas_index_container& btas_indices() const;
#endif

 private:
  //< A complex number that multiplies corresponding data-tensor ultimately
  constant_type scalar_{1};

  //< Represents the binary operation eg. sum or product
  Operation operation_{Operation::Eval};

  ///
  /// @brief The non-contracting indices.
  ///
  /// for a product type evaluation.
  /// Or, a set of indices from the two identical sets of of left and right
  /// tensors in case of a sum type evaluation.
  ///
  label_container_type indices_{};

  //< The unique identifier of the evaluation resulting into this object
  hash_type hash_value_{0};

  //< Left tensor for the binary evaluation
  EvTensorPtr left_tensor_{nullptr};

  //< Right tensor for the binary evaluation
  EvTensorPtr right_tensor_{nullptr};

  //< BTAS works with ordinal as index label rather than wstring
#ifdef SEQUANT_HAS_BTAS
  btas_index_container btas_indices_{};
#endif
};

/// Build binary evaluations from EvalTensor and Tensor
class BinaryOpTypeBuilder {
 public:
  BinaryOpTypeBuilder(EvalTensor::Operation opr = EvalTensor::Operation::Eval)
      : op{opr} {}

  ///
  /// Build EvalTensor from another EvalTensor object and 
  /// a sequant Expr.
  /// 
  /// Intended to be used with std::accumulate
  /// @return std::shared_ptr<EvalTensor>
  //
  EvTensorPtr operator()(EvTensorPtr&, const ExprPtr&) const;

 private:
  /// The type of operation between the parameters that gave the result.
  ///
  /// Could either be EvalTensor::Operation::Sum
  /// or EvalTensor::Operation::Product. Conceivably,
  /// EvalTensor::Operation::Eval should not be
  /// assigned to @param op from here.
  /// 
  EvalTensor::Operation op;
};

}  // namespace sequant::evaluate

#endif  // SEQUANT_EVAL_TENSOR_HPP
