#ifndef SEQUANT_FACTORIZER_INTERMEDIATE_EVAL_TENSOR_HPP
#define SEQUANT_FACTORIZER_INTERMEDIATE_EVAL_TENSOR_HPP

#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

namespace sequant::factorize {

///
/// Non-leaf type EvalTensor node in a binary evaluation of the
/// sequant::Expr at some context.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
class IntermediateEvalTensor : public EvalTensor {
 private:
  /// Evaluation node to the left.
  std::shared_ptr<EvalTensor> left_tensor_ptr_{nullptr};

  /// Evaluation node to the right
  std::shared_ptr<EvalTensor> right_tensor_ptr_{nullptr};

  /// The type of the binary evaluation.
  Operation op_type_;

 public:
  /// @return false
  bool is_leaf() const override;

  /// Setter method for the left evaluation node.
  void set_left_ptr(const std::shared_ptr<EvalTensor>&);

  /// Getter method for the left evaluation node.
  const EvalTensorPtr& get_left_tensor() const;

  /// Setter method for the right evaluation node.
  void set_right_ptr(const std::shared_ptr<EvalTensor>&);

  /// Getter method for the left evaluation node.
  const EvalTensorPtr& get_right_tensor() const;

  /// Setter method for the operation type.
  /// @param op An Operation.
  void set_operation(Operation op);

  /// Getter method for the operation type.
  /// @return Operation for this intermediate tensor.
  const Operation get_operation() const;
};

}  // namespace sequant::factorize

#endif /* SEQUANT_FACTORIZER_INTERMEDIATE_EVAL_TENSOR_HPP  */
