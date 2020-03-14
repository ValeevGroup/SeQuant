#ifndef SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP
#define SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP

#include <SeQuant/core/expr_fwd.hpp>

#include "eval_tensor_fwd.hpp"
// #include "path_tree.hpp"

///
/// Evaluation tree builder from sequant Expr.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
namespace sequant::evaluate {
class EvalTensorBuilder {
 private:
  /// By default assume the tensor data is real -- not complex.
  const bool complex_tensor_data{false};

 public:
  /// Constructor.
  /// @param complex_tensor_data false by default. Set to true when working with
  /// complex-valued tensors.
  explicit EvalTensorBuilder(bool complex_tensor_data = false);

  /// Build evaluation tree from a sequant expression.
  /// @note It is assumed that the expr is a sum and each summand is a product
  /// such that the first factor of the product is either an antisymmetrization
  /// tensor or symmetrization tensor.
  /// @param expr sequant ExprPtr.
  /// @return Pointer to the Evaluation Tensor
  /// @throws std::domain_error when the expr is neither Sum, Product or Tensor.
  EvalTensorPtr build_tree(const ExprPtr& expr) const;

 private:
  /// Build EvalTensor from a sequant Product.
  /// @param expr sequant ExprPtr to sequant Product.
  ///
  /// @return EvalTensor pointer.
  EvalTensorPtr build_product(const ExprPtr& expr) const;

  /// Build EvalTensor from a sequant Sum.
  /// @param expr sequant ExprPtr to sequant Sum.
  ///
  /// @return EvalTensor pointer.
  EvalTensorPtr build_sum(const ExprPtr& expr) const;

  /// Build leaf EvalTensor from sequant tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant Tensor.
  /// @return EvalTensorPtr to a EvalTensor.
  EvalTensorPtr build_leaf(const ExprPtr& expr) const;

  /// Build binary evaluation intermediate tensor from two evaluation tensors.
  ///
  /// @param left_eval_expr The left evaluation tensor.
  /// @param right_eval_expr The right evaluation tensor.
  /// @param op The binary evaluation type eg. sum, product, antisymmetrization
  /// @return EvalTensorPtr to a EvalTensor.
  /// @throw domain_error when the indices do not match for sum type evaluations
  /// of two tensors.
  EvalTensorPtr build_intermediate(const EvalTensorPtr& left_eval_expr,
                                   const EvalTensorPtr& right_eval_expr,
                                   Operation op) const;

  /// Hash leaf tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @param swap_bra_ket_labels Hash kets before bras if true. Should be set
  /// true while hashing real-valued tensors.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_leaf(const ExprPtr& expr, bool swap_bra_ket_labels) const;

  /// Hash intermediate tensor
  ///
  /// @param eval_expr sequant ExprPtr to a sequant expression.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_intermediate(const EvalTensorPtr& eval_expr) const;

  /// Should we swap the whole bra indices with the whole ket indices?
  ///
  /// When a sequant Tensor represents a real valued data tensor extra
  /// canonicalization is done.
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return True if swapping is necessary.
  /// @throw domain_error if bra rank and ket rank do not match.
  bool need_bra_ket_swap(const ExprPtr& expr) const;
};

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_EVAL_TENSOR_BUILDER_HPP */
