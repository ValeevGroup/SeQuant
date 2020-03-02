#ifndef SEQUANT_FACTORIZE_EVAL_TENSOR_BUILDER_HPP
#define SEQUANT_FACTORIZE_EVAL_TENSOR_BUILDER_HPP

#include <SeQuant/core/expr.hpp>

#include "eval_tensor_fwd.hpp"
#include "path_tree.hpp"

///
/// Evaluation tree builder from sequant Expr.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
namespace sequant::factorize {
class EvalTreeBuilder {
 private:
  /// The evaluation tree built by this builder.
  EvalTensorPtr eval_tree_{nullptr};

  /// By default assume the tensor data is real -- not complex.
  const bool complex_tensor_data{false};

 public:
  /// Constructor.
  /// @param complex_tensor_data false by default. Set to true when working with
  /// complex-valued tensors.
  EvalTreeBuilder(bool complex_tensor_data = false);

  /// Getter for the built EvalTree.
  /// @return EvalTensorPtr to the evaluation tree.
  const EvalTensorPtr& get_eval_tree() const;

  /// Build EvalTensor from a sequant Product.
  /// @param prod ProductPtr to sequant Product.
  /// @param path shared_ptr to PathTree that determines the sequence of
  /// evaluation.
  ///
  /// @return Evaluation tree pointer.
  const EvalTensorPtr build_from_product(const ProductPtr& prod,
                                         const PathTreePtr& path) const;

  /// Build leaf EvalTensor from sequant tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant Tensor.
  /// @return EvalTensorPtr to a EvalTensor.
  const EvalTensorPtr build_leaf(const ExprPtr& expr) const;

  /// Build intermediate tensor from sequant expr.
  ///
  /// @param ltensor The left evaluation tensor.
  /// @param rtensor The right evaluation tensor.
  /// @return EvalTensorPtr to a EvalTensor.
  /// @throw domain_error when the indices do not match for sum type evaluations
  /// of two tensors.
  const EvalTensorPtr build_intermediate(const EvalTensorPtr& ltensor,
                                         const EvalTensorPtr& rtensor,
                                         Operation op) const;

  /// Hash leaf tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @param swap_bra_ket_labels hash kets before bras if true. Should be set
  /// true while hashing real-valued tensors.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_leaf(const ExprPtr& expr, bool swap_bra_ket_labels) const;

  /// Hash intermediate tensor
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  HashType hash_intermediate(const EvalTensorPtr& expr) const;

  /// Should we swap the whole bra indices with the whole ket indices?
  ///
  /// When a sequant Tensor represents a real valued data tensor extra
  /// canonicalization is done.
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return True if swapping is necessary.
  /// @throw domain_error if bra rank and ket rank do not match.
  bool swap_bra_ket_labels(const ExprPtr& expr) const;
};

}  // namespace sequant::factorize

#endif /* ifndef SEQUANT_FACTORIZE_EVAL_TENSOR_BUILDER_HPP */
