#ifndef SEQUANT_FACTORIZE_EVAL_TREE_BUILDER_HPP
#define SEQUANT_FACTORIZE_EVAL_TREE_BUILDER_HPP

#include <SeQuant/core/expr_fwd.hpp>

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

 public:
  /// Construct EvalTree from sequant Expr.
  /// @param expr sequant ExprPtr to a sequant expression.
  EvalTreeBuilder(const ExprPtr& expr);

  /// Getter for the built EvalTree.
  /// @return EvalTensorPtr to the evaluation tree.
  const EvalTensorPtr& get_eval_tree() const;

  /// Build EvalTensor from a sequant Product.
  /// @param expr ExprPtr to sequant Product.
  /// @param path PathTreePtr that determines the sequence of evaluation.
  /// @return Evaluation tree pointer.
  static EvalTensorPtr& build_from_product(const ExprPtr& expr,
                                           const PathTreePtr& path);

  /// Build leaf EvalTensor from sequant tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant Tensor.
  /// @param real_valued Whether to assume the data-tensor is real valued. True
  /// by default, set false if the the data tensor should be assumed to be
  /// complex valued.
  ///
  /// @return EvalTensorPtr to a EvalTensor.
  static EvalTensorPtr build_leaf(const ExprPtr& expr, bool real_valued = true);

  /// Build intermediate tensor from sequant expr.
  ///
  /// @param ltensor The left evaluation tensor.
  /// @param rtensor The right evaluation tensor.
  /// @return EvalTensorPtr to a EvalTensor.
  static EvalTensorPtr build_intermediate(const EvalTensorPtr& ltensor,
                                          const EvalTensorPtr& rtensor,
                                          Operation op);

  /// Hash leaf tensor.
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @param swap_bra_ket_labels hash kets before bras if true. Should be set
  /// true while hashing real-valued tensors.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  static HashType hash_leaf(const ExprPtr& expr, bool swap_bra_ket_labels);

  /// Hash intermediate tensor
  ///
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return Hash value of the tensor based on its kind and the index space of
  /// its bra and ket indices.
  static HashType hash_intermediate(const EvalTensorPtr& expr);

  /// Should we swap the whole bra indices with the whole ket indices?
  ///
  /// When a sequant Tensor represents a real valued data tensor extra
  /// canonicalization is done.
  /// @param expr sequant ExprPtr to a sequant expression.
  /// @return True if swapping is necessary.
  /// @throw domain_error if bra rank and ket rank do not match.
  static bool swap_bra_ket_labels(const ExprPtr& expr);
};

}  // namespace sequant::factorize

#endif /* ifndef SEQUANT_FACTORIZE_EVAL_TREE_BUILDER_HPP */
