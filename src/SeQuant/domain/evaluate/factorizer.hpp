#ifndef SEQUANT_EVALUATE_FACTORIZER_HPP
#define SEQUANT_EVALUATE_FACTORIZER_HPP

///
/// Find optimal evaluation sequences of EvalTensor evaluation tree that result
/// into lower ops count.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///

#include "eval_tensor.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <tuple>

namespace sequant::evaluate {

/// Generate the corresponding EvalTensor's of @c lexpr and @c rexpr. Such that
/// the resulting EvalTensor's have the maximum number of common sub-expressions
/// between them.
///
/// @param lexpr left sequant Expr to be fused.
/// @param rexpr right sequant Expr to be fused.
/// @return A tuple of shared pointers to EvalTensor's.
template <typename DataTensorType>
std::tuple<EvalTensorPtr<DataTensorType>, EvalTensorPtr<DataTensorType>>
fuse_optimally(const ExprPtr& lexpr, const ExprPtr& rexpr);

namespace detail {

///
/// A functor for eval_tensor visit method.
///
/// It collects unique hash values and the ops counts of the corresponding
/// evaluation tensors. The intention is to initialize such an object once and
/// pass it as a parameter for visit method of EvalTensor objects to get the
/// ops count when such objects are fused.
///
class FusionOpsCounter {
 private:
  container::set<HashType> m_hash_values;

  OpsCount m_ops_count{0};

 public:
  template <typename DataTensorType>
  void operator()(const EvalTensorPtr<DataTensorType>& tree);

  const container::set<HashType>& get_hash_values() const;

  OpsCount get_ops_count() const;
};

/// @brief Translate a PathTree to a TensorNetwork object.
/// @param path A PathTree object that encodes how the product will be
/// factorized.
/// @param expr A expr to seuant product to be factorized.
/// @return ExprPtr to a sequant Product.
ExprPtr path_to_product(const std::shared_ptr<PathTree>& path,
                        const ExprPtr& expr);

}  // namespace detail

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_FACTORIZER_HPP */
