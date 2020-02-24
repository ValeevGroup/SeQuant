#ifndef SEQUANT_FACTORIZER_LEAF_EVAL_TENSOR_HPP
#define SEQUANT_FACTORIZER_LEAF_EVAL_TENSOR_HPP

#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

#include <boost/any.hpp>

#include <memory>

namespace sequant::factorize {

///
/// The leaf node of the EvalTensor tree.
///
/// Leaf type evaluations have an access to the appropriate data-tensor object
/// provided during evaluation. Consequently, evaluation should simply return
/// the data pointed by this LeafEvalTensor.
///
/// @author Bimal Gaudel
/// @version Feb 2020
///
class LeafEvalTensor : public EvalTensor {
 private:
  /// The data-tensor pointer assigned during the evaluation time based on this
  /// LeafEvalTensor's hash value. Eg. std::shared_ptr<btas::Tensor<double>>
  /// when using BTAS backend. See: https://github.com/BTAS/BTAS
  std::shared_ptr<boost::any> dtensor_ptr_{nullptr};

 public:
  /// Setter method for the dtensor_ptr_;
  /// @param dtensor_ptr A shared_ptr to the data-tensor associated with this
  /// leaf tensor.
  void set_dtensor_ptr(const std::shared_ptr<boost::any>& dtensor_ptr);


  /// @return True.
  bool is_leaf() const override;
};

}  // namespace sequant::factorize

#endif /* SEQUANT_FACTORIZER_LEAF_EVAL_TENSOR_HPP */
