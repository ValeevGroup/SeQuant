#ifndef SEQUANT_EVALUATE_EVAL_CONTEXT
#define SEQUANT_EVALUATE_EVAL_CONTEXT

#include "eval_tree.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/space.hpp>

#include <memory>
#include <string>

///
/// Context for evaluation of an EvalTensor.
///
/// @author Bimal Gaudel
/// @version Mar 14, 2020
///
namespace sequant::evaluate {

template <typename DataTensorType>
class EvalContext {
 private:
  using ContextMapType =
      container::map<HashType, std::shared_ptr<DataTensorType>>;
  ContextMapType hash_to_tensor_{};

 public:
  EvalContext(container::map<ExprPtr, std::shared_ptr<DataTensorType>>&
                  expr_to_dtensor_map) {
    for (auto& item : expr_to_dtensor_map) {
      auto hash_value = EvalTree<DataTensorType>(item.first).hash_value();
      this->hash_to_tensor_.insert(std::make_pair(hash_value, item.second));
    }
  }

  const auto& get_map() { return hash_to_tensor_; }
};

}  // namespace sequant::evaluate

#endif /* ifndef SEQUANT_EVALUATE_EVAL_CONTEXT */
