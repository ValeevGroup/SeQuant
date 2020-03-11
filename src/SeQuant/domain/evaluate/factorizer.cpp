#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"
#include "path_tree.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <tuple>

namespace sequant::evaluate {

/* To-Do
std::tuple<EvalTensorPtr, EvalTensorPtr> fuse_optimally(const ExprPtr& lexpr,
                                                    const ExprPtr& rexpr) {
// cast expressions to Product (TensorNetwork?)
const auto& lprod = std::dynamic_pointer_cast<Product>(lexpr);
const auto& rprod = std::dynamic_pointer_cast<Product>(rexpr);

//
}
*/

namespace detail {

container::set<HashType> get_hash_values(const EvalTensorPtr& tensor) {
  if (tensor->is_leaf())
    return container::set<HashType>{tensor->get_hash_value()};

  const auto& imed = std::dynamic_pointer_cast<EvalTensorIntermediate>(tensor);

  // get hash values of the left eval tensor
  auto lresult = get_hash_values(imed->get_left_tensor());

  // get hash values of the right eval tensor
  auto rresult = get_hash_values(imed->get_right_tensor());

  // merget them together and return
  lresult.insert(rresult.begin(), rresult.end());
  return lresult;
}

OpsCount get_unique_ops_count(const EvalTensorPtr& tensor,
                              const container::set<HashType>& hash_values) {
  // if a hash value exists, don't count it.
  if (hash_values.contains(tensor->get_hash_value())) return 0;

  if (tensor->is_leaf()) return tensor->get_ops_count();

  const auto& imed = std::dynamic_pointer_cast<EvalTensorIntermediate>(tensor);
  return imed->get_hash_value() +
         get_unique_ops_count(imed->get_left_tensor(), hash_values) +
         get_unique_ops_count(imed->get_right_tensor(), hash_values);
}

}  // namespace detail

}  // namespace sequant::evaluate
