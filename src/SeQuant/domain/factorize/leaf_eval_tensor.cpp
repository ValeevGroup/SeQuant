#include "leaf_eval_tensor.hpp"
#include "eval_tensor.hpp"
#include "eval_tensor_fwd.hpp"

#include <boost/any.hpp>

#include <memory>

namespace sequant::factorize {

void LeafEvalTensor::set_dtensor_ptr(
    const std::shared_ptr<boost::any>& dtensor_ptr) {
  dtensor_ptr_ = dtensor_ptr;
}

bool LeafEvalTensor::is_leaf() const { return true; }

}  // namespace sequant::factorize
