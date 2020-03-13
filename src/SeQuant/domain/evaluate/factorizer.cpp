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

class FusionOpsCounter {
 private:
  container::set<HashType> m_hash_values;
  OpsCount m_ops_count;

 public:
  void operator()(const EvalTensorPtr& tree) {
    HashType hvalue = tree->get_hash_value();

    if (m_hash_values.contains(hvalue)) return;

    m_hash_values.insert(hvalue);

    m_ops_count += tree->get_ops_count();
  }

  const container::set<HashType>& get_hash_values() const {
    return m_hash_values;
  }

  OpsCount get_ops_count() const { return m_ops_count; }
};

ExprPtr path_to_product(const std::shared_ptr<PathTree>& path,
                        const ExprPtr& expr) {
  auto product = std::dynamic_pointer_cast<Product>(expr);
  ProductPtr result(new Product{});
  if (path->is_leaf()) {
    result->append(product->factors()[path->get_label()]);
    return result;
  }

  result->append(product->factors()[path->get_label()]);
  for (const auto& i : path->get_children()) {
    auto res = path_to_product(i, product);
    auto& res_product = res->as<Product>();
    // product of single factors is flattened
    if (res_product.factors().size() == 1)
      result->append(1.0, std::move(res_product.factors()[0]));
    // else append as it is
    else
      result->append(std::move(res));
  }
  return result;
}

}  // namespace detail

}  // namespace sequant::evaluate
