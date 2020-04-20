#include "eval_tree.hpp"
#include "eval_tree_node.hpp"

#include <SeQuant/core/tensor.hpp>

#include <memory>

namespace sequant::evaluate {

HashType EvalTree::hash_value() { return root->hash_value(); }

OpsCount EvalTree::ops_count(
    const container::map<IndexSpace::TypeAttr, size_t>& ispace_size_map) {
  return _ops_count(root, ispace_size_map);
}

OpsCount EvalTree::_ops_count(
    const EvalNodePtr& node,
    const container::map<IndexSpace::TypeAttr, size_t>& ispace_size_map) {
  if (node->is_leaf()) return 0;

  auto intrnl_node = std::dynamic_pointer_cast<EvalTreeInternalNode>(node);

  OpsCount count = 0;
  auto op_type = intrnl_node->operation();

  if (op_type == Operation::PRODUCT) {
    auto& left_indices = intrnl_node->left()->indices();
    auto& right_indices = intrnl_node->right()->indices();

    container::set<Index> unique_indices;

    for (const auto& idx : left_indices) unique_indices.insert(idx);
    for (const auto& idx : right_indices) unique_indices.insert(idx);

    auto contraction_ops = 1;
    for (const auto& idx : unique_indices)
      contraction_ops *= (ispace_size_map.find(idx.space().type()))->second;
    count += contraction_ops;
  }

  return count + _ops_count(intrnl_node->left(), ispace_size_map) +
         _ops_count(intrnl_node->right(), ispace_size_map);
}

void EvalTree::visit(const std::function<void(const EvalNodePtr&)>& visitor) {
  visitor(root);
  if (root->is_leaf()) return;

  auto intrnl_node = std::dynamic_pointer_cast<EvalTreeInternalNode>(root);
  visitor(intrnl_node->left());
  visitor(intrnl_node->right());
}

EvalTree::EvalTree(const ExprPtr& expr, bool canonize_leaf_braket) {
  root = build_expr(expr, canonize_leaf_braket);
}

EvalNodePtr EvalTree::build_expr(const ExprPtr& expr,
                                 bool canonize_leaf_braket) {
  if (expr->is<Tensor>())
    return std::make_shared<EvalTreeLeafNode>(
        EvalTreeLeafNode(expr, canonize_leaf_braket));
  else if (expr->is<Sum>())
    return build_sum(expr, canonize_leaf_braket);
  else if (expr->is<Product>())
    return build_prod(expr, canonize_leaf_braket);
  else
    throw std::logic_error(
        "Only sum, product or tensor is allowed for eval node construction");
}

EvalNodePtr EvalTree::build_sum(const ExprPtr& expr,
                                bool canonize_leaf_braket) {
  auto sum_accumulator = [canonize_leaf_braket](const EvalNodePtr& lexpr,
                                                const ExprPtr& summand) {
    return std::make_shared<EvalTreeInternalNode>(EvalTreeInternalNode(
        lexpr, build_expr(summand, canonize_leaf_braket), Operation::SUM));
  };

  auto& sum = expr->as<Sum>();

  return std::accumulate(sum.begin() + 1, sum.end(),
                         build_expr(sum.summand(0), canonize_leaf_braket),
                         sum_accumulator);
}

EvalNodePtr EvalTree::build_prod(const ExprPtr& expr,
                                 bool canonize_leaf_braket) {
  //
  auto prod_accumulator = [canonize_leaf_braket](const EvalNodePtr& lexpr,
                                                 const ExprPtr& factor) {
    return std::make_shared<EvalTreeInternalNode>(EvalTreeInternalNode(
        lexpr, build_expr(factor, canonize_leaf_braket), Operation::PRODUCT));
  };

  auto& prod = expr->as<Product>();
  auto& fac0 = prod.factor(0);

  if (auto label = fac0->is<Tensor>() ? fac0->as<Tensor>().label() : L"";
      label == L"A" || label == L"P") {
    // (anti-)symmetrization tensor encountered
    auto right = std::make_shared<Product>(prod.begin() + 1, prod.end());
    right->scale(prod.scalar());
    return std::make_shared<EvalTreeInternalNode>(EvalTreeInternalNode(
        build_expr(fac0, canonize_leaf_braket),
        build_expr(right, canonize_leaf_braket),
        (label == L"A") ? Operation::ANTISYMMETRIZE : Operation::SYMMETRIZE));
  } else {
    auto init = build_expr(fac0, canonize_leaf_braket);
    init->scale(prod.scalar().real());
    return std::accumulate(prod.begin() + 1, prod.end(), init,
                           prod_accumulator);
  }
}

container::svector<std::tuple<int, container::svector<size_t>>>
EvalTree::_phase_perm(container::svector<size_t>& ords, size_t begin,
                      size_t swaps_count) {
  if (begin + 1 == ords.size()) {
    // found a new permutation
    // even permutation phase: +1 odd: -1.
    int phase = (swaps_count % 2 == 0) ? 1 : -1;

    return container::svector<std::tuple<int, container::svector<size_t>>>{
        std::make_tuple(phase, ords)};
  }

  // recursively call the function to compute the permutations
  // @note swaps_count is incremented only if the swapped indices are not
  // the same ie. to swap an element with itself is not counted.
  container::svector<std::tuple<int, container::svector<size_t>>> result;
  for (auto ii = begin; ii < ords.size(); ++ii) {
    std::swap(ords[begin], ords[ii]);

    auto recursed_result = _phase_perm(
        ords, begin + 1, ii == begin ? swaps_count : swaps_count + 1);

    for (auto& p : recursed_result) result.push_back(p);

    // undo the swap so the side effect is nullified for next call
    std::swap(ords[begin], ords[ii]);
  }
  return result;
}

}  // namespace sequant::evaluate
