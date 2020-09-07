#include "eval_tree.hpp"
#include "eval_tree_node.hpp"

#include <SeQuant/core/tensor.hpp>

#include <ctime>
#include <memory>
#include <random>
#include <sstream>

namespace sequant::evaluate {

EvalTree::EvalTree(const ExprPtr& expr, bool canonize_leaf_braket) {
  root = build_expr(expr, canonize_leaf_braket);
}

HashType EvalTree::hash_value() { return root->hash_value(); }

OpsCount EvalTree::ops_count(
    const container::map<IndexSpace::TypeAttr, size_t>& ispace_size_map) {
  return _ops_count(root, ispace_size_map);
}

void EvalTree::visit(
    std::function<void(const std::shared_ptr<const EvalTreeNode>&)> visitor)
    const {
  _visit(root, visitor);
}

auto hsv_to_rgb = [](double h, double s, double v) {
  // https://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/
  int h_i = (int)(h * 6);
  double f = h * 6 - h_i;
  double p = v * (1 - s);
  double q = v * (1 - f * s);
  double t = v * (1 - (1 - f) * s);

  double r, g, b;
  r = g = b = -1.;

  if (h_i == 0) {
    r = v;
    g = t;
    b = p;
  } else if (h_i == 1) {
    r = q;
    g = v;
    b = p;
  } else if (h_i == 2) {
    r = p;
    g = v;
    b = t;
  } else if (h_i == 3) {
    r = p;
    g = q;
    b = v;
  } else if (h_i == 4) {
    r = t;
    g = p;
    b = v;
  } else if (h_i == 5) {
    r = v;
    g = p;
    b = q;
  }

  std::wostringstream rgb;
  rgb << std::hex;
  for (auto c : {r, g, b}) rgb << (int)(256 * c);

  return rgb.str();
};

void EvalTree::digraph(std::wostream& stream) const {
  container::vector<NodeInfo> node_definitions{};
  size_t globalNodeCount = 0;
  node_definitions.push_back(
      NodeInfo{globalNodeCount, root->hash_value(), root->to_latex()});

  std::wostringstream oss;
  _digraph(root, 0, globalNodeCount, node_definitions, oss);
  oss.flush();

  container::map<size_t, size_t> hash_counts;
  for (auto& def : node_definitions) {
    if (auto found = hash_counts.find(def.hash); found != hash_counts.end())
      found->second += 1;
    else
      hash_counts.insert(decltype(hash_counts)::value_type(def.hash, 1));
  }

  //
  // Generate random RGB color code for repeating node labels
  //
  std::random_device seeder;
  const auto seed = seeder.entropy() ? seeder() : std::time(nullptr);
  std::mt19937_64 randEngine(static_cast<std::mt19937::result_type>(seed));
  std::uniform_real_distribution<double> dist;
  const double GOLDEN_RATIO_CONJ = 0.618033988749895;

  container::map<size_t, std::wstring> hash_to_color;
  for (const auto& cc : hash_counts) {
    if (cc.second > 1) {  // repeating hash
      auto hue = GOLDEN_RATIO_CONJ + dist(randEngine);
      hue = hue > 1 ? hue - 1 : hue;
      hash_to_color.insert(decltype(hash_to_color)::value_type(
          cc.first, hsv_to_rgb(hue, 0.5, 0.95)));
    }
  }

  stream << "digraph EvalTree {\n";

  for (const auto& ndef : node_definitions) {
    stream << "node" << ndef.id << " [texlbl = \"$" << ndef.label << "$\"";
    if (auto cc = hash_to_color.find(ndef.hash); cc != hash_to_color.end()) {
      stream << ", color=\"#" << cc->second << "\"";
      stream << ", style=filled";
    }
    stream << "];\n";
  }
  stream << std::endl << oss.str();
  stream << "}" << std::endl;
}

bool EvalTree::swap_labels(const ExprPtr& expr) {
  auto predicate = [&expr](const auto& leaf) {
    return leaf.expr()->to_latex() == expr->to_latex();
  };
  return _swap_braket_labels(root, predicate);
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

void EvalTree::_visit(
    const std::shared_ptr<const EvalTreeNode>& node,
    std::function<void(const std::shared_ptr<const EvalTreeNode>&)> visitor) {
  if (node->is_leaf()) {
    visitor(node);
    return;
  }

  visitor(node);
  auto intrnl_node =
      std::dynamic_pointer_cast<const EvalTreeInternalNode>(node);
  _visit(intrnl_node->left(), visitor);
  _visit(intrnl_node->right(), visitor);
}

void EvalTree::_digraph(std::shared_ptr<const EvalTreeNode> node,
                        size_t node_count, size_t& global_node_count,
                        container::vector<NodeInfo>& node_defs,
                        std::wostream& stream) {
  if (node->is_leaf()) return;
  const auto intrnl_node =
      std::dynamic_pointer_cast<const EvalTreeInternalNode>(node);

  ++global_node_count;
  node_defs.push_back(NodeInfo{global_node_count,
                               intrnl_node->left()->hash_value(),
                               intrnl_node->left()->to_latex()});
  ++global_node_count;
  node_defs.push_back(NodeInfo{global_node_count,
                               intrnl_node->right()->hash_value(),
                               intrnl_node->right()->to_latex()});

  auto leftNodeCount = global_node_count - 1;
  auto rightNodeCount = global_node_count;
  stream << "node" << node_count << " -> {"
         << "node" << leftNodeCount << ", node" << rightNodeCount << "};\n";

  _digraph(intrnl_node->left(), leftNodeCount, global_node_count, node_defs,
           stream);
  _digraph(intrnl_node->right(), rightNodeCount, global_node_count, node_defs,
           stream);
}

bool EvalTree::_swap_braket_labels(
    const EvalNodePtr& node,
    const std::function<bool(const EvalTreeLeafNode&)>& predicate) {
  if (node->is_leaf()) {
    auto leaf_node = std::dynamic_pointer_cast<EvalTreeLeafNode>(node);
    if (predicate(*leaf_node)) {
      leaf_node->swap_labels();
      return true;
    } else
      return false;
  }
  // internal node
  auto intrnl_node = std::dynamic_pointer_cast<EvalTreeInternalNode>(node);

  auto swap_success = _swap_braket_labels(intrnl_node->left(), predicate) ||
                      _swap_braket_labels(intrnl_node->right(), predicate);

  if (swap_success) intrnl_node->update_hash();
  return swap_success;
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
      label == L"A" || label == L"S") {
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
