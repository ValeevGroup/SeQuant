#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/optimize/sum.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/operations.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/tail.hpp>
#include <range/v3/view/transform.hpp>

#include <stack>

namespace sequant::opt {

bool has_only_single_atom(const ExprPtr& term) {
  if (term->is_atom()) {
    return true;
  }

  // Recursively check that all elements in the expression tree have only a
  // single element in them. At this point this means checking for Sum or
  // Product objects that only have a single addend or factor respectively.
  return term->size() == 1 && has_only_single_atom(*term->begin());
}

container::vector<container::vector<size_t>> clusters(
    Sum const& expr, container::vector<FullBinaryNode<EvalExpr>> const& nodes) {
  SEQUANT_ASSERT(nodes.size() == expr.size());
  using ranges::views::tail;
  using ranges::views::transform;
  using hash_t = size_t;
  using pos_t = size_t;
  using stack_t = std::stack<pos_t, container::vector<pos_t>>;

  container::map<hash_t, container::set<pos_t>> positions;
  {
    pos_t pos = 0;
    auto visitor = [&positions, &pos](auto const& n) {
      auto h = hash::value(*n);
      if (auto&& found = positions.find(h); found != positions.end()) {
        found->second.emplace(pos);
      } else {
        positions.emplace(h, decltype(positions)::mapped_type{pos});
      }
    };

    for (auto const& term : expr) {
      auto const& node = nodes[pos];
      if (has_only_single_atom(term)) {
        visitor(node);
      } else {
        node.visit_internal(visitor);
      }
      ++pos;
    }
  }

  container::map<pos_t, container::vector<pos_t>> connections;
  {
    for (auto const& [_, v] : positions) {
      auto const v0 = ranges::front(v);
      auto const v_ = ranges::views::tail(v) |
                      ranges::to<decltype(connections)::mapped_type>;
      if (auto&& found = connections.find(v0); found != connections.end())
        for (auto p : v_) found->second.push_back(p);
      else
        connections.emplace(v0, v_);
    }
  }
  positions.clear();

  container::vector<container::vector<pos_t>> result;
  {
    container::set<pos_t> visited;
    for (auto k : connections | ranges::views::keys)
      if (!visited.contains(k)) {
        stack_t dfs_stack;
        dfs_stack.push(k);
        container::vector<pos_t> clstr;
        while (!dfs_stack.empty()) {
          auto p = dfs_stack.top();
          dfs_stack.pop();
          if (!visited.contains(p)) {
            clstr.push_back(p);
            visited.emplace(p);
          }
          if (auto&& found = connections.find(p); found != connections.end())
            for (auto p_ : ranges::views::reverse(found->second))
              dfs_stack.push(p_);
        }
        result.emplace_back(std::move(clstr));
      }
  }
  return result;
}

container::vector<container::vector<size_t>> clusters(Sum const& expr) {
  container::vector<FullBinaryNode<EvalExpr>> nodes;
  nodes.reserve(expr.size());
  for (auto const& term : expr) nodes.push_back(binarize(term));
  return clusters(expr, nodes);
}

Sum reorder(Sum const& sum,
            container::vector<FullBinaryNode<EvalExpr>> const& nodes) {
  Sum result;

  for (auto const& clstr : clusters(sum, nodes))
    for (auto p : clstr) result.append(sum.at(p));

  SEQUANT_ASSERT(result.size() == sum.size());
  return result;
}

Sum reorder(Sum const& sum) {
  container::vector<FullBinaryNode<EvalExpr>> nodes;
  nodes.reserve(sum.size());
  for (auto const& term : sum) nodes.push_back(binarize(term));
  return reorder(sum, nodes);
}

}  // namespace sequant::opt
