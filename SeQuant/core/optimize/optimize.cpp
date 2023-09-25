#include <SeQuant/core/optimize.hpp>
#include <stack>

namespace sequant {

namespace opt {

ExprPtr tail_factor(ExprPtr const& expr) noexcept {
  if (expr->is<Tensor>())
    return expr->clone();

  else if (expr->is<Product>()) {
    auto scalar = expr->as<Product>().scalar();
    if (scalar == 1 && expr->size() == 2) {
      // product with
      //   -single factor that is a tensor
      //   -scalar is just 1
      //  will not be formed because of this block
      return expr->at(1);
    }
    auto facs = ranges::views::tail(*expr);
    return ex<Product>(Product{scalar, ranges::begin(facs), ranges::end(facs)});
  } else {
    // sum
    auto summands = *expr | ranges::views::transform(
                                [](auto const& x) { return tail_factor(x); });
    return ex<Sum>(Sum{ranges::begin(summands), ranges::end(summands)});
  }
}

void pull_scalar(ExprPtr expr) noexcept {
  if (!expr->is<Product>()) return;
  auto& prod = expr->as<Product>();

  auto scal = prod.scalar();
  for (auto&& x : *expr)
    if (x->is<Product>()) {
      auto& p = x->as<Product>();
      scal *= p.scalar();
      p.scale(1 / p.scalar());
    }

  prod.scale(1 / prod.scalar());
  prod.scale(scal);
}

container::vector<container::vector<size_t>> clusters(Sum const& expr) {
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
      auto const node = eval_node<EvalExpr>(term);
      node.visit_internal(visitor);
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

}  // namespace opt

ExprPtr reorder(ExprPtr const& expr) {
  if (!expr->is<Sum>()) return expr;

  Sum result;

  for (auto const& clstr : opt::clusters(expr->as<Sum>()))
    for (auto p : clstr) result.append(expr->at(p));

  assert(result.size() == expr->size());
  return ex<Sum>(std::move(result));
}

}  // namespace sequant