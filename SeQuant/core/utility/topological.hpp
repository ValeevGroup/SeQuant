#ifndef SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP
#define SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/enumerate.hpp>

#include <algorithm>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

namespace sequant {

template <typename Func, typename T>
concept dependency_query = requires(const Func &f, const T &val) {
  { f(val) } -> std::ranges::range;
  requires(std::same_as<
           std::remove_cvref_t<std::ranges::range_value_t<decltype(f(val))>>,
           T>);
};

/// @brief Determines the topological ordering (in the computer-science sense)
/// of the indices in [0, n), given their dependency edges directly as
/// indices into that same range
///
/// @param n The number of elements to order
/// @param get_dependencies A function that, given an index in [0, n), yields
/// a range of indices in [0, n) that it depends on
/// @param comp If provided, a strict weak order over indices used to pick
/// among candidates with no (remaining) dependencies; among candidates that
/// compare equivalent under it (or when no comp is given), the smallest
/// index is chosen, so that the result is always deterministic
/// @returns The topological ordering as a list of indices in [0, n)
template <typename IndexDepFunc, typename Comp = std::identity>
  requires((std::same_as<Comp, std::identity> ||
            std::relation<Comp, std::size_t, std::size_t>) &&
           dependency_query<IndexDepFunc, std::size_t>)
std::vector<std::size_t> topological_order_indexed(
    std::size_t n, const IndexDepFunc &get_dependencies, Comp comp = {}) {
  container::svector<std::size_t> indegree(n, 0);
  container::svector<container::svector<std::size_t>> dependents(n);

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t dep : get_dependencies(i)) {
      ++indegree.at(i);
      dependents.at(dep).push_back(i);
    }
  }

  auto tie_broken_less = [&](std::size_t a, std::size_t b) {
    if constexpr (std::same_as<Comp, std::identity>) {
      return a < b;
    } else {
      if (comp(a, b)) return true;
      if (comp(b, a)) return false;
      return a < b;
    }
  };

  container::svector<std::size_t> ready;
  for (std::size_t i = 0; i < n; ++i)
    if (indegree[i] == 0) ready.push_back(i);

  // Kahn's algorithm to select indices in (a) topological order
  std::vector<std::size_t> order;
  order.reserve(n);
  while (!ready.empty()) {
    auto best = std::ranges::min_element(ready, tie_broken_less);
    std::size_t idx = *best;
    ready.erase(best);

    order.push_back(idx);

    for (std::size_t dep : dependents[idx]) {
      SEQUANT_ASSERT(indegree[dep] > 0);
      if (--indegree[dep] == 0) ready.push_back(dep);
    }
  }

  if (order.size() != n) {
    throw Exception(
        "Impossible dependencies encountered in topological_order_indexed()");
  }

  return order;
}

/// @brief Partitions the indices in [0, n) into the strongly connected
/// components of their dependency graph: two indices share a component iff
/// each depends on the other, directly or through other indices. An index
/// with no cycle through it is a component of its own.
///
/// @param n The number of elements
/// @param get_dependencies A function that, given an index in [0, n), yields
/// a range of indices in [0, n) that it depends on
/// @returns per index, the id of its component, in [0, number of components);
/// components are numbered in the order they complete, so an index depends
/// (transitively) only on indices in components with ids at most its own
template <typename IndexDepFunc>
  requires dependency_query<IndexDepFunc, std::size_t>
std::vector<std::size_t> strongly_connected_components_indexed(
    std::size_t n, const IndexDepFunc &get_dependencies) {
  // Tarjan's algorithm, iteratively (the graphs this serves can be deep).
  constexpr std::size_t unseen = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> component(n, unseen);
  container::svector<std::size_t> index(n, unseen), low(n, 0), stack;
  container::svector<char> on_stack(n, 0);
  // the walk: each entry is a node and how many of its dependencies have
  // been visited
  container::svector<std::pair<std::size_t, std::size_t>> walk;
  std::size_t counter = 0, n_components = 0;
  auto const enter = [&](std::size_t x) {
    index[x] = low[x] = counter++;
    stack.push_back(x);
    on_stack[x] = 1;
    walk.push_back({x, 0});
  };
  for (std::size_t start = 0; start < n; ++start) {
    if (index[start] != unseen) continue;
    enter(start);
    while (!walk.empty()) {
      std::size_t const x = walk.back().first;
      auto const deps = get_dependencies(x);
      auto const n_deps = static_cast<std::size_t>(std::ranges::size(deps));
      if (walk.back().second < n_deps) {
        std::size_t const y =
            *std::ranges::next(std::ranges::begin(deps), walk.back().second++);
        if (index[y] == unseen)
          enter(y);
        else if (on_stack[y])
          low[x] = std::min(low[x], index[y]);
        continue;
      }
      walk.pop_back();
      if (!walk.empty())
        low[walk.back().first] = std::min(low[walk.back().first], low[x]);
      if (low[x] != index[x]) continue;
      for (std::size_t y = unseen; y != x;) {
        y = stack.back();
        stack.pop_back();
        on_stack[y] = 0;
        component[y] = n_components;
      }
      ++n_components;
    }
  }
  return component;
}

/// @brief Determines the topological ordering (in the computer-science sense)
/// of the provided elements
///
/// @param range The range of objects whose ordering shall be determined
/// @param get_dependencies A function that yields a range of dependencies
/// for the given object. All these dependencies must be (references to) objects
/// in range.
/// @param comp If provided, this is used to determine the order of elements
/// for which the topological ordering is not unique
/// @returns The topological ordering as a list of indices into range
template <std::ranges::random_access_range Range, typename DepFunc,
          typename Comp = std::identity>
  requires((std::same_as<Comp, std::identity> ||
            std::relation<Comp, std::ranges::range_value_t<Range>,
                          std::ranges::range_value_t<Range>>) &&
           dependency_query<DepFunc, std::ranges::range_value_t<Range>> &&
           std::equality_comparable<std::ranges::range_value_t<Range>>)
std::vector<std::size_t> topological_order(Range &&range,
                                           const DepFunc &get_dependencies,
                                           Comp comp = {}) {
  using std::ranges::begin;
  using std::ranges::end;
  using std::ranges::size;

  using Value = std::ranges::range_value_t<Range>;

  const std::size_t n = size(range);

  // Pre-compute dependencies between elements in range, expressed as indices
  // into range
  container::svector<container::svector<std::size_t>> deps_by_index(n);
  for (const auto &[i, current] : ranges::views::enumerate(range)) {
    for (const Value &current_dep : get_dependencies(current)) {
      auto it = std::ranges::find(range, current_dep);
      SEQUANT_ASSERT(it != end(range));
      deps_by_index[i].push_back(std::ranges::distance(begin(range), it));
    }
  }

  auto get_dep_indices = [&](std::size_t i) -> decltype(auto) {
    return deps_by_index.at(i);
  };

  if constexpr (std::same_as<Comp, std::identity>) {
    return topological_order_indexed(n, get_dep_indices);
  } else {
    return topological_order_indexed(
        n, get_dep_indices, [&](std::size_t a, std::size_t b) {
          return comp(*(begin(range) + a), *(begin(range) + b));
        });
  }
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP
