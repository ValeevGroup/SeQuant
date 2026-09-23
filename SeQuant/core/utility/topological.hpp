#ifndef SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP
#define SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP

#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/enumerate.hpp>

#include <algorithm>
#include <ranges>
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
  std::vector<std::size_t> indegree(n, 0);
  std::vector<std::vector<std::size_t>> dependents(n);

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

  std::vector<std::size_t> ready;
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
  std::vector<std::vector<std::size_t>> deps_by_index(n);
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
