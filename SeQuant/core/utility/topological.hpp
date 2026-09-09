#ifndef SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP
#define SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP

#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view/enumerate.hpp>

#include <algorithm>
#include <limits>
#include <map>
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

  std::vector<std::size_t> order;
  order.reserve(size(range));

  std::map<std::size_t, std::size_t> num_deps;
  std::map<std::size_t, std::vector<std::size_t>> dependents;

  // Pre-compute dependencies between elements in range
  for (const auto &[i, current] : ranges::views::enumerate(range)) {
    auto deps = get_dependencies(current);

    num_deps.emplace(i, 0);

    for (const Value &current_dep : deps) {
      auto it = std::ranges::find(range, current_dep);
      SEQUANT_ASSERT(it != end(range));
      std::size_t dep_idx = std::ranges::distance(begin(range), it);
      ++num_deps[i];
      dependents[dep_idx].emplace_back(i);
    }
  }

  // Kahn's algorithm to select entries from range in (a) topological order
  do {
    auto candidates = num_deps | std::views::filter([](const auto &pair) {
                        return pair.second == 0;
                      });

    // If a comparator was specified, use it to determine which of the
    // candidates to select first, otherwise just use the first
    auto elem = [&]() {
      if constexpr (std::same_as<Comp, std::identity>) {
        return begin(candidates);
      } else {
        return std::ranges::min_element(candidates, comp,
                                        [&](const auto &pair) -> const Value & {
                                          return *(begin(range) + pair.first);
                                        });
      }
    }();

    if (elem == end(candidates)) {
      throw Exception(
          "Impossible dependencies encountered in topological_order()");
    }

    std::size_t idx = elem->first;
    // mark as used
    num_deps.erase(elem.base());

    order.emplace_back(idx);

    auto deps = dependents.find(idx);

    if (deps != end(dependents)) {
      for (std::size_t current : deps->second) {
        SEQUANT_ASSERT(num_deps.at(current) > 0);
        --num_deps.at(current);
      }
    }
  } while (order.size() != size(range));

  return order;
}

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_TOPOLOGICAL_HPP
