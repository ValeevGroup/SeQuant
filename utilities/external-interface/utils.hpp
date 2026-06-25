#ifndef SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP
#define SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/indices.hpp>

#include <concepts>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <vector>

class IndexSpaceMeta {
 public:
  struct Entry {
    std::string tag;
    std::string name;
  };

  IndexSpaceMeta() = default;

  std::size_t getSize(const sequant::IndexSpace &space) const;

  std::size_t getSize(const sequant::Index &index) const;

  std::string getName(const sequant::IndexSpace &space) const;

  std::string getTag(const sequant::IndexSpace &space) const;

  void registerSpace(sequant::IndexSpace space, Entry entry);

 private:
  std::map<sequant::IndexSpace, Entry> m_entries;
};

sequant::container::svector<sequant::container::svector<sequant::Index>>
getExternalIndexPairs(const sequant::ExprPtr &expression);

bool needsSymmetrization(const sequant::ExprPtr &expression);

sequant::ExprPtr generateResultSymmetrization(const sequant::ResultExpr &result,
                                              std::wstring_view precursorName);

sequant::ExprPtr generateResultSymmetrization(const sequant::Tensor &result,
                                              std::wstring_view precursorName);

sequant::ExprPtr generateResultSymmetrization(
    std::wstring_view precursorName,
    const sequant::IndexGroups<std::vector<sequant::Index>> &externals,
    const sequant::Tensor &ref);

std::optional<sequant::ExprPtr> pop_symmetrizer(sequant::ResultExpr &expr);

template <std::ranges::random_access_range Container,
          std::ranges::random_access_range Order>
  requires(std::integral<std::ranges::range_value_t<Order>>)
void reorder(Container &&vec, Order &&order) {
  SEQUANT_ASSERT(std::ranges::size(vec) == std::ranges::size(order));
  SEQUANT_ASSERT(std::set<std::size_t>(order.begin(), order.end()).size() ==
                 std::ranges::size(order));

  std::vector<std::ranges::range_value_t<Container>> tmp;
  tmp.reserve(std::ranges::size(vec));
  tmp.insert(tmp.begin(), std::make_move_iterator(vec.begin()),
             std::make_move_iterator(vec.end()));

  const auto offset = std::ranges::min(order);

  for (std::size_t i = 0; i < std::ranges::size(order); ++i) {
    vec[i] = std::move(tmp.at(order[i] - offset));
  }
}

#endif
