//
// Created by Eduard Valeyev on 2019-02-14.
//

#ifndef SEQUANT_HUGENHOLTZ_HPP
#define SEQUANT_HUGENHOLTZ_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/all.hpp>

#include <algorithm>
#include <utility>

namespace sequant {

/// @brief HugenholtzVertex represents a sequence of edges arranged into groups
/// by (topological) equivalence.
/// @tparam Edge the edge type; it must be EqualityComparable
/// @tparam EdgeEquality type of equality tester for Edge objects;
///         `EdgeEquality(const Edge&, const Edge&)` must be implicitly
///         convertible to bool and evaluate to true for equivalent arguments.
template <typename Edge, typename EdgeEquality>
class HugenholtzVertex {
  using Group = std::pair<Edge, container::set<size_t>>;

 public:
  HugenholtzVertex() = default;
  HugenholtzVertex(const HugenholtzVertex&) = default;
  HugenholtzVertex(HugenholtzVertex&&) = default;

  template <typename EdgeRange, typename = std::enable_if_t<!std::is_same_v<
                                    HugenholtzVertex, std::decay_t<EdgeRange>>>>
  HugenholtzVertex(EdgeRange&& edge_range, EdgeEquality equals = EdgeEquality{})
      : equals_(std::move(equals)) {
    const auto num_edges = ranges::size(edge_range);
    edge_to_group_.reserve(num_edges);

    // compute the groups
    size_t edge_idx = 0;
    ranges::for_each(edge_range, [this, &edge_idx](const Edge& edge) {
      const auto grp_it = std::find_if(
          begin(groups_), end(groups_),
          [&edge, this](const Group& grp) { return equals_(grp.first, edge); });
      if (grp_it == end(groups_)) {  // no group yet? create
        groups_.emplace_back(
            std::make_pair(edge, typename Group::second_type{edge_idx}));
        edge_to_group_.push_back(groups_.size() - 1);
      } else {
        [[maybe_unused]] auto result = grp_it->second.insert(edge_idx);
        SEQUANT_ASSERT(result.second);
        edge_to_group_.push_back(grp_it - begin(groups_));
      }
      ++edge_idx;
    });
  }

  /// reports the total number of edges
  /// @return the total number of edges
  auto num_edges() const { return edge_to_group_.size(); }

  /// reports the total number of groups
  /// @return the total number of groups
  /// @note since @c erase() never erases groups, this returns the total number
  /// of groups ever constructed,
  ///       i.e. some of the groups may be empty. Use @c num_nonempty_groups()
  ///       to get the number of groups with edges assigned to them.
  /// @sa HugenholtzVertex::num_nonempty_groups()
  auto num_groups() const { return groups_.size(); }

  /// reports the number of nonempty groups
  /// @return the number of nonempty groups
  /// @sa HugenholtzVertex::num_groups()
  auto num_nonempty_groups() const {
    return std::accumulate(begin(groups_), end(groups_), size_t(0),
                           [](const size_t& left, const Group& right) {
                             return left + (right.second.empty() ? 0 : 1);
                           });
  }

  /// Group accessor
  /// @param edge_idx the edge index
  /// @return the group to which this edge belongs
  const auto& group(size_t edge_idx) const {
    const auto grp_idx = edge_to_group_.at(edge_idx);
    return groups_.at(grp_idx);
  }

  /// Group accessor
  /// @param grp_idx the group ordinal index
  /// @return the group whose ordinal index is @p group_idx
  const auto& group_at(size_t grp_idx) const { return groups_.at(grp_idx); }

  /// Group size accessor
  /// @param edge_idx the edge index
  /// @return the size of the group to which this edge belongs
  /// @note this is equivalent to @c group(edge_idx).size()
  size_t group_size(size_t edge_idx) const {
    const auto result = group(edge_idx).second.size();
    SEQUANT_ASSERT(result > 0);
    return result;
  }

  /// @param edge_idx the ordinal index of the edge to be removed
  /// @param edge the edge descriptor, only used to assert logic when NDEBUG is
  /// not defined
  void erase(size_t edge_idx, [[maybe_unused]] const Edge& edge) {
    // preconditions
    const auto grp_idx = edge_to_group_.at(edge_idx);
    Group& grp = groups_.at(grp_idx);
    SEQUANT_ASSERT(std::find_if(begin(groups_), end(groups_),
                                [this, &edge](const Group& grp) {
                                  return equals_(grp.first, edge);
                                }) != end(groups_));
    SEQUANT_ASSERT(grp.second.find(edge_idx) != end(grp.second));
    SEQUANT_ASSERT(equals_(grp.first, edge));

    // remove the edge from the map and update the groups
    edge_to_group_.erase(edge_to_group_.begin() + edge_idx);
    grp.second.erase(edge_idx);
    for (auto&& g : groups_) {
      ranges::for_each(g.second, [&edge_idx](size_t& e) {
        if (e >= edge_idx) {
          SEQUANT_ASSERT(e > 0);
          --e;
        }
      });
    }
    // groups are never erased!
  }

  /// @param edge_idx the ordinal index of the edge to be inserted
  /// @param edge the edge descriptor
  void insert(const size_t edge_idx, const Edge& edge) {
    // preconditions
    if (edge_idx > edge_to_group_.size()) {
      throw std::out_of_range(
          "HugenholtzVertex::insert : can only insert or append");
    }

    auto grp_it = std::find_if(
        begin(groups_), end(groups_),
        [this, &edge](const Group& grp) { return equals_(grp.first, edge); });

    // update existing edge indices if inserting in the middle
    if (edge_idx < edge_to_group_.size()) {
      for (auto&& g : groups_) {
        ranges::for_each(g.second, [&edge_idx](size_t& e) {
          if (e >= edge_idx) {
            ++e;
          }
        });
      }
    }

    // now update the group with the new index and update the edge->grp map
    if (grp_it == end(groups_)) {  // no group yet? create
      groups_.emplace_back(
          std::make_pair(edge, typename Group::second_type{edge_idx}));
      edge_to_group_.insert(edge_to_group_.begin() + edge_idx,
                            groups_.size() - 1);
    } else {  // just insert but keep the group indices sorted
      [[maybe_unused]] auto result = grp_it->second.insert(edge_idx);
      SEQUANT_ASSERT(result.second);
      edge_to_group_.insert(edge_to_group_.begin() + edge_idx,
                            grp_it - groups_.begin());
    }
  }

 private:
  container::svector<size_t> edge_to_group_;  // maps edge index to group idx
  container::svector<Group>
      groups_;  // group = Edge + (sorted) indices of the equivalent edges
  EdgeEquality equals_;
};

}  // namespace sequant

#endif  // SEQUANT_HUGENHOLTZ_HPP
