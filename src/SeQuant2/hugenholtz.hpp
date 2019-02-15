//
// Created by Eduard Valeyev on 2019-02-14.
//

#ifndef SEQUANT2_HUGENHOLTZ_HPP
#define SEQUANT2_HUGENHOLTZ_HPP

namespace sequant2 {

template <typename Edge, typename Compare>
class HugenholtzVertex {
  using Group = std::pair<Edge, container::set<size_t>>;

 public:
  HugenholtzVertex() = default;
  HugenholtzVertex(const HugenholtzVertex&) = default;
  HugenholtzVertex(HugenholtzVertex&&) = default;

  template <typename EdgeRange, typename = std::enable_if_t<!std::is_same_v<
                                    HugenholtzVertex, std::decay_t<EdgeRange>>>>
  HugenholtzVertex(EdgeRange&& edge_range, Compare compare = Compare{})
      : compare_(std::move(compare)) {
    const auto num_edges = ranges::size(edge_range);
    edge_to_group_.reserve(num_edges);

    // compute the groups
    size_t edge_counter = 0;
    ranges::for_each(edge_range, [this, &edge_counter](const Edge& edge) {
      const auto grp_it = std::find_if(begin(groups_), end(groups_),
                                       [&edge, this](const Group& grp) {
                                         return compare_(grp.first, edge);
                                       });
      if (grp_it == end(groups_)) {  // no group yet? create
        groups_.emplace_back(edge, {edge_counter});
      } else {
        grp_it->second.push_back(edge_counter);
      }
      ++edge_counter;
    });
  }

  /// Group accessor
  /// @param edge_idx the edge index
  /// @return the group to which this edge belongs
  const auto& group(size_t edge_idx) const {
    const auto grp_idx = edge_to_group_.at(edge_idx);
    return groups_.at(grp_idx).second;
  }

  /// Group size accessor
  /// @param edge_idx the edge index
  /// @return the size of the group to which this edge belongs
  /// @note this is equivalent to @c group(edge_idx).size()
  const size_t group_size(size_t edge_idx) const {
    const auto result = group(edge_idx).size();
    assert(result > 0);
    return result;
  }

  /// @param edge_idx the ordinal index of the edge to be removed
  /// @param edge the edge descriptor, only used to assert logic when NDEBUG is
  /// not defined
  void erase(size_t edge_idx, const Edge& edge) {
    // preconditions
    const auto grp_idx = edge_to_group_.at(edge_idx);
    Group& grp = groups_.at(grp_idx);
    assert(std::find_if(begin(groups_), end(groups_),
                        [this, &edge](const Group& grp) {
                          return compare_(grp.first, edge);
                        }) != end(groups_));
    assert(grp.second.find(edge_idx) != end(grp.second));
    assert(compare_(grp.first, edge));

    // remove the edge from the map and update the groups
    edge_to_group_.erase(edge_to_group_.begin() + edge_idx);
    grp.second.erase(edge_idx);
    for (auto&& g : groups_) {
      ranges::for_each(g.second, [&edge_idx](size_t& e) {
        if (e >= edge_idx) {
          assert(e > 0);
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
    assert(edge_idx <=
           edge_to_group_.size());  // can insert into the set, or at the end

    auto grp_it = std::find_if(
        begin(groups_), end(groups_),
        [this, &edge](const Group& grp) { return compare_(grp.first, edge); });

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
      auto result = grp_it->second.insert(edge_idx);
      assert(result.second);
    }
  }

 private:
  container::svector<size_t> edge_to_group_;  // maps edge index to group idx
  container::vector<Group>
      groups_;  // group = Edge + (sorted) indices of the equivalent edges
  Compare compare_;
};

}  // namespace sequant2

#endif  // SEQUANT2_HUGENHOLTZ_HPP
