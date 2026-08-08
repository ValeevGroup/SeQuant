#ifndef SEQUANT_EVAL_PLACEMENT_ROUTER_HPP
#define SEQUANT_EVAL_PLACEMENT_ROUTER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace sequant::eval {

/// \brief Override home descriptor.
///
/// \details A home has two parts, resolved to a depth by \c home_depth against
/// the live batch context:
/// - \c slot_positions -- CANONICAL result-slot positions (indices into
///   \c canon_indices()) of the value's home modes that lie ON its own slots.
///   DAG-GLOBAL: identical for every occurrence, and resolved PER USE by
///   mapping each position to THAT occurrence's physical index (\c
///   canon_indices()[p]), so one overlay places two occurrences whose batched
///   slot binds different labels (the g.C legs, i_3 vs i_4) at DIFFERENT
///   depths. This is the case real remat produces (\c shrink_candidates is
///   enclosing INTERSECT carried).
/// - \c free_modes -- home modes NOT on the value's own slots (the value is
///   invariant to them but homed within their loop). These cannot be a slot
///   position, so they are matched by Index identity directly. Empty for real
///   remat moves; used when a value is relocated into a loop it does not carry.
/// Both empty => invariant to the whole nest (the chain root). \c split_index
/// discriminates two cells of one value at the SAME resolved depth
/// (register-allocation split; 0 unless remat splits).
struct HomeTarget {
  container::svector<std::size_t> slot_positions;
  container::svector<Index> free_modes;
  std::size_t split_index = 0;
};

/// \brief A read+store override seam keyed on an occurrence (see
///        \c occurrence_key.hpp): maps a \c TensorNetwork::
///        SlotCanonicalizationMetadata to a \c HomeTarget that overrides
///        where the value is placed/read, in place of the default
///        (\c place_at_this_level's rl+1 walk, \c eval.hpp).
///
/// \details The eval Enter-stage read and the \c place_at_this_level store both
/// consult \c route() (see \c eval.hpp); an absent override (the default, and
/// the only state in Phase 2 -- overrides are populated in Phase 3/4) leaves
/// both sides deriving home exactly as before, so an empty router is inert.
///
/// \tparam TreeNode the evaluation tree node type (same as the owning \c
///         CacheManager's key type).
template <typename TreeNode>
class PlacementRouter {
 public:
  using BatchContext = typename CacheManager<TreeNode, false>::BatchContext;
  using Key = TensorNetwork::SlotCanonicalizationMetadata;

  [[nodiscard]] bool empty() const noexcept {
    return overrides_.empty() && moved_hashes_.empty();
  }

  void set_override(Key key, HomeTarget home) {
    overrides_.insert_or_assign(std::move(key), std::move(home));
  }

  /// \return the override for \p key, or \c nullptr if none is registered --
  ///         the Phase 2 default, meaning the caller derives home as today.
  [[nodiscard]] HomeTarget const* route(Key const& key) const {
    auto it = overrides_.find(key);
    return it == overrides_.end() ? nullptr : &it->second;
  }

  /// Record that a VALUE (by its \c EvalExpr::hash_value canonical identity) is
  /// moved by remat -- i.e. it has an override at ONE OR MORE of its occurrence
  /// keys. Complements \c route() (keyed by the CONTEXT-DEPENDENT occurrence
  /// key): the hoist consults this to learn, from an OUTER scope whose
  /// occurrence-key query would miss, that a value is destined for a DEEPER
  /// home and so must NOT be built full here (see place_at_this_level in
  /// eval.hpp). Populated by \c remat_to_router alongside the per-occurrence
  /// overrides.
  void mark_moved(std::size_t value_hash) { moved_hashes_.insert(value_hash); }

  /// \return true iff \p value_hash was marked moved (see \c mark_moved). A
  ///         context-invariant "is this value demoted anywhere?" query.
  [[nodiscard]] bool moved(std::size_t value_hash) const {
    return moved_hashes_.find(value_hash) != moved_hashes_.end();
  }

  /// \return \c 1 + (deepest index \c i in \p ctx whose mode equals \p node's
  ///         \c canon_indices()[p] for some home slot position \c p), or \c 0
  ///         if none (invariant to the whole nest, i.e. the chain root).
  ///         Resolves the DAG-global \c home.slot_positions to THIS
  ///         occurrence's physical indices via \p node, so two occurrences of
  ///         one value (sharing one overlay) whose batched slot binds different
  ///         labels resolve to DIFFERENT depths. Static mirror of \c
  ///         place_at_this_level's \c rl+1 walk (in \c eval.hpp).
  [[nodiscard]] std::size_t home_depth(HomeTarget const& home,
                                       TreeNode const& node,
                                       BatchContext const& ctx) const {
    auto const& canon = node->canon_indices();
    for (int i = static_cast<int>(ctx.size()) - 1; i >= 0; --i) {
      for (auto const& p : home.slot_positions)
        if (p < canon.size() && ctx[i].first == canon[p])
          return static_cast<std::size_t>(i) + 1;
      for (auto const& m : home.free_modes)
        if (ctx[i].first == m) return static_cast<std::size_t>(i) + 1;
    }
    return 0;
  }

 private:
  std::unordered_map<Key, HomeTarget, RouterKeyHash, RouterKeyEqual> overrides_;
  std::unordered_set<std::size_t> moved_hashes_;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PLACEMENT_ROUTER_HPP
