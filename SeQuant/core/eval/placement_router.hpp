#ifndef SEQUANT_EVAL_PLACEMENT_ROUTER_HPP
#define SEQUANT_EVAL_PLACEMENT_ROUTER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <cstddef>
#include <unordered_map>
#include <utility>

namespace sequant::eval {

/// \brief Override home descriptor.
///
/// \details \c residency is the enclosing batch-loop modes the value should
/// be variant to at its overridden home; it resolves to a depth via the
/// \c home_depth() rl+1 walk against the live batch_context (empty =>
/// invariant to the whole nest, i.e. the chain root). \c split_index is
/// ALWAYS 0 in Phase 2 (Phase 4 introduces peak splits; it is carried here so
/// that landing splits later does not need a key refactor).
struct HomeTarget {
  container::svector<Index> residency;
  std::size_t split_index = 0;
};

/// \brief A read+store override seam keyed on an occurrence (see
///        \c occurrence_key.hpp): maps a \c TensorNetwork::
///        SlotCanonicalizationMetadata to a \c HomeTarget that overrides
///        where the value is placed/read, in place of the default
///        (\c place_at_this_level's rl+1 walk, \c eval.hpp).
///
/// \details Phase 2 builds this container and the \c CacheManager plumbing it
/// needs (\c cache_manager.hpp), but does not yet wire it into the eval
/// read/store paths -- \c route() is not consulted anywhere yet.
///
/// \tparam TreeNode the evaluation tree node type (same as the owning \c
///         CacheManager's key type).
template <typename TreeNode>
class PlacementRouter {
 public:
  using BatchContext = typename CacheManager<TreeNode, false>::BatchContext;
  using Key = TensorNetwork::SlotCanonicalizationMetadata;

  [[nodiscard]] bool empty() const noexcept { return overrides_.empty(); }

  void set_override(Key key, HomeTarget home) {
    overrides_.insert_or_assign(std::move(key), std::move(home));
  }

  /// \return the override for \p key, or \c nullptr if none is registered --
  ///         the Phase 2 default, meaning the caller derives home as today.
  [[nodiscard]] HomeTarget const* route(Key const& key) const {
    auto it = overrides_.find(key);
    return it == overrides_.end() ? nullptr : &it->second;
  }

  /// \return \c 1 + (deepest index \c i in \p ctx whose mode is in \p
  ///         home.residency), or \c 0 if none (invariant to the whole nest,
  ///         i.e. the chain root). Static mirror of \c place_at_this_level's
  ///         \c rl+1 walk (\c eval.hpp:1730-1739).
  [[nodiscard]] std::size_t home_depth(HomeTarget const& home,
                                       BatchContext const& ctx) const {
    for (int i = static_cast<int>(ctx.size()) - 1; i >= 0; --i)
      for (auto const& ix : home.residency)
        if (ctx[i].first == ix) return static_cast<std::size_t>(i) + 1;
    return 0;
  }

 private:
  std::unordered_map<Key, HomeTarget, RouterKeyHash, RouterKeyEqual> overrides_;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PLACEMENT_ROUTER_HPP
