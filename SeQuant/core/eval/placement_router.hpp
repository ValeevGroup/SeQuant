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
/// \details A home is a label-free DAG-SCOPE: an ordered sequence of SPACES
/// describing a batch-loop NEST PREFIX (e.g. \c [occ, occ, aux]), independent
/// of any value or tree. It names a level by its POSITION in the sequence, so a
/// space repeated ALONG the nest (the g.C \c i_3-outer / \c i_4-inner chain) is
/// disambiguated by position without a separate instance index. A DAG-scope is
/// DAG-GLOBAL: identical for every occurrence of the value. It is resolved PER
/// USE by \c home_depth, which maps each scope space to THAT occurrence's
/// physical batched index (of the same space) via the occurrence key, so one
/// overlay places two occurrences whose batched slot binds different labels at
/// DIFFERENT depths -- realizing the split with a single overlay.
///
/// Every \c dag_scope entry is a space the value CARRIES: \c home_depth
/// resolves it to a physical loop via \c key.get_indices() (the occurrence's
/// own batched slots), and \c remat_to_router only ever emits carried-mode
/// spaces (a cell's \c home_modes are a subset of its \c carried). A level the
/// value is homed WITHIN but does NOT carry (the value is invariant to that
/// loop but rebuilt inside it -- leg B's \c i_3) is therefore NOT a \c
/// dag_scope entry: it has no per-occurrence physical anchor, so it cannot be
/// resolved here. It is instead derived from the enclosing NEST (levels below
/// the home that the cell does not carry) where it is priced as recompute --
/// the replication factor (see the design spec, Design section 3). This is why
/// the old \c free_modes special case is gone rather than folded in: a
/// homed-within level is a nest fact, not a home coordinate. An empty \c
/// dag_scope => invariant to the whole nest (the chain root). Split cells of
/// one value differ in \c dag_scope; the DAG-scope alone identifies the home,
/// so no separate split index is needed.
///
/// caveat: the space-sequence is a unique coordinate only when no nest level
/// has same-space SIBLINGS (see the design spec section "DAG-scopes form a
/// tree", Coordinate caveat). Our construction only ever produces same-space
/// NESTING, never un-collapsed same-space siblings, so the sequence suffices.
struct HomeTarget {
  container::svector<IndexSpace> dag_scope;
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

  /// \return \c 1 + (the DEEPEST index \c i in \p ctx that this occurrence
  ///         binds to some space of \p home's DAG-scope), or \c 0 if none
  ///         (invariant to the whole nest, i.e. the chain root). Resolves the
  ///         DAG-global \c home.dag_scope to THIS occurrence's physical loop
  ///         via \p key (its batched-slot canonicalization): each scope space
  ///         \c S is matched to the occurrence's in-scope batched physical
  ///         index(es) of space \c S (\c key.get_indices()), each of which is
  ///         mapped to its live-loop depth in \p ctx; the deepest such depth
  ///         wins. So two occurrences of one value (sharing one overlay) whose
  ///         batched slot binds different labels (the g.C legs, \c i_3 vs
  ///         \c i_4) resolve to DIFFERENT depths. Static mirror of \c
  ///         place_at_this_level's \c rl+1 walk (in \c eval.hpp).
  ///
  /// \note Taking the DEEPEST matched depth is insensitive to same-space
  ///       ORDINAL ordering within \c ctx (which \c i_3 vs \c i_4 sits outer),
  ///       and correct even when a space appears \c k times in \c dag_scope: at
  ///       most \c k physical indices of that space participate, and the
  ///       deepest is what the nest prefix resolves to. caveat: this collapses
  ///       to one depth if ONE occurrence bound MULTIPLE same-space batched
  ///       indices that a nest with same-space SIBLINGS would need to keep
  ///       distinct -- out of scope here (see the design spec "DAG-scopes form
  ///       a tree").
  [[nodiscard]] std::size_t home_depth(HomeTarget const& home,
                                       BatchContext const& ctx,
                                       Key const& key) const {
    auto const phys = key.template get_indices<container::svector<Index>>();
    std::size_t deepest = 0;  // 0 => no match (invariant to the whole nest)
    for (auto const& space : home.dag_scope)
      for (auto const& ix : phys) {
        if (!(ix.space() == space)) continue;
        for (std::size_t i = 0; i < ctx.size(); ++i)
          if (ctx[i].first == ix && i + 1 > deepest) deepest = i + 1;
      }
    return deepest;
  }

  /// \return true iff a router-directed fetch resolved to depth \p hd for
  ///         occurrence \p key is CONSISTENT with that occurrence: \p hd == 0
  ///         (the chain root -- nothing sliced) OR the live loop at \p hd
  ///         (\c ctx[hd-1]) is one of THIS occurrence's own in-scope batched
  ///         physical indices (\c key.get_indices()) whose space \p home's
  ///         DAG-scope names.
  ///
  /// \details The RELEASE-SAFE defense-in-depth criterion the Enter-stage
  /// router read gates on (design spec sec.4). By construction \c home_depth
  /// only ever returns such an \p hd, so on every correct schedule this HOLDS
  /// and the read proceeds byte-identically. It is not a no-op assert: if a
  /// FUTURE \c home_depth / overlay regression resolved an occurrence to a
  /// scope it does not bind (e.g. collapsing two divergently-relabeled
  /// occurrences onto ONE shared cache entry -- the exact hole the DAG-scope
  /// per-use resolution closes), this returns false and the read REFUSES the
  /// share, recomputing this occurrence's own value instead of serving a
  /// wrong-slice entry -- in RELEASE as well as debug. Pass the \p hd that \c
  /// home_depth returned and the SAME \p ctx / \p key.
  [[nodiscard]] bool home_resolution_consistent(HomeTarget const& home,
                                                BatchContext const& ctx,
                                                std::size_t hd,
                                                Key const& key) const {
    if (hd == 0) return true;
    if (hd > ctx.size()) return false;
    Index const& hm = ctx[hd - 1].first;
    auto const phys = key.template get_indices<container::svector<Index>>();
    if (std::find(phys.begin(), phys.end(), hm) == phys.end()) return false;
    for (auto const& space : home.dag_scope)
      if (hm.space() == space) return true;
    return false;
  }

 private:
  std::unordered_map<Key, HomeTarget, RouterKeyHash, RouterKeyEqual> overrides_;
  std::unordered_set<std::size_t> moved_hashes_;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_PLACEMENT_ROUTER_HPP
