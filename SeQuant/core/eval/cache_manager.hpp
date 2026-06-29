#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/expr.hpp>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>

namespace sequant {

///
/// This class implements a cache manager useful for the cases when the number
/// of times the cached objects will be accessed is known.
///
/// \tparam TreeNode The evaluation tree node type used as the cache key.
/// \tparam force_hash_collisions If true, forces all hash values to 0 (for
///         testing collision safety).
///
template <typename TreeNode, bool force_hash_collisions>
class CacheManager {
 public:
  using key_type = TreeNode;

  /// A custom evaluator type. `evaluate()` consults the cache's custom
  /// evaluator (if set) before applying its standard recursive scheme to each
  /// *non-leaf* node, invoking it as `custom_evaluator(node, cache)`:
  ///   - a non-null `ResultPtr` means "I evaluated this subtree myself" (e.g.
  ///     blocked over a contracted index to bound peak memory) and is used (and
  ///     cached) as-is;
  ///   - a null result means "decline", and the standard scheme evaluates the
  ///     node (its children are in turn consulted via the threaded cache).
  /// The leaf evaluator is captured by the callable. A custom evaluator that
  /// (re)evaluates the subtree by the standard scheme on a transformed operand
  /// set should do so on a *scratch* cache (e.g. `CacheManager::empty()`) to
  /// avoid both re-interception and polluting this cache with partial results.
  using custom_evaluator_type =
      std::function<ResultPtr(key_type const&, CacheManager&)>;

 private:
  using hasher_type = TreeNodeHasher<TreeNode, force_hash_collisions>;
  using comparator_type = TreeNodeEqualityComparator<TreeNode>;

  class entry {
   private:
    size_t max_life;

    size_t life_c;

    ResultPtr data_p;

    /// Lazily-computed, memoized data_p->size_in_bytes(): disengaged until the
    /// size is first queried, then computed once and cached until data_p is
    /// replaced or released. size_in_bytes() walks the DistArray's tiles, which
    /// is expensive; it is only queried to populate the eval trace's memory
    /// diagnostics, so for held data we compute it at most once and never when
    /// it is not asked for. mutable so size_in_bytes() can stay const (the
    /// classic lazy-memoization pattern).
    mutable std::optional<size_t> size_bytes_;

    /// Persistent (P) entries are never drained on access and survive reset(),
    /// so their data lives across multiple evaluations (e.g. CC iterations);
    /// non-persistent (NP) entries are released after their last use and
    /// cleared by reset() (the default, and historical, behavior).
    bool persistent_;

   public:
    explicit entry(size_t count, bool persistent = false) noexcept
        : max_life{count},
          life_c{count},
          data_p{nullptr},
          persistent_{persistent} {}

    [[nodiscard]] ResultPtr access() noexcept {
      if (!data_p) return nullptr;
      if (persistent_) return data_p;  // never drain a persistent entry
      if (decay() == 0) {              // last use: release the data
        size_bytes_.reset();
        return std::move(data_p);
      }
      return data_p;
    }

    void store(ResultPtr&& data) noexcept {
      data_p = std::move(data);
      size_bytes_
          .reset();  // (re)computed lazily on demand; see size_in_bytes()
    }

    void reset() noexcept {
      life_c = max_life;
      if (!persistent_) {  // persistent data (and its size) survives reset()
        data_p = nullptr;
        size_bytes_.reset();
      }
    }

    [[nodiscard]] bool persistent() const noexcept { return persistent_; }

    [[nodiscard]] size_t life_count() const noexcept { return life_c; }

    [[nodiscard]] size_t max_life_count() const noexcept { return max_life; }

    [[nodiscard]] size_t size_in_bytes() const noexcept {
      if (!data_p) return 0;
      if (!size_bytes_)
        size_bytes_ = data_p->size_in_bytes();  // lazy, memoized
      return *size_bytes_;
    }

    [[nodiscard]] bool alive() const noexcept { return data_p ? true : false; }

   private:
    [[nodiscard]] int decay() noexcept {
      return life_c > 0 ? static_cast<int>(--life_c) : 0;
    }

  };  // entry

  static ResultPtr store(entry& ent, ResultPtr&& data) noexcept {
    ent.store(std::move(data));
    return ent.access();
  }

  std::unordered_map<TreeNode, entry, hasher_type, comparator_type> cache_map_;

  /// Running high-water mark (bytes) of the eval engine's live working set,
  /// updated by note_working_set() and cleared by reset(). Held here rather
  /// than in the recursive evaluate() so it persists across the whole
  /// evaluation of one term and is naturally reset between terms.
  size_t working_set_hwmark_ = 0;

  /// Optional custom evaluator consulted by evaluate() (see
  /// custom_evaluator_type). Empty => always defer to the standard scheme.
  custom_evaluator_type custom_evaluator_{};

 public:
  /// Sets the custom evaluator (see custom_evaluator_type). Pass an empty
  /// std::function to clear it.
  void set_custom_evaluator(custom_evaluator_type fn) noexcept {
    custom_evaluator_ = std::move(fn);
  }

  /// \return the custom evaluator (empty if none is set).
  [[nodiscard]] custom_evaluator_type const& custom_evaluator() const noexcept {
    return custom_evaluator_;
  }
  /// Default persistence classifier: every entry is non-persistent (NP).
  struct all_non_persistent {
    bool operator()(key_type const&) const noexcept { return false; }
  };

  /// \param decaying iterable of (node, use-count) pairs to register for
  ///        caching.
  /// \param is_persistent predicate classifying each node as persistent (P,
  ///        never released on access, survives reset()) or non-persistent (NP,
  ///        released after its last use and by reset()). Defaults to all-NP.
  /// \note P nodes should be registered with whatever use-count is convenient
  ///       (it is not consulted for P entries) and may have a count of 1 even
  ///       though NP caching only registers nodes repeated min_repeats times.
  template <typename Iterable, typename PersistencePred = all_non_persistent>
    requires(!std::same_as<std::remove_cvref_t<Iterable>, CacheManager>)
  explicit CacheManager(Iterable&& decaying,
                        PersistencePred is_persistent = {}) noexcept {
    for (auto&& [k, c] : decaying)
      cache_map_.try_emplace(k, entry{c, is_persistent(k)});
  }

  ///
  /// Resets all cached data.
  ///
  void reset() noexcept {
    for (auto&& [k, v] : cache_map_) v.reset();
    working_set_hwmark_ = 0;
  }

  /// Fold the per-op live working set @p current_bytes into the running
  /// high-water mark and return the updated mark. Reported as `hw=` in the
  /// per-op eval trace; monotonically non-decreasing until reset().
  size_t note_working_set(size_t current_bytes) noexcept {
    working_set_hwmark_ = std::max(working_set_hwmark_, current_bytes);
    return working_set_hwmark_;
  }

  /// Current running high-water mark (bytes) of the live working set.
  [[nodiscard]] size_t working_set_hwmark() const noexcept {
    return working_set_hwmark_;
  }

  ///
  /// @brief Access cached data.
  ///
  /// @param key The key that identifies the cached data.
  /// @return ResultPtr to Result
  ResultPtr access(key_type const& key) noexcept {
    if (auto found = cache_map_.find(key); found != cache_map_.end())
      return found->second.access();
    return nullptr;
  }

  ///
  /// @param key The key to identify the cached data.
  /// @param data The data to be cached.
  /// \return Pointer to the stored data. Implictly accesses the stored data,
  ///         hence, decays the lifetime if the key accesses a decaying cache
  ///         entry. Passing @c key that was not present during construction of
  ///         this CacheManager object, stores nothing, but still returns a
  ///         valid pointer to @c data.
  [[nodiscard]] ResultPtr store(key_type const& key, ResultPtr data) noexcept {
    if (auto found = cache_map_.find(key); found != cache_map_.end())
      return store(found->second, std::move(data));
    return data;
  }

  ///
  /// \brief Check if the key exists in the database: does not check if cache
  ///        exists
  ///
  [[nodiscard]] bool exists(key_type const& key) const noexcept {
    return cache_map_.find(key) != cache_map_.end();
  }

  ///
  /// \brief Invokes \p fn with a const reference to every registered key.
  ///
  /// Keys are the canonical evaluation-tree nodes registered at construction;
  /// iteration order is unspecified. Use together with persistent()/alive() to
  /// enumerate, e.g., persistent entries that have not been populated yet.
  ///
  template <typename F>
    requires std::invocable<F&, key_type const&>
  void for_each_key(F&& fn) const {
    for (auto const& [k, v] : cache_map_) fn(k);
  }

  /// if the key exists in the database, return the current lifetime count of
  /// the cached data otherwise return -1
  [[nodiscard]] int life(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.life_count());
  }

  /// if the key exists in the database, return the maximum lifetime count of
  /// the cached data that implies the maximum number of accesses allowed for
  /// this key before the cache is released. This value was set by the c'tor.
  [[nodiscard]] int max_life(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.max_life_count());
  }

  /// \return true iff the key is registered for caching and currently holds
  ///         stored data (i.e. has been stored and not yet drained by its
  ///         final access).
  [[nodiscard]] bool alive(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && iter->second.alive();
  }

  /// \return true iff the key is registered for caching and classified
  ///         persistent (P: never released on access, survives reset()).
  [[nodiscard]] bool persistent(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && iter->second.persistent();
  }

  /// \return size in bytes of the data currently held for @p key, or 0 if
  ///         the key is not registered or no data is currently stored.
  [[nodiscard]] size_t entry_size_in_bytes(key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter == cache_map_.end() ? 0 : iter->second.size_in_bytes();
  }

  ///
  /// \return The number of entries with life_count greater than zero.
  ///
  [[nodiscard]] size_t alive_count() const noexcept {
    using ranges::views::filter;
    using ranges::views::transform;
    using ranges::views::values;
    return ranges::accumulate(cache_map_                            //
                                  | values                          //
                                  | filter(&entry::alive)           //
                                  | transform(&entry::life_count),  //
                              size_t{0});
  }

  ///
  /// \return Returns the sum of `Result::size_in_bytes` of alive entries.
  ///
  [[nodiscard]] size_t size_in_bytes() const noexcept {
    using ranges::views::transform;
    using ranges::views::values;
    return ranges::accumulate(
        cache_map_ | values | transform(&entry::size_in_bytes), size_t{0});
  }

  ///
  /// Get an empty cache manager.
  ///
  static CacheManager empty() noexcept {
    using map_type =
        std::unordered_map<TreeNode, size_t, hasher_type, comparator_type>;
    return CacheManager{map_type{}};
  }

  // for unit testing
  template <typename T>
  struct access_by;
  template <typename T>
  friend struct access_by;

};  // CacheManager

///
/// \brief Make a cache manager from an iterable of evaluable nodes.
///
/// \param nodes An iterable of eval nodes.
///
/// \param min_repeats Minimum number of repeats for a node to be cached. By
///                    default anything repeated twice or more is cached.
///
/// \return A cache manager.
///
/// \see CacheManager
///
template <bool force_hash_collisions = false>
auto cache_manager(meta::eval_node_range auto const& nodes,
                   size_t min_repeats = 2) noexcept {
  using TreeNode =
      std::ranges::range_value_t<std::remove_cvref_t<decltype(nodes)>>;
  using Hasher = TreeNodeHasher<TreeNode, force_hash_collisions>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;

  // Phase 1: Scan with pointer-based map (low memory)
  std::unordered_map<const TreeNode*, size_t, Hasher, Comp> imed_counts;

  auto imed_visitor = [&imed_counts](auto&& n) -> bool {
    if (auto found = imed_counts.find(&n); found != imed_counts.end()) {
      ++found->second;
      return false;
    }
    imed_counts.emplace(&n, 1);
    return true;
  };

  ranges::for_each(nodes, [&imed_visitor](auto&& tree) {
    tree.visit_internal(imed_visitor);
  });

  // Phase 2: Copy repeated entries (node by value)
  std::unordered_map<TreeNode, size_t, Hasher, Comp> filtered;
  for (auto&& [ptr, count] : imed_counts) {
    if (count >= min_repeats) filtered.emplace(*ptr, count);
  }

  return CacheManager<TreeNode, force_hash_collisions>{std::move(filtered)};
}

///
/// \brief Make a cache manager that distinguishes persistent (P) from
///        non-persistent (NP) intermediates, deriving persistence from a
///        solver-supplied volatility predicate and the evaluation DAG.
///
/// A node is *volatile* (V) if its value changes between evaluations (e.g. it
/// depends on the amplitudes being solved). \p is_volatile flags intrinsically
/// volatile nodes (typically the amplitude leaves); volatility is then
/// propagated up the DAG (a node is V iff it is intrinsically volatile or any
/// child is V). Persistence is derived from volatility and the consumer
/// (parent) relationship:
///
///   - V node                       -> NP (released after last use)
///   - NV node with >=1 V consumer   -> P  (the NV/V frontier: constant data
///                                          feeding per-iteration work; kept
///                                          across evaluations)
///   - NV node with no V consumer    -> NP (only feeds other NV nodes, so it is
///                                          absorbed into them and not needed
///                                          across evaluations)
///
/// Only *internal* (non-leaf) nodes are cached: NP nodes that repeat at least
/// \p min_repeats times (the usual CSE rule), plus *all* P nodes regardless of
/// repeat count (a P node is reused across evaluations even if used once each).
///
/// Default footprint accessor for cache_manager: reports zero footprint for
/// every node, so the footprint gate is inert (combined with max_footprint==0).
struct zero_footprint {
  double operator()(auto const&) const noexcept { return 0.; }
};

/// Default batchability predicate for cache_manager: no index is batchable, so
/// the free-batchable-axis caching veto is inert (preserves the pre-batch
/// behavior for callers that do not pass a predicate).
struct never_batchable {
  bool operator()(auto const&) const noexcept { return false; }
};

/// \param nodes the evaluation forest.
/// \param is_volatile `bool(TreeNode const&)`: true if the node is
///        intrinsically volatile. Only its value on leaves matters in practice
///        (volatility propagates up), but it is consulted on every node.
/// \param min_repeats minimum NP repeats to cache (default 2).
/// \param footprint_of `double(TreeNode const&)`: the materialized storage
///        footprint of a node's result (e.g. its element count or byte size).
///        Consulted only when \p max_footprint > 0.
/// \param max_footprint footprint gate: any node whose \p footprint_of exceeds
///        this is NOT cached (neither as an NP repeat nor as a P frontier
///        node), so it is recomputed by each consumer instead of being
///        materialized whole and held. This bounds the peak/sustained footprint
///        of huge intermediates that carry a free large-space index (e.g. a
///        half-transformed DF integral with a free projected-AO index), at the
///        cost of recomputation. 0 (default) disables the gate.
/// \param is_batchable_index `bool(Index const&)`: an index the runtime batched
///        evaluator slices over (typically the DF/RI auxiliary). A node whose
///        *result* (canonical) indices contain such an index carries a
///        batchable axis FREE: the evaluator slices it per batch and the
///        single-term optimizer prices it sliced, so caching it whole would
///        hold an intermediate both other components mean to slice. Such nodes
///        are NOT cached (neither NP repeat nor P frontier) -- recomputed
///        (sliced under each consumer's batch trigger) instead of materialized
///        whole and held. This is the structural counterpart of \p
///        max_footprint: the batch axis, not a byte threshold, identifies the
///        free-large-index intermediates. The default never_batchable accepts
///        nothing, leaving the veto inert.
/// \see CacheManager, cache_manager
template <bool force_hash_collisions = false,
          typename FootprintOf = zero_footprint,
          typename IsBatchableIndex = never_batchable>
auto cache_manager(meta::eval_node_range auto const& nodes, auto&& is_volatile,
                   size_t min_repeats = 2, FootprintOf footprint_of = {},
                   double max_footprint = 0.,
                   IsBatchableIndex is_batchable_index = {})
  requires requires(
      std::ranges::range_value_t<std::remove_cvref_t<decltype(nodes)>> const&
          n) {
    { is_volatile(n) } -> std::convertible_to<bool>;
    { footprint_of(n) } -> std::convertible_to<double>;
  }
{
  using TreeNode =
      std::ranges::range_value_t<std::remove_cvref_t<decltype(nodes)>>;
  using Hasher = TreeNodeHasher<TreeNode, force_hash_collisions>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;

  std::unordered_map<TreeNode, size_t, Hasher, Comp> counts;  // internal uses
  std::unordered_map<TreeNode, bool, Hasher, Comp> volatile_of;  // memoized
  std::unordered_set<TreeNode, Hasher, Comp> persistent;  // NV/V frontier

  // Single DAG walk: count internal-node uses (CSE), memoize volatility
  // bottom-up, and mark the NV/V frontier. Every (parent, child) edge is
  // visited exactly once (children are recursed only on a node's first visit),
  // so a child is marked persistent iff some volatile parent consumes it.
  auto visit = [&](auto&& self, TreeNode const& n) -> bool {
    bool const first = !volatile_of.contains(n);
    if (!n.leaf()) ++counts[n];  // count this use of an internal node
    if (!first) return volatile_of.at(n);
    bool v;
    if (n.leaf()) {
      v = is_volatile(n);
    } else {
      bool const vl = self(self, n.left());
      bool const vr = self(self, n.right());
      v = is_volatile(n) || vl || vr;
      if (v) {  // n is a volatile consumer => its NV internal children are P
        if (!vl && !n.left().leaf()) persistent.insert(n.left());
        if (!vr && !n.right().leaf()) persistent.insert(n.right());
      }
    }
    volatile_of.emplace(n, v);
    return v;
  };
  for (auto&& tree : nodes) visit(visit, tree);

  // Cache NP repeats + every P node; persistence = membership in `persistent`.
  // Footprint gate: a node whose result is larger than max_footprint is never
  // cached (so it is recomputed by each consumer rather than materialized whole
  // and held), bounding the footprint of huge free-large-index intermediates.
  // Free-batchable-axis veto: a node whose result carries an index the runtime
  // slices over (is_batchable_index) is, by construction, a free-large-index
  // intermediate the evaluator builds one batch-slice at a time and the
  // optimizer prices sliced. Caching it -- as an NP repeat or an NV/V-frontier
  // P node -- would materialize and hold it whole, contradicting both. Veto its
  // caching (the structural form of the max_footprint gate) so each consumer
  // recomputes it sliced under its own batch trigger.
  std::unordered_map<TreeNode, size_t, Hasher, Comp> filtered;
  for (auto&& [n, c] : counts) {
    if (!(c >= min_repeats || persistent.contains(n))) continue;
    bool free_batchable_axis = false;
    for (auto const& ix : n->canon_indices())
      if (is_batchable_index(ix)) {
        free_batchable_axis = true;
        break;
      }
    if (free_batchable_axis ||
        (max_footprint > 0. && footprint_of(n) > max_footprint)) {
      persistent.erase(n);  // keep is_persistent consistent with what is cached
      continue;
    }
    filtered.emplace(n, c);
  }

  auto is_persistent = [persistent = std::move(persistent)](TreeNode const& n) {
    return persistent.contains(n);
  };
  return CacheManager<TreeNode, force_hash_collisions>{
      std::move(filtered), std::move(is_persistent)};
}

///
/// \brief Estimates the peak memory required to hold the intermediates that
///        repeat when a Sum is evaluated term by term.
/// \note Reordering the terms in a Sum affects the peak cache memory.
///
/// \param expr A Sum whose terms will be evaluated by reusing intermediates.
/// \param min_repeats Minimum number of repeats for a node to be cached. If not
/// provided, will use the default of \c cache_manager().
/// \return AsyCost object representing the peak memory as a polynomial in
///         the index-space sizes of the stored tensors.
///
AsyCost peak_cache(Sum const& expr,
                   std::optional<size_t> min_repeats = std::nullopt);

}  // namespace sequant

#endif  // SEQUANT_EVAL_CACHE_MANAGER_HPP
