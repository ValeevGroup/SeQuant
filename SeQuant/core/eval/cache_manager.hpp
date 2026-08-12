#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/expr.hpp>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <any>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sequant::eval {

// Forward declaration only (not a full include of placement_router.hpp):
// PlacementRouter's BatchContext alias needs CacheManager's full definition,
// so placement_router.hpp includes this header; CacheManager only needs a
// non-owning pointer to the (incomplete) PlacementRouter type, so the
// forward declaration here avoids the include cycle.
template <typename TreeNode>
class PlacementRouter;

/// \brief Destination for the structured schedule dump (SCHEDULE_RUN_EVENT
/// records) emitted by `evaluate()`. A caller sets one on the (root) cache to
/// capture one evaluation's batched schedule to \c os; \c fired is a fire-once
/// latch the caller raises after the first capture so a later re-entry (e.g. a
/// subsequent CC iteration) does not re-dump. Null sink / null \c os => no dump
/// (the default, byte-identical to before this seam). Non-owning: \c os must
/// outlive the cache.
struct ScheduleSink {
  std::ostream* os = nullptr;
  bool fired = false;
};

}  // namespace sequant::eval

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

  /// A shaped-product hook type. `evaluate()` consults the cache's shaped-
  /// product hook (if set) at each binary-Product node *before* the product is
  /// computed.  It receives the node (wrapped in a std::any as a
  /// std::reference_wrapper<key_type const>) plus the already-evaluated
  /// left/right operands and the [left, right, result] annotations.  It returns
  /// a non-null ResultPtr to *replace* the normal product (e.g. a shape-
  /// constrained emission of it), or a null ResultPtr to decline (the standard
  /// prod() then runs).  Empty (default) => never consulted; existing behavior
  /// is byte-identical.
  ///
  /// All backend-specific types (TA shapes, tranges, set_shape) stay inside the
  /// hook's closure (built by the backend, e.g. TAEvalContext::make_hook());
  /// the generic CacheManager and eval see only Result/ResultPtr/std::any.
  using shaped_product_hook_type = std::function<ResultPtr(
      std::any const& node, Result const& left, Result const& right,
      std::array<std::any, 3> const& annot)>;

  /// A whole-scope driver type. When set, the forest-range
  /// `evaluate(Nodes const&, layout, leaf, cache)` entrypoint (eval.hpp) routes
  /// the WHOLE forest through this driver instead of its per-tree descent,
  /// passing the forest and the result layout. The driver type-erases the
  /// `sequant::evaluate(forest, BatchPolicy, layout, leaf, cache, mode_order,
  /// make_scope_guard)` overload (scope_executor.hpp) -- which builds the
  /// batched-DAG schedule and calls `eval::evaluate_whole_scope` -- so that
  /// eval.hpp need not include scope_executor.hpp (that would be a cycle:
  /// scope_executor.hpp includes eval.hpp). The captured leaf evaluator,
  /// BatchPolicy, mode_order and scope-guard factory live inside the closure,
  /// which is built at the call site that owns those concrete types (e.g.
  /// MPQC's batched CSV-CCk residual install). The forest is a
  /// `std::vector<key_type>` and the layout a `std::string` (the residual
  /// annotation); other evaluate() instantiations never match and stay on the
  /// standard scheme. Empty (default) => no whole-scope routing; behavior is
  /// byte-identical.
  using whole_scope_driver_type =
      std::function<ResultPtr(std::vector<key_type> const& forest,
                              std::string const& layout, CacheManager&)>;

  /// The batch context: an ordered stack (outermost-first) of the enclosing
  /// realized batch loops, one entry per loop, `{axis K, {block_lo, block_hi}}`
  /// (element range). Set on the per-block scratch by the batched evaluator
  /// before it re-enters evaluate(); read by the Enter-stage slice-on-use so a
  /// cached intermediate fetched from an ancestor scope is sliced to the modes
  /// of the loops the fetch crossed (see eval.hpp). Empty (default) => no
  /// enclosing batch loop, so slice-on-use is inert and behavior is
  /// byte-identical to the pre-slice-on-use path.
  using BatchContext =
      container::svector<std::pair<Index, std::pair<std::size_t, std::size_t>>>;

  /// Result of access_at(): the fetched pointer plus the hop distance (number
  /// of parent links crossed) to the scope that held it. hops == 0 means a
  /// local hit; a null ptr carries hops == 0.
  struct AccessResult {
    ResultPtr ptr;
    std::size_t hops;
  };

  /// DIAGNOSTIC (dry-run costing): per-DISTINCT-value build tally for the
  /// avoidable-recompute rollup, keyed by the SAME node identity the cache
  /// dedups on (TreeNodeHasher + TreeNodeEqualityComparator = topological hash
  /// bin + Bliss connectivity 3-way cmp + recursive child compare), so two
  /// topologically-distinct nodes sharing a 64-bit hash are NOT folded and
  /// per-block / alpha-renamed builds of ONE value ARE folded.
  ///
  /// Recompute is measured with ACTUAL replay FLOPs, deduped at the (value,
  /// SLICE) granularity -- NOT against a build-once "full extent" denominator,
  /// which is ill-defined when slicing is non-uniform. \c slices maps a SLICE
  /// signature -- the enclosing batch context PROJECTED onto the modes THIS
  /// value actually carries (empty for a value invariant to every live loop) --
  /// to that slice's {build count, one build's actual cost}. Then:
  ///   total     = sum over slices of builds*cost   (== the replay's dryrun
  ///   sum) build-once = sum over slices of cost         (each DISTINCT slice
  ///   once) avoidable = sum over slices of (builds-1)*cost.
  /// A value tiled over DISTINCT slices (different blocks) has builds==1 per
  /// slice -> 0 avoidable (tiling is not recompute, even if the blocks are
  /// unequal). A value rebuilt at the SAME slice -- e.g. a node invariant to an
  /// enclosing loop, whose projected signature is identical every block -- has
  /// builds>1 at one slice -> (builds-1)*cost avoidable. Costs need not be
  /// uniform across slices; each slice carries its own realized cost.
  struct BuildRecord {
    std::size_t count = 0;  // number of builds of this exact (value, slice)
    double flops = 0;       // this slice's actual realized-extent FLOPs
    double exec = 0;        // this slice's actual roofline exec-cost estimate
  };
  struct BuildTally {
    std::unordered_map<std::string, BuildRecord> slices;
  };

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

    /// \return true iff this entry currently holds the SAME buffer (pointer
    ///         identity) as @p other. A sliced/permuted/phase-shifted read of
    ///         this entry is a DISTINCT buffer, so it compares unequal. Used by
    ///         the peak trace to detect an operand that aliases a cached buffer
    ///         (whose bytes are then already counted, and must not be added
    ///         again).
    [[nodiscard]] bool holds(ResultPtr const& other) const noexcept {
      return data_p && data_p.get() == other.get();
    }

    /// Upgrade this entry to an unbounded (resident-until-reset) non-persistent
    /// life, preserving any currently stored data and its cached size. Used to
    /// re-home an existing finite-life CSE entry at a scope where the value
    /// must survive ALL of its (possibly per-block-repeated) reads within one
    /// evaluation, so a partially drained entry becomes resident again rather
    /// than being freed and rebuilt by a later consumer. See \c
    /// CacheManager::ensure_home_slot.
    void make_resident() noexcept {
      max_life = std::numeric_limits<size_t>::max();
      life_c = std::numeric_limits<size_t>::max();
    }

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

  /// DIAGNOSTIC: per-DISTINCT-value build tally (see BuildTally), keyed by the
  /// same node identity as cache_map_. Populated by tally_build() from the eval
  /// loop's build choke point (eval.hpp finish_phase_b) for EVERY product
  /// build, whether that value is a cache entry, a footprint-gated recompute,
  /// or a per-batch rebuild -- so the rollup is complete. Held only on the
  /// scope- chain ROOT (tally_build routes there); scratch caches never
  /// populate it, and reset() does NOT clear it (the tally spans the whole
  /// forest replay).
  std::unordered_map<TreeNode, BuildTally, hasher_type, comparator_type>
      recompute_tally_;

  /// Gate for tally_build(): false (default) => tally_build is a no-op, so the
  /// wet (TA) eval path never populates recompute_tally_ and stays byte-
  /// identical. The dry-run costing replay (cost_profile) sets this true on the
  /// root cache before the replay. Held on the root only (tally_build routes
  /// there and checks it there).
  bool recompute_tally_enabled_ = false;

  /// Parent cache for the scope chain (loop-nest visibility). A batch scratch
  /// sets this to the cache one level up; access() delegates on a local miss
  /// so a loop-invariant node stored once at an ancestor level is found by
  /// every inner body without copy-down. Null (default) => standalone cache,
  /// byte-identical to pre-scope-chain behavior. Non-owning; the parent must
  /// outlive this cache.
  CacheManager* parent_ = nullptr;

  /// Running high-water mark (bytes) of the eval engine's live working set,
  /// updated by note_working_set() and cleared by reset(). Held here rather
  /// than in the recursive evaluate() so it persists across the whole
  /// evaluation of one term and is naturally reset between terms.
  size_t working_set_hwmark_ = 0;

  /// Optional custom evaluator consulted by evaluate() (see
  /// custom_evaluator_type). Empty => always defer to the standard scheme.
  custom_evaluator_type custom_evaluator_{};

  shaped_product_hook_type shaped_product_hook_{};

  /// Optional whole-scope driver consulted by the forest-range evaluate() (see
  /// whole_scope_driver_type). Empty => no whole-scope routing.
  whole_scope_driver_type whole_scope_driver_{};

  /// Enclosing realized batch loops for slice-on-use (see BatchContext). Empty
  /// (default) => no enclosing batch loop; the batched evaluator sets it on the
  /// per-block scratch before each re-entry. Not cleared by reset() (it is
  /// per-loop-iteration structural, re-set each block by the evaluator).
  BatchContext batch_context_{};

  /// Non-owning placement router (see \c placement_router.hpp). Null
  /// (default) => no override wired; \c placement_router() falls through to
  /// \c parent_ (only the root cache is wired in practice). The pointee must
  /// outlive this cache.
  eval::PlacementRouter<TreeNode> const* placement_router_ = nullptr;

  /// Non-owning hierarchy-wide co-resident high-water tracker (see
  /// \c eval::PeakMonitor). Null (default) => \c note_working_set() only
  /// updates this cache's own \c working_set_hwmark_; \c peak_monitor() falls
  /// through to \c parent_ (only the root cache is wired in practice). The
  /// pointee must outlive this cache.
  eval::PeakMonitor* peak_monitor_ = nullptr;

  /// Optional OWNING backing for \c placement_router_. A router built by a
  /// pre-pass (e.g. the remat placement pass) is a local at the build site; a
  /// CacheManager returned BY VALUE from such a builder must carry the router
  /// alive with it. \c adopt_placement_router stores it here (shared, so
  /// CacheManager stays copyable) and points \c placement_router_ at it. The
  /// forward-declared \c PlacementRouter is fine in a \c shared_ptr member (its
  /// deleter is type-erased at construction, where the type is complete).
  std::shared_ptr<eval::PlacementRouter<TreeNode> const> owned_router_{};

  /// Non-owning schedule-dump sink (see \c eval::ScheduleSink). Null (default)
  /// => `evaluate()` emits no SCHEDULE_RUN_EVENT records; falls through to
  /// \c parent_ (only the root cache is wired in practice). The pointee must
  /// outlive this cache.
  eval::ScheduleSink* schedule_sink_ = nullptr;

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

  /// Sets the shaped-product hook (see shaped_product_hook_).  Pass an empty
  /// std::function to clear it.
  void set_shaped_product_hook(shaped_product_hook_type fn) noexcept {
    shaped_product_hook_ = std::move(fn);
  }

  /// \return the shaped-product hook (empty if none is set).
  [[nodiscard]] shaped_product_hook_type const& shaped_product_hook()
      const noexcept {
    return shaped_product_hook_;
  }

  /// Sets the whole-scope driver (see whole_scope_driver_type). Pass an empty
  /// std::function to clear it.
  void set_whole_scope_driver(whole_scope_driver_type fn) noexcept {
    whole_scope_driver_ = std::move(fn);
  }

  /// \return the whole-scope driver (empty if none is set).
  [[nodiscard]] whole_scope_driver_type const& whole_scope_driver()
      const noexcept {
    return whole_scope_driver_;
  }

  /// Sets the batch context (see batch_context_). Pass an empty context to
  /// clear it.
  void set_batch_context(BatchContext c) noexcept {
    batch_context_ = std::move(c);
  }

  /// \return the batch context (empty if none is set).
  [[nodiscard]] BatchContext const& batch_context() const noexcept {
    return batch_context_;
  }

  /// Sets the scope-chain parent (see parent_). Pass nullptr to detach.
  void set_parent(CacheManager* p) noexcept { parent_ = p; }

  /// \return the scope-chain parent (see parent_), or nullptr if this is a
  ///         standalone / chain-root cache. Used by the batched evaluator to
  ///         walk up to a target ancestor level when hoisting an invariant.
  [[nodiscard]] CacheManager* parent() const noexcept { return parent_; }

  /// Sets the local placement router (see placement_router_). Pass nullptr
  /// to detach. Non-owning; the pointee must outlive this cache.
  void set_placement_router(eval::PlacementRouter<TreeNode> const* r) noexcept {
    placement_router_ = r;
  }

  /// Takes OWNERSHIP of a placement router (see owned_router_) and wires it as
  /// the local router. Use this when the router is a build-site local and the
  /// CacheManager is returned by value: the shared_ptr keeps the router alive
  /// for the cache's whole lifetime. Pass an empty shared_ptr to detach both.
  void adopt_placement_router(
      std::shared_ptr<eval::PlacementRouter<TreeNode> const> r) noexcept {
    owned_router_ = std::move(r);
    placement_router_ = owned_router_.get();
  }

  /// \return the local router if set, else the one inherited from parent_
  ///         (only the root cache is wired in practice); nullptr if none is
  ///         wired anywhere along the chain. Non-owning.
  [[nodiscard]] eval::PlacementRouter<TreeNode> const* placement_router()
      const noexcept {
    return placement_router_ ? placement_router_
           : parent_         ? parent_->placement_router()
                             : nullptr;
  }

  /// Sets the local peak monitor (see peak_monitor_). Pass nullptr to detach.
  /// Non-owning; the pointee must outlive this cache.
  void set_peak_monitor(eval::PeakMonitor* m) noexcept { peak_monitor_ = m; }

  /// \return the local peak monitor if set, else the one inherited from
  ///         \c parent_ (only the root cache is wired in practice); nullptr
  ///         if none is wired anywhere along the chain. Non-owning.
  [[nodiscard]] eval::PeakMonitor* peak_monitor() const noexcept {
    return peak_monitor_ ? peak_monitor_
           : parent_     ? parent_->peak_monitor()
                         : nullptr;
  }

  /// Sets the schedule-dump sink (see schedule_sink_). Pass nullptr to detach.
  /// Non-owning; the pointee (and its \c os) must outlive this cache.
  void set_schedule_sink(eval::ScheduleSink* s) noexcept { schedule_sink_ = s; }

  /// \return the local schedule sink if set, else the one inherited from
  ///         \c parent_ (only the root cache is wired in practice); nullptr if
  ///         none is wired anywhere along the chain. Non-owning.
  [[nodiscard]] eval::ScheduleSink* schedule_sink() const noexcept {
    return schedule_sink_ ? schedule_sink_
           : parent_      ? parent_->schedule_sink()
                          : nullptr;
  }

  /// Ensure a scope-hoist slot exists for @p key so a loop-invariant
  /// intermediate can be stored here (store() is a no-op for an unregistered
  /// key). The slot is NON-persistent with an effectively unbounded life, so it
  /// is never drained by access() and lives until the next reset() -- per-batch
  /// for a batch scratch (rebuilt for the next batch of the loop it is scoped
  /// to), per-term for the real cache (rebuilt for the next term). Idempotent:
  /// an existing entry (with any stored data) is left untouched. The unbounded
  /// life -- rather than the emitted effective_count -- is deliberate: a
  /// whole-nest invariant's escaped-outer set is empty, so its emitted
  /// effective_count is 1, which as a life would drain the entry on first use;
  /// reset() is the correct lifetime boundary for a hoisted invariant.
  void ensure_hoist_slot(key_type const& key) {
    cache_map_.try_emplace(
        key, entry{std::numeric_limits<size_t>::max(), /*persistent=*/false});
  }

  /// Ensure @p key has a RESIDENT home slot at THIS cache: an unbounded
  /// non-persistent life (like \c ensure_hoist_slot -- lives until the next
  /// \c reset()), so the value, once stored, is read by every consumer rather
  /// than drained and rebuilt. Unlike \c ensure_hoist_slot (which is a no-op
  /// on an existing entry), this UPGRADES an existing finite-life CSE entry to
  /// the same unbounded life via \c entry::make_resident, preserving any
  /// stored data. Used by the ordered executor (\c ordered_executor.hpp) to
  /// home a root-scope \c BuildStep value (its \c CellLegality::home_floor is
  /// empty -- a whole-nest invariant) at the root cache so it is built ONCE
  /// and read by every consumer, including a block-internal consumer reaching
  /// it through the scope chain -- the "de-alias the composite to root"
  /// property the ordered schedule's per-value homing exists to provide, and
  /// which the plain per-forest CSE life (drained after its unbatched use
  /// count) does not, since a realized batch loop reads a root-homed invariant
  /// once per block.
  void ensure_home_slot(key_type const& key) {
    auto [it, inserted] = cache_map_.try_emplace(
        key, entry{std::numeric_limits<size_t>::max(), /*persistent=*/false});
    if (!inserted) it->second.make_resident();
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
  /// per-op eval trace; monotonically non-decreasing until reset(). @p op_hash
  /// (default 0) identifies the op node being evaluated at the call site (0
  /// when no node is in scope there); forwarded to \c peak_monitor()'s
  /// \c observe() so a wired \c PeakMonitor can report WHERE its hierarchy-
  /// wide high-water was observed.
  size_t note_working_set(size_t current_bytes, size_t op_hash = 0) noexcept {
    working_set_hwmark_ = std::max(working_set_hwmark_, current_bytes);
    if (auto* m = peak_monitor()) {
      // DIAGNOSTIC (analysis-only): if a live-set capture hook is installed,
      // enumerate the chain's alive entries BEFORE observe() advances the mark,
      // on each real high-water advance. Gated on on_peak_liveset being set, so
      // the default path is byte-identical (no enumeration).
      if (m->on_peak_liveset && current_bytes > m->hwmark_bytes) {
        std::vector<eval::PeakLiveEntry> live;
        for (CacheManager const* c = this; c; c = c->parent_)
          for (auto const& [k, e] : c->cache_map_)
            if (e.alive()) live.push_back({k->hash_value(), e.size_in_bytes()});
        m->on_peak_liveset(current_bytes, live);
      }
      m->observe(current_bytes, op_hash);
    }
    // DIAGNOSTIC (SEQUANT_UT_PEAK_COMPOSE): on each new GLOBAL max working set,
    // print what composes it -- the co-resident cache chain vs the single
    // transient result/scratch being formed (current_bytes - chain_residency),
    // plus the largest single alive entry. Answers whether the realized peak is
    // cache-co-residency-bound or transient-working-set-bound. Env-gated, off
    // by default; harmless (a fprintf on monotone maxima only).
    static bool const compose =
        std::getenv("SEQUANT_UT_PEAK_COMPOSE") != nullptr;
    if (compose) {
      static size_t g_max = 0;
      if (current_bytes > g_max) {
        g_max = current_bytes;
        size_t const chain = chain_residency();
        size_t const transient =
            current_bytes > chain ? current_bytes - chain : 0;
        size_t max_entry = 0, n_alive = 0;
        for (CacheManager const* c = this; c; c = c->parent_)
          for (auto const& [k, e] : c->cache_map_)
            if (e.alive()) {
              ++n_alive;
              max_entry = std::max(max_entry, e.size_in_bytes());
            }
        std::fprintf(stderr,
                     "[peak-compose] max=%.1f GB = cache_chain %.1f + "
                     "transient(result) %.1f | n_alive_chain=%zu "
                     "max_single_alive=%.1f GB\n",
                     current_bytes / 1e9, chain / 1e9, transient / 1e9, n_alive,
                     max_entry / 1e9);
      }
    }
    return working_set_hwmark_;
  }

  /// Current running high-water mark (bytes) of the live working set.
  [[nodiscard]] size_t working_set_hwmark() const noexcept {
    return working_set_hwmark_;
  }

  /// DIAGNOSTIC: record one product build of @p key at slice @p slice_sig
  /// costing @p flops (this build's actual, realized-extent cost). @p slice_sig
  /// is the enclosing batch context projected onto the modes @p key carries
  /// (empty when the value is invariant to every live loop), so repeats of ONE
  /// slice fold (recompute) while distinct slices stay separate (tiling).
  /// Routes to the scope-chain ROOT so every build -- from any per-batch
  /// scratch -- accumulates in ONE map keyed by node identity (see
  /// recompute_tally_). Called only in the dry-run costing replay.
  void tally_build(key_type const& key, std::string const& slice_sig,
                   double flops, double exec) noexcept {
    if (parent_) {
      parent_->tally_build(key, slice_sig, flops, exec);
      return;
    }
    if (!recompute_tally_enabled_) return;  // wet path: no-op
    auto& slice = recompute_tally_[key].slices[slice_sig];
    slice.count += 1;     // one more build of this exact (value, slice)
    slice.flops = flops;  // this slice's actual cost (same for repeats)
    slice.exec = exec;    // this slice's actual exec-cost (same for repeats)
  }

  /// Enable/disable the per-node recompute tally (see
  /// recompute_tally_enabled_). Set on the root cache by the dry-run costing
  /// replay; left false everywhere else so tally_build() is a no-op on the wet
  /// eval path.
  void set_recompute_tally_enabled(bool on) noexcept {
    recompute_tally_enabled_ = on;
  }

  /// \return the per-DISTINCT-value build tally accumulated by tally_build()
  ///         on this (root) cache (see recompute_tally_). Read after the replay
  ///         to roll up avoidable recompute per node identity.
  [[nodiscard]] std::unordered_map<TreeNode, BuildTally, hasher_type,
                                   comparator_type> const&
  recompute_tally() const noexcept {
    return recompute_tally_;
  }

  /// Sum over ALIVE entries of this cache's own residency (bytes). Unlike
  /// working_set_hwmark() (a high-water MAX over time), this is the CURRENT
  /// live residency at the instant of the call.
  [[nodiscard]] size_t current_residency() const noexcept {
    size_t s = 0;
    for (auto const& [k, e] : cache_map_)
      if (e.alive()) s += e.size_in_bytes();
    return s;
  }
  /// current_residency() of this cache plus every ancestor along the scope
  /// chain (parent_): the total live residency visible at this scope at one
  /// instant.
  [[nodiscard]] size_t chain_residency() const noexcept {
    return current_residency() + (parent_ ? parent_->chain_residency() : 0);
  }

  /// \return true iff some ALIVE entry on this cache or any ancestor along the
  ///         scope chain physically holds @p value (pointer identity).
  ///         Read-only: unlike access_at() it decays no lifetime. The peak
  ///         trace uses it to skip an operand whose bytes are already counted
  ///         -- locally in \c bytes(cache,...) or up-chain in \c
  ///         chain_residency()
  ///         -- because the operand aliases that resident buffer. A sliced (or
  ///         permuted, or phase-shifted) read of a resident value is a DISTINCT
  ///         buffer with a different pointer, so it is correctly NOT skipped.
  [[nodiscard]] bool chain_holds(ResultPtr const& value) const noexcept {
    if (!value) return false;
    for (auto const& [k, e] : cache_map_)
      if (e.holds(value)) return true;
    return parent_ ? parent_->chain_holds(value) : false;
  }

  ///
  /// @brief Access cached data.
  ///
  /// @param key The key that identifies the cached data.
  /// @return the fetched pointer plus the hop distance (number of parent links
  ///         crossed) to the scope that held it; {nullptr, 0} on a total miss.
  ///
  /// A local entry only "hits" if it is currently holding data; a key
  /// registered locally but never (yet) stored here -- e.g. a hoisted
  /// loop-invariant node whose value lives only at an ancestor level -- is a
  /// local miss just like an unregistered key, and must fall through the
  /// same way. Standalone (parent_ == nullptr) behavior is unchanged: a total
  /// miss returns {nullptr, 0}. The hop distance surfaces the value's lifetime
  /// scope so the caller (Enter-stage slice-on-use) can slice it to exactly the
  /// batch loops the fetch crossed.
  [[nodiscard]] AccessResult access_at(key_type const& key) noexcept {
    if (auto found = cache_map_.find(key); found != cache_map_.end())
      if (auto data = found->second.access(); data) return {data, 0};
    if (!parent_) return {nullptr, 0};
    auto up = parent_->access_at(key);
    return {up.ptr, up.hops + 1};  // count the link we just crossed
  }

  /// @param key The key that identifies the cached data.
  /// @return ResultPtr to Result. Thin forwarder to access_at() that drops the
  ///         hop distance, for the non-batched callers that do not slice.
  ResultPtr access(key_type const& key) noexcept { return access_at(key).ptr; }

  /// Fetch @p key from EXACTLY @p hops scopes up the chain (walk @p hops
  /// parent links, then one LOCAL entry::access() there), rather than
  /// searching the chain like access_at() does. Used by a router-directed
  /// read that already knows the target scope (e.g. from a HomeTarget's
  /// home_depth) and wants to enforce that home rather than fall through to
  /// whatever scope happens to hold the value.
  ///
  /// @param key The key that identifies the cached data.
  /// @param hops The exact number of parent links to walk before accessing.
  /// @return the fetched pointer, decaying that scope's entry the same single
  ///         lifetime step as access_at() (entry::access()); nullptr if the
  ///         walk runs off the root before @p hops links, or if the target
  ///         scope does not currently hold @p key.
  [[nodiscard]] ResultPtr access_at_hops(key_type const& key,
                                         std::size_t hops) noexcept {
    CacheManager* c = this;
    for (std::size_t i = 0; i < hops && c; ++i) c = c->parent_;
    if (!c) return nullptr;
    if (auto found = c->cache_map_.find(key); found != c->cache_map_.end())
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
///
/// A node is also refused run-scope residence when it is BATCH-VARIANT: its
/// cross-occurrence lifetime mask is non-empty (\c !EvalExpr::mask_all_full();
/// this builder itself calls \c stamp_lifetime_masks over \p nodes before the
/// DAG walk below, so the mask is always current here regardless of caller --
/// see \c lifetime_mask.hpp). Such a node is sliced by some enclosing External
/// batch mode in every occurrence, so its value differs per batch of that mode
/// -- caching it whole at run scope would serve a wrong-batch value to a deeper
/// consumer on cache fall-through (the F1 hazard). Only an all-full node (empty
/// mask, including every node on the OFF path) is admitted.
/// \see CacheManager, cache_manager
template <bool force_hash_collisions = false,
          typename FootprintOf = zero_footprint>
auto cache_manager(meta::eval_node_range auto const& nodes, auto&& is_volatile,
                   size_t min_repeats = 2, FootprintOf footprint_of = {},
                   double max_footprint = 0.)
  requires requires(
      std::ranges::range_value_t<std::remove_cvref_t<decltype(nodes)>> const&
          n) {
    { is_volatile(n) } -> std::convertible_to<bool>;
    { footprint_of(n) } -> std::convertible_to<double>;
  }
{
  // Stamp the cross-occurrence lifetime mask on this SAME forest before the
  // DAG walk / veto below reads it (the batch-variant veto reads
  // EvalExpr::mask_all_full()). Doing this here -- rather than leaving it to
  // each caller -- makes "mask is current for the veto" an invariant of this
  // builder instead of a per-caller obligation: every caller of this overload
  // (SeQuant's build_dryrun_cache, mpqc's build_cache_manager) is covered
  // uniformly. Unconditional and idempotent; a no-op when the forest carries
  // no External batched_here() stamps (every mask stays empty/all-full), so
  // this never changes behavior on the OFF path.
  sequant::stamp_lifetime_masks(nodes);

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
  // Batch-variant veto ("a batched node cannot be run-scope"): this builder
  // populates the outermost / persistent (run-scope) cache, so it must refuse
  // any node that is batch-VARIANT -- one whose cached value would depend on
  // which batch is live. Such a node is refused for two reasons: caching it
  // whole contradicts the runtime slicing it (and the optimizer pricing it
  // sliced), and -- the F1 safety invariant -- a child batch scratch that
  // misses locally falls through to this cache for ANY key, so a batch-variant
  // final left here could be served full (wrong-batch) to an inner body. A node
  // is batch-variant iff its cross-occurrence lifetime mask is non-empty (\c
  // !n->mask_all_full(), \c lifetime_mask.hpp) -- some External batch mode (of
  // this node or an enclosing ancestor, over ALL its occurrences under the
  // canonical meet) slices it, so its value differs per batch of that mode even
  // if it slices nothing itself.
  // A node that is invariant to every batched mode is NOT vetoed and stays
  // cacheable at run scope -- this is where a hoisted loop-invariant
  // intermediate (all-full mask; or an External-only / no batched_here entry,
  // e.g. gC) lands. OFF path (no order-aware annotations, hence no \c
  // stamp_lifetime_masks External stamps): every mask is empty (all-full,
  // \c EvalExpr::sliced_modes_ default-constructed), so the veto never fires
  // and admits exactly what it did before -- byte-identical.
  std::unordered_map<TreeNode, size_t, Hasher, Comp> filtered;
  for (auto&& [n, c] : counts) {
    if (!(c >= min_repeats || persistent.contains(n))) continue;
    // Batch-variant: a node whose cross-occurrence lifetime mask is non-empty
    // is sliced by some enclosing external mode in every occurrence => its
    // value differs per batch => refused run-scope residence. all-full (empty
    // mask; incl. the OFF path) is admitted.
    bool const batch_variant = !n->mask_all_full();
    if (batch_variant ||
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
