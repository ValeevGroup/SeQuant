#ifndef SEQUANT_EVAL_CACHE_MANAGER_HPP
#define SEQUANT_EVAL_CACHE_MANAGER_HPP

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/peak_monitor.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/value_id.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <any>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
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

/// \brief DIAGNOSTIC (analysis-only, OFF by default): a global monotonic
///        "access clock" stamped on every genuine cache READ, plus a
///        per-value (keyed by canonical node hash) record of the LAST clock at
///        which that value was read.
///
/// Home (root-homed) cache entries are pinned (life_c == SIZE_MAX), so their
/// lifetime counter never drains and cannot be used to infer their genuine last
/// use. This clock records the real thing: \c CacheManager::access_at and
/// \c access_at_hops (the ONLY genuine consumer-read paths -- the store-return
/// \c entry::access() bypasses both) call \c tick() and stamp the read value's
/// hash into \c last_access_map on every hit when \c enabled(). The record is a
/// GLOBAL map keyed by node hash (not a per-entry field) deliberately: a
/// batch-loop tier-B value lives on a per-block scratch cache that is destroyed
/// at block close, so a per-entry field would be lost; the global map keeps the
/// value's FINAL last-read clock across the whole run regardless of which
/// (root or transient scratch) scope held it. Reset by the harness before a
/// measured run. Single-threaded dry-run only. When \c enabled() is false every
/// stamp site is a no-op and the eval path stays byte-identical.
struct AccessClock {
  /// One-shot env gate (SEQUANT_UT_ACCESS_CLOCK). Read once; when unset every
  /// stamp site below is inert.
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_ACCESS_CLOCK") != nullptr;
    return on;
  }
  static std::size_t& counter() noexcept {
    static std::size_t c = 0;
    return c;
  }
  /// hash -> final (max) clock at which a value with that hash was read.
  static std::unordered_map<std::size_t, std::size_t>&
  last_access_map() noexcept {
    static std::unordered_map<std::size_t, std::size_t> m;
    return m;
  }
  /// Advance and return the clock (one genuine read == one tick).
  static std::size_t tick() noexcept { return ++counter(); }
  /// Current clock value WITHOUT advancing (used to timestamp the peak).
  static std::size_t now() noexcept { return counter(); }
  /// Record a genuine read of the value with hash @p h at a fresh clock tick.
  static void stamp(std::size_t h) noexcept {
    if (!enabled()) return;
    last_access_map()[h] = tick();
  }
  /// Clear the clock and the per-value record before a measured run.
  static void reset() noexcept {
    counter() = 0;
    last_access_map().clear();
  }
};

/// \brief DIAGNOSTIC (analysis-only, OFF by default): counts and times every
///        `cache_map_.find(key)` performed by CacheManager. Each such find runs
///        TreeNodeEqualityComparator on any bucket match -- a recursive
///        structural compare of the whole subtree (memoized-hash O(1) per node
///        + linear bliss ConstGraphCmp on the connectivity graph, recursed into
///        left/right). This meter isolates the aggregate wall time that lookup
///        costs, so a schedule that issues more read-from-home lookups
///        (ordered) can be compared against one that inlines (forest). Env gate
///        SEQUANT_UT_LOOKUP_METER; when unset every site is a no-op passthrough
///        and the eval path stays byte-identical. Single-threaded runs only.
struct LookupMeter {
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_LOOKUP_METER") != nullptr;
    return on;
  }
  static std::size_t& calls() noexcept {
    static std::size_t c = 0;
    return c;
  }
  static std::uint64_t& nanos() noexcept {
    static std::uint64_t n = 0;
    return n;
  }
  /// Time ONLY the find() call in @p f (a nullary functor returning the
  /// iterator); accumulate count + elapsed ns; return f()'s result. Works for
  /// const and non-const maps alike (return type is deduced).
  template <typename F>
  static auto timed(F&& f) {
    if (!enabled()) return f();
    auto const t0 = std::chrono::steady_clock::now();
    auto r = f();
    auto const t1 = std::chrono::steady_clock::now();
    ++calls();
    nanos() += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    return r;
  }
  static void reset() noexcept {
    calls() = 0;
    nanos() = 0;
  }
  /// Prints the accumulated totals to stderr at program teardown when enabled.
  struct Reporter {
    ~Reporter() {
      if (!enabled()) return;
      std::fprintf(stderr,
                   "[lookup-meter] cache_map_.find calls=%zu total_ms=%.3f "
                   "mean_ns=%.1f\n",
                   calls(), static_cast<double>(nanos()) / 1e6,
                   calls() ? static_cast<double>(nanos()) / calls() : 0.0);
    }
  };
  inline static Reporter reporter_{};
};

/// \brief DIAGNOSTIC (analysis-only, OFF by default): scoped wall-clock timer
///        that accumulates elapsed time into named buckets, for locating where
///        the (single-threaded) ordered-executor driver spends time. Times are
///        INCLUSIVE (a region's time includes any nested regions and any TA
///        dispatch inside it), so compare siblings and drill into the hot one.
///        Env gate SEQUANT_UT_PHASE; inert when unset. Single-threaded only.
struct PhaseTimer {
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_PHASE") != nullptr;
    return on;
  }
  // Backend-supplied barrier that drains all pending async work (set by the TA
  // backend to a world fence). When present, each Scope fences at BOTH
  // boundaries so async work cannot leak across timer regions and every
  // region's time reflects the work it actually dispatched-and-completed.
  static std::function<void()>& fence_hook() noexcept {
    static std::function<void()> h;
    return h;
  }
  static std::size_t& barrier_count() noexcept {
    static std::size_t c = 0;
    return c;
  }
  static std::uint64_t& barrier_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  static void barrier() noexcept {
    if (auto const& h = fence_hook()) {
      auto const t0 = std::chrono::steady_clock::now();
      h();
      barrier_ns() += static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(
              std::chrono::steady_clock::now() - t0)
              .count());
      ++barrier_count();
    }
  }
  // name -> {total ns, call count}. Leaked so teardown Reporter can read it.
  static std::map<std::string, std::pair<std::uint64_t, std::uint64_t>>&
  acc() noexcept {
    static auto* const m =
        new std::map<std::string, std::pair<std::uint64_t, std::uint64_t>>();
    return *m;
  }
  struct Scope {
    char const* name;
    std::chrono::steady_clock::time_point t0;
    bool on;
    explicit Scope(char const* n) noexcept : name(n), on(enabled()) {
      if (on) t0 = std::chrono::steady_clock::now();
    }
    ~Scope() {
      if (!on) return;
      barrier();  // dtor-only: drain the async WE dispatched before stopping
      auto const d = std::chrono::duration_cast<std::chrono::nanoseconds>(
                         std::chrono::steady_clock::now() - t0)
                         .count();
      auto& e = acc()[name];
      e.first += static_cast<std::uint64_t>(d);
      ++e.second;
    }
  };
  struct Reporter {
    ~Reporter() {
      if (!enabled()) return;
      std::fprintf(stderr, "[phase-timer] (inclusive wall, name: s / calls)\n");
      for (auto const& [n, e] : acc())
        std::fprintf(stderr, "  %-28s %8.3f s  %8llu\n", n.c_str(),
                     e.first / 1e9, static_cast<unsigned long long>(e.second));
    }
  };
  inline static Reporter reporter_{};
};

/// \brief DIAGNOSTIC (analysis-only, OFF by default): splits CC-eval wall into
///        (a) time INSIDE top-level evaluate_impl (the "body") and (b) the GAP
///        between successive top-level evaluate_impl calls (the ordered
///        executor's schedule-walk / block-structure / call-site machinery).
///        Fences at body exit (via PhaseTimer::barrier) so each call owns its
///        async and the gap is pure between-call executor time. Depth-guarded:
///        nested evaluate_impl re-entries (forest's custom evaluator) fold into
///        the enclosing top-level call. Env gate SEQUANT_UT_EVALIMPL.
struct EvalImplTimeline {
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_EVALIMPL") != nullptr;
    return on;
  }
  static int& depth() noexcept {
    static thread_local int d = 0;
    return d;
  }
  static std::uint64_t& body_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  static std::uint64_t& gap_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  static std::size_t& calls() noexcept {
    static std::size_t v = 0;
    return v;
  }
  static std::chrono::steady_clock::time_point& last_exit() noexcept {
    static std::chrono::steady_clock::time_point t{};
    return t;
  }
  static bool& have_last() noexcept {
    static bool b = false;
    return b;
  }
  static std::uint64_t& max_gap_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  // Node-evaluation (prod / TA contraction) time INSIDE the body, fenced so it
  // captures the contraction's full sync+async cost. body - prod == "other
  // inside evaluate_impl" (stack machine, cache access, apply_phase, store).
  static std::uint64_t& prod_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  static void note_prod(std::chrono::steady_clock::time_point t0) noexcept {
    if (!enabled()) return;
    PhaseTimer::barrier();  // drain this contraction's async before stopping
    prod_ns() += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0)
            .count());
  }
  // apply_phase (canonicalization-phase tensor multiply) time inside the body:
  // node-eval-adjacent TA work that is NOT the contraction.
  static std::uint64_t& phase_ns() noexcept {
    static std::uint64_t v = 0;
    return v;
  }
  static void note_phase(std::chrono::steady_clock::time_point t0) noexcept {
    if (!enabled()) return;
    PhaseTimer::barrier();
    phase_ns() += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0)
            .count());
  }
  struct Scope {
    std::chrono::steady_clock::time_point entry;
    bool on;
    bool top;
    Scope() noexcept : on(enabled()), top(false) {
      if (!on) return;
      top = (depth() == 0);
      ++depth();
      if (top) {
        entry = std::chrono::steady_clock::now();
        if (have_last()) {
          auto const g = std::chrono::duration_cast<std::chrono::nanoseconds>(
                             entry - last_exit())
                             .count();
          gap_ns() += static_cast<std::uint64_t>(g);
          if (static_cast<std::uint64_t>(g) > max_gap_ns())
            max_gap_ns() = static_cast<std::uint64_t>(g);
        }
      }
    }
    ~Scope() {
      if (!on) return;
      --depth();
      if (!top) return;
      PhaseTimer::barrier();  // drain this body's async so the gap is pure
      auto const e = std::chrono::steady_clock::now();
      body_ns() += static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(e - entry)
              .count());
      last_exit() = e;
      have_last() = true;
      ++calls();
    }
  };
  struct Reporter {
    ~Reporter() {
      if (!enabled()) return;
      std::fprintf(stderr,
                   "[evalimpl-timeline] top_calls=%zu  body=%.3f s  "
                   "gap(between-calls)=%.3f s  max_gap=%.1f ms\n",
                   calls(), body_ns() / 1e9, gap_ns() / 1e9,
                   max_gap_ns() / 1e6);
      auto const accounted = prod_ns() + phase_ns();
      std::fprintf(
          stderr,
          "[evalimpl-timeline]   of body: node-eval(prod)=%.3f s  "
          "apply_phase=%.3f s  other-inside=%.3f s\n",
          prod_ns() / 1e9, phase_ns() / 1e9,
          (body_ns() > accounted ? (body_ns() - accounted) / 1e9 : 0.0));
      std::fprintf(stderr,
                   "[evalimpl-timeline]   instrumentation barrier(): calls=%zu "
                   "total=%.3f s (this is fence overhead, NOT real work)\n",
                   PhaseTimer::barrier_count(), PhaseTimer::barrier_ns() / 1e9);
    }
  };
  inline static Reporter reporter_{};
};

/// \brief DIAGNOSTIC (analysis-only, OFF by default): counts ACTUAL builds at
///        the single build chokepoint (finish_phase_b) -- every freshly
///        computed non-leaf node, cached or not. Unlike DefUseMeter (which
///        counts cache.store() calls, and so also counts re-home/placement
///        stores of an already-built value), this counts real contraction
///        executions, so build_count > iterations for a node == genuine
///        recompute (re-run contraction), the ground truth the FLOP count
///        reflects. Env gate SEQUANT_UT_BUILD_METER. Single-threaded only.
struct BuildMeter {
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_BUILD_METER") != nullptr;
    return on;
  }
  static std::unordered_map<std::size_t, std::size_t>& build_count() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
    return *m;
  }
  static std::size_t& total() noexcept {
    static std::size_t t = 0;
    return t;
  }
  static std::unordered_map<std::size_t, std::string>& label_of() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::string>();
    return *m;
  }
  static std::unordered_map<std::size_t, std::size_t>& read_count() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
    return *m;
  }
  static std::unordered_map<std::size_t, std::size_t>& max_life_of() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
    return *m;
  }
  /// Record one genuine cache read (access_at hit) of @p hash and its entry's
  /// @p max_life -- ground truth for "is the lifetime >= the real read count?".
  static void on_read(std::size_t hash, std::size_t max_life) noexcept {
    if (!enabled()) return;
    ++read_count()[hash];
    max_life_of()[hash] = max_life;
  }
  /// Record one fresh build of the non-leaf node @p hash (optional human
  /// @p label, kept from the first sighting for the dump).
  static void on_build(std::size_t hash,
                       std::string label = std::string{}) noexcept {
    if (!enabled()) return;
    ++total();
    ++build_count()[hash];
    if (!label.empty()) label_of().try_emplace(hash, std::move(label));
  }
  struct Reporter {
    ~Reporter() {
      if (!enabled()) return;
      std::map<std::size_t, std::size_t> hist;
      for (auto const& [h, c] : build_count()) ++hist[c];
      std::fprintf(stderr, "[build-meter] builds=%zu distinct=%zu", total(),
                   build_count().size());
      std::fprintf(stderr, " histogram(k:nodes):");
      for (auto const& [k, n] : hist) std::fprintf(stderr, " %zu:%zu", k, n);
      std::fprintf(stderr, "\n");
      if (char const* dump = std::getenv("SEQUANT_UT_BUILD_DUMP")) {
        if (std::FILE* fp = std::fopen(dump, "w")) {
          for (auto const& [h, c] : build_count()) {
            auto const it = label_of().find(h);
            auto const rit = read_count().find(h);
            auto const mit = max_life_of().find(h);
            std::fprintf(fp, "%zu\t%zu\t%zu\t%zu\t%s\n", h, c,
                         rit != read_count().end() ? rit->second : 0,
                         mit != max_life_of().end() ? mit->second : 0,
                         it != label_of().end() ? it->second.c_str() : "");
          }
          std::fclose(fp);
        }
      }
    }
  };
  inline static Reporter reporter_{};
};

/// \brief DIAGNOSTIC (analysis-only, OFF by default): measures def-to-use
///        distance -- how many production "ops" elapse between a value being
///        stored (produced by a contraction) and each time it is consumed. A
///        global op-clock ticks once per genuine \c store(key,data); each value
///        records its production tick, and every \c access_at hit accumulates
///        (current_clock - production_tick), weighted by the operand's bytes.
///        The byte-weighted aggregate proxies CPU-cache reload pressure: a
///        large operand that waited many ops between production and use is
///        likely evicted from cache and must be reloaded from DRAM when finally
///        read. Uniform across schedulers (both route produced values through
///        store() and reads through access_at). Robust to forest's
///        immediate-consume bypass: a distance-0 consumption contributes ~0 to
///        the byte-weighted sum. Env gate SEQUANT_UT_DEFUSE_METER; inert when
///        unset. Single-threaded runs only.
struct DefUseMeter {
  static bool enabled() noexcept {
    static bool const on = std::getenv("SEQUANT_UT_DEFUSE_METER") != nullptr;
    return on;
  }
  static std::size_t& op_clock() noexcept {
    static std::size_t c = 0;
    return c;
  }
  // Leaked (never destroyed) so the teardown Reporter can safely read them
  // regardless of static-destruction order.
  static std::unordered_map<std::size_t, std::size_t>& produce_at() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
    return *m;
  }
  /// distinct-node -> number of store (production) events, exposing recompute:
  /// a node stored more than (iterations) times was rebuilt while cacheable.
  static std::unordered_map<std::size_t, std::size_t>& store_count() noexcept {
    static auto* const m = new std::unordered_map<std::size_t, std::size_t>();
    return *m;
  }
  static std::size_t& reads() noexcept {
    static std::size_t r = 0;
    return r;
  }
  static std::uint64_t& sum_dist() noexcept {
    static std::uint64_t s = 0;
    return s;
  }
  static long double& sum_dist_bytes() noexcept {
    static long double s = 0;
    return s;
  }
  static long double& sum_bytes() noexcept {
    static long double s = 0;
    return s;
  }
  static std::size_t& max_dist() noexcept {
    static std::size_t m = 0;
    return m;
  }
  /// Tick the op-clock and record this value's production time.
  static void on_store(std::size_t hash) noexcept {
    if (!enabled()) return;
    produce_at()[hash] = ++op_clock();
    ++store_count()[hash];
  }
  /// Accumulate the def-to-use distance for a consumption of @p hash whose
  /// current stored size is @p bytes. Values never stored (leaves/inputs) are
  /// skipped -- their production op is undefined.
  static void on_read(std::size_t hash, std::size_t bytes) noexcept {
    if (!enabled()) return;
    auto const it = produce_at().find(hash);
    if (it == produce_at().end()) return;
    std::size_t const d = op_clock() - it->second;
    ++reads();
    sum_dist() += d;
    sum_dist_bytes() += static_cast<long double>(d) * bytes;
    sum_bytes() += bytes;
    if (d > max_dist()) max_dist() = d;
  }
  struct Reporter {
    ~Reporter() {
      if (!enabled()) return;
      double const mean =
          reads() ? static_cast<double>(sum_dist()) / reads() : 0.0;
      double const meanw =
          sum_bytes() > 0 ? static_cast<double>(sum_dist_bytes() / sum_bytes())
                          : 0.0;
      std::fprintf(stderr,
                   "[defuse-meter] ops=%zu distinct=%zu reads=%zu "
                   "mean_dist=%.2f byteweighted_mean_dist=%.2f max_dist=%zu "
                   "sum_dist=%llu sum_dist_bytes=%.4e\n",
                   op_clock(), produce_at().size(), reads(), mean, meanw,
                   max_dist(), static_cast<unsigned long long>(sum_dist()),
                   static_cast<double>(sum_dist_bytes()));
      // Store-count histogram: how many distinct nodes were stored k times.
      // k > iterations => recompute (a cacheable node rebuilt while resident).
      std::map<std::size_t, std::size_t> hist;
      for (auto const& [h, c] : store_count()) ++hist[c];
      std::fprintf(stderr, "[defuse-meter] store-count histogram (k:nodes):");
      for (auto const& [k, n] : hist) std::fprintf(stderr, " %zu:%zu", k, n);
      std::fprintf(stderr, "\n");
      // Optional: dump the distinct stored-node hash set (hash store_count),
      // one per line, so two runs' cache-membership sets can be diffed.
      if (char const* const path = std::getenv("SEQUANT_UT_DEFUSE_DUMP")) {
        if (std::FILE* f = std::fopen(path, "w")) {
          for (auto const& [h, c] : store_count())
            std::fprintf(f, "%zu %zu\n", h, c);
          std::fclose(f);
        }
      }
    }
  };
  inline static Reporter reporter_{};
};

}  // namespace sequant::eval

namespace sequant {

/// \brief One entry of a \c CacheManager::BatchContext: the enclosing realized
///        batch loop's axis, DAG-scope nest position, and element range.
///
/// \details A free (non-template) type -- like \c DagScopeLevel / \c
/// ModeToLevel (\c dag_scope.hpp), which it embeds -- since none of its
/// fields depend on the owning \c CacheManager's \c TreeNode /
/// \c force_hash_collisions template parameters.
///
/// \c axis is the KEPT field: the OLD resolution path (\c evaluate_impl's
/// Enter-stage \c slice_to_use, \c PlacementRouter::home_depth /
/// \c home_resolution_consistent) still looks a fetched node's mode up by
/// this exact \c Index (via \c index_position(nd, axis), with a
/// space-mapped fallback for the ordered executor) -- unchanged by this
/// struct's introduction. \c level is the NEW DAG-scope nest position a later
/// task (Task 7) will resolve against instead. \c range is the loop's current
/// element block (renamed from the old pair's \c second, same meaning).
/// \c exact_axis is NEW: the ordered executor (which co-evaluates a whole
/// type-bucketed block under one canonical \c block.axis) leaves it
/// \c nullopt, while the forest / whole-scope evaluator (which pushes each
/// member's own physical axis) fills it with that same axis -- redundant with
/// \c axis today, but keeping the two fields separate lets Task 7 tell the two
/// resolution strategies apart without re-deriving which one produced this
/// entry.
///
/// Purely ADDITIVE relative to the old `std::pair<Index,
/// std::pair<std::size_t, std::size_t>>`: every existing reader that used to
/// destructure `{axis, {lo, hi}}` continues to compile once updated to name
/// `.axis` / `.range` instead of `.first` / `.second` (a mechanical rename;
/// see the callers this task updates), and every writer still populates
/// exactly the same `{axis, {lo, hi}}` pair, now alongside `level` and
/// `exact_axis`.
struct BatchContextEntry {
  Index axis;
  DagScopeLevel level{};
  std::pair<std::size_t, std::size_t> range{};
  std::optional<Index> exact_axis{};
};

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
  /// The NODE type. Used for everything node-facing: the custom evaluator, the
  /// whole-scope / multiroot drivers, the persistence predicate, and
  /// \c for_each_key -- all of which see forest nodes, unchanged by Pillar 1.
  using key_type = TreeNode;

  /// Pillar 1 (slice-colored value identity): the CACHE-MAP key. A
  /// \c CachedValue pairs a node with its recorded home-slice coloring, so two
  /// values of one node sliced on different slots occupy distinct entries; an
  /// unsliced value (empty coloring) is byte-identical to keying by the node.
  /// \c CachedValue is IMPLICITLY constructible from a node, so every existing
  /// call site that passes a bare node keeps compiling and keys the map exactly
  /// as before; only the ordered executor passes a genuinely colored key.
  using cache_key_type = eval::CachedValue<TreeNode>;

  /// Pillar 1 / Task 7: a per-build coloring context -- node hash -> the
  /// home-slice coloring of the VALUE that node represents in the current
  /// build. The ordered executor sets one (RAII) before each `evaluate_impl`
  /// build, from `operand_vids` + `home_mode_depth`; the map-keying methods
  /// re-color an EMPTY-coloring key from it so `evaluate_impl`'s bare-node
  /// store/access build the right home-colored `CachedValue`. Per-build =>
  /// node->coloring is 1:1.
  using ValueColoringCtx =
      std::unordered_map<std::size_t, eval::ValueIdColoring>;

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

  /// A multi-root driver type. When set, the free function
  /// `evaluate_multiroot(roots, layouts, leaf, cache)` (eval.hpp) routes the
  /// WHOLE set of independent \p roots through this driver instead of
  /// evaluating each root as its own separate forest. Unlike
  /// `whole_scope_driver_type`, which SUMS its forest's roots into one
  /// `ResultPtr`, this driver returns a MAP: one result per root, in input
  /// order, with NO cross-root summation -- the roots need not even be
  /// commensurate in shape (e.g. independent CC residual equations). \p
  /// layouts holds one layout string per root, in the same order as \p
  /// roots, so heterogeneous roots (e.g. distinct CC residual annotations
  /// like R1 `{a;i}` vs R2 `{ab;ij}`) can each be permuted to their own
  /// result layout. The intended installer (ordered_executor.hpp)
  /// concatenates \p roots into ONE schedule so a subexpression shared
  /// across roots is built once (the same CSE `compute_dag_boulevard`
  /// already gives a concatenated node list), then returns each root's own
  /// (unsummed) `value_result`. The captured leaf evaluator and any
  /// batching policy live inside the closure, exactly as
  /// `whole_scope_driver_type` documents. Empty (default) =>
  /// `evaluate_multiroot` throws (there is no per-root fallback: unlike
  /// whole-scope routing, which falls back to the ordinary per-tree descent
  /// when unset, a multi-root caller MUST explicitly wire a driver that
  /// understands the cross-root CSE contract).
  using multiroot_driver_type = std::function<container::svector<ResultPtr>(
      container::svector<key_type> const& roots,
      container::svector<std::string> const& layouts, CacheManager&)>;

  /// The batch context: an ordered stack (outermost-first) of the enclosing
  /// realized batch loops, one entry per loop (see \c BatchContextEntry: axis
  /// K, its DAG-scope level, and the current `{block_lo, block_hi}` element
  /// range). Set on the per-block scratch by the batched evaluator before it
  /// re-enters evaluate(); read by the Enter-stage slice-on-use so a cached
  /// intermediate fetched from an ancestor scope is sliced to the modes of the
  /// loops the fetch crossed (see eval.hpp). Empty (default) => no enclosing
  /// batch loop, so slice-on-use is inert and behavior is byte-identical to
  /// the pre-slice-on-use path.
  using BatchContext = container::svector<BatchContextEntry>;

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
  using hasher_type = eval::CachedValueHasher<TreeNode, force_hash_collisions>;
  using comparator_type = eval::CachedValueEqual<TreeNode>;

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

    /// Regression tripwire (see \c store()): true iff a NON-persistent entry
    /// has already been \c store()'d since its last \c reset(). After the
    /// combined single-DAG evaluation (every value built exactly once per
    /// evaluation), a second \c store() into the same NP entry with no
    /// intervening \c reset() means a duplicate producer survived -- a bug.
    /// Not consulted (or set) for PERSISTENT entries, which legitimately
    /// re-store across batch replays. Cleared by \c reset().
    bool stored_this_eval_ = false;

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
        // The value has now been fully consumed and released, so a later
        // store() of this key is a fresh RE-PRODUCTION (e.g. the CCk energy
        // observable, evaluated on the residual's cache, recomputing a volatile
        // whose residual-side lifetime just drained), not a duplicate producer
        // within one eval. Clear the re-store tripwire so that legitimate
        // re-production is allowed; a true duplicate (two store()s with no
        // access draining the value in between) leaves the flag set and still
        // trips.
        stored_this_eval_ = false;
        return std::move(data_p);
      }
      return data_p;
    }

    // NOT noexcept: SEQUANT_ASSERT below can throw (SEQUANT_ASSERT_BEHAVIOR
    // = THROW) -- see the tripwire comment on stored_this_eval_. In every
    // default build config (IGNORE, or Debug's ABORT) this never actually
    // throws, so callers written against the historical noexcept contract
    // are unaffected in practice.
    void store(ResultPtr&& data) {
      // Regression tripwire: a NON-persistent entry re-stored without an
      // intervening reset() means the same value was produced twice within
      // one evaluation (a duplicate producer). PERSISTENT entries are
      // excluded: they legitimately re-store across batch replays.
      if (!persistent_) {
        SEQUANT_ASSERT(!stored_this_eval_);
        stored_this_eval_ = true;
      }
      data_p = std::move(data);
      size_bytes_
          .reset();  // (re)computed lazily on demand; see size_in_bytes()
    }

    void reset() noexcept {
      life_c = max_life;
      stored_this_eval_ = false;
      if (!persistent_) {  // persistent data (and its size) survives reset()
        data_p = nullptr;
        size_bytes_.reset();
      }
    }

    [[nodiscard]] bool persistent() const noexcept { return persistent_; }

    /// \return whether this NON-persistent entry has been \c store()'d since
    ///         its last \c reset() -- the re-store tripwire's flag state, used
    ///         by tests to observe the guard without depending on assert
    ///         behavior (which is elided unless \c SEQUANT_ASSERT_ENABLED).
    ///         Always \c false for a persistent entry (never set).
    [[nodiscard]] bool stored_this_eval() const noexcept {
      return stored_this_eval_;
    }

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

    /// Set (or reset) this entry's bounded life to exactly @p count uses,
    /// preserving any currently stored data. Used by \c
    /// CacheManager::ensure_home_slot(key, use_count, persistent) to home a
    /// value with a genuine use-count-bounded lifetime -- released at its
    /// count-th access -- rather than \c make_resident's unconditional
    /// unbounded pin.
    void set_life(size_t count) noexcept {
      max_life = count;
      life_c = count;
    }

    /// Upgrade this entry to persistent (never drained on access, survives
    /// reset()), preserving any currently stored data. Used by \c
    /// CacheManager::ensure_home_slot(key, use_count, persistent) to promote
    /// an existing entry when a later caller discovers the key is actually
    /// iteration-invariant.
    void make_persistent() noexcept { persistent_ = true; }

   private:
    [[nodiscard]] int decay() noexcept {
      return life_c > 0 ? static_cast<int>(--life_c) : 0;
    }

  };  // entry

  // NOT noexcept: forwards to entry::store(), which is not noexcept (see
  // there).
  static ResultPtr store(entry& ent, ResultPtr&& data) {
    ent.store(std::move(data));
    return ent.access();
  }

  std::unordered_map<cache_key_type, entry, hasher_type, comparator_type>
      cache_map_;

  /// DIAGNOSTIC: per-DISTINCT-value build tally (see BuildTally), keyed by the
  /// same node identity as cache_map_. Populated by tally_build() from the eval
  /// loop's build choke point (eval.hpp finish_phase_b) for EVERY product
  /// build, whether that value is a cache entry, a footprint-gated recompute,
  /// or a per-batch rebuild -- so the rollup is complete. Held only on the
  /// scope- chain ROOT (tally_build routes there); scratch caches never
  /// populate it, and reset() does NOT clear it (the tally spans the whole
  /// forest replay).
  // Keyed by NODE identity (not the cache's value-id): a per-DISTINCT-value
  // build diagnostic the dry-run costing rolls up, consumed node-by-node by the
  // eval tests. Pillar 1's slice-colored value identity lives in cache_map_
  // (correctness); this rollup stays node-keyed (byte-identical to before).
  std::unordered_map<TreeNode, BuildTally,
                     TreeNodeHasher<TreeNode, force_hash_collisions>,
                     TreeNodeEqualityComparator<TreeNode>>
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

  /// Backend realizations of external-axis batching's array ops (zero
  /// destination + axis chunking); see \c BackendArrayOps. Inherited from
  /// \c parent_ (only the root cache is wired in practice). Non-owning; the
  /// pointee must outlive this cache.
  BackendArrayOps const* array_ops_ = nullptr;

  /// Non-owning loop-colored slice seam (see \c LoopColoredSliceSeam,
  /// dag_scope.hpp): the per-value (hash-keyed) sliced-mode -> loop assignment
  /// \c slice_to_use (eval.hpp) reads to resolve a fetched value's physical
  /// slice mode off the loop-colored canonical layout -- the ordered path's
  /// sole slice-mode source. Populated by the ordered executor's shared core
  /// (\c run_ordered_schedule_pre_results) from the \c OrderedSchedule's
  /// \c compute_sliced_mode_assignment, and inherited from \c parent_ (only
  /// the root cache is wired in practice), mirroring \c array_ops_ / \c
  /// placement_router_. Null (default) => no seam wired => the ordered arm
  /// leaves the fetch unsliced (byte-identical). Non-owning; the pointee
  /// must outlive this cache.
  LoopColoredSliceSeam const* loop_colored_slice_seam_ = nullptr;

  /// The CONSUMER identity (an eval-node hash) currently fetching values under
  /// this cache (sliced-value canonical-layout / loop-coloring design, PILLAR
  /// 2). The ordered executor sets this to the hash of the member-root / output
  /// it is about to \c evaluate_impl (ordered_executor.hpp), bracketing each
  /// such call, so \c slice_to_use (eval.hpp) can disambiguate WHICH use-site
  /// is fetching a shared symmetric value -- the datum \c
  /// LoopColoredSliceSeam::by_hash_consumer needs to bind each occurrence to
  /// its own free mode. Nullopt (default) => no consumer tracked => the seam
  /// falls back to its consumer-blind first-match (byte-identical). Inherited
  /// from \c parent_ so a per-block child scratch sees the enclosing consumer.
  std::optional<std::size_t> current_consumer_{};

  /// TEST-ONLY observer of every slice_to_use decision (sliced-value
  /// canonical-layout / loop-coloring design, PILLAR 2): invoked -- only when
  /// set -- with (fetched value hash, current consumer, resolved physical slice
  /// position, loop level ordinal) just before the slice is applied. Lets a
  /// unit test witness that two consumers of one symmetric shared value slice
  /// DIFFERENT physical modes (the w8 fix), which is otherwise internal to the
  /// per-batch fetch and invisible in the final accumulated result. Unset
  /// (default) => zero cost, byte-identical. Inherited from \c parent_ so the
  /// per-block scratch reports through the root cache's observer.
  std::function<void(std::size_t, std::optional<std::size_t>, std::size_t, int)>
      slice_observer_{};

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

  /// Optional multi-root driver consulted by evaluate_multiroot() (see
  /// multiroot_driver_type). Empty => evaluate_multiroot() throws (no
  /// per-root fallback).
  multiroot_driver_type multiroot_driver_{};

  /// Enclosing realized batch loops for slice-on-use (see BatchContext). Empty
  /// (default) => no enclosing batch loop; the batched evaluator sets it on the
  /// per-block scratch before each re-entry. Not cleared by reset() (it is
  /// per-loop-iteration structural, re-set each block by the evaluator).
  BatchContext batch_context_{};

  /// When true, \c evaluate_impl treats a cache MISS on a non-leaf, non-top
  /// node as a hard error (\c sequant::Exception) instead of silently
  /// recomputing or serving an empty/unfilled array. The ordered read-from-home
  /// discipline statically pre-schedules every value, so a vanished operand
  /// (premature eviction -- e.g. an under-predicted use count in \c
  /// ordered_home_reads) is a real defect that must surface loudly, never hang
  /// a downstream contraction waiting on tiles that will never be produced.
  /// Set on every read-from-home scratch (see \c make_batched_scratch);
  /// parent-inheriting so nested scratches enforce the same invariant. Default
  /// false => forest/recursive evaluation (miss => compute) is unchanged.
  bool require_resident_reads_ = false;

  /// Non-owning placement router (see \c placement_router.hpp). Null
  /// (default) => no override wired; \c placement_router() falls through to
  /// \c parent_ (only the root cache is wired in practice). The pointee must
  /// outlive this cache.
  eval::PlacementRouter<TreeNode> const* placement_router_ = nullptr;

  /// Pillar 1 / Task 7: the LOCAL per-build coloring context (see \c
  /// ValueColoringCtx). Non-owning, set by the ordered executor around an
  /// `evaluate_impl` build. LOCAL only (NOT inherited through the parent
  /// chain): a map-keying method builds its key using THIS cache's context, and
  /// access_at then walks parents with that already-built key. Null (default)
  /// => no re-color => byte-identical.
  ValueColoringCtx const* value_coloring_ctx_ = nullptr;

  /// Pillar 1 / Task 7: re-color an EMPTY-coloring key from the local per-build
  /// coloring context, so a bare-node key (`evaluate_impl` passes nodes, which
  /// implicit-convert to an empty-coloring `CachedValue`) picks up the home
  /// coloring of the value that node represents in this build. A key that
  /// already carries a coloring (the executor's own `value_of(vid)` calls) is
  /// respected as-is; a node absent from the context, or a null context,
  /// returns the key unchanged => byte-identical without a context.
  [[nodiscard]] cache_key_type recolor(cache_key_type const& key) const {
    if (value_coloring_ctx_ && key.coloring.empty()) {
      auto const it = value_coloring_ctx_->find(key.node->hash_value());
      if (it != value_coloring_ctx_->end() && !it->second.empty())
        return cache_key_type{key.node, it->second};
    }
    return key;
  }

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

  /// Sets the multi-root driver (see multiroot_driver_type). Pass an empty
  /// std::function to clear it.
  void set_multiroot_driver(multiroot_driver_type fn) noexcept {
    multiroot_driver_ = std::move(fn);
  }

  /// \return the multi-root driver (empty if none is set).
  [[nodiscard]] multiroot_driver_type const& multiroot_driver() const noexcept {
    return multiroot_driver_;
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

  /// Enable/disable the resident-reads invariant (see \c
  /// require_resident_reads_). Set true on read-from-home scratches.
  void set_require_resident_reads(bool v) noexcept {
    require_resident_reads_ = v;
  }

  /// \return whether a non-leaf, non-top cache miss must be a hard error --
  ///         this cache's own flag, else inherited from the parent chain.
  [[nodiscard]] bool require_resident_reads() const noexcept {
    return require_resident_reads_
               ? true
               : (parent_ ? parent_->require_resident_reads() : false);
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

  /// Sets the local per-build coloring context (see value_coloring_ctx_). Pass
  /// nullptr to detach. Non-owning; the pointee must outlive the build. LOCAL
  /// (not inherited): the RAII guard is the ordered executor's.
  void set_value_coloring_ctx(ValueColoringCtx const* ctx) noexcept {
    value_coloring_ctx_ = ctx;
  }

  /// \return the LOCAL per-build coloring context (or null). No parent
  /// fall-through: keys are colored using the cache that builds them.
  [[nodiscard]] ValueColoringCtx const* value_coloring_ctx() const noexcept {
    return value_coloring_ctx_;
  }

  /// Pillar 1 / Task 7: re-key every registered entry through the current
  /// coloring context, so a scratch whose members were registered by bare
  /// (uncolored) node becomes genuinely VALUE-keyed -- each member entry keyed
  /// by its home-slice-colored value-id, matching the colored store/access.
  /// Call once, right after construction, with the per-scope context set. A
  /// null context (or a scope with nothing sliced) leaves every key unchanged
  /// => byte-identical.
  void recolor_registered_entries() {
    if (!value_coloring_ctx_) return;
    std::unordered_map<cache_key_type, entry, hasher_type, comparator_type>
        rekeyed;
    rekeyed.reserve(cache_map_.size());
    for (auto& [k, e] : cache_map_) rekeyed.emplace(recolor(k), std::move(e));
    cache_map_ = std::move(rekeyed);
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

  /// Sets the backend array-ops (see BackendArrayOps). Non-owning; the pointee
  /// must outlive this cache. Absent (nullptr) anywhere along the chain means a
  /// batched external-axis scatter has no way to build its destination.
  void set_array_ops(BackendArrayOps const* a) noexcept { array_ops_ = a; }

  /// \return the local array-ops if set, else the one inherited from
  ///         \c parent_ (only the root cache is wired in practice); nullptr if
  ///         none is wired anywhere along the chain. Non-owning.
  [[nodiscard]] BackendArrayOps const* array_ops() const noexcept {
    return array_ops_ ? array_ops_ : parent_ ? parent_->array_ops() : nullptr;
  }

  /// Sets the loop-colored slice seam (see loop_colored_slice_seam_). Pass
  /// nullptr to detach. Non-owning; the pointee must outlive this cache.
  void set_loop_colored_slice_seam(LoopColoredSliceSeam const* s) noexcept {
    loop_colored_slice_seam_ = s;
  }

  /// \return the local loop-colored slice seam if set, else the one inherited
  ///         from \c parent_ (only the root cache is wired in practice);
  ///         nullptr if none is wired anywhere along the chain. Non-owning.
  [[nodiscard]] LoopColoredSliceSeam const* loop_colored_slice_seam()
      const noexcept {
    return loop_colored_slice_seam_ ? loop_colored_slice_seam_
           : parent_                ? parent_->loop_colored_slice_seam()
                                    : nullptr;
  }

  /// Sets the current consumer identity (see current_consumer_). Pass nullopt
  /// to clear. Cheap value set; the ordered executor RAII-restores it around
  /// each member-root evaluate_impl.
  void set_current_consumer(std::optional<std::size_t> consumer) noexcept {
    current_consumer_ = consumer;
  }

  /// \return the consumer identity set on this cache, else the one inherited
  ///         from \c parent_ (a per-block child scratch defers to its enclosing
  ///         scope's consumer); nullopt if none is set anywhere along the
  ///         chain.
  [[nodiscard]] std::optional<std::size_t> current_consumer() const noexcept {
    return current_consumer_ ? current_consumer_
           : parent_         ? parent_->current_consumer()
                             : std::nullopt;
  }

  /// Sets the TEST-ONLY slice observer (see slice_observer_). Pass an empty
  /// function to clear.
  void set_slice_observer(
      std::function<void(std::size_t, std::optional<std::size_t>, std::size_t,
                         int)>
          obs) {
    slice_observer_ = std::move(obs);
  }

  /// \return the slice observer set on this cache, else the one inherited from
  ///         \c parent_; an empty function if none is set along the chain.
  [[nodiscard]] std::function<void(std::size_t, std::optional<std::size_t>,
                                   std::size_t, int)> const&
  slice_observer() const noexcept {
    if (slice_observer_) return slice_observer_;
    if (parent_) return parent_->slice_observer();
    return slice_observer_;  // empty
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
  void ensure_hoist_slot(cache_key_type const& key) {
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
  void ensure_home_slot(cache_key_type const& key) {
    auto [it, inserted] = cache_map_.try_emplace(
        key, entry{std::numeric_limits<size_t>::max(), /*persistent=*/false});
    if (!inserted) it->second.make_resident();
  }

  /// Home @p key with a bounded (use-count) or persistent lifetime, instead
  /// of \c ensure_home_slot(key)'s unconditional \c make_resident pin. A
  /// persistent slot survives \c reset() (iteration-invariant); a
  /// non-persistent slot is released at its @p use_count-th access (its
  /// genuine last use) rather than living unbounded until the next reset().
  /// Idempotent like \c ensure_home_slot(key): an already-present entry is
  /// upgraded in place (to persistent, or to the new bounded life) rather
  /// than replaced, preserving any stored data.
  void ensure_home_slot(cache_key_type const& key, std::size_t use_count,
                        bool persistent) {
    auto [it, inserted] = cache_map_.try_emplace(
        key, entry{persistent ? std::numeric_limits<size_t>::max() : use_count,
                   persistent});
    if (!inserted) {
      if (persistent)
        it->second.make_persistent();
      else
        it->second.set_life(use_count);
    }
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

  /// \return whether @p key's cache entry (this cache or an ancestor) is
  ///         classified PERSISTENT. The shared cache builder
  ///         (\c sequant::cache_manager) computes the correct persistence
  ///         predicate -- non-volatile AND has a volatile DIRECT consumer --
  ///         and stamps it here at construction. The ordered executor's
  ///         home-slot seeding consults this instead of re-deriving
  ///         persistence, so it cannot over-enroll a non-volatile value that
  ///         has no volatile consumer. \c false when the key has no entry (the
  ///         shared cache chose not to register it, i.e. it is not persistent).
  [[nodiscard]] bool entry_is_persistent(
      cache_key_type const& key) const noexcept {
    for (CacheManager const* c = this; c; c = c->parent_) {
      auto const it = c->cache_map_.find(key);
      if (it != c->cache_map_.end()) return it->second.persistent();
    }
    return false;
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
  void tally_build(cache_key_type const& key, std::string const& slice_sig,
                   double flops, double exec) noexcept {
    if (parent_) {
      parent_->tally_build(key, slice_sig, flops, exec);
      return;
    }
    if (!recompute_tally_enabled_) return;  // wet path: no-op
    auto& slice = recompute_tally_[key.node].slices[slice_sig];
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
  [[nodiscard]] std::unordered_map<
      TreeNode, BuildTally, TreeNodeHasher<TreeNode, force_hash_collisions>,
      TreeNodeEqualityComparator<TreeNode>> const&
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

  /// \return true iff some ALIVE entry on this cache or any ancestor along the
  ///         scope chain physically holds @p value (pointer identity) AND that
  ///         entry is either PERSISTENT (\c persistent(), never drained --
  ///         survives across evaluations) or has MORE THAN ONE consumer (\c
  ///         max_life_count() > 1) -- i.e. @p value is a SHARED buffer with at
  ///         least one read still pending (or a resident/persistent home read
  ///         by every consumer), so mutating it in place would corrupt those
  ///         other reads.
  ///
  /// This is the runtime safety gate the in-place \c Sum accumulation
  /// (eval.hpp) needs, and is strictly the "the accumulator is shared" test:
  /// a PRIVATE single-use buffer is NOT reported here, so in-place still fires
  /// for it. Two cases yield a private (unshared) accumulator, both correctly
  /// returning false:
  ///   - a TRANSIENT running total (freshly allocated by a prior \c sum() /
  ///     \c prod(), never \c store'd) is held by no entry at all; and
  ///   - a single-use CSE entry (\c max_life == 1) MOVES its buffer out on its
  ///     sole \c entry::access() (\c decay() reaches 0 -> \c
  ///     std::move(data_p)), so once it has been read as an operand it no
  ///     longer \c holds() it -- and even were it somehow still alive, \c
  ///     max_life > 1 excludes it.
  /// A value homed RESIDENT (\c ensure_home_slot, \c max_life == SIZE_MAX) or
  /// with a genuine multi-use count (\c max_life > 1, e.g. a subexpression
  /// shared across two roots, or a per-batch-reread home) is never drained by
  /// a single read, so it stays held and IS reported here -- the case the
  /// elided \c SEQUANT_ASSERT could not catch at runtime. A PERSISTENT entry
  /// (registered via the \c CacheManager(Iterable&&, PersistencePred) ctor,
  /// e.g. the \c make_batched_scratch path) is ALSO reported even when its
  /// \c max_life == 1: \c entry::access() never drains a persistent entry, so
  /// it stays resident-forever and is shared across the batched replays --
  /// mutating it in place would corrupt every later read. Hence the guard is
  /// \c persistent() OR \c max_life > 1, not \c max_life alone.
  /// A held NON-persistent single-use entry (\c max_life == 1) is NOT reported
  /// (it drained its buffer out on its sole read, so it no longer holds it),
  /// keeping in-place enabled for the private common case.
  [[nodiscard]] bool chain_holds_shared(ResultPtr const& value) const noexcept {
    if (!value) return false;
    for (auto const& [k, e] : cache_map_)
      if (e.holds(value) && (e.persistent() || e.max_life_count() > 1))
        return true;
    return parent_ ? parent_->chain_holds_shared(value) : false;
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
  [[nodiscard]] AccessResult access_at(
      cache_key_type const& key,
      std::optional<std::size_t> origin_consumer = std::nullopt) noexcept {
    // Task 7: color the key from THIS cache's per-build context (no-op without
    // one), then look up AND walk parents with the ALREADY-colored key. The
    // cache genuinely keys by VALUE: two values of one node coexist as distinct
    // colored entries (a home value found by its own home identity, never a
    // same-node sibling), which is the whole point of the value-keyed cache.
    // DIAGNOSTIC: capture the consuming value at the ORIGIN cache (where the
    // fetch starts) and thread it through the parent recursion so a root-homed
    // hit can name WHO read it.
    if (!origin_consumer) origin_consumer = current_consumer();
    cache_key_type const rk = recolor(key);
    if (auto found =
            eval::LookupMeter::timed([&] { return cache_map_.find(rk); });
        found != cache_map_.end()) {
      static long const _ax_target = [] {
        char const* dh = std::getenv("SEQUANT_COUNT_ACCESS");
        return dh ? std::strtol(dh, nullptr, 10) : -1L;
      }();
      bool const _ax_match =
          _ax_target >= 0 && (found->first->hash_value() % 100000u) ==
                                 static_cast<unsigned>(_ax_target);
      if (auto data = found->second.access(); data) {
        if (_ax_match) {
          static std::size_t _ax_n = 0;
          std::cerr << "[axcount] hash="
                    << (found->first->hash_value() % 100000u) << " read#"
                    << (++_ax_n)
                    << " remaining_life=" << found->second.life_count()
                    << " max_life=" << found->second.max_life_count()
                    << " consumer="
                    << (origin_consumer ? *origin_consumer % 100000u : 0u)
                    << std::endl;
        }
        // DIAGNOSTIC (SEQUANT_UT_ACCESS_CLOCK): stamp this genuine local-hit
        // read into the global access clock. No-op when the gate is off.
        eval::AccessClock::stamp(found->first->hash_value());
        eval::BuildMeter::on_read(found->first->hash_value(),
                                  found->second.max_life_count());
        // NB: size_in_bytes() walks the array; only pay it when metering.
        if (eval::DefUseMeter::enabled())
          eval::DefUseMeter::on_read(found->first->hash_value(),
                                     found->second.size_in_bytes());
        return {data, 0};
      } else if (_ax_match) {
        std::cerr << "[axcount] hash=" << (found->first->hash_value() % 100000u)
                  << " MISS-DRAINED (data absent) remaining_life="
                  << found->second.life_count()
                  << " max_life=" << found->second.max_life_count()
                  << " consumer="
                  << (origin_consumer ? *origin_consumer % 100000u : 0u)
                  << std::endl;
      }
    }
    if (!parent_) return {nullptr, 0};
    auto up = parent_->access_at(rk, origin_consumer);
    return {up.ptr, up.hops + 1};  // count the link we just crossed
  }

  /// @param key The key that identifies the cached data.
  /// @return ResultPtr to Result. Thin forwarder to access_at() that drops the
  ///         hop distance, for the non-batched callers that do not slice.
  ResultPtr access(cache_key_type const& key) noexcept {
    return access_at(key).ptr;
  }

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
  [[nodiscard]] ResultPtr access_at_hops(cache_key_type const& key,
                                         std::size_t hops) noexcept {
    cache_key_type const rk = recolor(key);  // Task 7: color at THIS cache
    CacheManager* c = this;
    for (std::size_t i = 0; i < hops && c; ++i) c = c->parent_;
    if (!c) return nullptr;
    if (auto found =
            eval::LookupMeter::timed([&] { return c->cache_map_.find(rk); });
        found != c->cache_map_.end()) {
      auto data = found->second.access();
      static long const _axh_target = [] {
        char const* dh = std::getenv("SEQUANT_COUNT_ACCESS");
        return dh ? std::strtol(dh, nullptr, 10) : -1L;
      }();
      if (_axh_target >= 0 && (found->first->hash_value() % 100000u) ==
                                  static_cast<unsigned>(_axh_target)) {
        std::cerr << "[axcount-hops] hash="
                  << (found->first->hash_value() % 100000u) << " hops=" << hops
                  << (data ? " HIT" : " MISS-DRAINED")
                  << " remaining_life=" << found->second.life_count()
                  << " max_life=" << found->second.max_life_count()
                  << std::endl;
      }
      // DIAGNOSTIC (SEQUANT_UT_ACCESS_CLOCK): stamp this genuine
      // router-directed read into the global access clock. No-op when the gate
      // is off.
      if (data) {
        eval::AccessClock::stamp(found->first->hash_value());
        if (eval::DefUseMeter::enabled())
          eval::DefUseMeter::on_read(found->first->hash_value(),
                                     found->second.size_in_bytes());
      }
      return data;
    }
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
  // NOT noexcept: forwards to entry::store(), which is not noexcept (see
  // there).
  [[nodiscard]] ResultPtr store(cache_key_type const& key, ResultPtr data) {
    cache_key_type const rk = recolor(key);  // Task 7: color at this cache
    if (auto found =
            eval::LookupMeter::timed([&] { return cache_map_.find(rk); });
        found != cache_map_.end()) {
      // INSTRUMENTATION (SEQUANT_REBUILD_TRACE, analysis-only): a store onto an
      // entry that is ALREADY alive means the value was still resident in cache
      // yet got rebuilt anyway -- a hash-bypassing recompute. Emit the hash and
      // the persistent flag so a persistent (should-never-recompute) value that
      // is nonetheless rebuilt-while-resident is caught. No-op when the gate is
      // off.
      static bool const rebuild_trace =
          std::getenv("SEQUANT_REBUILD_TRACE") != nullptr;
      if (rebuild_trace && found->second.alive())
        std::cerr << "REBUILD-RESIDENT persist=" << found->second.persistent()
                  << " hash=" << key->hash_value() << "\n";
      eval::DefUseMeter::on_store(key->hash_value());
      return store(found->second, std::move(data));
    }
    return data;
  }

  ///
  /// \brief Check if the key exists in the database: does not check if cache
  ///        exists
  ///
  [[nodiscard]] bool exists(cache_key_type const& key) const noexcept {
    return cache_map_.find(recolor(key)) != cache_map_.end();
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
    // The map key is a CachedValue; callers enumerate NODES (key_type), so hand
    // them the node. The coloring is a cache-internal identity detail.
    for (auto const& [k, v] : cache_map_) fn(k.node);
  }

  /// if the key exists in the database, return the current lifetime count of
  /// the cached data otherwise return -1
  [[nodiscard]] int life(cache_key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.life_count());
  }

  /// if the key exists in the database, return the maximum lifetime count of
  /// the cached data that implies the maximum number of accesses allowed for
  /// this key before the cache is released. This value was set by the c'tor.
  [[nodiscard]] int max_life(cache_key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    auto end = cache_map_.end();
    return iter == end ? -1 : static_cast<int>(iter->second.max_life_count());
  }

  /// \return true iff the key is registered for caching and currently holds
  ///         stored data (i.e. has been stored and not yet drained by its
  ///         final access).
  [[nodiscard]] bool alive(cache_key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && iter->second.alive();
  }

  /// \return true iff @p key is alive (holding data) at THIS cache or ANY
  ///         ancestor scope up the parent chain. Non-decrementing (unlike \c
  ///         access_at): a pure residency probe. Used by \c
  ///         make_batched_scratch to decide that a batch-invariant value
  ///         already resident at its home is read from there each batch (the
  ///         parent-chain fall-through), so it is neither registered nor
  ///         rebuilt in the per-batch scratch.
  [[nodiscard]] bool resident_in_chain(
      cache_key_type const& key) const noexcept {
    if (auto iter =
            eval::LookupMeter::timed([&] { return cache_map_.find(key); });
        iter != cache_map_.end() && iter->second.alive())
      return true;
    return parent_ ? parent_->resident_in_chain(key) : false;
  }

  /// \return true iff the key is registered for caching and classified
  ///         persistent (P: never released on access, survives reset()).
  [[nodiscard]] bool persistent(cache_key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && iter->second.persistent();
  }

  /// \return true iff the key is registered, NON-persistent, and has been
  ///         \c store()'d since its last \c reset() -- the re-store
  ///         tripwire's flag state (see \c entry::store()). Test-facing
  ///         accessor; always \c false for a persistent key or an
  ///         unregistered one.
  [[nodiscard]] bool stored_this_eval(
      cache_key_type const& key) const noexcept {
    auto iter = cache_map_.find(key);
    return iter != cache_map_.end() && iter->second.stored_this_eval();
  }

  /// \return size in bytes of the data currently held for @p key, or 0 if
  ///         the key is not registered or no data is currently stored.
  [[nodiscard]] size_t entry_size_in_bytes(
      cache_key_type const& key) const noexcept {
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
///                    default (1) everything is cached, so use-count tracking
///                    is exact.
///
/// \return A cache manager.
///
/// \see CacheManager
///
template <bool force_hash_collisions = false>
auto cache_manager(meta::eval_node_range auto const& nodes,
                   size_t min_repeats = 1) noexcept {
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
/// \param min_repeats minimum NP repeats to cache (default 1).
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
                   size_t min_repeats = 1, FootprintOf footprint_of = {},
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
  // no External node_slice_mask() stamps (every mask stays empty/all-full), so
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
  // intermediate (all-full mask; or an External-only / no node_slice_mask
  // entry, e.g. gC) lands. OFF path (no order-aware annotations, hence no \c
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
