#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/eval/fwd.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/member_axis.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/placement_router.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/schedule_dump.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/range/operations.hpp>

#include <algorithm>
#include <any>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>

// Headers for process_rss_bytes() — see log::process_rss_bytes() below.
#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#include <fstream>
#endif

namespace sequant {

namespace log {

using Duration = std::chrono::nanoseconds;

struct Bytes {
  size_t value;
};

/// \return whether the eval trace is being printed (logger level > 0). The
/// eval engine's per-op memory diagnostics (the hwmark working set and the
/// cache total) feed only the trace output, so the standalone functions that
/// compute them -- the cache-aware bytes() overload below and eval()/cache()
/// -- short-circuit when this is false, and are thus computed only when a
/// trace line will actually be emitted.
[[nodiscard]] inline bool printing() noexcept {
  return Logger::instance().eval.level > 0;
}

template <typename T, typename... Ts>
[[nodiscard]] inline auto bytes(T const& arg, Ts const&... args) {
  auto one = [](auto const& a) -> size_t {
    if constexpr (requires {
                    static_cast<bool>(a);
                    a->size_in_bytes();
                  }) {
      // Smart-pointer-like operand: tolerate null so callers (e.g. the
      // EvalOp::Adjoint dispatcher, which leaves `right` unevaluated) can
      // pass an empty ResultPtr without an external guard.
      return a ? a->size_in_bytes() : size_t{0};
    } else if constexpr (requires { a->size_in_bytes(); })
      return a->size_in_bytes();
    else
      return a.size_in_bytes();
  };
  return Bytes{(one(arg) + ... + one(args))};
}

/// Cache-aware bytes(): bytes(cache) + bytes(args...), where bytes(cache) sums
/// the (lazily memoized) sizes of the cache's alive entries. This is the
/// working-set/total walk used only to populate the trace's hwmark/total
/// fields, so it short-circuits to 0 when no trace line will be printed --
/// avoiding the per-op walk of every alive entry, which with persistent
/// entries is otherwise paid on every op of every iteration.
template <typename N, bool F, typename... Ts>
[[nodiscard]] inline Bytes bytes(CacheManager<N, F> const& cache,
                                 Ts const&... args) {
  if (!printing()) return Bytes{0};
  return Bytes{cache.size_in_bytes() + (size_t{0} + ... + bytes(args).value)};
}

[[nodiscard]] inline auto to_string(Bytes bs) noexcept {
  return std::format("{}B", bs.value);
}

/// \return the process physical-memory footprint, in bytes. On macOS this is
/// `phys_footprint` from `TASK_VM_INFO` — the same accounted footprint Activity
/// Monitor reports in its "Memory" column, and what jetsam limits act on. It
/// excludes shared/reclaimable pages (frameworks, the shared cache, file-backed
/// clean pages), so it is much smaller than the raw resident-set size
/// (`mach_task_basic_info::resident_size`), which double-counts shared text and
/// is what made this column read far larger than Activity Monitor. On Linux we
/// read resident pages from `/proc/self/statm`. Returns 0 on other platforms
/// and on read failure (no exception, so safe to call from logging paths).
/// Cheap (~µs) — intended to be called once per log record.
[[nodiscard]] inline std::size_t process_rss_bytes() noexcept {
#if defined(__APPLE__)
  ::task_vm_info_data_t vm_info{};
  ::mach_msg_type_number_t vm_count = TASK_VM_INFO_COUNT;
  if (::task_info(::mach_task_self(), TASK_VM_INFO,
                  reinterpret_cast<::task_info_t>(&vm_info),
                  &vm_count) == KERN_SUCCESS &&
      vm_count >= TASK_VM_INFO_COUNT) {
    return static_cast<std::size_t>(vm_info.phys_footprint);
  }
  // Fallback: raw resident-set size (larger; includes shared pages).
  ::mach_task_basic_info_data_t info{};
  ::mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (::task_info(::mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<::task_info_t>(&info),
                  &count) != KERN_SUCCESS) {
    return 0;
  }
  return static_cast<std::size_t>(info.resident_size);
#elif defined(__linux__)
  // /proc/self/statm columns are page counts:
  //   total resident shared text lib data dt
  std::ifstream f("/proc/self/statm");
  std::size_t pages_total = 0, pages_resident = 0;
  if (!(f >> pages_total >> pages_resident)) return 0;
  static const long page_size = ::sysconf(_SC_PAGESIZE);
  if (page_size <= 0) return 0;
  return pages_resident * static_cast<std::size_t>(page_size);
#else
  return 0;
#endif
}

/// Convenience wrapper around process_rss_bytes() returning a Bytes.
[[nodiscard]] inline Bytes rss() noexcept { return Bytes{process_rss_bytes()}; }

/// type of data or operation
enum struct EvalMode {
  Constant,
  Variable,
  Power,
  Tensor,
  Permute,
  Product,
  MultByPhase,
  Sum,
  SumInplace,
  Symmetrize,
  Antisymmetrize,
  Unknown
};

[[nodiscard]] EvalMode eval_mode(meta::eval_node auto const& node) {
  if (node.leaf()) {
    return node->is_constant()   ? EvalMode::Constant
           : node->is_variable() ? EvalMode::Variable
           : node->is_power()    ? EvalMode::Power
           : node->is_tensor()   ? EvalMode::Tensor
                                 : EvalMode::Unknown;
  } else {
    return node->is_product()   ? EvalMode::Product
           : node->is_sum()     ? EvalMode::Sum
           : node->is_adjoint() ? EvalMode::Permute
                                : EvalMode::Unknown;
  }
}

[[nodiscard]] constexpr auto to_string(EvalMode mode) noexcept {
  return (mode == EvalMode::Constant)         ? "Constant"
         : (mode == EvalMode::Variable)       ? "Variable"
         : (mode == EvalMode::Power)          ? "Power"
         : (mode == EvalMode::Tensor)         ? "Tensor"
         : (mode == EvalMode::Permute)        ? "Permute"
         : (mode == EvalMode::Product)        ? "Product"
         : (mode == EvalMode::MultByPhase)    ? "MultByPhase"
         : (mode == EvalMode::Sum)            ? "Sum"
         : (mode == EvalMode::SumInplace)     ? "SumInplace"
         : (mode == EvalMode::Symmetrize)     ? "Symmetrize"
         : (mode == EvalMode::Antisymmetrize) ? "Antisymmetrize"
                                              : "??";
}

enum struct CacheMode { Store, Access, Release };

[[nodiscard]] constexpr auto to_string(CacheMode mode) noexcept {
  return (mode == CacheMode::Store)    ? "Store"
         : (mode == CacheMode::Access) ? "Access"
                                       : "Release";
}

enum struct TermMode { Begin, End };

[[nodiscard]] constexpr auto to_string(TermMode mode) noexcept {
  return (mode == TermMode::Begin) ? "Begin" : "End";
}

/// One log record per eval op. Line format:
///
// clang-format off
/// Eval | <mode> | <time> | [left=L | right=R |] result=X | alloc=A | hw=H | rss=R | [<heap suffix> |] <label>
// clang-format on
///
/// Which fields are set depends on the op's arity:
///
///   mode                                          | left/right | alloc
///   ----------------------------------------------+------------+--------
///   Constant / Variable / Tensor (leaf)           | —          | result
///   Permute / MultByPhase /                       | —          | result
///     Symmetrize / Antisymmetrize                 |            |
///   SumInplace                                    | —          | 0B
///   Sum / Product                                 | set        | result
///
/// Only Sum and Product set left/right, since their operand sizes can
/// differ from the result. Other modes omit those fields rather than
/// zeroing them, so a logged 0B always means an empty buffer.
///
/// mem_result is the size of the buffer the op produces; for SumInplace
/// it's the size of the accumulator after the add. mem_alloc is what the
/// op allocated — equal to mem_result everywhere except SumInplace,
/// which writes into the accumulator and allocates nothing.
///
/// mem_hwmark is the eval engine's high-water mark: the running maximum,
/// over all ops since the cache was last reset, of the per-op live working
/// set
///
///   bytes(cache) + bytes(result) + bytes of each operand not aliased
///                                  to a cache entry
///
/// (aliasing is decided at each call site by CacheManager::chain_holds --
/// pointer identity against every alive entry on the scope chain -- so an
/// operand read full from a cache at ANY scope is counted once via the cache
/// residency, while a sliced/permuted/phase-shifted read is a distinct buffer
/// and is added). bytes(cache) here is this cache's own residency; the
/// ancestors' residency is added separately as CacheManager::chain_residency().
/// It is reported as a running max so it is
/// monotonically non-decreasing within one evaluation — the peak memory the
/// engine reaches — rather than the instantaneous per-op working set, which
/// oscillates as the cache fills and drains. The max is held by the
/// CacheManager and cleared by CacheManager::reset() (called per term), so
/// each term reports its own peak.
///
/// rss is the process physical-memory footprint measured immediately before
/// the record is emitted (`phys_footprint` via `TASK_VM_INFO` on macOS — the
/// value Activity Monitor's "Memory" column shows; resident pages from
/// `/proc/self/statm` on Linux; 0 on other platforms). Use it to triage
/// memory held outside the eval engine — long-lived tensors not in the
/// cache, runtime/library overhead, allocator fragmentation. mem_hwmark and
/// rss diverge by roughly that "everything else" component.
struct EvalStat {
  EvalMode mode;
  Duration time;
  Bytes mem_result{};
  Bytes mem_alloc{};
  Bytes mem_hwmark{};
  std::optional<Bytes> mem_left;
  std::optional<Bytes> mem_right;
};

struct CacheStat {
  CacheMode mode;
  size_t key;
  int curr_life, max_life;
  size_t num_alive;
  Bytes entry_memory;
  Bytes total_memory;
};

template <typename Arg, typename... Args>
void log(Arg const& arg, Args const&... args) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) write_log(l, arg, std::format(" | {}", args)..., '\n');
}

template <typename... Args>
auto eval(EvalStat const& stat, Args const&... args) {
  if (!printing()) return;  // nothing to format/emit; skip rss() and formatting
  auto const result_s = std::format("result={}", to_string(stat.mem_result));
  auto const alloc_s = std::format("alloc={}", to_string(stat.mem_alloc));
  auto const hw_s = std::format("hw={}", to_string(stat.mem_hwmark));
  // Reduce this rank's RSS to the value to report (e.g. sum over ranks = true
  // total app memory). This runs on every rank (printing() is level>0,
  // identical across ranks), so an injected collective reducer is matched.
  auto const rss_local = process_rss_bytes();
  auto const& rss_reduce = Logger::instance().eval.rss_reduce;
  auto const rss_s = std::format(
      "rss={}",
      to_string(Bytes{rss_reduce ? rss_reduce(rss_local) : rss_local}));
  // Optional backend-supplied suffix (already reduced across ranks); omitted
  // entirely when the hook is unset or returns empty.
  auto const& heap_stats = Logger::instance().eval.heap_stats;
  auto const heap_s = heap_stats ? heap_stats() : std::string{};
  auto emit = [&](auto const&... trailer) {
    if (stat.mem_left) {
      SEQUANT_ASSERT(stat.mem_right);
      log("Eval",                                               //
          to_string(stat.mode),                                 //
          stat.time,                                            //
          std::format("left={}", to_string(*stat.mem_left)),    //
          std::format("right={}", to_string(*stat.mem_right)),  //
          result_s, alloc_s, hw_s, rss_s,                       //
          trailer...);
    } else {
      log("Eval",                //
          to_string(stat.mode),  //
          stat.time,             //
          result_s, alloc_s, hw_s, rss_s, trailer...);
    }
  };
  if (heap_s.empty())
    emit(args...);
  else
    emit(heap_s, args...);
}

template <typename... Args>
auto cache(CacheStat const& stat, Args const&... args) {
  log("Cache",                                                   //
      to_string(stat.mode),                                      //
      std::format("key={}", stat.key),                           //
      std::format("life={}/{}", stat.curr_life, stat.max_life),  //
      std::format("alive={}", stat.num_alive),                   //
      std::format("entry={}", to_string(stat.entry_memory)),     //
      std::format("total={}", to_string(stat.total_memory)),     //
      args...);
}

template <typename N, bool F, typename... Args>
auto cache(N const& node, CacheManager<N, F>& cm, Args const&... args) {
  // Structured runtime-schedule event, emitted to the cache's ScheduleSink when
  // one is wired (independent of the trace level). Keyed by node->hash_value()
  // -- the SAME identity the IR emitter (schedule_dump.hpp) writes -- so the
  // schedule visualizer joins runtime lifetimes onto the IR DAG by hash. Mode:
  // Store (first build), Access (reuse), Release (last use). Store-count > 1
  // for a hash means the value was rebuilt (recompute).
  if (auto* const sink = cm.schedule_sink(); sink && sink->os && !sink->fired) {
    auto const cl = cm.life(node);
    auto const ml = cm.max_life(node);
    char const* const evm = (cl == 0)        ? "Release"
                            : (cl + 1 == ml) ? "Store"
                                             : "Access";
    *sink->os << "SCHEDULE_RUN_EVENT {\"hash\":\"" << node->hash_value()
              << "\",\"mode\":\"" << evm << "\",\"life\":" << cl
              << ",\"max_life\":" << ml << "}\n";
  }
  if (!printing()) return;  // skip the entry/total size walks and formatting
  using CacheMode::Access;
  using CacheMode::Release;
  using CacheMode::Store;
  auto const key = hash::value(*node);
  auto const cur_l = cm.life(node);
  auto const max_l = cm.max_life(node);
  bool const release = cur_l == 0;
  bool const store = cur_l + 1 == max_l;
  cache(CacheStat{.mode = store     ? Store
                          : release ? Release
                                    : Access,
                  .key = key,
                  .curr_life = cur_l,
                  .max_life = max_l,
                  .num_alive = cm.alive_count(),
                  .entry_memory = {cm.entry_size_in_bytes(node)},
                  .total_memory = {bytes(cm)}},
        args...);
}

inline auto term(TermMode mode, std::string_view term) {
  log("Term", to_string(mode), term);
}

/// Invoke the optional post-op memory-release hook
/// (Logger::eval.release_memory) if set. Called after each freshly evaluated op
/// regardless of trace level, so a large transient's freed pages can be
/// returned before the next op allocates. The injected hook self-throttles;
/// empty hook = no-op.
inline void release_after_op() {
  if (auto const& rel = Logger::instance().eval.release_memory; rel) rel();
}

[[nodiscard]] auto label(meta::eval_node auto const& node) {
  return node->is_primary()
             ? node->label()
             : std::format("{} {} {} -> {}", node.left()->label(),
                           (node->is_product() ? "*"
                            : node->is_sum()   ? "+"
                                               : "??"),  //
                           node.right()->label(), node->label());
}

/// Scope annotation for the batch-loop enter/leave markers (and reused by
/// slice_home_annot): the base_keys of the currently-open batch loops,
/// outermost-first, e.g. `scope={i,i,K,}`. \p active is a CacheManager
/// batch_context() -- a sequence of {Index mode, element-range} entries.
template <typename BatchContext>
std::string scope_annot(BatchContext const& active) {
  std::string scope;
  for (auto const& e : active) {
    if (!scope.empty()) scope += ",";
    scope += toUtf8(e.axis.space().base_key());
  }
  return std::format("scope={{{}}}", scope);
}

/// Per-op slice/home-scope metadata, reported IDENTICALLY for the forest
/// evaluator and the whole-scope (DAG) executor so their traces diff cleanly:
///   canon=[<canon_indices full_labels>] sliced=[<ordinal>:<base_key> ...]
///   scope={..}
/// - canon = node->canon_indices() full_labels: THIS occurrence's array layout
///   (labels are valid here -- they identify the layout the ordinals index).
/// - sliced = this node's OWN batched modes -- its \c sliced_modes() stamp --
///   as `p:base_key` at each stamped mode's canonical POSITION p.
///   Ordinal:base-space-key ONLY -- never a bare label (labels are tree-scoped;
///   ordinals are DAG-safe). This is the value's truthful batched set, NOT a
///   space-match against the open loops (which would misleadingly flag every
///   same-space mode, e.g. a spectator's contracted occ under an occ loop).
/// - scope = the open loops, outer..inner (see scope_annot).
template <meta::eval_node Node, typename BatchContext>
std::string slice_home_annot(Node const& node, BatchContext const& active) {
  std::string canon, sliced;
  auto const& idxs = node->canon_indices();
  for (auto const& ix : idxs) {
    canon += toUtf8(ix.full_label());
    canon += " ";
  }
  container::svector<std::size_t> spos;
  for (auto const& sm : node->sliced_modes())
    if (auto const p = index_position(node, sm);
        p && std::find(spos.begin(), spos.end(), *p) == spos.end())
      spos.push_back(*p);
  std::sort(spos.begin(), spos.end());
  for (auto const p : spos)
    sliced += std::format("{}:{} ", p,
                          toUtf8(std::wstring(idxs[p].space().base_key())));
  auto annot = std::format("canon=[{}] sliced=[{}] {}", canon, sliced,
                           scope_annot(active));
  // Append any schedule-derived per-node metadata (e.g. the value's remat
  // home + use scopes -- properties the running annotation cannot see). Empty
  // provider => nothing appended => byte-identical to the base annotation.
  if (auto const& nm = Logger::instance().eval.node_meta; nm)
    annot += " " + nm(node->hash_value());
  return annot;
}

/// Enriched node-info trailer for trace emission: the node label (see the
/// single-arg overload) followed by the per-op slice/home-scope metadata (see
/// slice_home_annot), so every compute and cache-read trace line carries which
/// of the node's modes are sliced by the active batch loops and the loop-nest
/// scope. Centralized here so the forest and DAG paths reach it uniformly
/// through the one node-info trailer. Only invoked under printing() (at the
/// trace emission sites), so the canon/sliced formatting cost is trace-only.
template <meta::eval_node Node, typename BatchContext>
[[nodiscard]] std::string label(Node const& node, BatchContext const& active) {
  return std::format("{} {}", label(node), slice_home_annot(node, active));
}

}  // namespace log

// implementation details of the eval engine; prefer sequant::detail over an
// unnamed namespace in a header (see CppCoreGuidelines SF.21)
namespace detail {

///
/// Invokes @c fun that returns void on the arguments @c args and returns the
/// time duration as @c std::chrono::duration<double>.
template <typename F, typename... Args>
[[nodiscard]] log::Duration timed_eval_inplace(F&& fun, Args&&... args)
  requires(std::is_invocable_r_v<void, F, Args...>)
{
  using Clock = std::chrono::high_resolution_clock;
  auto tstart = Clock::now();
  std::forward<F>(fun)(std::forward<Args>(args)...);
  auto tend = Clock::now();
  return {tend - tstart};
}

template <typename T>
constexpr bool is_cache_manager_v = false;

template <typename N, bool F>
constexpr bool is_cache_manager_v<CacheManager<N, F>> = true;

// True if ANY of Args is a CacheManager. The fold over `||` is well-defined
// (and false) for an empty pack, so there is no tuple_element underflow to
// guard against. Used to keep the variadic cache-appending
// evaluate()/evaluate_impl() forwarders below from matching a call that ALREADY
// carries a cache -- e.g. the scope-executor overload
// evaluate(forest, policy, layout, leaf, cache, mode_order, guard), whose
// CacheManager is NOT the last argument. A last-argument-only check let the
// forwarder fire on that overload, append a second cache, and fail to resolve.
template <typename... Args>
concept any_type_is_cache_manager =
    (... || is_cache_manager_v<std::remove_cvref_t<Args>>);

template <typename... Args>
auto&& arg0(Args&&... args) {
  return std::get<0>(std::forward_as_tuple(std::forward<Args>(args)...));
}

auto&& node0(auto&& val) { return std::forward<decltype(val)>(val); }
auto&& node0(std::ranges::range auto&& rng) {
  return ranges::front(std::forward<decltype(rng)>(rng));
}

enum struct CacheCheck { Checked, Unchecked };

}  // namespace detail

enum struct Trace {
  On,
  Off,
  Default =
#ifdef SEQUANT_EVAL_TRACE
      On
#else
      Off
#endif
};
static_assert(Trace::Default == Trace::On || Trace::Default == Trace::Off);

// implementation details of the eval engine; prefer sequant::detail over an
// unnamed namespace in a header (see CppCoreGuidelines SF.21)
namespace detail {
[[nodiscard]] consteval bool trace(Trace t) noexcept { return t == Trace::On; }
}  // namespace detail

/// \brief The indices contracted at a binary evaluation node.
///
/// These are the indices present in *both* children's (canonical) result
/// indices but absent from the node's own result indices -- i.e. the indices
/// summed over by this node's product. Empty for leaves, for sums, and for
/// products with no contracted index (e.g. a pure outer/Hadamard product).
///
/// Each such index `K` is a valid mode to evaluate the subtree rooted at
/// `node` in batches: the node computes `R = sum_K f(K)`, so
/// `R = sum_{blocks b} sum_{K in b} f(K)` -- evaluating per-block and summing
/// bounds the peak memory of `K`-carrying intermediates in the subtree. A
/// custom evaluator (see CacheManager::custom_evaluator_type) can use this to
/// implement batched evaluation.
[[nodiscard]] inline Index::index_vector contracted_indices(
    meta::eval_node auto const& node) {
  Index::index_vector result;
  if (node.leaf() || !node->is_product()) return result;
  auto const& l = node.left()->canon_indices();
  auto const& r = node.right()->canon_indices();
  auto const& c = node->canon_indices();
  auto contains = [](auto const& vec, Index const& ix) {
    return std::find(vec.begin(), vec.end(), ix) != vec.end();
  };
  for (Index const& ix : l)
    if (contains(r, ix) && !contains(c, ix)) result.push_back(ix);
  return result;
}

/// \brief A default mode to batch the subtree at \p node over: the contracted
/// index (see contracted_indices) that satisfies \p accept, choosing the one
/// with the largest IndexSpace approximate size -- typically the auxiliary/RI
/// index, whose elimination most reduces the peak intermediate.
///
/// \param accept a predicate `bool(Index const&)` selecting which contracted
///        indices are eligible to batch over (e.g. only those in a given
///        IndexSpace). This lets a caller scope batching to specific modes.
/// \return nullopt if no contracted index satisfies \p accept.
template <typename IndexPredicate>
[[nodiscard]] inline std::optional<Index> batch_axis(
    meta::eval_node auto const& node, IndexPredicate const& accept) {
  std::optional<Index> best;
  for (Index const& ix : contracted_indices(node)) {
    if (!accept(ix)) continue;
    if (!best ||
        best->space().approximate_size() < ix.space().approximate_size())
      best = ix;
  }
  return best;
}

/// \overload Batches over any contracted index (largest approximate size).
[[nodiscard]] inline std::optional<Index> batch_axis(
    meta::eval_node auto const& node) {
  return batch_axis(node, [](Index const&) { return true; });
}

// index_position() moved to SeQuant/core/eval/slicing_signature.hpp (shared
// with the hoist-path split and router-read guard); still available here via
// that include.

/// \return the first leaf in the subtree rooted at \p node whose canonical
///         indices contain \p ix, paired with the position of \p ix there; or
///         nullopt if no such leaf. Used to learn \p ix's tile structure from
///         a tensor that carries it.
template <typename Node>
[[nodiscard]] std::optional<std::pair<Node, std::size_t>> find_leaf_carrying(
    Node const& node, Index const& ix) {
  if (node.leaf()) {
    if (auto const p = index_position(node, ix)) return std::pair{node, *p};
    return std::nullopt;
  }
  if (auto found = find_leaf_carrying(node.left(), ix)) return found;
  return find_leaf_carrying(node.right(), ix);
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \tparam Cache If CacheCache::Checked (default) the root \p node is looked up
///               in \p cache before evaluating; a hit short-circuits it. Child
///               nodes are always evaluated Checked. Unchecked skips the lookup
///               for the root only.
///
/// \note The traversal is iterative: it maintains its own explicit work stack
///       (a `std::deque` of frames) rather than recursing, so evaluation depth
///       is bounded by the heap, not the C++ call stack. This keeps deep trees
///       (e.g. a Sum or product chain with thousands of operands) stack-safe.
///       Custom-evaluator interception is preserved exactly: it is consulted
///       when a frame is first visited and a non-null result short-circuits the
///       subtree (its children are never pushed), so subtree pruning -- the
///       mechanism batched eval relies on -- is unaffected.
///
/// \param node A node that can be evaluated using \p leaf_evaluator as the leaf
///             evaluator.
/// \param leaf_evaluator The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
// The recursive evaluation engine. NOT the top-level entry: the batched
// evaluator re-enters THIS (evaluate_impl), never `evaluate`, so `evaluate`
// reads as the single outermost call and the whole batched recursion is
// contained in evaluate_impl. External callers use the thin `evaluate`
// overloads (below); internal re-entries (the scatter/contraction per-block
// re-evaluations and the hoisted-invariant builds) call evaluate_impl directly.

/// \brief Pillar 1 / B-full stage A: the single-op COMPUTE kernel.
///
/// \details The raw op applied to already-evaluated operand results, dispatched
/// by \p node's op type, with the contraction annotation computed from \p node.
/// VALUE-in, VALUE-out: it takes the operands' \c Result buffers and returns
/// the op's raw result, touching NO cache and NO value identity. It is exactly
/// the innermost `left->adjoint/sum/prod(...)` compute that both the
/// tree-walking
/// \c evaluate_impl and the value/occurrence-driven ordered executor perform,
/// so extracting it is a byte-identical seam. NOT handled here (all
/// caller-side): a LEAF (a \c leaf_evaluator fetch, not an op), the
/// shaped-product hook (the caller calls this only when the hook declines), \c
/// apply_phase + store, the in-place-Sum fast path, and all tally / trace /
/// timing / \c last_op_flops sentinel / \c force_sync fence.
template <meta::can_evaluate Node>
[[nodiscard]] ResultPtr apply_one_op(Node const& node, ResultPtr const& left,
                                     ResultPtr const& right) {
  if (node->op_type() == EvalOp::Adjoint) {
    // Unary: only the left operand; the right child is the Constant(1)
    // sentinel.
    std::array<std::any, 2> const adj_ann{node.left()->annot(), node->annot()};
    return left->adjoint(adj_ann);
  }
  std::array<std::any, 3> const ann{node.left()->annot(), node.right()->annot(),
                                    node->annot()};
  if (node->op_type() == EvalOp::Sum) return left->sum(*right, ann);
  SEQUANT_ASSERT(node->op_type() == EvalOp::Product);
  bool const de_nest =
      node.left()->tot() && node.right()->tot() && !node->tot();
  return left->prod(*right, ann, de_nest ? DeNest::True : DeNest::False);
}

template <Trace EvalTrace = Trace::Default,
          detail::CacheCheck Cache = detail::CacheCheck::Checked,
          meta::can_evaluate Node, typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate_impl(Node const& node,         //
                        F const& leaf_evaluator,  //
                        CacheManager<N, FHC>& cache) {
  // DIAGNOSTIC (SEQUANT_UT_EVALIMPL): split CC-eval into time inside
  // top-level evaluate_impl (body) vs the gap between successive top-level
  // calls (executor/call-site machinery). See EvalImplTimeline.
  eval::EvalImplTimeline::Scope _tl_evalimpl;
  // Multiply a (possibly cached) result by its node's canonicalization phase.
  // Formerly the `mult_by_phase` lambda local to the Checked wrapper.
  auto apply_phase = [&cache](auto const& nd, ResultPtr res) -> ResultPtr {
    auto phase = nd->canon_phase();
    if (phase == 1) return res;

    ResultPtr post;
    auto const _ph0 = std::chrono::steady_clock::now();
    auto time =
        detail::timed_eval_inplace([&]() { post = res->mult_by_phase(phase); });
    eval::EvalImplTimeline::note_phase(_ph0);

    if constexpr (detail::trace(EvalTrace)) {
      size_t hwmark = log::bytes(cache, post).value;
      if (!cache.alive(nd)) hwmark += log::bytes(res).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{
          .mode = log::EvalMode::MultByPhase,
          .time = time,
          .mem_result = log::bytes(post),
          .mem_alloc = log::bytes(post),
          .mem_hwmark = {cache.note_working_set(hwmark, nd->hash_value())}};
      log::eval(stat, std::format("{} * {}", phase, nd->label()),
                log::slice_home_annot(nd, cache.batch_context()));
    }
    return post;
  };

  // Slice-on-use: slice a value fetched at the Enter stage to the current batch
  // block for the `hops` INNERMOST enclosing batch loops that `nd` carries and
  // that the value does not yet have baked in. `hops` == number of enclosing
  // loops crossed to reach the value's lifetime scope (0 for a local hit / a
  // freshly built value in this scope; d == batch_context().size() for a fresh
  // leaf whose lifetime is top). The slice set is exactly
  // (use scope MINUS lifetime scope) INTERSECT carried(nd): the `hops`
  // innermost batch_context entries, filtered by index_position(nd, axis).
  // slice_mode is non-mutating, so a cached full value is left undisturbed.
  // Empty batch_context (the OFF path) => d == hops == 0 => the loop is empty
  // and the value is returned unchanged, byte-identical to the pre-slice-on-use
  // path.
  auto slice_to_use = [&cache](ResultPtr value, auto const& nd,
                               std::size_t hops) -> ResultPtr {
    auto const& ctx = cache.batch_context();
    std::size_t const d = ctx.size();
    // hops (parent links access_at crossed) must not exceed d (batch_context
    // entries): each realized loop pushes exactly one entry and wires at most
    // one parent link, so hops <= d always. A violation would underflow
    // `d - hops` and silently UNDER-slice (oversized result); assert loudly.
    SEQUANT_ASSERT(hops <= d);
    bool const _home_diag = std::getenv("SEQUANT_UT_HOME_DIAG") != nullptr;
    if (_home_diag && d > 0) {
      std::cerr << "[HOME] node#" << (nd->hash_value() % 100000u)
                << " use_d=" << d << " hops=" << hops
                << " home_d=" << (d - hops) << " scope=[";
      for (std::size_t j = 0; j < d; ++j)
        std::cerr << toUtf8(std::wstring(ctx[j].axis.space().base_key()))
                  << "@o" << ctx[j].level.latitude_ordinal
                  << (j < d - hops ? "(home) " : "(x) ");
      std::cerr << "] canon=[";
      for (auto const& ix : nd->canon_indices())
        std::cerr << toUtf8(ix.full_label()) << " ";
      std::cerr << "]" << std::endl;
    }
    for (std::size_t i = d - hops; i < d; ++i) {
      auto const& axis = ctx[i].axis;
      auto const& blk = ctx[i].range;
      // Diagnostic-only exact match of this batch-context axis on nd (feeds the
      // [SLICE] trace below). The ACTUAL slice mode is p_new, resolved off the
      // loop-colored seam further down -- eval never deduces a slice mode from
      // an index's space, so a node that does not carry the exact axis simply
      // does not match here.
      auto p = index_position(nd, axis);
      bool const exact = p.has_value();
      if (p && std::getenv("SEQUANT_UT_SLICE_DIAG"))
        std::cerr << "[SLICE] node#" << (nd->hash_value() % 100000u)
                  << " space=" << toUtf8(axis.space().base_key())
                  << " pmode=" << *p
                  << " via=" << (exact ? "exact" : "fallback") << " idx@pmode="
                  << toUtf8(nd->canon_indices()[*p].full_label()) << " slice=["
                  << blk.first << "," << blk.second << ")" << std::endl;
      // p_new (Task 7, sliced-value canonical-layout / loop-coloring design):
      // the ACTUAL slice mode this entry resolves to, per path.
      //  - `exact_axis` present (forest / whole-scope, which push each
      //    member's OWN physical axis): the same intra-tree exact match as old
      //    `p` (exact_axis == axis on those paths), byte-identical to before.
      //  - else (the ordered executor, which pushes ONE canonical block.axis
      //    per type-bucketed loop): resolve OFF THE LOOP-COLORED CANONICAL
      //    LAYOUT (design sec.4) rather than the per-cell mode_to_level map --
      //    the loop-colored slice seam maps this loop (ctx[i].level) to the
      //    fetched value's OWN sliced-mode Index, whose physical slot on `nd`
      //    is then read directly (index_position). The seam returns the mode
      //    nd actually carries for this loop; there is NO space-deduction
      //    fallback -- eval never guesses a slice mode. A value with no seam
      //    entry for this loop (an unsliced value,
      //    or one built inside this loop -- the built-within participation
      //    gate) leaves p_new nullopt => full/unsliced. If the seam is unwired
      //    (a direct caller that did not set it) p_new likewise stays nullopt;
      //    the ordered executor always wires it in practice.
      std::optional<std::size_t> p_new;
      if (ctx[i].exact_axis) {
        p_new = index_position(nd, *ctx[i].exact_axis);
      } else if (auto const* seam = cache.loop_colored_slice_seam()) {
        if (std::optional<LoopId> const loop =
                seam->loop_of_level(ctx[i].level)) {
          // FRAME-CORRECT (2026-08-24 slot-slicing design): mode_of returns the
          // sliced mode's PHYSICAL POSITION, computed at schedule time in the
          // fetched value's OWN index-frame -- per-occurrence via the consumer
          // for a divergent (relabeled) or symmetric value. It is used directly
          // as the slice mode here: NO index_position, and NO re-matching a
          // canonical Index label against this (differently-framed) node.
          p_new =
              seam->mode_of(nd->hash_value(), *loop, cache.current_consumer());
          // COMPLETENESS GUARD (2026-08-25 loop-open design): the value
          // PARTICIPATES in this loop -- the seam has a sliced-mode fact for it
          // under `loop` for SOME consumer -- yet THIS fetch got none. That is
          // a scheduler gap: the operand would be served UNSLICED while its
          // contraction partner is sliced, and TA's tiled-range assert is
          // elided in Release, so the mismatched DistEval would DEADLOCK
          // silently instead of erroring. Fail loud. A value with NO fact under
          // `loop` (participates() == false) is genuinely invariant to it and
          // is correctly left unsliced -- the guard does not fire there.
          if (!p_new && seam->participates(nd->hash_value(), *loop))
            throw Exception(
                "slice_to_use: value participates in a batch loop but this "
                "fetch has no sliced-mode fact for the current consumer -- an "
                "incomplete sliced-mode assignment (occ_facts gap) that would "
                "leave the operand unsliced and mismatch its contraction "
                "partner. The occurrence's slice fact must be recorded at "
                "schedule time, never left to a silent unsliced fetch.");
        }
      }
      // Task 7-part-2 / Task 8: the transitional equivalence assert this used
      // to cross-check against the old per-cell mode_to_level map is RETIRED,
      // and that map and its populators are now deleted entirely (Task 8).
      // With per-consumer binding, the loop-colored resolution
      // INTENTIONALLY diverges from what the old consumer-blind map would
      // have given for the symmetric case -- that divergence IS the fix.
      if (_home_diag)
        std::cerr << "  [HOME-SLICE] node#" << (nd->hash_value() % 100000u)
                  << " crossed-axis="
                  << toUtf8(std::wstring(ctx[i].axis.space().base_key()))
                  << "@o" << ctx[i].level.latitude_ordinal
                  << " -> p_new=" << (p_new ? std::to_string(*p_new) : "none")
                  << " blk=[" << blk.first << "," << blk.second << ")"
                  << std::endl;
      if (p_new) {
        // Test-only witness of the consumer-aware slice decision (PILLAR 2):
        // no-op unless an observer is installed.
        if (auto const& obs = cache.slice_observer(); obs)
          obs(nd->hash_value(), cache.current_consumer(), *p_new,
              ctx[i].level.latitude_ordinal);
        value = value->slice_mode(*p_new, blk.first, blk.second);
      }
    }
    return value;
  };

  // One entry of the explicit evaluation stack. `stage` records how far a node
  // has progressed; `left`/`right` hold its evaluated operands; `store_after`
  // marks a Checked node that exists in the cache map but has not been stored
  // yet, so its computed result must be cached (this replaces the recursive
  // wrapper's `evaluate<..., Unchecked>` re-entry).
  enum class Stage { Enter, NeedLeft, NeedRight, NeedLeftAdj };
  struct Frame {
    Node node;
    bool checked;
    Stage stage = Stage::Enter;
    bool store_after = false;
    ResultPtr left, right;
  };

  // Finalize a freshly computed Phase-B result: if this Checked node needs
  // storing, cache it (phase-applied) and hand back the phase-applied cached
  // pointer -- exactly the recursive Checked wrapper's store path. Otherwise
  // pass the raw result through unchanged.
  auto finish_phase_b = [&cache, &apply_phase](Frame const& f,
                                               ResultPtr rb) -> ResultPtr {
    // DIAGNOSTIC (SEQUANT_UT_BUILD_METER): count actual builds at this single
    // chokepoint (non-leaf only = real contraction executions). No-op when off.
    if (!f.node.leaf())
      eval::BuildMeter::on_build(f.node->hash_value(),
                                 eval::BuildMeter::enabled()
                                     ? log::label(f.node, cache.batch_context())
                                     : std::string{});
    // Per-op BUILD event: finish_phase_b is the single choke point every
    // freshly computed node passes through -- leaves, custom-eval subtrees, and
    // standard contractions -- so this counts EVERY build (cached or not),
    // giving the schedule visualizer full per-node recompute coverage (the
    // cache Store/ Access/Release events above cover only cached nodes). Keyed
    // by hash_value() to join onto the IR DAG. Emitted to the cache's
    // ScheduleSink when one is wired (set_schedule_sink); no sink => no dump.
    if (auto* const sink = cache.schedule_sink();
        sink && sink->os && !sink->fired) {
      std::ostream& os = *sink->os;
      // ctx = the active batch loops (mode -> block offset) this build ran
      // under. The visualizer counts DISTINCT ctx projections onto the modes a
      // node depends on: builds at a repeated projected slice are avoidable
      // recompute (a value rebuilt where an enclosing loop it is invariant to
      // advanced) vs the inherent per-block work batching requires.
      // Leaves are ACCESSED, not computed (leaf_evaluator just hands back a
      // ref to a precomputed input), so tag them "Fetch" -- cheap, not
      // recompute. Only internal nodes (contractions) are real "Build" work.
      char const* const bmode = f.node.leaf() ? "Fetch" : "Build";
      os << "SCHEDULE_RUN_EVENT {\"hash\":\"" << f.node->hash_value()
         << "\",\"mode\":\"" << bmode << "\"";
      // sig = the avoidable-recompute join key, matching the IR node record and
      // cost_profile's per-node label (result + sorted operand pair). Internal
      // nodes only; leaves (Fetch) are not tallied. Lets the visualizer join
      // this hash to cost_profile's per-node avoidable without reconstructing
      // the signature in the renderer.
      // Per-build flops (cm->flops) keyed to this node by hash above; the
      // recompute rollup and the visualizer join on the topological hash, not
      // a dummy-/slice-dependent signature.
      if (!f.node.leaf()) os << ",\"flops\":" << eval::detail::last_op_flops();
      os << ",\"ctx\":[";
      bool first = true;
      for (auto const& entry : cache.batch_context()) {
        Index const& ix = entry.axis;
        auto const& blk = entry.range;
        // dep = does the node's subtree carry this loop mode (free or
        // contracted below)? If not, the node is INVARIANT to it and rebuilding
        // per block of it is avoidable recompute. find_leaf_carrying works in
        // the node's own label space, so no alpha-renaming reconciliation.
        bool const dep = find_leaf_carrying(f.node, ix).has_value();
        os << (first ? "" : ",") << "[\"" << toUtf8(ix.full_label()) << "\","
           << blk.first << "," << (dep ? 1 : 0) << "]";
        first = false;
      }
      os << "]}\n";
    }
    if (!f.store_after) return rb;
    auto ptr = cache.store(f.node, apply_phase(f.node, std::move(rb)));
    if constexpr (detail::trace(EvalTrace))
      log::cache(f.node, cache, log::label(f.node, cache.batch_context()));
    return apply_phase(f.node, ptr);
  };

  // A `std::deque` is used so that a reference to the top frame stays valid
  // across push_back (which reallocates a `std::vector`).
  std::deque<Frame> stk;
  stk.push_back(
      Frame{.node = node, .checked = (Cache == detail::CacheCheck::Checked)});

  ResultPtr ret;  // result handed up by the frame that most recently finalized

  // Deliver `r` to the parent frame and pop the just-completed frame.
  auto finalize = [&stk, &ret](ResultPtr r) {
    ret = std::move(r);
    stk.pop_back();
  };

  while (!stk.empty()) {
    Frame& f = stk.back();

    switch (f.stage) {
      case Stage::Enter: {
        // --- Checked cache wrapper: a hit returns directly; a miss on a node
        //     that exists in the map schedules a store once computed. ---
        if (f.checked) {
          // --- Router consult: an override seam ahead of the default
          //     access_at() below (see placement_router.hpp). The
          //     `router && !router->empty()` short-circuit is FIRST so an
          //     empty/null router (the Phase 2 default, and every production
          //     eval) never computes an occurrence_key -- zero hot-path cost
          //     and byte-identical behavior (the `routed` flag stays false,
          //     and the default path immediately below is untouched). ---
          bool routed = false;
          // Only a value the remat pass actually MOVED can have a router
          // override -- and the build keys occurrence_key ONLY for such nodes
          // (placement_remat.hpp). Gate on moved() first: it is a hash-set
          // lookup (no occurrence_key), so a non-moved node route()-misses
          // exactly as before (byte-identical), but we never compute an
          // occurrence_key for a node the router could not have keyed. This
          // also keeps occurrence_key off Sum (and other non-tensorial) nodes,
          // which are never moved and are not tensor networks.
          if (auto const* router = cache.placement_router();
              router && !router->empty() &&
              router->moved(f.node->hash_value())) {
            auto const& ctx = cache.batch_context();
            container::svector<Index> ctx_modes;
            for (auto const& e : ctx) ctx_modes.push_back(e.axis);
            auto const key = eval::occurrence_key(f.node, ctx_modes);
            if (auto const* home = router->route(key)) {
              std::size_t const use_depth = ctx.size();
              std::size_t const hd = router->home_depth(*home, ctx, key);
              SEQUANT_ASSERT(hd <= use_depth);
              // LIVE defense-in-depth guard (design sec.4, RELEASE-safe): the
              // router-directed fetch keys by CANONICAL hash at the resolved
              // scope, so serving is only sound if the resolution is CONSISTENT
              // with THIS occurrence (home_resolution_consistent -- the live
              // loop at hd is one of this occurrence's own batched indices
              // whose space the overlay names). By construction home_depth
              // returns only such an hd, so on every correct schedule this
              // holds and the fetch proceeds byte-identically. If a future
              // overlay/home_depth regression resolved an occurrence to a scope
              // it does not bind (collapsing two divergently-relabeled
              // occurrences onto ONE entry -- the hole the DAG-scope resolution
              // closes), the share is REFUSED here: `routed` stays false and
              // the default access_at path below serves this occurrence its OWN
              // value (recompute), never a wrong-slice entry -- in release, not
              // only under an assert. The hoist split keeps consistent
              // occurrences shareable.
              if (router->home_resolution_consistent(*home, ctx, hd, key)) {
                std::size_t const hops = use_depth - hd;
                if (ResultPtr ptr = cache.access_at_hops(f.node, hops); ptr) {
                  if constexpr (detail::trace(EvalTrace))
                    log::cache(f.node, cache,
                               log::label(f.node, cache.batch_context()));
#ifdef SEQUANT_ROUTER_SHADOW
                  // Dev-only correctness check (default OFF): for a NO-OP
                  // override (residency == the value's ACTUAL current home),
                  // the router-directed fetch must reproduce access_at()
                  // pointer-for-pointer. Not safe to enable in production: this
                  // extra access_at() call decays a non-persistent entry's
                  // lifetime a second time.
                  {
                    auto const shadow = cache.access_at(f.node);
                    SEQUANT_ASSERT(shadow.ptr.get() == ptr.get());
                  }
#endif
                  finalize(
                      slice_to_use(apply_phase(f.node, ptr), f.node, hops));
                  routed = true;
                }
              }
            }
          }
          if (routed) break;
          if (auto m = cache.access_at(f.node); m.ptr) {
            if constexpr (detail::trace(EvalTrace))
              log::cache(f.node, cache,
                         log::label(f.node, cache.batch_context()));
            // Slice-on-use: a value fetched `m.hops` scopes up does not have
            // this scope's (and any intervening) batch slices baked in, so
            // slice it to the current block for the loops the fetch crossed. A
            // local hit (hops == 0) or the OFF path (empty batch_context) is a
            // no-op, so this stays byte-identical to apply_phase() alone there.
            finalize(slice_to_use(apply_phase(f.node, m.ptr), f.node, m.hops));
            break;
          }
          // Resident-reads invariant (read-from-home ordered scratches): a
          // non-leaf value other than this call's own top node MUST already be
          // resident -- it was statically scheduled and built by a prior step.
          // A miss here means it vanished (premature eviction / under-predicted
          // use count in ordered_home_reads), which must be a hard error, never
          // a silent recompute or empty-array serve that hangs a downstream
          // contraction. Leaves (evaluated fresh) and the top node (being built
          // now) legitimately miss.
          if (cache.require_resident_reads() && !f.node.leaf() &&
              f.node->hash_value() != node->hash_value()) {
            std::string lbl;
            for (auto const& ix : f.node->canon_indices())
              lbl += toUtf8(ix.full_label()) + " ";
            throw Exception(
                "evaluate_impl: a read-from-home value vanished before use "
                "(premature eviction -- likely an under-predicted use count in "
                "ordered_home_reads); a missing value must never be served as "
                "an empty array. canon=[" +
                lbl + "]");
          }
          f.store_after = cache.exists(f.node);
        }

        // --- Custom-evaluator interception (non-leaf only): a non-null result
        //     short-circuits the subtree -- children are never pushed. This is
        //     the subtree pruning batched eval relies on; see the class note. A
        //     null return declines to the standard scheme below. ---
        if (!f.node.leaf()) {
          if (auto const& custom_eval = cache.custom_evaluator(); custom_eval) {
            ResultPtr intercepted;
            auto time = detail::timed_eval_inplace(
                [&]() { intercepted = custom_eval(f.node, cache); });
            if (intercepted) {
              if constexpr (detail::trace(EvalTrace)) {
                size_t hwmark = log::bytes(cache, intercepted).value;
                hwmark +=
                    cache.parent() ? cache.parent()->chain_residency() : 0;
                log::eval(log::EvalStat{.mode = log::eval_mode(f.node),
                                        .time = time,
                                        .mem_result = log::bytes(intercepted),
                                        .mem_alloc = log::bytes(intercepted),
                                        .mem_hwmark = {cache.note_working_set(
                                            hwmark, f.node->hash_value())}},
                          log::label(f.node, cache.batch_context()));
              }
              log::release_after_op();
              finalize(finish_phase_b(f, std::move(intercepted)));
              break;
            }
          }
        }

        // --- Leaf. ---
        if (f.node.leaf()) {
          ResultPtr result;  // the FULL leaf (traced and cached full)
          auto time = detail::timed_eval_inplace(
              [&]() { result = leaf_evaluator(f.node); });
          if constexpr (detail::trace(EvalTrace)) {
            size_t hwmark = log::bytes(cache, result).value;
            hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
            log::eval(log::EvalStat{.mode = log::eval_mode(f.node),
                                    .time = time,
                                    .mem_result = log::bytes(result),
                                    .mem_alloc = log::bytes(result),
                                    .mem_hwmark = {cache.note_working_set(
                                        hwmark, f.node->hash_value())}},
                      log::label(f.node, cache.batch_context()));
          }
          log::release_after_op();
          // Store the FULL leaf under its canonical key (a block slice would
          // corrupt the cache), then return it SLICED to the current block: a
          // freshly built leaf's lifetime is top, so every enclosing carried
          // batch loop is unbaked and must be sliced (hops == batch_context
          // size). This reproduces the old per-block leaf slicing (le_g) on the
          // main value path; the OFF path (empty batch_context) is a no-op.
          ResultPtr stored = finish_phase_b(f, std::move(result));
          finalize(slice_to_use(stored, f.node, cache.batch_context().size()));
          break;
        }

        // --- Internal node: request the left operand (always Checked). The
        //     stage must advance before the push (push may grow the deque). ---
        f.stage = (f.node->op_type() == EvalOp::Adjoint) ? Stage::NeedLeftAdj
                                                         : Stage::NeedLeft;
        stk.push_back(Frame{.node = f.node.left(), .checked = true});
        break;
      }

      case Stage::NeedLeftAdj: {
        // Unary IR op (Adjoint): only the left operand is evaluated; the right
        // child is the Constant(1) sentinel kept to preserve FullBinaryNode's
        // invariant, and is intentionally never pushed.
        f.left = std::move(ret);
        SEQUANT_ASSERT(f.left);
        ResultPtr result;
        auto time = detail::timed_eval_inplace(
            [&]() { result = apply_one_op(f.node, f.left, f.right); });

        if constexpr (detail::trace(EvalTrace)) {
          // `right` is null here (see log::bytes() null tolerance).
          size_t hwmark = log::bytes(cache, result).value;
          if (!cache.chain_holds(f.left)) hwmark += log::bytes(f.left).value;
          hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
          log::eval(log::EvalStat{.mode = log::eval_mode(f.node),
                                  .time = time,
                                  .mem_result = log::bytes(result),
                                  .mem_alloc = log::bytes(result),
                                  .mem_hwmark = {cache.note_working_set(
                                      hwmark, f.node->hash_value())},
                                  .mem_left = log::bytes(f.left),
                                  .mem_right = log::bytes(f.right)},
                    log::label(f.node, cache.batch_context()));
        }
        log::release_after_op();
        finalize(finish_phase_b(f, std::move(result)));
        break;
      }

      case Stage::NeedLeft: {
        f.left = std::move(ret);
        SEQUANT_ASSERT(f.left);
        f.stage = Stage::NeedRight;
        stk.push_back(Frame{.node = f.node.right(), .checked = true});
        break;
      }

      case Stage::NeedRight: {
        f.right = std::move(ret);
        SEQUANT_ASSERT(f.left);
        SEQUANT_ASSERT(f.right);

        std::array<std::any, 3> const ann{
            f.node.left()->annot(), f.node.right()->annot(), f.node->annot()};
        ResultPtr result;
        log::Duration time;
        // In-place accumulation is eligible only when f.left's PROVENANCE is
        // known to be an evaluation-local, exclusively-owned buffer:
        //  - !f.node.left().leaf(): a leaf's ResultPtr comes straight out of
        //    the caller-supplied leaf_evaluator, whose provenance this engine
        //    cannot see -- a memoizing evaluator (the norm for AO integrals /
        //    amplitudes reused across calls, e.g. rand_tensor_yield in the
        //    unit tests) can hand out the SAME buffer to unrelated callers,
        //    so mutating it here would silently corrupt those other reads.
        //    Only an INTERNAL node's result is guaranteed freshly built by
        //    THIS evaluation (prod()/sum()/permute() always allocate), so
        //    only internal nodes are safe candidates.
        //  - !cache.chain_holds_shared(f.left): even an internal node's result
        //    must not be a live, SHARED entry on the CacheManager's scope chain
        //    -- i.e. not referenced again elsewhere in THIS evaluation's
        //    tree/forest (multi-use, so some other read is pending). This is a
        //    RUNTIME condition, not merely an assert: under this Release build
        //    SEQUANT_ASSERT expands to a no-op (SEQUANT_ASSERT_ENABLED
        //    undefined), so a guard that only asserted would be ELIDED and the
        //    mutation would silently corrupt a value shared across roots (see
        //    the ordered-executor multi-root path, where a subexpression CSE'd
        //    across two independent roots is homed resident and read by both).
        //    chain_holds_shared() tests, by pointer identity, whether a LIVE
        //    entry with more than one consumer (max_life > 1 -- a resident home
        //    or a not-yet-drained multi-use CSE entry) still holds f.left; a
        //    private single-use accumulator is NOT reported (a transient
        //    running total is held by no entry, and a single-use CSE entry
        //    moves its buffer out on its sole read, so it no longer holds it),
        //    so in-place still fires for the private common case -- see
        //    CacheManager::chain_holds_shared's own doc comment.
        // A marked Sum whose left child IS a leaf (the common case for the
        // innermost Sum of a chain, whose left is the chain seed), OR whose
        // left operand is a shared cache-resident value, therefore falls back
        // to the allocating sum() below despite being marked -- one bounded
        // extra allocation, not per term -- and every OTHER (non-leaf-seeded,
        // private) Sum in the chain still accumulates in place from there on,
        // since each already-computed running total is a fresh,
        // evaluation-local buffer.
        bool const inplace_eligible =
            f.node->op_type() == EvalOp::Sum && f.node->accumulate_in_place() &&
            !f.node.left().leaf() && !cache.chain_holds_shared(f.left);
        if (inplace_eligible) {
          SEQUANT_ASSERT(!cache.chain_holds_shared(f.left));

          // The accumulator (f.left) is used AS IS, in its own layout --
          // never repermuted -- since binarize() pins a marked Sum's
          // canon_indices_ to its LEFT operand's (see EvalExpr::binarize
          // (Sum)'s make_sum lambda), so this node's own target layout IS
          // f.left's layout. The right addend generally is NOT already in
          // that layout (each summand of the original N-ary Sum
          // canonicalizes independently), so -- exactly as the allocating
          // sum() path below permutes BOTH operands into the result's
          // layout (see e.g. ResultTensorBTAS::sum) -- it must be permuted
          // into f.left's layout (ann[0]) before the raw elementwise
          // add_inplace; skipping this for a tensor Sum silently adds
          // mismatched layouts. Scalars carry no layout: ResultScalar::
          // permute() is unimplemented (throws) and none is needed, since
          // ann[0]/ann[1]/ann[2] are all trivially empty for a scalar Sum.
          bool const needs_permute =
              f.node->result_type() == ResultType::Tensor;
          ResultPtr right_aligned = f.right;
          log::Duration perm_time{};
          if (needs_permute) {
            perm_time = detail::timed_eval_inplace([&]() {
              right_aligned =
                  f.right->permute(std::array<std::any, 2>{ann[1], ann[0]});
            });
          }
          if constexpr (detail::trace(EvalTrace)) {
            if (needs_permute) {
              // Mirrors the top-level evaluate(node, layout, ...) Permute
              // event above: a genuine fresh allocation for the (typically
              // much smaller, single-term) right addend, logged separately
              // from the SumInplace event below.
              size_t hwmark = log::bytes(cache, right_aligned).value;
              if (!cache.chain_holds(f.right))
                hwmark += log::bytes(f.right).value;
              hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
              log::eval(log::EvalStat{.mode = log::EvalMode::Permute,
                                      .time = perm_time,
                                      .mem_result = log::bytes(right_aligned),
                                      .mem_alloc = log::bytes(right_aligned),
                                      .mem_hwmark = {cache.note_working_set(
                                          hwmark, f.node->hash_value())}},
                        log::label(f.node, cache.batch_context()));
            }
          }

          // Move the accumulator out of f.left (leaving it null; it is not
          // read again below) and add the (now aligned) right addend into
          // it: zero allocation here, unlike the fresh buffer sum() below
          // returns.
          time = detail::timed_eval_inplace([&]() {
            result = std::move(f.left);
            result->add_inplace(*right_aligned);
          });

          if constexpr (detail::trace(EvalTrace)) {
            // Mirrors the forest-level SumInplace accounting (the
            // Nodes-range evaluate() overload, below): bytes(cache, result)
            // already counts the (mutated) former-accumulator buffer, since
            // result now IS that buffer, so only the (aligned) right addend
            // is added, and only if it is not already resident on the scope
            // chain. mem_alloc is zero -- SumInplace itself allocates
            // nothing (the right-alignment allocation, if any, was already
            // logged above as its own Permute event) -- and mem_left/
            // mem_right are left unset, matching the "SumInplace | --- |
            // 0B" row of the EvalStat doc table above.
            size_t hwmark = log::bytes(cache, result).value;
            if (!cache.chain_holds(right_aligned))
              hwmark += log::bytes(right_aligned).value;
            hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
            log::eval(log::EvalStat{.mode = log::EvalMode::SumInplace,
                                    .time = time,
                                    .mem_result = log::bytes(result),
                                    .mem_alloc = {0},
                                    .mem_hwmark = {cache.note_working_set(
                                        hwmark, f.node->hash_value())}},
                      log::label(f.node, cache.batch_context()));
          }
          log::release_after_op();
          finalize(finish_phase_b(f, std::move(result)));
          break;
        }
        if (f.node->op_type() == EvalOp::Sum) {
          time = detail::timed_eval_inplace(
              [&]() { result = apply_one_op(f.node, f.left, f.right); });
        } else {
          SEQUANT_ASSERT(f.node->op_type() == EvalOp::Product);
          // Consult the shaped-product hook (if set) before evaluating the
          // product. The hook receives the node (wrapped in a std::any as a
          // std::reference_wrapper so the full IR node is inspectable) plus the
          // evaluated operands and annotations; a non-null return *replaces*
          // the normal product (e.g. a shape-constrained emission of it), a
          // null return declines and the standard prod() below runs. An empty
          // hook is never consulted; default-empty => byte-identical behavior.
          // DIAGNOSTIC (SEQUANT_UT_FORCE_SYNC): force the product's async
          // execution to complete INSIDE the timed region so the prod() timer
          // captures execution, not just dispatch. No-op when the gate is off.
          static bool const force_sync =
              std::getenv("SEQUANT_UT_FORCE_SYNC") != nullptr;
          // DIAGNOSTIC (SEQUANT_UT_PROD_TR): log operand tiled ranges BEFORE
          // the contraction, so a product that DEADLOCKS inside the backend
          // einsum (a sliced-vs-unsliced mode mismatch whose TA_ASSERT is
          // elided in Release) still surfaces its operands -- the post-prod
          // trace never fires for it because prod() never returns.
          if (std::getenv("SEQUANT_UT_PROD_TR"))
            std::cerr << "[PROD-PRE] " << toUtf8(f.node->label()) << "\n    L."
                      << f.left->trange_annot() << "\n    R."
                      << (f.right ? f.right->trange_annot() : std::string{})
                      << std::endl;
          auto const _tp0 =
              std::chrono::steady_clock::now();  // node-eval start
          if (auto const& hook = cache.shaped_product_hook(); hook) {
            time = detail::timed_eval_inplace([&]() {
              result =
                  hook(std::any{std::cref(f.node)}, *f.left, *f.right, ann);
              if (force_sync && result) result->fence();
            });
          }
          if (!result) {
            // Sentinel so the recompute tally below fires ONLY when a DryRun
            // prod actually computed fresh flops for THIS node. DryRunOps::prod
            // early-returns without setting last_op_flops for a scalar*tensor
            // product (and it is never set at all by the wet TA backend); a
            // stale value from a previous op must not be attributed here (it
            // would fold garbage into the identity-keyed tally). prod sets
            // last_op_flops >= 0 only on the real contraction path.
            eval::detail::last_op_flops() = -1.0;
            time = detail::timed_eval_inplace([&]() {
              result = apply_one_op(f.node, f.left, f.right);
              if (force_sync && result) result->fence();
            });
            // Record this product build against f.node's IDENTITY (keyed by the
            // exact cache identity: hash-bin + Bliss, so 64-bit hash collisions
            // are NOT folded) in the (root) cache's recompute tally. Deduped at
            // the (value, SLICE) granularity using ACTUAL replay FLOPs (DryRun
            // prod just stashed this build's realized cost in last_op_flops): a
            // value built over DISTINCT slices is tiling (not recompute); the
            // same value rebuilt at the SAME slice -- including a node
            // invariant to an enclosing loop, rebuilt every block -- is
            // recompute. The slice signature is the live batch context
            // PROJECTED onto the modes f.node carries (find_leaf_carrying), so
            // an invariant's projection is identical (empty for that loop)
            // every block and its rebuilds fold. A no-op unless the dry-run
            // replay enabled the tally on the root cache (the wet TA path
            // leaves it disabled), so this is byte-identical off the costing
            // path.
            if (eval::detail::last_op_flops() >= 0.0) {
              std::string slice_sig;
              for (auto const& entry : cache.batch_context()) {
                Index const& ix = entry.axis;
                auto const& blk = entry.range;
                if (find_leaf_carrying(f.node, ix).has_value()) {
                  slice_sig += toUtf8(ix.full_label());
                  slice_sig += ':';
                  slice_sig += std::to_string(blk.first);
                  slice_sig += ';';
                }
              }
              cache.tally_build(f.node, slice_sig,
                                eval::detail::last_op_flops(),
                                eval::detail::last_op_exec());
            }
          }
          // node-eval done: fence + accumulate the contraction's full cost.
          eval::EvalImplTimeline::note_prod(_tp0);
        }

        SEQUANT_ASSERT(result);

        if constexpr (detail::trace(EvalTrace)) {
          // Skip an operand's bytes only when it aliases a cache buffer that is
          // already counted (locally in bytes(cache,...) or up-chain in
          // chain_residency()). chain_holds() tests pointer identity against
          // every alive entry on the scope chain: a cached child fetched full
          // aliases and is skipped; an apply_phase / sliced / permuted child is
          // a distinct buffer and is added.
          size_t hwmark = log::bytes(cache, result).value;
          if (!cache.chain_holds(f.left)) hwmark += log::bytes(f.left).value;
          if (f.right && !cache.chain_holds(f.right))
            hwmark += log::bytes(f.right).value;
          hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
          log::eval(
              log::EvalStat{.mode = log::eval_mode(f.node),
                            .time = time,
                            .mem_result = log::bytes(result),
                            .mem_alloc = log::bytes(result),
                            .mem_hwmark = {cache.note_working_set(
                                hwmark, f.node->hash_value())},
                            .mem_left = log::bytes(f.left),
                            .mem_right = log::bytes(f.right)},
              log::label(f.node, cache.batch_context()) + " | L." +
                  f.left->trange_annot() +
                  (f.right ? " R." + f.right->trange_annot() : std::string{}) +
                  " O." + result->trange_annot());
        }
        log::release_after_op();
        finalize(finish_phase_b(f, std::move(result)));
        break;
      }
    }
  }

  return ret;
}

/// Top-level single-node evaluation entry: a thin redirect to the recursive
/// engine \c evaluate_impl. Kept distinct so that `evaluate` denotes the
/// outermost call while all internal recursion (including the batched
/// evaluator's per-block re-entries) is spelled \c evaluate_impl.
template <Trace EvalTrace = Trace::Default,
          detail::CacheCheck Cache = detail::CacheCheck::Checked,
          meta::can_evaluate Node, typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate(Node const& node, F const& leaf_evaluator,
                   CacheManager<N, FHC>& cache) {
  return evaluate_impl<EvalTrace, Cache>(node, leaf_evaluator, cache);
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param node A node that can be evaluated using \p leaf_evaluator as the leaf
///             evaluator.
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
/// \param leaf_evaluator The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate Node, typename F,
          typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>  //
ResultPtr evaluate(Node const& node,           //
                   auto const& layout,         //
                   F const& leaf_evaluator,    //
                   CacheManager<N, FHC>& cache) {
  // if the layout is not the default constructed value need to permute
  bool const perm = layout != decltype(layout){};

  std::string xpr;
  if constexpr (detail::trace(EvalTrace)) {
    xpr = toUtf8(io::serialization::to_string(to_expr(node)));
    log::term(log::TermMode::Begin, xpr);
  }

  struct {
    ResultPtr pre, post;
  } result;

  result.pre = evaluate_impl<EvalTrace>(node, leaf_evaluator, cache);

  auto time = detail::timed_eval_inplace([&]() {
    result.post = perm ? result.pre->permute(
                             std::array<std::any, 2>{node->annot(), layout})
                       : result.pre;
  });

  SEQUANT_ASSERT(result.post);

  // logging
  if constexpr (detail::trace(EvalTrace)) {
    if (perm) {
      // result.pre aliases a cache buffer only when the inner evaluate returned
      // it unchanged (node cached at some scope, no mult_by_phase fresh alloc);
      // chain_holds() tests that by pointer identity across the scope chain. A
      // permuted/phase-shifted pre is a distinct buffer and is added.
      size_t hwmark = log::bytes(cache, result.post).value;
      if (!cache.chain_holds(result.pre))
        hwmark += log::bytes(result.pre).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{
          .mode = log::EvalMode::Permute,
          .time = time,
          .mem_result = log::bytes(result.post),
          .mem_alloc = log::bytes(result.post),
          .mem_hwmark = {cache.note_working_set(hwmark, node->hash_value())}};
      log::eval(stat, node->label(),
                log::slice_home_annot(node, cache.batch_context()));
    }
    log::term(log::TermMode::End, xpr);
  }
  return result.post;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param nodes A range of node that can be evaluated using \p leaf_evaluator
/// as the
///              leaf evaluator. The evaluation result of the elements of
///              \p nodes will be summed up.
///
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
///               The results of each element from \p nodes will be permuted
///               to this layout before being summed.
///
/// \param leaf_evaluator The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& nodes,  //
                   auto const& layout,  //
                   F const& leaf_evaluator, CacheManager<N, FHC>& cache) {
  // Whole-scope dispatch (Task 9): if the cache carries a whole-scope driver
  // (installed by a caller that requested batch:whole_scope), route the WHOLE
  // forest through it instead of the per-tree descent below. The driver
  // type-erases the scope_executor.hpp overload -- which lives downstream of
  // this header and cannot be included here without a cycle. Guarded by
  // `if constexpr` so it is a no-op (and byte-identical) for every
  // instantiation except the concrete residual-forest shape the driver is
  // built for: a `std::vector<N>` forest with a `std::string` layout. With no
  // driver set (the default), this is a single null-function test and the
  // standard scheme runs unchanged.
  if constexpr (std::is_same_v<std::remove_cvref_t<Nodes>, std::vector<N>> &&
                std::is_same_v<std::remove_cvref_t<decltype(layout)>,
                               std::string>) {
    if (auto const& drv = cache.whole_scope_driver(); drv)
      return drv(nodes, layout, cache);
  }

  ResultPtr result;

  for (auto&& n : nodes) {
    if (!result) {
      result = evaluate<EvalTrace>(n, layout, leaf_evaluator, cache);
      continue;
    }

    ResultPtr pre = evaluate<EvalTrace>(n, layout, leaf_evaluator, cache);
    auto time =
        detail::timed_eval_inplace([&]() { result->add_inplace(*pre); });

    // logging
    if constexpr (detail::trace(EvalTrace)) {
      // SumInplace allocates nothing: it writes into the accumulator. hwmark
      // counts the cache plus both operands live at this moment; skip pre's
      // bytes only when pre aliases a chain-resident cache buffer (fetched
      // full). pre comes back from the permute-wrapping evaluate, so a
      // permuted/phase-shifted read is a distinct buffer with its own pointer
      // and is added; chain_holds() decides by pointer identity.
      size_t hwmark = log::bytes(cache, result).value;
      if (!cache.chain_holds(pre)) hwmark += log::bytes(pre).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{
          .mode = log::EvalMode::SumInplace,
          .time = time,
          .mem_result = log::bytes(result),
          .mem_alloc = {0},
          .mem_hwmark = {cache.note_working_set(hwmark, n->hash_value())}};
      log::eval(stat, n->label(),
                log::slice_home_annot(n, cache.batch_context()));
    }
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param nodes A range of node that can be evaluated using \p leaf_evaluator
/// as the
///              leaf evaluator. The evaluation result of the elements of
///              \p nodes will be summed up.
///
/// \param leaf_evaluator The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
/// \note Because this function takes no layout argument, it is only useful
///       to evaluate summations of the elements in the \p nodes when they
///       are scalar results.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& nodes,  //
                   F const& leaf_evaluator, CacheManager<N, FHC>& cache) {
  using annot_type = decltype([](std::ranges::range_value_t<Nodes> const& n) {
    return n->annot();
  });

  static_assert(std::is_default_constructible_v<annot_type>);
  return evaluate(nodes, annot_type{}, leaf_evaluator, cache);
}

///
/// \brief Task 4 (multi-root single-DAG eval): the multi-root dispatch entry.
/// Mirrors the forest-range `evaluate(Nodes const&, layout, leaf_evaluator,
/// cache)` overload's whole-scope dispatch above, but for the cache's
/// `multiroot_driver()` (cache_manager.hpp) instead of its
/// `whole_scope_driver()`: consults the driver and routes the WHOLE set of
/// independent \p roots through it, returning one result PER ROOT with NO
/// cross-root summation (a map, not a sum) -- unlike the whole-scope driver,
/// which SUMS its forest into one `ResultPtr` (so a "not installed" caller
/// can safely fall back to the ordinary per-tree descent), a multi-root
/// caller has no equally-cheap per-tree fallback that would still give the
/// cross-root CSE this entry point exists for (see `multiroot_driver_type`'s
/// own doc comment) -- so with no driver installed this throws rather than
/// silently degrading to N independent (non-CSE'd) evaluations. The intended
/// driver (installed via `cache.set_multiroot_driver(...)`) is the ordered
/// executor's `evaluate_ordered_multiroot` (ordered_executor.hpp), closed
/// over the caller's own `OrderedSchedule`/`RichSchedule`/`target`/batching
/// policy exactly as `whole_scope_driver_type`'s doc comment describes for
/// the summed driver.
///
/// \param roots The independent root trees to evaluate; a subexpression
///        shared across two or more of them is built exactly once by a
///        driver that concatenates them into one schedule (the ordered
///        executor's own contract -- see `evaluate_ordered_multiroot`).
/// \param layouts The layout each root's own result is permuted to, one
///        entry per element of \p roots in the same order; same meaning as
///        the forest-range `evaluate`'s \p layout, applied per-root (not
///        once across a cross-root sum) -- this is what lets heterogeneous
///        roots (e.g. distinct CC residual annotations) each land in their
///        own layout.
/// \param leaf_evaluator Unused when a driver is installed (its own leaf
///        evaluator lives inside the driver's closure, exactly as
///        `whole_scope_driver_type` documents) -- present only so this
///        entry's own call signature matches the rest of the `evaluate`
///        family; kept for interface symmetry, not consulted directly here
///        (there being no per-tree fallback path to consult it on).
/// \param cache The cache whose `multiroot_driver()` is consulted.
/// \return One `ResultPtr` per element of \p roots, in \p roots' own order.
/// \throws std::logic_error if `cache.multiroot_driver()` is unset -- there
///         is no per-root fallback; a multi-root caller must explicitly
///         install a driver that understands the cross-root CSE contract.
/// \throws std::logic_error if `layouts.size() != roots.size()`.
///
template <Trace EvalTrace = Trace::Default, typename node_t, typename F,
          typename N, bool FHC>
  requires meta::leaf_node_evaluator<node_t, F>
container::svector<ResultPtr> evaluate_multiroot(
    container::svector<node_t> const& roots,
    container::svector<std::string> const& layouts,
    [[maybe_unused]] F const& leaf_evaluator, CacheManager<N, FHC>& cache) {
  static_assert(
      std::is_same_v<node_t, N>,
      "evaluate_multiroot: the roots' node type must match the cache's key "
      "type");
  if (layouts.size() != roots.size())
    throw std::logic_error(
        "evaluate_multiroot: layouts.size() must equal roots.size() -- one "
        "layout per root is required");
  auto const& drv = cache.multiroot_driver();
  if (!drv)
    throw std::logic_error(
        "evaluate_multiroot: no multiroot driver installed on the cache -- "
        "there is no per-root fallback; install one via "
        "cache.set_multiroot_driver(...) (e.g. ordered_executor.hpp's "
        "evaluate_ordered_multiroot, closed over an OrderedSchedule built "
        "from the SAME roots)");
  return drv(roots, layouts, cache);
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Evaluate given node (or a range of nodes) using an empty cache
///        manager. Calls the other sequant::evaluate function overloads.
/// \see evaluate.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
  requires(!detail::any_type_is_cache_manager<Args...>)
ResultPtr evaluate(Args&&... args) {
  using Node = std::remove_cvref_t<decltype(detail::node0(
      detail::arg0(std::forward<Args>(args)...)))>;
  auto cache = CacheManager<Node>::empty();
  return evaluate<EvalTrace>(std::forward<Args>(args)..., cache);
}

/// Empty-cache overload of the recursive engine (see \c evaluate_impl): builds
/// on a fresh CacheManager. Used by the batched evaluator's hoisted-invariant
/// build so that re-entry stays spelled evaluate_impl rather than the top-level
/// evaluate.
template <Trace EvalTrace = Trace::Default, typename... Args>
  requires(!detail::any_type_is_cache_manager<Args...>)
ResultPtr evaluate_impl(Args&&... args) {
  using Node = std::remove_cvref_t<decltype(detail::node0(
      detail::arg0(std::forward<Args>(args)...)))>;
  auto cache = CacheManager<Node>::empty();
  return evaluate_impl<EvalTrace>(std::forward<Args>(args)..., cache);
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Calls sequant::evaluate followed by the particle-symmetrization
///        function.
///        The number of particles is inferred by the tensor present in the
///        evaluation node(s). Presence of odd-ranked tensors in the evaluation
///        node(s) is an error.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
ResultPtr evaluate_symm(Args&&... args) {
  ResultPtr pre = evaluate<EvalTrace>(std::forward<Args>(args)...);
  SEQUANT_ASSERT(pre);
  ResultPtr result;
  auto time = detail::timed_eval_inplace([&]() { result = pre->symmetrize(); });

  // logging
  if constexpr (detail::trace(EvalTrace)) {
    // cache is owned by the inner evaluate call and out of scope here;
    // hwmark reflects only the local working set (pre + freshly allocated
    // result both live during the symmetrize op).
    auto stat = log::EvalStat{.mode = log::EvalMode::Symmetrize,
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = log::bytes(pre, result)};
    log::eval(
        stat,
        detail::node0(detail::arg0(std::forward<Args>(args)...))->label());
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Calls sequant::evaluate followed by the anti-symmetrization function
/// on
///        the bra indices and the ket indices. The bra and ket indices are
///        inferred from the evaluation node(s).
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
ResultPtr evaluate_antisymm(Args&&... args) {
  ResultPtr pre = evaluate<EvalTrace>(std::forward<Args>(args)...);
  SEQUANT_ASSERT(pre);

  auto const& n0 = detail::node0(detail::arg0(std::forward<Args>(args)...));

  ResultPtr result;
  auto time = detail::timed_eval_inplace(
      [&]() { result = pre->antisymmetrize(n0->as_tensor().bra_rank()); });

  // logging
  if constexpr (detail::trace(EvalTrace)) {
    // See Symmetrize for the rationale on hwmark.
    auto stat = log::EvalStat{.mode = log::EvalMode::Antisymmetrize,
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = log::bytes(pre, result)};
    log::eval(stat, n0->label());
  }
  return result;
}

/// \brief Builds a custom evaluator (see CacheManager::custom_evaluator_type)
/// that evaluates a subtree in batches over a contracted index, to bound the
/// peak memory of intermediates that carry that index.
///
/// For each node it is consulted on, the returned evaluator chooses a batch
/// mode \c K from the modes the optimizer annotated at that node
/// (\c EvalExpr::node_slice_mask; see the \c pick_sliceable lambda); it
/// declines if the node carries no accepted annotation. Annotations are
/// authoritative -- there is no heuristic fallback, so a node the
/// peak-constrained optimizer did not ask to batch is never batched. It asks
/// the backend to partition \c K into contiguous, whole-tile element-range
/// batches of at most
/// \p target_batch_size(K) elements each -- the target is an upper bound, not a
/// goal (Result::mode_batches). Mode selection is sliceability-aware: it takes
/// the first accepted mode that actually partitions into more than one batch in
/// the current (possibly already-outer-sliced) context, so a mode already
/// sliced by an outer re-entry is skipped and the node advances to its next
/// annotated mode; if no candidate is sliceable it declines (leaving small /
/// unselected / already-fully-sliced indices to the standard scheme).
///
/// A node the optimizer priced with *several* annotated modes is sliced on each
/// of them, and annotated modes may also sit at *different* nodes of one tree;
/// either way the batching nests. The per-batch scratch cache carries a
/// reinstalled copy of this evaluator, so when the standard-scheme replay of an
/// outer batch reaches an annotated node -- an inner one, or the SAME node with
/// a still-unsliced mode -- the evaluator fires again and slices the next mode
/// WITHIN the outer batch -- `for outer-batch: for inner-batch: replay`.
/// The reinstalled evaluator closes over the outer-sliced leaf evaluator, so
/// inner slicing composes on top of the outer slice; nesting is exact by the
/// same `sum_K = sum_{batches} sum_{K in batch}` identity applied per mode, and
/// a captured depth counter backstops runaway re-entry.
///
/// Otherwise it *replays the build of every compatible persistent final* in
/// the same batch passes: the group is the trigger node plus every key of
/// \p cache that is registered persistent, not yet alive, and batches over an
/// mode with the identical realized partition. Per batch, each group member is
/// evaluated by the standard scheme -- with every leaf carrying the member's
/// batch mode sliced to the batch's element range -- on a shared *registered*
/// scratch cache (see detail::make_batched_scratch), so sub-intermediates
/// repeated within a member (canonically-equal siblings) or shared between
/// members are evaluated once per batch, exactly as the real cache would share
/// them; the per-member partials are summed across batches. This is exact
/// because `sum_K = sum_{batches} sum_{K in batch}`, and never materializes
/// the whole batch-mode extent of any intermediate at once. Completed members
/// are stored into \p cache (canonical-phase convention); the trigger's result
/// is returned for evaluate() to cache as usual. Members nested inside other
/// members evaluate in earlier passes and are then seeded (slice-free w.r.t.
/// the outer batch mode) or re-derived sliced in the outer pass. Considering a
/// group candidate costs one leaf evaluation (the mode_batches probe); with an
/// unregistered (empty) cache the group is just the trigger.
///
/// Why a *group* of trees rather than the trigger alone: sub-intermediates are
/// shared between separately-intercepted finals, and a scratch scoped to one
/// final cannot see the other consumers. Concretely, in DF-based PNO-CCSD the
/// half-transformed DF factor gC = g.C (g the 3-index DF factor carrying the
/// aux index K, C the PNO coefficients) feeds both canonically-equal gCC
/// children of the particle-particle-ladder intermediate W = gCC.gCC *and* the
/// triply-transformed final gCCC. Unbatched, the real cache builds gC once and
/// serves all three uses (its keys are canonical, max_life = 3). Batching each
/// final in isolation rebuilds gC n_batches times *per final* -- the shared
/// scratch of a single pass dedups W's two gCC children within each batch, but
/// cross-final sharing with gCCC is restored only by streaming both finals
/// over the same batch partition in the same passes, which brings gC back to
/// one evaluation per batch (work parity with the unbatched path, at sliced
/// rather than full intermediate peak memory).
///
/// \param leaf_evaluator the leaf evaluator (captured).
/// \param target_batch_size per-index function
///        `std::function<std::size_t(Index const&)>` returning the per-batch
///        slice size (in elements) for a given (batch-mode) index -- an UPPER
///        BOUND, not a goal. Backend-neutral: a tiled backend rounds batch
///        boundaries to tile boundaries (down to a tile multiple), so realized
///        batches are uneven and each covers at most this many elements, except
///        the one-tile floor (a lone tile larger than the target).
/// \param accept predicate selecting which contracted indices may be batched
///        (e.g. only those in the auxiliary/RI IndexSpace). Defaults to any.
/// \param make_scope_guard factory, called with the batch count, returning an
///        RAII object held for the duration of the batched partial
///        contractions; a backend may use it to relax block-sparse screening
///        (scaled by the batch count) so per-batch screening does not drop
///        small contributions that are significant once summed over the full
///        batch mode. Defaults to a no-op (make_no_scope_guard).
///        NESTED levels: when annotated modes sit at different nodes of one
///        tree (see the nesting note above the class), the re-entrant inner
///        evaluator is built with this SAME factory (unchanged, along with
///        \p accept, \p is_volatile, \p persistent_only and \p
///        target_batch_size -- see the reinstall below), so the inner level
///        constructs its own guard from `make_scope_guard(inner_batches)`.
///        The outer guard is held for the outer level's entire batch loop,
///        which includes the per-batch evaluate() calls that re-enter and
///        construct the inner guard -- so both guards are alive
///        simultaneously while the innermost contractions run. A backend
///        factory that relaxes screening scaled by ITS OWN level's batch
///        count therefore composes multiplicatively across nesting depth:
///        net relaxation = product of batch counts over all nesting levels
///        (outer_batches * inner_batches * ...), matching the invariant that
///        a contribution significant over the full product of batch modes
///        must not be screened away in any individual (outer-cell,
///        inner-cell, ...) combination. No extra bookkeeping is needed to
///        achieve this -- it falls out of RAII scoping plus unchanged
///        threading of the factory through re-entry; see
///        "nested scope guards compose multiplicatively" in test_eval_ta.cpp
///        for a structural proof (dense TensorD has no real screening to
///        relax, so it spies on guard construction/destruction instead).
/// \param is_volatile predicate flagging a volatile leaf node (e.g. an
///        amplitude tensor); the evaluator declines to batch any node whose
///        subtree contains such a leaf, so only persistent (build-once)
///        subtrees are streamed. Defaults to never_volatile (no persistence
///        gate). Same classification as the eval cache's volatility predicate.
///        Kept last so the prior 4-argument form (…, accept, make_scope_guard)
///        still compiles unchanged.
struct accept_any_index {
  bool operator()(Index const&) const noexcept { return true; }
};

/// Default scope-guard factory for make_batched_custom_evaluator: produces a
/// no-op guard. A backend may supply a factory whose returned RAII object
/// relaxes block-sparse screening for the duration of the batched partial
/// contractions, so that a result block whose norm clears the screening
/// threshold over the *full* batch mode is not dropped in every individual
/// batch (which would lose its contribution to the sum). The factory is called
/// with the batch count, so the backend can scale the relaxation accordingly
/// (e.g. divide a Cauchy-Schwarz norm-product screening threshold by n_batches:
/// the bound for a sub-sum over 1/n of the batch mode is ~1/n of the full
/// bound). See make_batched_custom_evaluator's \p make_scope_guard parameter.
///
/// Nesting: the factory is reused unchanged at every re-entrant (nested) level
/// (see \p make_scope_guard's doc on make_batched_custom_evaluator), so nested
/// levels each construct their own guard, scaled by their own level's batch
/// count, and the outer guard(s) remain alive while the inner one is
/// constructed and used. A per-level relaxation therefore composes
/// MULTIPLICATIVELY: net relaxation = product of batch counts over all alive
/// nesting levels, which is exactly the factor needed for a contribution
/// significant over the full product of batch modes to survive per-cell
/// screening at every nesting depth.
struct no_scope_guard {};
struct make_no_scope_guard {
  no_scope_guard operator()(std::size_t /*n_batches*/) const noexcept {
    return {};
  }
};

/// Default node-volatility predicate for make_batched_custom_evaluator: no node
/// is volatile, so batching is gated only by the index predicate. A caller may
/// instead supply a predicate flagging volatile (e.g. amplitude-dependent) leaf
/// nodes; the evaluator then declines to batch any node whose subtree contains
/// such a leaf, so only persistent (build-once) subtrees are streamed over the
/// batch mode (a volatile subtree is rebuilt every evaluation, so batching it
/// would pay the partition + relaxed-screening cost on every pass for no
/// lasting memory benefit).
struct never_volatile {
  template <typename Node>
  bool operator()(Node const&) const noexcept {
    return false;
  }
};

/// \return whether any node in the subtree rooted at \p n satisfies \p pred.
template <typename Node, typename Pred>
[[nodiscard]] bool subtree_any(Node const& n, Pred const& pred) {
  if (pred(n)) return true;
  if (n.leaf()) return false;
  return subtree_any(n.left(), pred) || subtree_any(n.right(), pred);
}

namespace detail {

/// The scratch cache for one batched replay pass, plus the alive persistent
/// real-cache entries to pre-seed it with. In the DEFAULT (seeding) mode those
/// seeds are registered persistent in the scratch so they survive the per-batch
/// reset(); the caller copies their values in before the batch loop. In the
/// \c read_from_home mode (ordered executor) there is no seed set -- a
/// batch-invariant home-resident subnode is read straight from the parent chain
/// each batch instead (see make_batched_scratch's \p read_from_home).
template <typename TreeNode, bool FHC>
struct BatchedScratch {
  CacheManager<TreeNode, FHC> cache;
  std::vector<TreeNode const*> seeds;
};

/// \brief Builds the scratch CacheManager for one batched replay pass over
/// \p members (each a subtree root paired with its batch mode).
///
/// Walks every member subtree with the same pruned counting walk as
/// cache_manager() (descend on first visit of a canonical-equal node, so
/// counts match access counts under caching -- and also on a re-encounter
/// whose slicing signature differs from the first visit's, so that
/// descendants' signatures under an inconsistently-sliced occurrence are
/// recorded rather than hidden by the prune) and registers every internal
/// subnode that repeats AND has a consistent slicing signature -- the position
/// of the containing member's batch mode in the subnode's canon_indices(), or
/// its absence -- across all occurrences. Signature consistency is what makes
/// a scratch hit exact across members: canonical equality maps the index at
/// canonical position p to the index at position p, so equal signatures plus
/// equal realized element ranges (guaranteed by the caller's grouping) imply
/// identical slices. Inconsistently-sliced subnodes are not registered and so
/// are evaluated per occurrence, unshared. Count inexactness arising from the
/// pruned walk is benign: an undercount makes evaluate() recompute a drained
/// entry, an overcount keeps an entry until the per-batch reset().
///
/// Subnodes whose signature is consistently 'absent' (no leaf below carries
/// the mode -- the mode is contracted at the member's root, so a subtree
/// containing a mode-carrying leaf carries the mode free in its
/// canon_indices()) have batch-invariant full values; in the DEFAULT mode those
/// that are alive persistent entries of \p real are returned as seeds and the
/// caller copies their values into the scratch before the batch loop.
///
/// \param read_from_home ORDERED-executor discipline (default off --
/// whole-scope and forest-descent keep the seeding behavior verbatim, so MPQC's
/// current batched CC is byte-unchanged). When ON: (1) a batch-invariant
/// subnode ALREADY RESIDENT anywhere up \p real's chain is NEITHER registered
/// here NOR seeded -- \c evaluate_impl reads it straight from the parent chain
/// each batch (access_at falls through an empty local entry), the single
/// read-from-home access discipline; and (2) member ROOTS are registered too
/// (not just their subnodes), so a member consuming another member reads it
/// from the scratch rather than re-evaluating it. Together these remove every
/// non-cached node class, making a homed value's home read count exactly its
/// direct-DAG-parent count over the ordered scopes (see ordered_schedule.hpp
/// ordered_home_reads).
template <typename TreeNode, bool FHC, typename Members>
[[nodiscard]] BatchedScratch<TreeNode, FHC> make_batched_scratch(
    Members const& members, CacheManager<TreeNode, FHC> const& real,
    bool read_from_home = false,
    typename CacheManager<TreeNode, FHC>::ValueColoringCtx const* coloring =
        nullptr) {
  using Hasher = TreeNodeHasher<TreeNode, FHC>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;

  // The batch's EXTERNAL modes: obtained exactly as the evaluator obtains them
  // (partition each member root's node_slice_mask() by BatchModeType). An
  // External mode is an external that survives FREE onto a node's result, so a
  // node carrying one is NOT batch-invariant under that mode -- its value
  // depends on the external slice. When the caller nests an External mode
  // OUTSIDE a Contracted one (the External mode is sliced by an outer re-entry,
  // then this scratch batches an inner Contracted mode), a persistent
  // intermediate that carries the External mode but not the Contracted `mode`
  // would look seedable under `mode` alone -- yet seeding its full
  // (unsliced-external) value would be wrong under the outer slice. Tracking
  // the External modes in the signature (below) forbids seeding/sharing such
  // nodes. When there is no External mode this list is empty and every
  // External-derived test is a no-op, keeping the Contracted-only behavior
  // byte-identical.
  container::svector<Index> ext_axes;
  for (auto const& [root, mode] : members) {
    if (root->leaf()) continue;
    for (auto const& [ix, knd] : (*root)->node_slice_mask())
      if (knd == BatchModeType::External &&
          std::find(ext_axes.begin(), ext_axes.end(), ix) == ext_axes.end())
        ext_axes.push_back(ix);
  }

  struct Meta {
    std::size_t count = 0;
    std::optional<std::size_t> sig;
    // Positions of each External mode (in ext_axes order) in this node's
    // canonical result indices, or nullopt if the node does not carry it. This
    // is a function of the (canonical) node alone, so it is identical across
    // all occurrences of a canonically-equal node.
    container::svector<std::optional<std::size_t>> ext_sig;
    bool consistent = true;
  };
  std::unordered_map<TreeNode const*, Meta, Hasher, Comp> meta;

  // The external-mode signature (positions of each ext axis on n's result, or
  // absent) is exactly slicing_signature(n, ext_axes); the batch-mode `sig`
  // below is the single-mode case, index_position(n, mode).
  auto ext_sig_of = [&ext_axes](TreeNode const& n) {
    return slicing_signature(n, ext_axes);
  };

  auto visit = [&meta, &ext_sig_of](auto&& self, TreeNode const& n,
                                    Index const& mode) -> void {
    if (n.leaf()) return;
    auto const sig = index_position(n, mode);
    auto const esig = ext_sig_of(n);
    auto const [it, first] = meta.try_emplace(&n);
    auto& e = it->second;
    if (first) {
      e.sig = sig;
      e.ext_sig = esig;
    } else if (e.sig != sig || e.ext_sig != esig) {
      e.consistent = false;
    }
    ++e.count;
    // Prune a re-encounter only when its signature matches the first one:
    // canonical equality maps canonical position p to position p, so an equal
    // signature here implies the descendants' signatures equal those already
    // recorded on the first walk (deeper accesses shared and counted). A
    // differing signature gives no such guarantee -- descend so descendants'
    // signatures under this occurrence are recorded too; otherwise a
    // descendant sliced differently only under this (unshared, pruned)
    // occurrence could pass the guard and serve wrong slices. The extra
    // descendant counts are real accesses: an inconsistently-sliced occurrence
    // is evaluated per occurrence, not served from the scratch at n. The
    // External signature is invariant across occurrences (a function of the
    // canonical node), so folding it into the match only tightens the guard.
    if (!first && e.sig == sig && e.ext_sig == esig) return;
    self(self, n.left(), mode);
    self(self, n.right(), mode);
  };
  for (auto const& [root, mode] : members) {
    if (root->leaf()) continue;
    if (read_from_home) {
      // ordered: register the member ROOT too, so a member consuming another
      // member reads it from the scratch instead of re-evaluating it.
      visit(visit, *root, mode);
    } else {
      // default: member roots are accumulated by the caller, not cached here.
      visit(visit, root->left(), mode);
      visit(visit, root->right(), mode);
    }
  }

  std::unordered_map<TreeNode, std::size_t, Hasher, Comp> reg;
  std::unordered_set<TreeNode, Hasher, Comp> seed_keys;
  std::vector<TreeNode const*> seeds;
  for (auto const& [ptr, e] : meta) {
    if (!e.consistent) continue;  // ambiguous slicing: never share
    // A node carrying ANY batched External mode has an external slice a
    // seeded/home-read full value would ignore -- so it is never
    // shareable-full.
    bool const carries_ext =
        std::any_of(e.ext_sig.begin(), e.ext_sig.end(),
                    [](auto const& p) { return p.has_value(); });
    if (read_from_home) {
      // ORDERED: a batch-invariant subnode already resident anywhere up the
      // chain is read straight from home each batch -- neither registered here
      // (a registered-but-reset local entry would rebuild every batch) nor
      // copied in. Single access discipline, no seeds. See the \p
      // read_from_home doc above.
      if (!e.sig && !carries_ext && real.resident_in_chain(*ptr)) continue;
      if (e.count >= 2) reg.emplace(*ptr, e.count);
    } else {
      // DEFAULT (whole-scope / forest-descent): seed an alive PERSISTENT
      // batch-invariant real entry into the scratch (persistent so it survives
      // reset()), else register a repeated subnode. Byte-identical to the
      // pre-read-from-home behavior MPQC's batched CC relies on.
      bool const seedable =
          !e.sig && !carries_ext && real.persistent(*ptr) && real.alive(*ptr);
      if (seedable) {
        seeds.push_back(ptr);
        seed_keys.insert(*ptr);
        reg.emplace(*ptr, e.count);  // count ignored for persistent entries
      } else if (e.count >= 2) {
        reg.emplace(*ptr, e.count);
      }
    }
  }
  auto is_persistent = [seed_keys = std::move(seed_keys)](TreeNode const& n) {
    return seed_keys.contains(n);
  };
  CacheManager<TreeNode, FHC> scratch{std::move(reg), std::move(is_persistent)};
  // Pillar 1 / Task 7: this sub-top scratch is VALUE-keyed. `reg` above
  // registered its members by bare node; re-key them through the per-scope
  // coloring context so each becomes its home-slice-colored value-id --
  // matching the colored store/access. Set the context so runtime store/access
  // recolor through it too. Null context => unchanged (byte-identical top-level
  // path).
  if (coloring) {
    scratch.set_value_coloring_ctx(coloring);
    scratch.recolor_registered_entries();
  }
  // Read-from-home scratches statically pre-schedule every value and read
  // batch-invariant operands from home each batch (no seeding); a miss on such
  // an operand is a real defect (premature eviction / under-predicted use
  // count), so require it to surface as a hard error rather than an empty-array
  // hang. Seeded scratches (read_from_home=false) keep miss=>compute.
  if (read_from_home) scratch.set_require_resident_reads(true);
  return {std::move(scratch), std::move(seeds)};
}

}  // namespace detail

/// Opt-in sink for the batched evaluator's per-batch scratch high-watermarks.
/// The batched replay runs each aux/mu~ batch in a SEPARATE scratch
/// CacheManager whose transients never enter the outer cache, so the outer
/// cache's working_set_hwmark() misses the batched-inner peak. When a non-null
/// PeakSink is threaded through make_batched_custom_evaluator (and its nested
/// re-instantiations), each scratch's high-watermark folds (max) into this one
/// global accumulator, yielding the true batched-replay peak. A null sink
/// (the default) leaves all existing behavior byte-identical.
using PeakSink = std::atomic<double>*;

template <typename F, typename IndexPredicate = accept_any_index,
          typename ScopeGuardFactory = make_no_scope_guard,
          typename IsVolatile = never_volatile>
[[nodiscard]] auto make_batched_custom_evaluator(
    F leaf_evaluator,
    std::function<std::size_t(Index const&)> target_batch_size,
    IndexPredicate accept = {}, ScopeGuardFactory make_scope_guard = {},
    IsVolatile is_volatile = {}, bool persistent_only = false,
    std::size_t depth = 0, PeakSink peak = nullptr) {
  return [leaf_evaluator = std::move(leaf_evaluator),
          target_batch_size = std::move(target_batch_size), accept, is_volatile,
          persistent_only, depth, peak,
          make_scope_guard](auto const& node, auto& cache) -> ResultPtr {
    // Runaway backstop: nesting re-enters this evaluator on the per-batch
    // scratch (see the reinstall below), incrementing depth once per nested
    // mode. Real trees nest a handful of modes deep; a large depth signals a
    // non-terminating re-entry (e.g. a mis-annotated mode that never shrinks).
    SEQUANT_ASSERT(depth < 8);
    // Slice a value to the current batch block for the `hops` innermost
    // enclosing batch loops that `nd` carries (see the Enter-stage
    // slice_to_use in evaluate() -- this is the same primitive, reachable from
    // the closure-internal probes below that bypass the Enter stage). `cache`
    // is the cache the closure fired on, so cache.batch_context() is THIS
    // node's ENCLOSING context.
    auto slice_to_use = [&cache](ResultPtr value, auto const& nd,
                                 std::size_t hops) -> ResultPtr {
      auto const& ctx = cache.batch_context();
      std::size_t const d = ctx.size();
      // See the assert on the evaluate() copy of this lambda: hops <= d always
      // (one batch_context push + <=1 parent link per level); a violation would
      // underflow `d - hops` and silently UNDER-slice.
      SEQUANT_ASSERT(hops <= d);
      for (std::size_t i = d - hops; i < d; ++i) {
        auto const& axis = ctx[i].axis;
        auto const& blk = ctx[i].range;
        if (auto const p = index_position(nd, axis))
          value = value->slice_mode(*p, blk.first, blk.second);
      }
      return value;
    };
    // A leaf evaluator that slices each fetched leaf to the enclosing batch
    // blocks: the REPLACEMENT for the old per-block `le_g`, generalized from
    // one mode to the whole enclosing nest. Used at the three closure-internal
    // sites (pick_sliceable probe, carrier_full pre-size, hoist build) that
    // consume leaf slicing but bypass the Enter stage; feeding the RAW
    // leaf_evaluator there would un-slice the enclosing loops and break the
    // "K is not re-picked" invariant (a re-entered probe would see K's FULL
    // extent and re-batch it). The MAIN value path does NOT use this: it re-
    // enters evaluate() with the raw leaf_evaluator and lets the Enter-stage
    // slice_to_use slice, which reproduces le_g exactly on that path.
    auto sliced_leaf = [&](auto const& ln) -> ResultPtr {
      return slice_to_use(leaf_evaluator(ln), ln, cache.batch_context().size());
    };
    // Synthesized DAG-scope level for a forest-evaluator push: this firing
    // realizes exactly ONE loop over `cache`'s own enclosing context, so every
    // push site below shares the same depth (`cache.batch_context().size() +
    // 1`, matching build_ordered_schedule's `d + 1` convention -- see
    // DagScopeLevel's doc comment) and differs only in the pushed axis's
    // space. Plumbing only (Task 6): nothing resolves by `level` yet on this
    // path -- the forest evaluator's OLD resolution is exact-axis
    // (`exact_axis`, filled at each push site below), which Task 7 will
    // consult instead.
    auto const synth_level = [&cache](Index const& ax) -> DagScopeLevel {
      return DagScopeLevel{.depth = cache.batch_context().size() + 1,
                           .space = std::wstring(ax.space().base_key())};
    };
    // Mode selection is SLICEABILITY-AWARE and realizes the optimizer's
    // multi-mode nesting one mode per depth level. candidate_axes lists this
    // node's batch modes in the optimizer's annotated order (see
    // EvalExpr::node_slice_mask), keeping the accepted annotations.
    //
    // The optimizer's annotations are AUTHORITATIVE at every depth: a node
    // carrying no accepted annotation means "do not batch this node", and is
    // left unbatched. There is deliberately NO heuristic fallback -- batching
    // is only ever realized where the peak-constrained optimizer asked for it,
    // so every realized batch loop is one the cost model priced. (Callers
    // therefore cannot batch without a peak budget: no budget => the optimizer
    // emits no annotations => nothing batches. See BatchPolicy::peak_threshold
    // and, on the MPQC side, validate_batch_config.)
    auto candidate_axes =
        [&accept](auto const& n) -> container::svector<Index> {
      container::svector<Index> out;
      for (auto const& entry : n->node_slice_mask())
        if (accept(entry.first)) out.push_back(entry.first);
      return out;
    };
    // Pick the FIRST candidate mode that is actually sliceable (partitions into
    // > 1 batch) in THIS (possibly already-outer-sliced) context, returning it
    // together with its realized partition. A mode already sliced by an outer
    // re-entry yields a single batch on the sliced leaf and is skipped, so a
    // nested re-entry on the SAME node advances to the node's next annotated
    // mode -- realizing `for K-batch: for mu1-batch: replay` at one multi-mode
    // node. The recursive reinstall below walks one mode per depth level (the
    // depth < 8 backstop bounds the re-entry).
    auto pick_sliceable = [&](auto const& n)
        -> std::optional<std::pair<
            Index, container::svector<std::pair<std::size_t, std::size_t>>>> {
      BackendArrayOps const* const aops = cache.array_ops();
      auto const& ectx = cache.batch_context();
      // Is ix's mode already sliced by an enclosing block, so it must not be
      // re-picked? Mirrors slice_to_use exactly: an EXACT context axis slices
      // ix. This replaces the old "slice the carrier leaf and see if it yields
      // a single batch" probe, which needed mode_batches.
      auto already_sliced = [&](Index const& ix) {
        for (auto const& e : ectx)
          if (e.axis == ix) return true;
        return false;
      };
      for (Index const& ix : candidate_axes(n)) {
        if (already_sliced(ix)) continue;
        SEQUANT_ASSERT(aops &&
                       "batched forest eval requires backend array-ops "
                       "(CacheManager::set_array_ops)");
        auto b = aops->axis_batches(ix, target_batch_size(ix));
        if (b.size() > 1) return std::make_pair(ix, std::move(b));
      }
      return std::nullopt;
    };

    // Persistence gate (opt-in via persistent_only): when set, decline to batch
    // any subtree containing a volatile leaf -- such a subtree is rebuilt every
    // evaluation, so batching pays the partition + relaxed-screening cost each
    // pass to amortize over nothing. By default (persistent_only == false) we
    // batch ACROSS THE BOARD: slicing the batch mode reduces the footprint of
    // any mode-carrying intermediate regardless of volatility, and the cost
    // model credits it accordingly, so the runtime must realize it too. (When
    // is_volatile is never_volatile the gate is moot either way.)
    if (persistent_only && subtree_any(node, is_volatile)) {
      return nullptr;
    }

    auto picked = pick_sliceable(node);
    if (!picked) {
      return nullptr;  // no accepted, sliceable mode (nothing to gain)
    }
    Index const K = std::move(picked->first);
    auto const batches = std::move(picked->second);

    using node_t = std::remove_cvref_t<decltype(node)>;
    using member_t = std::pair<node_t const*, Index>;
    TreeNodeEqualityComparator<node_t> const eq;

    // Classify the picked mode by BatchModeType. K is a batch mode of `node`;
    // the optimizer stamps it CONTRACTED (summed away -> block partials
    // ACCUMULATE) or EXTERNAL (an external index free on the node's result ->
    // block partials are DISJOINT slices, SCATTERED into a pre-sized result).
    // The depth-0 heuristic fallback only ever yields a contracted index, so an
    // mode absent from node_slice_mask() is Contracted -- keeping the
    // Contracted-only path (no External entry) byte-identical to before this
    // branch existed.
    BatchModeType picked_kind = BatchModeType::Contracted;
    for (auto const& [ix, knd] : node->node_slice_mask())
      if (ix == K) {
        picked_kind = knd;
        break;
      }

    // Per-level placement (order-aware only). Replaces hoist_invariants. This
    // firing realizes a batch loop over `K` at runtime `depth`. A
    // member-subtree node INVARIANT to this loop (it does not carry `K` on its
    // result) is built ONCE at its home level and served to every batch body
    // through the scope chain, rather than rebuilt per batch. A node's
    // residency (the batch modes it is variant to) is its \c sliced_modes():
    // the cross-occurrence lifetime-mask meet of ALL batched modes -- External
    // (occ) AND Contracted (aux) alike -- that live on the node's own result
    // slots (consistent placement across occurrences -> CSE). A node variant to
    // an outer aux loop carries that aux free on a result slot, so the aux mode
    // survives the meet into sliced_modes; the former per-occurrence
    // contracted_modes bolt-on is thus subsumed and gone.
    // A node is invariant to THIS loop iff `K` is NOT in its residency. Its
    // HOME level is the deepest ENCLOSING batch_context entry whose mode is in
    // the residency (the innermost enclosing loop it is variant to); -1 (the
    // chain root / run-term cache) if it is invariant to the whole nest. The
    // node is built once at that level (walk-up), sliced to its home blocks,
    // and reused across this loop's batches. A node carrying `K` (K in its
    // residency) is loop-LOCAL: it is left to inline evaluation (descend),
    // which finds any deeper-hoisted invariants through the chain -- so
    // descending never rebuilds them.
    //
    // This reproduces hoist_invariants' walk-up structure exactly, with the
    // former per-node scalar placement level replaced by the residency-derived
    // home level and the `sl < depth && !ext_loop_local` predicate replaced by
    // `K not in sliced_modes` (which subsumes ext_loop_local: an External
    // carrier has K in sliced_modes -> loop-local -> descended, never hoisted).
    // The order-aware GATE is the emitted `batch_order_aware()` bit (true for
    // every node the order-aware cost model emitted, including a whole-nest
    // invariant whose residency is empty): on the OFF path every node is
    // order-blind, so `targets` is empty, set_parent is NOT wired, and the
    // per-batch replay runs exactly as before. The bit is a positive signal an
    // empty residency cannot provide -- it is what distinguishes an OFF-path
    // all-full node (do not hoist) from an order-aware whole-nest invariant
    // (hoist to the root).
    auto place_at_this_level =
        [&](auto& scratch_cache, auto& parent_cache,
            std::vector<node_t const*> const& member_roots) {
          // The ENCLOSING batch loops (strictly OUTER to this firing); this
          // level's mode K and any INNER loop are NOT in it. A node is
          // hoistable here iff EVERY residency (sliced_modes) mode is one of
          // these outer loops: then it is invariant to this loop AND to every
          // inner loop, so its home is its deepest enclosing residency level
          // and it is built once there. If a residency mode is K (loop-local)
          // or an INNER mode (its home is a deeper loop), it is not all-outer
          // -> descend, so the deeper level handles it sliced. This is exactly
          // hoist_invariants' `scope_level < depth`: the deepest residency
          // being outer <=> ALL residency outer (deepest is the max), and a
          // node carrying an inner mode is (correctly) not hoisted at this
          // outer level -- the bug an over-eager `K not in residency` predicate
          // caused (holding an aux-carrier full over aux at the outer occ
          // level).
          auto const& ectx = parent_cache.batch_context();  // enclosing loops
          auto in_ectx = [&ectx](Index const& m) -> bool {
            for (auto const& e : ectx)
              if (e.axis == m) return true;
            return false;
          };
          auto residency_all_outer = [&in_ectx](node_t const& n) -> bool {
            for (auto const& ix : n->sliced_modes())
              if (!in_ectx(ix)) return false;
            return true;
          };
          // Phase 4b-3 runtime cutover: the demotion veto that used to force a
          // MEET-DEMOTED external carrier (a node carrying an EXTERNAL
          // node_slice_mask() stamp absent from its sliced_modes) to DESCEND is
          // gone. Placement is now purely router-override-or-seed: an
          // order-aware, residency-all-outer node is hoisted to its seed home
          // (FULL on any demoted mode) UNLESS the remat router (populated by
          // the MPQC pre-pass) overrides its home lower. The remat pass shrinks
          // such a giant and re-homes it via the router; with an empty router
          // (all SeQuant unit tests, no pre-pass) and the inf default there is
          // no batching, so the veto was inert anyway (design section 5). A
          // value cached at its seed home is the SAME value the descended path
          // produced -- the batched Enter-stage slice-on-use slices it to the
          // block when a nested external loop consumes it -- so this changes
          // PLACEMENT only, never the result.
          // Split gate (CSE-aware hoist): under a LIVE router, two canonically-
          // equal occurrences may SHARE one hoist materialization only if they
          // are slicing-signature CONSISTENT over the enclosing batch modes.
          // When inconsistent they bind different physical labels to the sliced
          // slot (the g.C legs, i_3 vs i_4) and must be built PER OCCURRENCE at
          // their own router-resolved depths -- otherwise one leg reads the
          // other's slice (see slicing_signature.hpp and the design spec
          // doc/dev/specs/2026-08-07-remat-cse-aware-split-design.md). With no
          // router the meet homes a divergent value FULL at the root, where
          // sharing is safe, so the split is gated on a live router to keep
          // every no-router path byte-identical.
          container::svector<Index> ectx_modes_sig;
          for (auto const& e : ectx) ectx_modes_sig.push_back(e.axis);
          bool const split_aware = parent_cache.placement_router() &&
                                   !parent_cache.placement_router()->empty();
          std::vector<node_t const*> targets;
          auto collect = [&](auto&& self, node_t const& n) -> void {
            if (n.leaf()) return;
            if (n->batch_order_aware() && residency_all_outer(n) &&
                !subtree_any(n, is_volatile)) {
              auto shares = [&](node_t const* p) {
                if (!eq(*p, n)) return false;
                if (!split_aware) return true;  // no router: dedup as before
                return slicing_signature(*p, ectx_modes_sig) ==
                       slicing_signature(n, ectx_modes_sig);
              };
              if (std::none_of(targets.begin(), targets.end(), shares))
                targets.push_back(&n);
              return;  // built as a unit -- do not descend into it
            }
            self(self, n.left());
            self(self, n.right());
          };
          for (node_t const* m : member_roots) {
            if (m->leaf()) continue;
            collect(collect, m->left());
            collect(collect, m->right());
          }
          if (targets.empty()) return;
          // Wire the scope chain only when there is something to hoist (matches
          // hoist_invariants; keeps the OFF path unwired -> byte-identical).
          scratch_cache.set_parent(&parent_cache);
          auto in_residency = [](node_t const& n, Index const& m) -> bool {
            auto const& sm = n->sliced_modes();
            return std::find(sm.begin(), sm.end(), m) != sm.end();
          };
          for (node_t const* dptr : targets) {
            node_t const& d = *dptr;
            // Home level = deepest enclosing-context entry whose mode is in d's
            // residency (sliced_modes); -1 => invariant to the whole nest
            // (chain root).
            int rl = -1;
            for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
              if (in_residency(d, ectx[i].axis)) {
                rl = i;
                break;
              }
            // Router override (see placement_router.hpp): a registered
            // occurrence override replaces the LEVEL only -- the collect/hoist
            // decision above (residency_all_outer alone, post Phase 4b-3) is
            // untouched. Empty/null router (Phase 2 default) => branch never
            // taken => rl is exactly the value computed above =>
            // byte-identical.
            // Gate on moved() FIRST (a hash-set lookup, no occurrence_key):
            // only a moved value can have an override here, and the build keys
            // occurrence_key only for moved nodes -- so this never computes an
            // occurrence_key for a node the router could not have keyed (incl.
            // non-tensorial Sum nodes, which are never moved and are not tensor
            // networks). Byte-identical: a non-moved d route()-missed and fell
            // through before, and it falls through (skips this block) now.
            if (auto const* router = parent_cache.placement_router();
                router && !router->empty() && router->moved(d->hash_value())) {
              container::svector<Index> ectx_modes;
              for (auto const& e : ectx) ectx_modes.push_back(e.axis);
              auto const key = eval::occurrence_key(d, ectx_modes);
              auto const* home = router->route(key);
              if (home) {
                rl = static_cast<int>(router->home_depth(*home, ectx, key)) - 1;
              } else {
                // d is moved but route() MISSED here: it is remat-demoted to a
                // DEEPER context than this hoist is inside, so its
                // per-occurrence override is keyed at that context. Do NOT
                // build it full at this (outer) home -- defer to the deeper
                // hoist where route() hits and it is built per-block.
                // occurrence_key is context- dependent, so the outer query
                // cannot see the deeper override; the context-invariant moved()
                // flag is what lets us skip. Without this, the eager outer full
                // build wins and the alive()-skip below blocks the intended
                // per-block rebuild.
                continue;
              }
            }
            // Locate the level-rl cache by walking UP from parent_cache (the
            // level depth-1 cache): rl == -1 => the chain root (the real/term
            // cache); rl >= 0 => the scratch (depth-1 - rl) hops up. Runtime
            // nest depth aligns with the scope-chain position (each realized
            // loop = one context entry = one scratch level).
            auto* target = &parent_cache;
            if (rl == -1) {
              while (target->parent()) target = target->parent();
            } else {
              // Release-safe guard (SEQUANT_ASSERT elides in release): never
              // walk the chain off its end and dereference a null parent.
              if (rl > static_cast<int>(depth) - 1)
                throw std::runtime_error(
                    "hoist home level not strictly outer to this loop");
              for (int lvl = static_cast<int>(depth) - 1; lvl > rl; --lvl) {
                auto* const p = target->parent();
                if (!p)
                  throw std::runtime_error(
                      "hoist walk-up exceeded the scope chain (a single-batch "
                      "sliced mode may have shifted the runtime nest depth)");
                target = p;
              }
            }
            target->ensure_hoist_slot(d);
            if (target->alive(d)) continue;  // built already in a broader scope
            // Build the whole invariant ONCE on a FRESH cache via the variadic
            // evaluate(n, sliced_leaf) (empty cache, NO custom evaluator, so no
            // re-entry into this batched evaluator). `sliced_leaf` slices the
            // enclosing loops d carries (up to its home level); the loops it
            // does not carry pass through unsliced (built full over its deeper
            // / invariant modes). Store under the same canonical-phase
            // convention the batched member store uses.
            ResultPtr built = evaluate_impl(d, sliced_leaf);
            if (auto const ph = d->canon_phase(); ph != 1)
              built = built->mult_by_phase(ph);
            (void)target->store(d, std::move(built));
          }
        };

    // DEBUG (behavior-neutral): log the trigger's depth, picked mode + kind,
    // and its full node_slice_mask() annotation + result indices, to diagnose
    // nested re-batching of a single aux mode. Emitted only when tracing is on.
    if (log::printing()) {
      std::string annot;
      for (auto const& [ix, knd] : node->node_slice_mask()) {
        annot += toUtf8(ix.full_label());
        annot += (knd == BatchModeType::External ? ":ext " : ":con ");
      }
      std::string res;
      for (auto const& ix : node->canon_indices()) {
        res += toUtf8(ix.full_label());
        res += " ";
      }
      auto scope_ctx = cache.batch_context();
      scope_ctx.push_back({K, synth_level(K), {0, 0}, K});
      log::log(
          "BatchAxes",
          std::format("depth={} picked={}:{} nbatches={} annot=[{}] "
                      "result=[{}] {}",
                      depth, toUtf8(K.full_label()),
                      picked_kind == BatchModeType::External ? "ext" : "con",
                      batches.size(), annot, res, log::scope_annot(scope_ctx)));
    }

    if (picked_kind == BatchModeType::External) {
      // SCATTER branch. K survives to node's result as a free external mode,
      // so the per-block partials are disjoint slices of one result (not
      // summands of a contraction): they are write_into_slice()d into a
      // pre-sized result, never add_inplace()d. Inner batch modes -- of either
      // kind, at this node or a descendant -- still nest through the same
      // per-block reinstall the contracted path uses, so External composes with
      // Contracted: within each external block the reinstalled evaluator slices
      // any remaining mode (a contracted mode accumulates within the block, an
      // inner external mode scatters into its own block-local result). K itself
      // is not re-picked on the re-entry: the pushed batch_context makes the
      // re-entry's sliced_leaf slice K's carrier to this block, so the block
      // yields a single batch on the sliced leaf and is skipped (same invariant
      // as the contracted nesting).
      auto const dest_mode = index_position(node, K);
      SEQUANT_ASSERT(dest_mode &&
                     "external batch mode is not free on the node's result");
      // Backend array-ops from the cache chain -- the SAME source the ordered
      // (DAG) executor reads, so forest and DAG build identical scatter
      // destinations from the node's own (unsliced) index list.
      BackendArrayOps const* const aops = cache.array_ops();
      SEQUANT_ASSERT(aops &&
                     "batched external-mode scatter requires backend array-ops "
                     "(CacheManager::set_array_ops)");

      // A single-node scratch: an external mode is not a
      // persistent-final sharing mode, so the group/replay machinery (which
      // co-batches cross-term contracted finals) does not apply -- scatter just
      // this node. The scratch still dedups repeats WITHIN the node's subtree.
      std::vector<member_t> solo{{&node, K}};
      auto bs = detail::make_batched_scratch(solo, cache);
      bs.cache.set_array_ops(cache.array_ops());  // inherit backend ops
      for (auto const* s : bs.seeds) (void)bs.cache.store(*s, cache.access(*s));
      place_at_this_level(bs.cache, cache, std::vector<node_t const*>{&node});

      auto const scope_guard = make_scope_guard(batches.size());
      (void)scope_guard;

      if (log::printing())
        log::log("BatchScatter", "Begin",
                 std::format("external mode over {} blocks", batches.size()));

      ResultPtr dest;
      for (auto const& [e_lo, e_hi] : batches) {
        if (e_lo == e_hi) continue;
        // Per-block slice marker for trace post-processing (avoidable-recompute
        // accounting), mirroring the contracted BatchGroup path: tags every op
        // replayed in this external block with the enclosing scatter loop's
        // (mode, element-range-low). Without it, the same expression scattered
        // into disjoint external blocks would carry an identical touched-slice
        // signature and be miscounted as a duplicate rebuild, when each block
        // is in fact legitimate (1/nblocks of the work). Gated on printing() so
        // it is inert unless a trace is being emitted.
        if (log::printing())
          log::log("BatchIter", toUtf8(std::wstring(K.full_label())), e_lo);
        bs.cache.reset();
        // Extend the enclosing batch context by THIS block and set it on the
        // scratch, so the re-entry's Enter-stage slice-on-use (and its own
        // sliced_leaf) slices every leaf carrying K to this block and composes
        // inner slices on top -- exactly what the old per-block `le_g` did on
        // the leaf path, plus it now also slices a cached intermediate fetched
        // from an ancestor scope (the slice-on-use fix). The raw leaf_evaluator
        // is threaded down; the Enter stage does the slicing.
        auto ctx = cache.batch_context();
        ctx.push_back({K, synth_level(K), {e_lo, e_hi}, K});
        bs.cache.set_batch_context(std::move(ctx));
        bs.cache.set_custom_evaluator(make_batched_custom_evaluator(
            std::function<ResultPtr(node_t const&)>{leaf_evaluator},
            target_batch_size, accept, make_scope_guard, is_volatile,
            persistent_only, depth + 1, peak));
        ResultPtr part = evaluate_impl(node, leaf_evaluator, bs.cache);
        // Pre-size the full-extent zero destination from the node's own
        // (unsliced) index list on the first block; the backend realizes it
        // (flat or nested) with no array in the DAG consulted.
        if (!dest) dest = aops->make_zeros(node->canon_indices());
        dest->write_into_slice(*part, *dest_mode, e_lo, e_hi);

        if (peak) {
          const double cand =
              static_cast<double>(bs.cache.working_set_hwmark());
          double cur = peak->load(std::memory_order_relaxed);
          while (cand > cur && !peak->compare_exchange_weak(
                                   cur, cand, std::memory_order_relaxed)) {
          }
        }
      }
      if (log::printing()) log::log("BatchScatter", "End");
      SEQUANT_ASSERT(dest);
      return dest;
    }

    // The replay group: the trigger plus every registered persistent key that
    // is not yet alive and batches over a mode with the identical realized
    // partition. All compatible persistent finals stream over the batch mode
    // in the same passes, so sub-intermediates shared between them (wherever
    // the scratch's slicing-signature guard admits sharing -- equal canonical
    // positions of the batch mode plus equal element ranges imply identical
    // slices) are evaluated once per batch instead of once per consumer.
    // The cost of considering a candidate is one leaf evaluation (the
    // mode_batches probe). With an unregistered (empty) real cache the group
    // is just the trigger.
    std::vector<member_t> group{{&node, K}};
    cache.for_each_key([&](node_t const& k) {
      if (!cache.persistent(k) || cache.alive(k)) return;
      if (eq(k, node)) return;  // the trigger occupies its own slot
      if (subtree_any(k, is_volatile)) return;  // defensive: P implies NV
      auto const pk = pick_sliceable(k);
      if (!pk) return;
      // Join iff this member's first sliceable mode realizes the identical
      // partition as the trigger (so all members stream over the same batches).
      if (pk->second != batches) return;
      group.emplace_back(&k, pk->first);
    });

    // Layer by nesting: a member whose subtree contains another member
    // evaluates in a later layer, with the inner result by then alive in the
    // real cache -- seeded into the outer pass when slice-free w.r.t. the
    // outer batch mode, re-derived sliced (correct, unshared) otherwise.
    auto contains = [&eq](node_t const& outer, node_t const& inner) -> bool {
      auto rec = [&eq, &inner](auto&& self, node_t const& n) -> bool {
        if (eq(n, inner)) return true;
        if (n.leaf()) return false;
        return self(self, n.left()) || self(self, n.right());
      };
      if (outer.leaf()) return false;
      return rec(rec, outer.left()) || rec(rec, outer.right());
    };
    std::vector<std::vector<member_t>> layers;
    {
      std::vector<member_t> remaining = std::move(group);
      while (!remaining.empty()) {
        std::vector<member_t> layer, rest;
        for (auto const& m : remaining) {
          bool const outer = std::any_of(
              remaining.begin(), remaining.end(), [&](member_t const& o) {
                return m.first != o.first && contains(*m.first, *o.first);
              });
          (outer ? rest : layer).push_back(m);
        }
        SEQUANT_ASSERT(!layer.empty());  // containment is a strict order
        layers.push_back(std::move(layer));
        remaining = std::move(rest);
      }
    }

    // Trace: the batched path co-evaluates a GROUP -- the trigger plus any
    // cross-term persistent finals that slice over the same aux partition --
    // streaming them together over the aux batches in one pass (so a
    // sub-intermediate shared between members is computed once per batch, not
    // once per consumer). The members are SIBLINGS computed alongside each
    // other, NOT a term hierarchy; the per-op Eval lines below interleave
    // across members and batches. Bracket the group and list its members so
    // those ops can be attributed. Distinct "BatchGroup"/"BatchMember" labels
    // (not Term|Begin/End) to avoid implying nesting; the top-level evaluate
    // still emits the enclosing per-term Term markers.
    if (log::printing()) {
      std::size_t n_members = 0;
      for (auto const& layer : layers) n_members += layer.size();
      auto scope_ctx = cache.batch_context();
      scope_ctx.push_back({K, synth_level(K), {0, 0}, K});
      log::log(
          "BatchGroup", "Begin",
          std::format("{} members co-evaluated over {} aux batches {}",
                      n_members, batches.size(), log::scope_annot(scope_ctx)));
      for (auto const& layer : layers)
        for (auto const& mk : layer)
          log::log("BatchMember",
                   toUtf8(io::serialization::to_string(to_expr(*mk.first))));
    }
    {
      // Structured BatchGroup for the visualizer: the co-evaluation unit (its
      // batch mode, block count, and member node hashes) so the DAG can draw a
      // subgraph enclosing the siblings streamed together over K -- the runtime
      // batching structure the IR forest cannot show. Gated by
      // SEQUANT_SCHED_DUMP.
      static bool const sched_dump =
          std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
      if (sched_dump) {
        std::cerr << "SCHEDULE_RUN_GROUP {\"kind\":\""
                  << (picked_kind == BatchModeType::External ? "external"
                                                             : "contracted")
                  << "\",\"mode\":\"" << toUtf8(K.full_label())
                  << "\",\"blocks\":" << batches.size() << ",\"members\":[";
        bool gfirst = true;
        for (auto const& layer : layers)
          for (auto const& mk : layer) {
            std::cerr << (gfirst ? "" : ",") << '"' << (*mk.first)->hash_value()
                      << '"';
            gfirst = false;
          }
        std::cerr << "]}\n";
      }
    }

    // RAII scope for the batched partial contractions; a backend-supplied
    // factory may relax block-sparse screening here (scaled by the batch count)
    // so per-batch screening does not drop contributions that survive over the
    // full batch mode. Held for the ENTIRE loop below, including the per-batch
    // evaluate() calls that may re-enter this evaluator on an inner annotated
    // node (see the reinstall's `make_scope_guard` argument): the inner
    // level's own guard is then constructed and destroyed while this (outer)
    // guard is still alive, so a backend that relaxes screening scaled by its
    // own level's batch count composes multiplicatively across nesting depth
    // (net relaxation = product of batch counts over all alive levels).
    auto const scope_guard = make_scope_guard(batches.size());
    (void)scope_guard;

    ResultPtr trigger_result;
    for (auto const& layer : layers) {
      // The layer's scratch cache: registered from the member subtrees (same
      // canonical-equality counting as the real cache), so repeated subtrees
      // -- canonically-equal siblings within a member as well as
      // sub-intermediates shared between members -- are evaluated once per
      // batch. Carries no custom evaluator (no re-interception) and keeps the
      // partial, sliced intermediates out of the real cache; reset() between
      // batches drops the previous batch's partials, while pre-seeded alive
      // persistent entries (registered persistent in the scratch) survive.
      auto bs = detail::make_batched_scratch(layer, cache);
      bs.cache.set_array_ops(cache.array_ops());  // inherit backend ops
      for (auto const* s : bs.seeds) (void)bs.cache.store(*s, cache.access(*s));
      {
        std::vector<node_t const*> roots;
        roots.reserve(layer.size());
        for (auto const& mk : layer) roots.push_back(mk.first);
        place_at_this_level(bs.cache, cache, roots);
      }

      std::vector<ResultPtr> acc(layer.size());
      for (auto const& [e_lo, e_hi] : batches) {
        if (e_lo == e_hi) continue;
        // Per-slice marker for trace post-processing (avoidable-recompute
        // accounting): tags every op replayed in this iteration with the
        // enclosing batch loop's (mode, element-range-low). The mode label K
        // and the range low bound uniquely identify this slice within the
        // (single-mode-per-level) nesting. Gated on printing() so it is inert
        // unless a trace is being emitted.
        if (log::printing())
          log::log("BatchIter", toUtf8(std::wstring(K.full_label())), e_lo);
        bs.cache.reset();
        for (std::size_t m = 0; m != layer.size(); ++m) {
          auto const& [mem, Km] = layer[m];
          // Extend the enclosing batch context by THIS member's block and set
          // it on the scratch, so the re-entry's Enter-stage slice-on-use (and
          // its own sliced_leaf) slices every leaf carrying Km to this block
          // and composes inner slices on top -- exactly what the old per-member
          // `le_g` did on the leaf path, plus it now also slices a cached
          // intermediate fetched from an ancestor scope (the slice-on-use fix).
          // Rebuilt from `cache.batch_context()` (the enclosing context) each
          // member so contexts do not accumulate across members. The raw
          // leaf_evaluator is threaded down (type-erased into a std::function
          // so this template's self-instantiation is finite: the leaf-evaluator
          // type is std::function at every deeper level); the Enter stage does
          // the slicing.
          auto ctx = cache.batch_context();
          ctx.push_back({Km, synth_level(Km), {e_lo, e_hi}, Km});
          bs.cache.set_batch_context(std::move(ctx));
          bs.cache.set_custom_evaluator(make_batched_custom_evaluator(
              std::function<ResultPtr(node_t const&)>{leaf_evaluator},
              target_batch_size, accept, make_scope_guard, is_volatile,
              persistent_only, depth + 1, peak));
          ResultPtr part = evaluate_impl(*mem, leaf_evaluator, bs.cache);
          if (!acc[m])
            acc[m] = std::move(part);
          else
            acc[m]->add_inplace(*part);
        }
        // Fold this batch's scratch high-watermark into the global sink. The
        // next iteration calls bs.cache.reset(), which ZEROES the scratch
        // hwmark, so the fold MUST happen here (per batch), not after the
        // batches loop. The loop is serial (the nested evaluate() re-entry is
        // serial too), but the sink is atomic; a relaxed fetch-max CAS keeps it
        // correct regardless. A null sink skips the fold entirely, leaving
        // existing (sink-less) callers byte-unchanged.
        if (peak) {
          const double cand =
              static_cast<double>(bs.cache.working_set_hwmark());
          double cur = peak->load(std::memory_order_relaxed);
          while (cand > cur && !peak->compare_exchange_weak(
                                   cur, cand, std::memory_order_relaxed)) {
          }
        }
      }

      // Store the members into the real cache under the canonical-phase
      // convention (mirroring evaluate()'s Checked store), eagerly per layer
      // so later layers can seed them. The trigger is returned instead: its
      // Checked wrapper stores it (a direct store here would double-decay a
      // non-persistent trigger's life count).
      for (std::size_t m = 0; m != layer.size(); ++m) {
        auto const* mem = layer[m].first;
        if (mem == &node) {
          trigger_result = std::move(acc[m]);
          continue;
        }
        ResultPtr v = std::move(acc[m]);
        if (auto const ph = (*mem)->canon_phase(); ph != 1)
          v = v->mult_by_phase(ph);
        (void)cache.store(*mem, std::move(v));
      }
    }
    if (log::printing()) {
      auto scope_ctx = cache.batch_context();
      scope_ctx.push_back({K, synth_level(K), {0, 0}, K});
      log::log("BatchGroup", "End", log::scope_annot(scope_ctx));
    }
    SEQUANT_ASSERT(trigger_result);
    return trigger_result;
  };
}

/// \brief Builds a batched custom evaluator (see make_batched_custom_evaluator)
/// from a \p policy object, lifting the policy's canonical Tensor-based
/// volatile-leaf predicate to the EvalNode predicate the batched evaluator
/// expects.
///
/// Exactly equivalent to calling make_batched_custom_evaluator with:
///   - target_batch_size = policy.batch_target_size
///   - accept           = policy.is_batchable_index()  (derived role union)
///   - make_scope_guard = make_scope_guard (forwarded)
///   - is_volatile      = EvalNode lift of policy.is_volatile_leaf:
///       n.leaf() && n->is_tensor() && policy.is_volatile_leaf(n->as_tensor())
///     (when policy.is_volatile_leaf is empty, no node is volatile)
///
/// \param policy       BatchPolicy carrying the three batchability predicates.
/// \param yielder      The leaf evaluator (captured and forwarded).
/// \param make_scope_guard  Optional scope-guard factory (same semantics as in
///        make_batched_custom_evaluator; defaults to make_no_scope_guard).
/// \param peak      Optional PeakSink (same semantics as in
///        make_batched_custom_evaluator); defaults to null (no folding).
///        NOTE: \p peak is the 4th positional argument, AFTER \p
///        make_scope_guard -- a caller who wants the sink but not a custom
///        scope guard must still pass the scope-guard factory explicitly
///        (e.g. `make_evaluator(policy, leaf, make_no_scope_guard{},
///        &sink)`); passing `&sink` in the 3rd slot silently binds it to
///        \p make_scope_guard (via template deduction) and leaves \p peak
///        null.
template <class F, class ScopeGuardFactory = make_no_scope_guard>
[[nodiscard]] auto make_evaluator(BatchPolicy const& policy, F yielder,
                                  ScopeGuardFactory make_scope_guard = {},
                                  PeakSink peak = nullptr) {
  auto is_volatile_node = [p = policy.is_volatile_leaf](auto const& n) -> bool {
    if (!n.leaf() || !n->is_tensor()) return false;
    return p && p(n->as_tensor());
  };
  // BatchPolicy docs: an empty is_batchable_index or batch_target_size means
  // "no batching". Forwarding an empty std::function would instead throw
  // std::bad_function_call from batch_axis()/target_batch_size() at evaluation
  // time, so when either is unset, substitute predicates that decline batching
  // (accept nothing => batch_axis returns nullopt => target_batch_size is never
  // called) rather than partially-filled ones.
  // Runtime accept = the DERIVED union of both batchability roles: a mode is
  // accepted at runtime if it is batchable in EITHER the contracted or the
  // external role (see BatchPolicy::is_batchable_index()).
  std::function<bool(Index const&)> accept = policy.is_batchable_index();
  std::function<std::size_t(Index const&)> target = policy.batch_target_size;
  if (!accept || !target) {
    accept = [](Index const&) { return false; };
    target = [](Index const&) -> std::size_t { return 0; };
  }
  return make_batched_custom_evaluator(
      std::move(yielder), std::move(target), std::move(accept),
      std::move(make_scope_guard), std::move(is_volatile_node),
      policy.persistent_only, /*depth=*/0, peak);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
