#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/eval/fwd.hpp>

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/schedule_dump.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
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
/// (aliasing is evaluated at each call site using cache.alive, canon_phase,
/// and the requested layout). It is reported as a running max so it is
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
  // Structured runtime-schedule event, gated by SEQUANT_SCHED_DUMP (independent
  // of the trace level so it works without full tracing). Keyed by
  // node->hash_value() -- the SAME identity the IR emitter (schedule_dump.hpp)
  // writes -- so the schedule visualizer joins runtime lifetimes onto the IR
  // DAG by hash. Mode: Store (first build), Access (reuse), Release (last use).
  // Store-count > 1 for a hash means the value was rebuilt (recompute).
  static bool const sched_dump = std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
  if (sched_dump) {
    auto const cl = cm.life(node);
    auto const ml = cm.max_life(node);
    char const* const evm = (cl == 0)        ? "Release"
                            : (cl + 1 == ml) ? "Store"
                                             : "Access";
    std::cerr << "SCHEDULE_RUN_EVENT {\"hash\":\"" << node->hash_value()
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

template <typename... Args>
concept last_type_is_cache_manager = is_cache_manager_v<std::remove_cvref_t<
    std::tuple_element_t<sizeof...(Args) - 1, std::tuple<Args...>>>>;

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

/// \return the position of index \p ix in \p node's canonical result indices
///         (i.e. the corresponding tensor mode), or nullopt if absent.
[[nodiscard]] inline std::optional<std::size_t> index_position(
    meta::eval_node auto const& node, Index const& ix) {
  auto const& idxs = node->canon_indices();
  for (std::size_t p = 0; p < idxs.size(); ++p)
    if (idxs[p] == ix) return p;
  return std::nullopt;
}

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
template <Trace EvalTrace = Trace::Default,
          detail::CacheCheck Cache = detail::CacheCheck::Checked,
          meta::can_evaluate Node, typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate(Node const& node,         //
                   F const& leaf_evaluator,  //
                   CacheManager<N, FHC>& cache) {
  // Multiply a (possibly cached) result by its node's canonicalization phase.
  // Formerly the `mult_by_phase` lambda local to the Checked wrapper.
  auto apply_phase = [&cache](auto const& nd, ResultPtr res) -> ResultPtr {
    auto phase = nd->canon_phase();
    if (phase == 1) return res;

    ResultPtr post;
    auto time =
        detail::timed_eval_inplace([&]() { post = res->mult_by_phase(phase); });

    if constexpr (detail::trace(EvalTrace)) {
      size_t hwmark = log::bytes(cache, post).value;
      if (!cache.alive(nd)) hwmark += log::bytes(res).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{.mode = log::EvalMode::MultByPhase,
                                .time = time,
                                .mem_result = log::bytes(post),
                                .mem_alloc = log::bytes(post),
                                .mem_hwmark = {cache.note_working_set(hwmark)}};
      log::eval(stat, std::format("{} * {}", phase, nd->label()));
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
    for (std::size_t i = d - hops; i < d; ++i) {
      auto const& [axis, blk] = ctx[i];
      if (auto const p = index_position(nd, axis))
        value = value->slice_mode(*p, blk.first, blk.second);
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
    // Per-op BUILD event: finish_phase_b is the single choke point every
    // freshly computed node passes through -- leaves, custom-eval subtrees, and
    // standard contractions -- so this counts EVERY build (cached or not),
    // giving the schedule visualizer full per-node recompute coverage (the
    // cache Store/ Access/Release events above cover only cached nodes). Keyed
    // by hash_value() to join onto the IR DAG. Gated by SEQUANT_SCHED_DUMP.
    static bool const sched_dump = std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
    if (sched_dump) {
      // ctx = the active batch loops (mode -> block offset) this build ran
      // under. The visualizer counts DISTINCT ctx projections onto the modes a
      // node depends on: builds at a repeated projected slice are avoidable
      // recompute (a value rebuilt where an enclosing loop it is invariant to
      // advanced) vs the inherent per-block work batching requires.
      // Leaves are ACCESSED, not computed (leaf_evaluator just hands back a
      // ref to a precomputed input), so tag them "Fetch" -- cheap, not
      // recompute. Only internal nodes (contractions) are real "Build" work.
      char const* const bmode = f.node.leaf() ? "Fetch" : "Build";
      std::cerr << "SCHEDULE_RUN_EVENT {\"hash\":\"" << f.node->hash_value()
                << "\",\"mode\":\"" << bmode << "\"";
      // sig = the avoidable-recompute join key, matching the IR node record and
      // cost_profile's per-node label (result + sorted operand pair). Internal
      // nodes only; leaves (Fetch) are not tallied. Lets the visualizer join
      // this hash to cost_profile's per-node avoidable without reconstructing
      // the signature in the renderer.
      if (!f.node.leaf())
        std::cerr << ",\"sig\":\""
                  << eval::detail::sched_json_escape(eval::cost_op_signature(
                         f.node->canon_indices(),
                         f.node.left()->canon_indices(),
                         f.node.right()->canon_indices()))
                  << "\"";
      std::cerr << ",\"ctx\":[";
      bool first = true;
      for (auto const& [ix, blk] : cache.batch_context()) {
        // dep = does the node's subtree carry this loop mode (free or
        // contracted below)? If not, the node is INVARIANT to it and rebuilding
        // per block of it is avoidable recompute. find_leaf_carrying works in
        // the node's own label space, so no alpha-renaming reconciliation.
        bool const dep = find_leaf_carrying(f.node, ix).has_value();
        std::cerr << (first ? "" : ",") << "[\"" << toUtf8(ix.full_label())
                  << "\"," << blk.first << "," << (dep ? 1 : 0) << "]";
        first = false;
      }
      std::cerr << "]}\n";
    }
    if (!f.store_after) return rb;
    auto ptr = cache.store(f.node, apply_phase(f.node, std::move(rb)));
    if constexpr (detail::trace(EvalTrace))
      log::cache(f.node, cache, log::label(f.node));
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
          if (auto m = cache.access_at(f.node); m.ptr) {
            if constexpr (detail::trace(EvalTrace))
              log::cache(f.node, cache, log::label(f.node));
            // Slice-on-use: a value fetched `m.hops` scopes up does not have
            // this scope's (and any intervening) batch slices baked in, so
            // slice it to the current block for the loops the fetch crossed. A
            // local hit (hops == 0) or the OFF path (empty batch_context) is a
            // no-op, so this stays byte-identical to apply_phase() alone there.
            finalize(slice_to_use(apply_phase(f.node, m.ptr), f.node, m.hops));
            break;
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
                log::eval(
                    log::EvalStat{
                        .mode = log::eval_mode(f.node),
                        .time = time,
                        .mem_result = log::bytes(intercepted),
                        .mem_alloc = log::bytes(intercepted),
                        .mem_hwmark = {cache.note_working_set(hwmark)}},
                    log::label(f.node));
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
            log::eval(
                log::EvalStat{.mode = log::eval_mode(f.node),
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = {cache.note_working_set(hwmark)}},
                log::label(f.node));
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
        std::array<std::any, 2> const adj_ann{f.node.left()->annot(),
                                              f.node->annot()};
        ResultPtr result;
        auto time = detail::timed_eval_inplace(
            [&]() { result = f.left->adjoint(adj_ann); });

        if constexpr (detail::trace(EvalTrace)) {
          // `right` is null here (see log::bytes() null tolerance).
          size_t hwmark = log::bytes(cache, result).value;
          if (!cache.alive(f.node.left()) || f.node.left()->canon_phase() != 1)
            hwmark += log::bytes(f.left).value;
          hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
          log::eval(
              log::EvalStat{.mode = log::eval_mode(f.node),
                            .time = time,
                            .mem_result = log::bytes(result),
                            .mem_alloc = log::bytes(result),
                            .mem_hwmark = {cache.note_working_set(hwmark)},
                            .mem_left = log::bytes(f.left),
                            .mem_right = log::bytes(f.right)},
              log::label(f.node));
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
        if (f.node->op_type() == EvalOp::Sum) {
          time = detail::timed_eval_inplace(
              [&]() { result = f.left->sum(*f.right, ann); });
        } else {
          SEQUANT_ASSERT(f.node->op_type() == EvalOp::Product);
          // Consult the shaped-product hook (if set) before evaluating the
          // product. The hook receives the node (wrapped in a std::any as a
          // std::reference_wrapper so the full IR node is inspectable) plus the
          // evaluated operands and annotations; a non-null return *replaces*
          // the normal product (e.g. a shape-constrained emission of it), a
          // null return declines and the standard prod() below runs. An empty
          // hook is never consulted; default-empty => byte-identical behavior.
          auto const de_nest =
              f.node.left()->tot() && f.node.right()->tot() && !f.node->tot();
          if (auto const& hook = cache.shaped_product_hook(); hook) {
            time = detail::timed_eval_inplace([&]() {
              result =
                  hook(std::any{std::cref(f.node)}, *f.left, *f.right, ann);
            });
          }
          if (!result) {
            time = detail::timed_eval_inplace([&]() {
              result = f.left->prod(*f.right, ann,
                                    de_nest ? DeNest::True : DeNest::False);
            });
          }
        }

        SEQUANT_ASSERT(result);

        if constexpr (detail::trace(EvalTrace)) {
          // A cached child is *distinct* from the local left/right when its
          // canon_phase != 1, because apply_phase allocates a fresh buffer
          // while the cache still holds the pre-phase data. So only skip the
          // local's bytes when the cache aliases the same buffer (phase == 1).
          size_t hwmark = log::bytes(cache, result).value;
          if (!cache.alive(f.node.left()) || f.node.left()->canon_phase() != 1)
            hwmark += log::bytes(f.left).value;
          if (f.right && (!cache.alive(f.node.right()) ||
                          f.node.right()->canon_phase() != 1))
            hwmark += log::bytes(f.right).value;
          hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
          log::eval(
              log::EvalStat{.mode = log::eval_mode(f.node),
                            .time = time,
                            .mem_result = log::bytes(result),
                            .mem_alloc = log::bytes(result),
                            .mem_hwmark = {cache.note_working_set(hwmark)},
                            .mem_left = log::bytes(f.left),
                            .mem_right = log::bytes(f.right)},
              log::label(f.node));
        }
        log::release_after_op();
        finalize(finish_phase_b(f, std::move(result)));
        break;
      }
    }
  }

  return ret;
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

  result.pre = evaluate<EvalTrace>(node, leaf_evaluator, cache);

  auto time = detail::timed_eval_inplace([&]() {
    result.post = perm ? result.pre->permute(
                             std::array<std::any, 2>{node->annot(), layout})
                       : result.pre;
  });

  SEQUANT_ASSERT(result.post);

  // logging
  if constexpr (detail::trace(EvalTrace)) {
    if (perm) {
      // result.pre aliases the cache only when the inner evaluate returned
      // the cached buffer unchanged — i.e. the node is cached AND no
      // mult_by_phase fresh allocation happened (phase == 1).
      size_t hwmark = log::bytes(cache, result.post).value;
      if (!cache.alive(node) || node->canon_phase() != 1)
        hwmark += log::bytes(result.pre).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{.mode = log::EvalMode::Permute,
                                .time = time,
                                .mem_result = log::bytes(result.post),
                                .mem_alloc = log::bytes(result.post),
                                .mem_hwmark = {cache.note_working_set(hwmark)}};
      log::eval(stat, node->label());
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
  ResultPtr result;

  // pre comes back from the permute-wrapping evaluate; it aliases the
  // cache only when the inner evaluate returned the cached buffer
  // unchanged — i.e. node cached, phase == 1, AND no permute happened.
  bool const layout_is_default = (layout == decltype(layout){});

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
      // SumInplace allocates nothing: it writes into the accumulator.
      // hwmark counts the cache plus both operands live at this moment;
      // skip pre's bytes only when pre is the cached buffer itself.
      size_t hwmark = log::bytes(cache, result).value;
      if (!cache.alive(n) || n->canon_phase() != 1 || !layout_is_default)
        hwmark += log::bytes(pre).value;
      hwmark += cache.parent() ? cache.parent()->chain_residency() : 0;
      auto stat = log::EvalStat{.mode = log::EvalMode::SumInplace,
                                .time = time,
                                .mem_result = log::bytes(result),
                                .mem_alloc = {0},
                                .mem_hwmark = {cache.note_working_set(hwmark)}};
      log::eval(stat, n->label());
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
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Evaluate given node (or a range of nodes) using an empty cache
///        manager. Calls the other sequant::evaluate function overloads.
/// \see evaluate.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
  requires(!detail::last_type_is_cache_manager<Args...>)
ResultPtr evaluate(Args&&... args) {
  using Node = std::remove_cvref_t<decltype(detail::node0(
      detail::arg0(std::forward<Args>(args)...)))>;
  auto cache = CacheManager<Node>::empty();
  return evaluate<EvalTrace>(std::forward<Args>(args)..., cache);
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
/// (\c EvalExpr::batched_here; see the \c pick_sliceable lambda); it declines
/// if the node carries no accepted annotation. Annotations are authoritative --
/// there is no heuristic fallback, so a node the peak-constrained optimizer did
/// not ask to batch is never batched. It asks the backend to
/// partition \c K into contiguous, whole-tile element-range batches of at most
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
/// real-cache entries to pre-seed it with (registered persistent in the
/// scratch, so they survive the per-batch reset()).
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
/// canon_indices()) have batch-invariant full values; those that are alive
/// persistent entries of \p real are returned as seeds, and the caller copies
/// their values into the scratch before the batch loop.
template <typename TreeNode, bool FHC, typename Members>
[[nodiscard]] BatchedScratch<TreeNode, FHC> make_batched_scratch(
    Members const& members, CacheManager<TreeNode, FHC> const& real) {
  using Hasher = TreeNodeHasher<TreeNode, FHC>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;

  // The batch's EXTERNAL modes: obtained exactly as the evaluator obtains them
  // (partition each member root's batched_here() by BatchModeType). An External
  // mode is an external that survives FREE onto a node's result, so a node
  // carrying one is NOT batch-invariant under that mode -- its value depends on
  // the external slice. When the caller nests an External mode OUTSIDE a
  // Contracted one (the External mode is sliced by an outer re-entry, then this
  // scratch batches an inner Contracted mode), a persistent intermediate that
  // carries the External mode but not the Contracted `mode` would look seedable
  // under `mode` alone -- yet seeding its full (unsliced-external) value would
  // be wrong under the outer slice. Tracking the External modes in the
  // signature (below) forbids seeding/sharing such nodes. When there is no
  // External mode this list is empty and every External-derived test is a
  // no-op, keeping the Contracted-only behavior byte-identical.
  container::svector<Index> ext_axes;
  for (auto const& [root, mode] : members) {
    if (root->leaf()) continue;
    for (auto const& [ix, knd] : (*root)->batched_here())
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

  auto ext_sig_of = [&ext_axes](TreeNode const& n) {
    container::svector<std::optional<std::size_t>> s;
    s.reserve(ext_axes.size());
    for (auto const& e : ext_axes) s.push_back(index_position(n, e));
    return s;
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
    // member roots themselves are accumulated by the caller, not cached here
    if (root->leaf()) continue;
    visit(visit, root->left(), mode);
    visit(visit, root->right(), mode);
  }

  std::unordered_map<TreeNode, std::size_t, Hasher, Comp> reg;
  std::unordered_set<TreeNode, Hasher, Comp> seed_keys;
  std::vector<TreeNode const*> seeds;
  for (auto const& [ptr, e] : meta) {
    if (!e.consistent) continue;  // ambiguous slicing: never share
    // A node carrying ANY batched External mode has an external slice its
    // seeded-full real-cache value would ignore -- so it is never seedable.
    bool const carries_ext =
        std::any_of(e.ext_sig.begin(), e.ext_sig.end(),
                    [](auto const& p) { return p.has_value(); });
    bool const seedable =
        !e.sig && !carries_ext && real.persistent(*ptr) && real.alive(*ptr);
    if (seedable) {
      seeds.push_back(ptr);
      seed_keys.insert(*ptr);
      reg.emplace(*ptr, e.count);  // count is ignored for persistent entries
    } else if (e.count >= 2) {
      reg.emplace(*ptr, e.count);
    }
  }
  auto is_persistent = [seed_keys = std::move(seed_keys)](TreeNode const& n) {
    return seed_keys.contains(n);
  };
  CacheManager<TreeNode, FHC> scratch{std::move(reg), std::move(is_persistent)};
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
        auto const& [axis, blk] = ctx[i];
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
    // Mode selection is SLICEABILITY-AWARE and realizes the optimizer's
    // multi-mode nesting one mode per depth level. candidate_axes lists this
    // node's batch modes in the optimizer's annotated order (see
    // EvalExpr::batched_here), keeping the accepted annotations.
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
      for (auto const& entry : n->batched_here())
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
      for (Index const& ix : candidate_axes(n)) {
        auto const lf = find_leaf_carrying(n, ix);
        if (!lf) continue;
        // sliced_leaf (not the raw evaluator): at depth > 0 the carrier leaf
        // must be sliced to the enclosing blocks so an already-outer-sliced
        // mode yields a SINGLE batch and is skipped -- the "K is not re-picked"
        // invariant. At depth 0 (empty enclosing context) this is the raw leaf.
        auto b = sliced_leaf(lf->first)->mode_batches(lf->second,
                                                      target_batch_size(ix));
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
    // mode absent from batched_here() is Contracted -- keeping the
    // Contracted-only path (no External entry) byte-identical to before this
    // branch existed.
    BatchModeType picked_kind = BatchModeType::Contracted;
    for (auto const& [ix, knd] : node->batched_here())
      if (ix == K) {
        picked_kind = knd;
        break;
      }

    // Per-level placement (order-aware only). Replaces hoist_invariants. This
    // firing realizes a batch loop over `K` at runtime `depth`. A
    // member-subtree node INVARIANT to this loop (it does not carry `K` on its
    // result) is built ONCE at its home level and served to every batch body
    // through the scope chain, rather than rebuilt per batch. A node's
    // residency (the batch modes it is variant to) is the UNION of two per-node
    // signals:
    //   - sliced_modes() : the EXTERNAL (occ) modes, from the cross-occurrence
    //     lifetime-mask meet (consistent placement across occurrences -> CSE);
    //   - contracted_modes() : the enclosing CONTRACTED (aux) modes the node
    //     carries open, emitted per-occurrence by the cost model (the piece the
    //     external-only mask cannot express -- a node is variant to an outer
    //     aux loop by carrying that aux free on its result).
    // A node is invariant to THIS loop iff `K` is NOT in its union. Its HOME
    // level is the deepest ENCLOSING batch_context entry whose mode is in the
    // union (the innermost enclosing loop it is variant to); -1 (the chain root
    // / run-term cache) if it is invariant to the whole nest. The node is built
    // once at that level (walk-up), sliced to its home blocks, and reused
    // across this loop's batches. A node carrying `K` (K in its union) is
    // loop-LOCAL: it is left to inline evaluation (descend), which finds any
    // deeper-hoisted invariants through the chain -- so descending never
    // rebuilds them.
    //
    // This reproduces hoist_invariants' walk-up structure exactly, with the
    // former per-node scalar placement level replaced by the union-derived
    // home level and the `sl < depth && !ext_loop_local` predicate replaced by
    // `K not in union` (which subsumes ext_loop_local: an External carrier has
    // K in sliced_modes -> in union -> loop-local -> descended, never
    // hoisted). The order-aware GATE is the emitted `batch_order_aware()` bit
    // (true for every node the order-aware cost model emitted, including a
    // whole-nest invariant whose union is empty): on the OFF path every node
    // is order-blind, so `targets` is empty, set_parent is NOT wired, and the
    // per-batch replay runs exactly as before. The bit is a positive signal an
    // empty union cannot provide -- it is what distinguishes an OFF-path
    // all-full node (do not hoist) from an order-aware whole-nest invariant
    // (hoist to the root).
    auto place_at_this_level =
        [&](auto& scratch_cache, auto& parent_cache,
            std::vector<node_t const*> const& member_roots) {
          // The ENCLOSING batch loops (strictly OUTER to this firing); this
          // level's mode K and any INNER loop are NOT in it. A node is
          // hoistable here iff EVERY residency (union) mode is one of these
          // outer loops: then it is invariant to this loop AND to every inner
          // loop, so its home is its deepest enclosing residency level and it
          // is built once there. If a union mode is K (loop-local) or an INNER
          // mode (its home is a deeper loop), it is not all-outer -> descend,
          // so the deeper level handles it sliced. This is exactly
          // hoist_invariants' `scope_level < depth`: the deepest residency
          // being outer <=> ALL residency outer (deepest is the max), and a
          // node carrying an inner mode is (correctly) not hoisted at this
          // outer level -- the bug an over-eager `K not in union` predicate
          // caused (holding an aux-carrier full over aux at the outer occ
          // level).
          auto const& ectx = parent_cache.batch_context();  // enclosing loops
          auto in_ectx = [&ectx](Index const& m) -> bool {
            for (auto const& e : ectx)
              if (e.first == m) return true;
            return false;
          };
          auto residency_all_outer = [&in_ectx](node_t const& n) -> bool {
            for (auto const& ix : n->sliced_modes())
              if (!in_ectx(ix)) return false;
            for (auto const& ix : n->contracted_modes())
              if (!in_ectx(ix)) return false;
            return true;
          };
          // A node carrying an EXTERNAL batched_here() stamp absent from its
          // sliced_modes is a MEET-DEMOTED external carrier: the
          // cross-occurrence meet intersected that external slot to empty
          // because occurrences bind it to DIFFERENT (proto-incompatible)
          // blocks (the cross-pair two-PNO giants), yet the node still carries
          // that external mode FREE in THIS occurrence -- so its per-occurrence
          // value is a per-external-BLOCK (scattered/sliced) result. Hoisting
          // it would build it at its home level with the demoted external mode
          // UNSLICED (not in any enclosing block), i.e. materialize the FULL
          // scattered giant. Never hoist such a node: descend so the batched
          // evaluator slices its external mode per occurrence (exactly as
          // hoist_invariants did via per-occurrence scope_level). A node with
          // only Contracted stamps (genuinely external- invariant, e.g. a
          // summed-K intermediate) has no such stamp and is still hoisted once;
          // a consistently-sliced external node has the mode in sliced_modes
          // (not demoted) and is placed/sliced at its loop.
          auto has_demoted_external = [](node_t const& n) -> bool {
            auto const& sm = n->sliced_modes();
            for (auto const& [ix, knd] : n->batched_here())
              if (knd == BatchModeType::External &&
                  std::find(sm.begin(), sm.end(), ix) == sm.end())
                return true;
            return false;
          };
          std::vector<node_t const*> targets;
          auto collect = [&](auto&& self, node_t const& n) -> void {
            if (n.leaf()) return;
            if (n->batch_order_aware() && residency_all_outer(n) &&
                !has_demoted_external(n) && !subtree_any(n, is_volatile)) {
              if (std::none_of(targets.begin(), targets.end(),
                               [&](node_t const* p) { return eq(*p, n); }))
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
          auto in_union = [](node_t const& n, Index const& m) -> bool {
            auto const& sm = n->sliced_modes();
            if (std::find(sm.begin(), sm.end(), m) != sm.end()) return true;
            auto const& cm = n->contracted_modes();
            return std::find(cm.begin(), cm.end(), m) != cm.end();
          };
          for (node_t const* dptr : targets) {
            node_t const& d = *dptr;
            // Home level = deepest enclosing-context entry whose mode is in d's
            // union; -1 => invariant to the whole nest (chain root).
            int rl = -1;
            for (int i = static_cast<int>(ectx.size()) - 1; i >= 0; --i)
              if (in_union(d, ectx[i].first)) {
                rl = i;
                break;
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
            ResultPtr built = evaluate(d, sliced_leaf);
            if (auto const ph = d->canon_phase(); ph != 1)
              built = built->mult_by_phase(ph);
            (void)target->store(d, std::move(built));
          }
        };

    // DEBUG (behavior-neutral): log the trigger's depth, picked mode + kind,
    // and its full batched_here() annotation + result indices, to diagnose
    // nested re-batching of a single aux mode. Emitted only when tracing is on.
    if (log::printing()) {
      std::string annot;
      for (auto const& [ix, knd] : node->batched_here()) {
        annot += toUtf8(ix.full_label());
        annot += (knd == BatchModeType::External ? ":ext " : ":con ");
      }
      std::string res;
      for (auto const& ix : node->canon_indices()) {
        res += toUtf8(ix.full_label());
        res += " ";
      }
      log::log("BatchAxes",
               std::format(
                   "depth={} picked={}:{} nbatches={} annot=[{}] result=[{}]",
                   depth, toUtf8(K.full_label()),
                   picked_kind == BatchModeType::External ? "ext" : "con",
                   batches.size(), annot, res));
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
      auto const carrier = find_leaf_carrying(node, K);
      SEQUANT_ASSERT(carrier && "no leaf carries the external batch mode");
      // The carrier leaf supplies K's tiling for pre-sizing. sliced_leaf (not
      // the raw evaluator): at depth 0 this is the FULL carrier (empty
      // enclosing context), but at depth > 0 it must be sliced to the OUTER
      // enclosing blocks so the scatter dest is pre-sized within the outer
      // block, not at the full outer extent.
      ResultPtr const carrier_full = sliced_leaf(carrier->first);

      // A single-node scratch: an external mode is not a
      // persistent-final sharing mode, so the group/replay machinery (which
      // co-batches cross-term contracted finals) does not apply -- scatter just
      // this node. The scratch still dedups repeats WITHIN the node's subtree.
      std::vector<member_t> solo{{&node, K}};
      auto bs = detail::make_batched_scratch(solo, cache);
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
        ctx.push_back({K, {e_lo, e_hi}});
        bs.cache.set_batch_context(std::move(ctx));
        bs.cache.set_custom_evaluator(make_batched_custom_evaluator(
            std::function<ResultPtr(node_t const&)>{leaf_evaluator},
            target_batch_size, accept, make_scope_guard, is_volatile,
            persistent_only, depth + 1, peak));
        ResultPtr part = evaluate(node, leaf_evaluator, bs.cache);
        // Pre-size on the first block (learns the result's non-mode extents and
        // kind from the block partial; the external mode is widened to full).
        if (!dest)
          dest = part->pre_sized_zeros_over_mode(*dest_mode, *carrier_full,
                                                 carrier->second);
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
      log::log("BatchGroup", "Begin",
               std::format("{} members co-evaluated over {} aux batches",
                           n_members, batches.size()));
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
          ctx.push_back({Km, {e_lo, e_hi}});
          bs.cache.set_batch_context(std::move(ctx));
          bs.cache.set_custom_evaluator(make_batched_custom_evaluator(
              std::function<ResultPtr(node_t const&)>{leaf_evaluator},
              target_batch_size, accept, make_scope_guard, is_volatile,
              persistent_only, depth + 1, peak));
          ResultPtr part = evaluate(*mem, leaf_evaluator, bs.cache);
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
    if (log::printing()) log::log("BatchGroup", "End");
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
