#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/eval/fwd.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_node.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/range/operations.hpp>

#include <algorithm>
#include <any>
#include <chrono>
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
/// Eval | <mode> | <time> | [left=L | right=R |] result=X | alloc=A | hw=H | rss=R | <label>
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
  auto const rss_s = std::format("rss={}", to_string(rss()));
  if (stat.mem_left) {
    SEQUANT_ASSERT(stat.mem_right);
    log("Eval",                                               //
        to_string(stat.mode),                                 //
        stat.time,                                            //
        std::format("left={}", to_string(*stat.mem_left)),    //
        std::format("right={}", to_string(*stat.mem_right)),  //
        result_s, alloc_s, hw_s, rss_s,                       //
        args...);
  } else {
    log("Eval",                //
        to_string(stat.mode),  //
        stat.time,             //
        result_s, alloc_s, hw_s, rss_s, args...);
  }
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
/// Each such index `K` is a valid axis to evaluate the subtree rooted at
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

/// \brief A default axis to batch the subtree at \p node over: the contracted
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
/// \tparam Cache If CacheCache::Checked (default) the \p cache will be
///               checked before evaluating. It is used to detect the base case
///               for recursion to prevent infinite recursion.
/// \param node A node that can be evaluated using \p le as the leaf
///             evaluator.
/// \param le The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default,
          detail::CacheCheck Cache = detail::CacheCheck::Checked,
          meta::can_evaluate Node, typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate(Node const& node,  //
                   F const& le,       //
                   CacheManager<N, FHC>& cache) {
  if constexpr (Cache == detail::CacheCheck::Checked) {  // return from cache if
                                                         // found

    auto mult_by_phase = [&node, &cache](ResultPtr res) {
      auto phase = node->canon_phase();
      if (phase == 1) return res;

      ResultPtr post;
      auto time = detail::timed_eval_inplace(
          [&]() { post = res->mult_by_phase(phase); });

      if constexpr (detail::trace(EvalTrace)) {
        size_t hwmark = log::bytes(cache, post).value;
        if (!cache.alive(node)) hwmark += log::bytes(res).value;
        auto stat =
            log::EvalStat{.mode = log::EvalMode::MultByPhase,
                          .time = time,
                          .mem_result = log::bytes(post),
                          .mem_alloc = log::bytes(post),
                          .mem_hwmark = {cache.note_working_set(hwmark)}};
        log::eval(stat, std::format("{} * {}", phase, node->label()));
      }
      return post;
    };

    if (auto ptr = cache.access(node); ptr) {
      if constexpr (detail::trace(EvalTrace))
        log::cache(node, cache, log::label(node));

      return mult_by_phase(ptr);
    } else if (cache.exists(node)) {
      auto ptr = cache.store(
          node,
          mult_by_phase(evaluate<EvalTrace, detail::CacheCheck::Unchecked>(
              node, le, cache)));
      if constexpr (detail::trace(EvalTrace))
        log::cache(node, cache, log::label(node));

      return mult_by_phase(ptr);
    } else {
      // do nothing
    }
  }

  ResultPtr result;
  ResultPtr left;
  ResultPtr right;

  log::Duration time;

  // Custom-evaluator interception: before the standard scheme, a non-leaf node
  // may be evaluated by the cache's custom evaluator (e.g. blocked over a
  // contracted index to bound peak memory). A non-null result is used (and
  // cached by the Checked wrapper) as-is; null declines to the standard scheme
  // below. See CacheManager::custom_evaluator_type.
  if (!node.leaf()) {
    if (auto const& custom_eval = cache.custom_evaluator(); custom_eval) {
      ResultPtr intercepted;
      time = detail::timed_eval_inplace(
          [&]() { intercepted = custom_eval(node, cache); });
      if (intercepted) {
        if constexpr (detail::trace(EvalTrace)) {
          log::eval(log::EvalStat{.mode = log::eval_mode(node),
                                  .time = time,
                                  .mem_result = log::bytes(intercepted),
                                  .mem_alloc = log::bytes(intercepted),
                                  .mem_hwmark = {cache.note_working_set(
                                      log::bytes(cache, intercepted).value)}},
                    log::label(node));
        }
        return intercepted;
      }
    }
  }

  if (node.leaf()) {
    time = detail::timed_eval_inplace([&]() { result = le(node); });
  } else if (node->op_type() == EvalOp::Adjoint) {
    // Unary IR op: dispatch on left operand only; right is the Constant(1)
    // sentinel kept around to preserve FullBinaryNode's invariant. We
    // intentionally skip evaluating the sentinel — leaf evaluators that
    // can't manufacture scalar constants (rare in practice but possible)
    // would otherwise be invoked needlessly.
    left = evaluate<EvalTrace>(node.left(), le, cache);
    SEQUANT_ASSERT(left);
    std::array<std::any, 2> const adj_ann{node.left()->annot(), node->annot()};
    time =
        detail::timed_eval_inplace([&]() { result = left->adjoint(adj_ann); });
  } else {
    left = evaluate<EvalTrace>(node.left(), le, cache);
    right = evaluate<EvalTrace>(node.right(), le, cache);
    SEQUANT_ASSERT(left);
    SEQUANT_ASSERT(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};
    if (node->op_type() == EvalOp::Sum) {
      time = detail::timed_eval_inplace(
          [&]() { result = left->sum(*right, ann); });
    } else {
      SEQUANT_ASSERT(node->op_type() == EvalOp::Product);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();
      time = detail::timed_eval_inplace([&]() {
        result =
            left->prod(*right, ann, de_nest ? DeNest::True : DeNest::False);
      });
    }
  }

  SEQUANT_ASSERT(result);

  // logging
  if constexpr (detail::trace(EvalTrace)) {
    if (node.leaf()) {
      log::eval(log::EvalStat{.mode = log::eval_mode(node),
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = {cache.note_working_set(
                                  log::bytes(cache, result).value)}},
                log::label(node));
    } else {
      // A cached child is *distinct* from the local left/right when its
      // canon_phase != 1, because mult_by_phase allocates a fresh buffer
      // while the cache still holds the pre-phase data. So only skip the
      // local's bytes when the cache aliases the same buffer (phase == 1).
      // Adjoint nodes evaluate only the left operand (the right child is the
      // sentinel Constant(1) — see the Adjoint branch above), so `right` is
      // null; log::bytes() tolerates a null shared_ptr for that reason.
      size_t hwmark = log::bytes(cache, result).value;
      if (!cache.alive(node.left()) || node.left()->canon_phase() != 1)
        hwmark += log::bytes(left).value;
      if (right &&
          (!cache.alive(node.right()) || node.right()->canon_phase() != 1))
        hwmark += log::bytes(right).value;
      log::eval(log::EvalStat{.mode = log::eval_mode(node),
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = {cache.note_working_set(hwmark)},
                              .mem_left = log::bytes(left),
                              .mem_right = log::bytes(right)},
                log::label(node));
    }
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param node A node that can be evaluated using \p le as the leaf
///             evaluator.
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
/// \param le The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate Node, typename F,
          typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>  //
ResultPtr evaluate(Node const& node,           //
                   auto const& layout,         //
                   F const& le,                //
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

  result.pre = evaluate<EvalTrace>(node, le, cache);

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
/// \param nodes A range of node that can be evaluated using \p le as the
///              leaf evaluator. The evaluation result of the elements of
///              \p nodes will be summed up.
///
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
///               The results of each element from \p nodes will be permuted
///               to this layout before being summed.
///
/// \param le The leaf evaluator that satisfies
///           `meta::leaf_node_evaluator<Node, F>`.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& nodes,  //
                   auto const& layout,  //
                   F const& le, CacheManager<N, FHC>& cache) {
  ResultPtr result;

  // pre comes back from the permute-wrapping evaluate; it aliases the
  // cache only when the inner evaluate returned the cached buffer
  // unchanged — i.e. node cached, phase == 1, AND no permute happened.
  bool const layout_is_default = (layout == decltype(layout){});

  for (auto&& n : nodes) {
    if (!result) {
      result = evaluate<EvalTrace>(n, layout, le, cache);
      continue;
    }

    ResultPtr pre = evaluate<EvalTrace>(n, layout, le, cache);
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
/// \param nodes A range of node that can be evaluated using \p le as the
///              leaf evaluator. The evaluation result of the elements of
///              \p nodes will be summed up.
///
/// \param le The leaf evaluator that satisfies
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
                   F const& le, CacheManager<N, FHC>& cache) {
  using annot_type = decltype([](std::ranges::range_value_t<Nodes> const& n) {
    return n->annot();
  });

  static_assert(std::is_default_constructible_v<annot_type>);
  return evaluate(nodes, annot_type{}, le, cache);
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
/// axis \c K via `batch_axis(node, accept)` (declining if none). It asks the
/// backend to partition \c K into contiguous element-range batches of about
/// \p target_batch_size elements each (Result::mode_batches); if that yields at
/// most one batch it declines (so small / unselected indices are left to the
/// standard scheme).
///
/// Otherwise it *replays the build of every compatible persistent final* in
/// the same batch passes: the group is the trigger node plus every key of
/// \p cache that is registered persistent, not yet alive, and batches over an
/// axis with the identical realized partition. Per batch, each group member is
/// evaluated by the standard scheme -- with every leaf carrying the member's
/// batch axis sliced to the batch's element range -- on a shared *registered*
/// scratch cache (see detail::make_batched_scratch), so sub-intermediates
/// repeated within a member (canonically-equal siblings) or shared between
/// members are evaluated once per batch, exactly as the real cache would share
/// them; the per-member partials are summed across batches. This is exact
/// because `sum_K = sum_{batches} sum_{K in batch}`, and never materializes
/// the whole batch-axis extent of any intermediate at once. Completed members
/// are stored into \p cache (canonical-phase convention); the trigger's result
/// is returned for evaluate() to cache as usual. Members nested inside other
/// members evaluate in earlier passes and are then seeded (slice-free w.r.t.
/// the outer batch axis) or re-derived sliced in the outer pass. Considering a
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
/// \param le the leaf evaluator (captured).
/// \param target_batch_size the desired size of each batch *in elements* (a
///        user knob; no memory model is assumed). Backend-neutral: a tiled
///        backend rounds batch boundaries to tile boundaries, so realized
///        batches are uneven and each covers at least this many elements where
///        possible.
/// \param accept predicate selecting which contracted indices may be batched
///        (e.g. only those in the auxiliary/RI IndexSpace). Defaults to any.
/// \param make_scope_guard factory, called with the batch count, returning an
///        RAII object held for the duration of the batched partial
///        contractions; a backend may use it to relax block-sparse screening
///        (scaled by the batch count) so per-batch screening does not drop
///        small contributions that are significant once summed over the full
///        batch axis. Defaults to a no-op (make_no_scope_guard).
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
/// threshold over the *full* batch axis is not dropped in every individual
/// batch (which would lose its contribution to the sum). The factory is called
/// with the batch count, so the backend can scale the relaxation accordingly
/// (e.g. divide a Cauchy-Schwarz norm-product screening threshold by n_batches:
/// the bound for a sub-sum over 1/n of the batch axis is ~1/n of the full
/// bound). See make_batched_custom_evaluator's \p make_scope_guard parameter.
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
/// batch axis (a volatile subtree is rebuilt every evaluation, so batching it
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
/// \p members (each a subtree root paired with its batch axis).
///
/// Walks every member subtree with the same pruned counting walk as
/// cache_manager() (descend on first visit of a canonical-equal node, so
/// counts match access counts under caching -- and also on a re-encounter
/// whose slicing signature differs from the first visit's, so that
/// descendants' signatures under an inconsistently-sliced occurrence are
/// recorded rather than hidden by the prune) and registers every internal
/// subnode that repeats AND has a consistent slicing signature -- the position
/// of the containing member's batch axis in the subnode's canon_indices(), or
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
/// the axis -- the axis is contracted at the member's root, so a subtree
/// containing an axis-carrying leaf carries the axis free in its
/// canon_indices()) have batch-invariant full values; those that are alive
/// persistent entries of \p real are returned as seeds, and the caller copies
/// their values into the scratch before the batch loop.
template <typename TreeNode, bool FHC, typename Members>
[[nodiscard]] BatchedScratch<TreeNode, FHC> make_batched_scratch(
    Members const& members, CacheManager<TreeNode, FHC> const& real) {
  using Hasher = TreeNodeHasher<TreeNode, FHC>;
  using Comp = TreeNodeEqualityComparator<TreeNode>;
  struct Meta {
    std::size_t count = 0;
    std::optional<std::size_t> sig;
    bool consistent = true;
  };
  std::unordered_map<TreeNode const*, Meta, Hasher, Comp> meta;

  auto visit = [&meta](auto&& self, TreeNode const& n,
                       Index const& axis) -> void {
    if (n.leaf()) return;
    auto const sig = index_position(n, axis);
    auto const [it, first] = meta.try_emplace(&n);
    auto& e = it->second;
    if (first)
      e.sig = sig;
    else if (e.sig != sig)
      e.consistent = false;
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
    // is evaluated per occurrence, not served from the scratch at n.
    if (!first && e.sig == sig) return;
    self(self, n.left(), axis);
    self(self, n.right(), axis);
  };
  for (auto const& [root, axis] : members) {
    // member roots themselves are accumulated by the caller, not cached here
    if (root->leaf()) continue;
    visit(visit, root->left(), axis);
    visit(visit, root->right(), axis);
  }

  std::unordered_map<TreeNode, std::size_t, Hasher, Comp> reg;
  std::unordered_set<TreeNode, Hasher, Comp> seed_keys;
  std::vector<TreeNode const*> seeds;
  for (auto const& [ptr, e] : meta) {
    if (!e.consistent) continue;  // ambiguous slicing: never share
    bool const seedable = !e.sig && real.persistent(*ptr) && real.alive(*ptr);
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
  return {CacheManager<TreeNode, FHC>{std::move(reg), std::move(is_persistent)},
          std::move(seeds)};
}

}  // namespace detail

template <typename F, typename IndexPredicate = accept_any_index,
          typename ScopeGuardFactory = make_no_scope_guard,
          typename IsVolatile = never_volatile>
[[nodiscard]] auto make_batched_custom_evaluator(
    F le, std::function<std::size_t(Index const&)> target_batch_size,
    IndexPredicate accept = {}, ScopeGuardFactory make_scope_guard = {},
    IsVolatile is_volatile = {}) {
  return [le = std::move(le), target_batch_size = std::move(target_batch_size),
          accept, is_volatile,
          make_scope_guard](auto const& node, auto& cache) -> ResultPtr {
    auto const K = batch_axis(node, accept);
    if (!K) return nullptr;

    // Persistence gate: only stream subtrees that are amplitude-independent
    // (built once). If the subtree contains a volatile leaf it is rebuilt on
    // every evaluation, so batching pays the partition + relaxed-screening cost
    // each pass for no lasting memory benefit -- decline to the standard
    // scheme. Default never_volatile => no gating (original behavior).
    if (subtree_any(node, is_volatile)) return nullptr;

    auto const leaf = find_leaf_carrying(node, *K);
    if (!leaf) return nullptr;
    auto const batches =
        le(leaf->first)->mode_batches(leaf->second, target_batch_size(*K));

    if (batches.size() <= 1)
      return nullptr;  // nothing to gain (or unbatchable)

    using node_t = std::remove_cvref_t<decltype(node)>;
    using member_t = std::pair<node_t const*, Index>;
    TreeNodeEqualityComparator<node_t> const eq;

    // The replay group: the trigger plus every registered persistent key that
    // is not yet alive and batches over an axis with the identical realized
    // partition. All compatible persistent finals stream over the batch axis
    // in the same passes, so sub-intermediates shared between them (wherever
    // the scratch's slicing-signature guard admits sharing -- equal canonical
    // positions of the batch axis plus equal element ranges imply identical
    // slices) are evaluated once per batch instead of once per consumer.
    // The cost of considering a candidate is one leaf evaluation (the
    // mode_batches probe). With an unregistered (empty) real cache the group
    // is just the trigger.
    std::vector<member_t> group{{&node, *K}};
    cache.for_each_key([&](node_t const& k) {
      if (!cache.persistent(k) || cache.alive(k)) return;
      if (eq(k, node)) return;  // the trigger occupies its own slot
      auto const Kk = batch_axis(k, accept);
      if (!Kk) return;
      if (subtree_any(k, is_volatile)) return;  // defensive: P implies NV
      auto const lk = find_leaf_carrying(k, *Kk);
      if (!lk) return;
      if (le(lk->first)->mode_batches(lk->second, target_batch_size(*Kk)) !=
          batches)
        return;
      group.emplace_back(&k, *Kk);
    });

    // Layer by nesting: a member whose subtree contains another member
    // evaluates in a later layer, with the inner result by then alive in the
    // real cache -- seeded into the outer pass when slice-free w.r.t. the
    // outer batch axis, re-derived sliced (correct, unshared) otherwise.
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

    // RAII scope for the batched partial contractions; a backend-supplied
    // factory may relax block-sparse screening here (scaled by the batch count)
    // so per-batch screening does not drop contributions that survive over the
    // full batch axis.
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

      std::vector<ResultPtr> acc(layer.size());
      for (auto const& [e_lo, e_hi] : batches) {
        if (e_lo == e_hi) continue;
        bs.cache.reset();
        for (std::size_t m = 0; m != layer.size(); ++m) {
          auto const& [mem, Km] = layer[m];
          // leaf evaluator that slices every leaf carrying the member's batch
          // axis to this element batch; others pass through unchanged.
          auto le_g = [&le, &Km, e_lo = e_lo,
                       e_hi = e_hi](auto const& leaf_node) -> ResultPtr {
            ResultPtr r = le(leaf_node);
            if (auto const p = index_position(leaf_node, Km))
              return r->slice_mode(*p, e_lo, e_hi);
            return r;
          };
          ResultPtr part = evaluate(*mem, le_g, bs.cache);
          if (!acc[m])
            acc[m] = std::move(part);
          else
            acc[m]->add_inplace(*part);
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
    SEQUANT_ASSERT(trigger_result);
    return trigger_result;
  };
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
