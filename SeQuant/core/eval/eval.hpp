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

#include <chrono>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <any>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace sequant {

namespace log {

using Duration = std::chrono::nanoseconds;

struct Bytes {
  size_t value;
};

template <typename T, typename... Ts>
[[nodiscard]] inline auto bytes(T const& arg, Ts const&... args) {
  auto one = [](auto const& a) -> size_t {
    if constexpr (requires { a->size_in_bytes(); })
      return a->size_in_bytes();
    else
      return a.size_in_bytes();
  };
  return Bytes{(one(arg) + ... + one(args))};
}

[[nodiscard]] inline auto to_string(Bytes bs) noexcept {
  return std::format("{}B", bs.value);
}

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
    return node->is_product() ? EvalMode::Product
           : node->is_sum()   ? EvalMode::Sum
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
/// Eval | <mode> | <time> | [left=L | right=R |] result=X | alloc=A | hw=H | <label>
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
/// which writes into the accumulator and allocates nothing. mem_hwmark
/// is the live working set during the op:
///
///   bytes(cache) + bytes(result) + bytes of each operand not aliased
///                                  to a cache entry
///
/// Aliasing is evaluated at each call site using cache.alive,
/// canon_phase, and the requested layout.
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
  auto const result_s = std::format("result={}", to_string(stat.mem_result));
  auto const alloc_s = std::format("alloc={}", to_string(stat.mem_alloc));
  auto const hw_s = std::format("hw={}", to_string(stat.mem_hwmark));
  if (stat.mem_left) {
    SEQUANT_ASSERT(stat.mem_right);
    log("Eval",                                               //
        to_string(stat.mode),                                 //
        stat.time,                                            //
        std::format("left={}", to_string(*stat.mem_left)),    //
        std::format("right={}", to_string(*stat.mem_right)),  //
        result_s, alloc_s, hw_s,                              //
        args...);
  } else {
    log("Eval",                //
        to_string(stat.mode),  //
        stat.time,             //
        result_s, alloc_s, hw_s, args...);
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

namespace {

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

}  // namespace

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

namespace {
[[nodiscard]] consteval bool trace(Trace t) noexcept { return t == Trace::On; }
}  // namespace

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
          CacheCheck Cache = CacheCheck::Checked, meta::can_evaluate Node,
          typename F, typename N, bool FHC>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate(Node const& node,  //
                   F const& le,       //
                   CacheManager<N, FHC>& cache) {
  if constexpr (Cache == CacheCheck::Checked) {  // return from cache if found

    auto mult_by_phase = [&node, &cache](ResultPtr res) {
      auto phase = node->canon_phase();
      if (phase == 1) return res;

      ResultPtr post;
      auto time =
          timed_eval_inplace([&]() { post = res->mult_by_phase(phase); });

      if constexpr (trace(EvalTrace)) {
        size_t hwmark = log::bytes(cache, post).value;
        if (!cache.alive(node)) hwmark += log::bytes(res).value;
        auto stat = log::EvalStat{.mode = log::EvalMode::MultByPhase,
                                  .time = time,
                                  .mem_result = log::bytes(post),
                                  .mem_alloc = log::bytes(post),
                                  .mem_hwmark = {hwmark}};
        log::eval(stat, std::format("{} * {}", phase, node->label()));
      }
      return post;
    };

    if (auto ptr = cache.access(node); ptr) {
      if constexpr (trace(EvalTrace)) log::cache(node, cache, log::label(node));

      return mult_by_phase(ptr);
    } else if (cache.exists(node)) {
      auto ptr = cache.store(
          node, mult_by_phase(evaluate<EvalTrace, CacheCheck::Unchecked>(
                    node, le, cache)));
      if constexpr (trace(EvalTrace)) log::cache(node, cache, log::label(node));

      return mult_by_phase(ptr);
    } else {
      // do nothing
    }
  }

  ResultPtr result;
  ResultPtr left;
  ResultPtr right;

  log::Duration time;

  if (node.leaf()) {
    time = timed_eval_inplace([&]() { result = le(node); });
  } else {
    left = evaluate<EvalTrace>(node.left(), le, cache);
    right = evaluate<EvalTrace>(node.right(), le, cache);
    SEQUANT_ASSERT(left);
    SEQUANT_ASSERT(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};
    if (node->op_type() == EvalOp::Sum) {
      time = timed_eval_inplace([&]() { result = left->sum(*right, ann); });
    } else {
      SEQUANT_ASSERT(node->op_type() == EvalOp::Product);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();
      time = timed_eval_inplace([&]() {
        result =
            left->prod(*right, ann, de_nest ? DeNest::True : DeNest::False);
      });
    }
  }

  SEQUANT_ASSERT(result);

  // logging
  if constexpr (trace(EvalTrace)) {
    if (node.leaf()) {
      log::eval(log::EvalStat{.mode = log::eval_mode(node),
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = log::bytes(cache, result)},
                log::label(node));
    } else {
      // A cached child is *distinct* from the local left/right when its
      // canon_phase != 1, because mult_by_phase allocates a fresh buffer
      // while the cache still holds the pre-phase data. So only skip the
      // local's bytes when the cache aliases the same buffer (phase == 1).
      size_t hwmark = log::bytes(cache, result).value;
      if (!cache.alive(node.left()) || node.left()->canon_phase() != 1)
        hwmark += log::bytes(left).value;
      if (!cache.alive(node.right()) || node.right()->canon_phase() != 1)
        hwmark += log::bytes(right).value;
      log::eval(log::EvalStat{.mode = log::eval_mode(node),
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = {hwmark},
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
  if constexpr (trace(EvalTrace)) {
    xpr = toUtf8(io::serialization::to_string(to_expr(node)));
    log::term(log::TermMode::Begin, xpr);
  }

  struct {
    ResultPtr pre, post;
  } result;

  result.pre = evaluate<EvalTrace>(node, le, cache);

  auto time = timed_eval_inplace([&]() {
    result.post = perm ? result.pre->permute(
                             std::array<std::any, 2>{node->annot(), layout})
                       : result.pre;
  });

  SEQUANT_ASSERT(result.post);

  // logging
  if constexpr (trace(EvalTrace)) {
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
                                .mem_hwmark = {hwmark}};
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
    auto time = timed_eval_inplace([&]() { result->add_inplace(*pre); });

    // logging
    if constexpr (trace(EvalTrace)) {
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
                                .mem_hwmark = {hwmark}};
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
  requires(!last_type_is_cache_manager<Args...>)
ResultPtr evaluate(Args&&... args) {
  using Node =
      std::remove_cvref_t<decltype(node0(arg0(std::forward<Args>(args)...)))>;
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
  auto time = timed_eval_inplace([&]() { result = pre->symmetrize(); });

  // logging
  if constexpr (trace(EvalTrace)) {
    // cache is owned by the inner evaluate call and out of scope here;
    // hwmark reflects only the local working set (pre + freshly allocated
    // result both live during the symmetrize op).
    auto stat = log::EvalStat{.mode = log::EvalMode::Symmetrize,
                              .time = time,
                              .mem_result = log::bytes(result),
                              .mem_alloc = log::bytes(result),
                              .mem_hwmark = log::bytes(pre, result)};
    log::eval(stat, node0(arg0(std::forward<Args>(args)...))->label());
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

  auto const& n0 = node0(arg0(std::forward<Args>(args)...));

  ResultPtr result;
  auto time = timed_eval_inplace(
      [&]() { result = pre->antisymmetrize(n0->as_tensor().bra_rank()); });

  // logging
  if constexpr (trace(EvalTrace)) {
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

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
