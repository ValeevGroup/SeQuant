#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval_fwd.hpp>
#include <SeQuant/domain/eval/result.hpp>

#include <btas/btas.h>
#include <tiledarray.h>

#include <chrono>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <any>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace sequant {

namespace log {

using Duration = std::chrono::nanoseconds;

struct Bytes {
  size_t value;
};

template <typename... T>
  requires((std::same_as<ResultPtr, T> && ...))
[[nodiscard]] auto bytes(T const&... args) {
  return Bytes{(args->size_in_bytes() + ...)};
}

[[nodiscard]] inline auto bytes(CacheManager const& cman) {
  return cman.size_in_bytes();
}

[[nodiscard]] inline auto to_string(Bytes bs) noexcept {
  return std::format("{}B", bs.value);
}

enum struct EvalMode {
  Constant,
  Variable,
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

struct EvalStat {
  EvalMode mode;
  Duration time;
  Bytes memory;
};

struct CacheStat {
  CacheMode mode;
  size_t key;
  int curr_life, max_life;
  size_t num_alive;
  Bytes memory;
};

template <typename Arg, typename... Args>
void log(Arg const& arg, Args const&... args) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) write_log(l, arg, std::format(" | {}", args)..., '\n');
}

template <typename... Args>
auto eval(EvalStat const& stat, Args const&... args) {
  log("Eval",                  //
      to_string(stat.mode),    //
      stat.time,               //
      to_string(stat.memory),  //
      args...);
}

template <typename... Args>
auto cache(CacheStat const& stat, Args const&... args) {
  log("Cache",                                              //
      to_string(stat.mode),                                 //
      stat.key,                                             //
      std::format("{}/{}", stat.curr_life, stat.max_life),  //
      stat.num_alive,                                       //
      to_string(stat.memory),                               //
      args...);
}

template <typename... Args>
auto cache(size_t key, CacheManager const& cm, Args const&... args) {
  using CacheMode::Access;
  using CacheMode::Release;
  using CacheMode::Store;
  auto const cur_l = cm.life(key);
  auto const max_l = cm.max_life(key);
  bool const release = cur_l == 0;
  bool const store = cur_l + 1 == max_l;
  cache(CacheStat{.mode = store     ? Store
                          : release ? Release
                                    : Access,
                  .key = key,
                  .curr_life = cur_l,
                  .max_life = max_l,
                  .num_alive = cm.alive_count(),
                  .memory = {bytes(cm)}},
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

template <typename... Args>
concept last_type_is_cache_manager =
    std::same_as<CacheManager, std::remove_cvref_t<std::tuple_element_t<
                                   sizeof...(Args) - 1, std::tuple<Args...>>>>;

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
/// \tparam Cache If CacheCache::Checked (default) the @param cache will be
///               checked before evaluating. It is used to detect the base case
///               for recursion to prevent infinite recursion.
/// \param node A node that can be evaluated using @param le as the leaf
///             evaluator.
/// \param le The leaf evaluator that satisfies
///           @code meta::leaf_node_evaluator<Node, F>.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default,
          CacheCheck Cache = CacheCheck::Checked, meta::can_evaluate Node,
          typename F>
  requires meta::leaf_node_evaluator<Node, F>
ResultPtr evaluate(Node const& node,  //
                   F const& le,       //
                   CacheManager& cache) {
  if constexpr (Cache == CacheCheck::Checked) {  // return from cache if found

    auto mult_by_phase = [&node](ResultPtr res) {
      auto phase = node->canon_phase();
      if (phase == 1) return res;

      ResultPtr post;
      auto time =
          timed_eval_inplace([&]() { post = res->mult_by_phase(phase); });

      if constexpr (trace(EvalTrace)) {
        auto stat = log::EvalStat{.mode = log::EvalMode::MultByPhase,
                                  .time = time,
                                  .memory = log::bytes(res, post)};
        log::eval(stat, std::format("{} * {}", phase, node->label()));
      }
      return post;
    };

    auto const h = hash::value(*node);
    if (auto ptr = cache.access(h); ptr) {
      if constexpr (trace(EvalTrace)) log::cache(h, cache);

      return mult_by_phase(ptr);
    } else if (cache.exists(h)) {
      auto ptr = cache.store(
          h, mult_by_phase(
                 evaluate<EvalTrace, CacheCheck::Unchecked>(node, le, cache)));
      if constexpr (trace(EvalTrace)) log::cache(h, cache);

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
    assert(left);
    assert(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};
    if (node->op_type() == EvalOp::Sum) {
      time = timed_eval_inplace([&]() { result = left->sum(*right, ann); });
    } else {
      assert(node->op_type() == EvalOp::Product);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();
      time = timed_eval_inplace([&]() {
        result = left->prod(*right, ann,
                            de_nest ? TA::DeNest::True : TA::DeNest::False);
      });
    }
  }

  assert(result);

  // logging
  if constexpr (trace(EvalTrace)) {
    auto stat =
        log::EvalStat{.mode = log::eval_mode(node),
                      .time = time,
                      .memory = node.leaf() ? log::bytes(result)
                                            : log::bytes(left, right, result)};
    log::eval(stat, log::label(node));
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param node A node that can be evaluated using @param le as the leaf
///             evaluator.
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
/// \param le The leaf evaluator that satisfies
///           @code meta::leaf_node_evaluator<Node, F>.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate Node, typename F>
  requires meta::leaf_node_evaluator<Node, F>  //
ResultPtr evaluate(Node const& node,           //
                   auto const& layout,         //
                   F const& le,                //
                   CacheManager& cache) {
  // if the layout is not the default constructed value need to permute
  bool const perm = layout != decltype(layout){};

  std::string xpr;
  if constexpr (trace(EvalTrace)) {
    xpr = to_string(deparse(to_expr(node)));
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

  assert(result.post);

  // logging
  if constexpr (trace(EvalTrace)) {
    if (perm) {
      auto stat = log::EvalStat{.mode = log::EvalMode::Permute,
                                .time = time,
                                .memory = log::bytes(result.pre, result.post)};
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
/// \param nodes A range of node that can be evaluated using @param le as the
///              leaf evaluator. The evaluation result of the elements of
///              @param nodes will be summed up.
///
/// \param layout The layout of the final result. Only meaningful if the result
///               has a layout (or supports permutation) eg. a tensor.
///               The results of each element from @param nodes will be permuted
///               to this layout before being summed.
///
/// \param le The leaf evaluator that satisfies
///           @code meta::leaf_node_evaluator<Node, F>.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& nodes,  //
                   auto const& layout,  //
                   F const& le, CacheManager& cache) {
  ResultPtr result;

  for (auto&& n : nodes) {
    if (!result) {
      result = evaluate<EvalTrace>(n, layout, le, cache);
      continue;
    }

    ResultPtr pre = evaluate<EvalTrace>(n, layout, le, cache);
    auto time = timed_eval_inplace([&]() { result->add_inplace(*pre); });

    // logging
    if constexpr (trace(EvalTrace)) {
      auto stat = log::EvalStat{.mode = log::EvalMode::SumInplace,
                                .time = time,
                                .memory = log::bytes(result, pre)};
      log::eval(stat, n->label());
    }
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \param nodes A range of node that can be evaluated using @param le as the
///              leaf evaluator. The evaluation result of the elements of
///              @param nodes will be summed up.
///
/// \param le The leaf evaluator that satisfies
///           @code meta::leaf_node_evaluator<Node, F>.
/// \param cache The cache for common sub-expression elimination.
/// \return Evaluated result as ResultPtr.
/// \note Because this function takes no layout argument, it is only useful
///       to evaluate summations of the elements in the @param nodes when they
///       are scalar results.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& nodes,  //
                   F const& le, CacheManager& cache) {
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
///        manager. Calls the other @code evalaute function overloads.
/// \see evaluate.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
  requires(!last_type_is_cache_manager<Args...>)
ResultPtr evaluate(Args&&... args) {
  auto cache = CacheManager::empty();
  return evaluate<EvalTrace>(std::forward<Args>(args)..., cache);
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Calls @code evaluate followed by the particle-symmetrization
///        function.
///        The number of particles is inferred by the tensor present in the
///        evaluation node(s). Presence of odd-ranked tensors in the evaluation
///        node(s) is an error.
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
ResultPtr evaluate_symm(Args&&... args) {
  ResultPtr pre = evaluate<EvalTrace>(std::forward<Args>(args)...);
  assert(pre);
  ResultPtr result;
  auto time = timed_eval_inplace([&]() { result = pre->symmetrize(); });

  // logging
  if constexpr (trace(EvalTrace)) {
    auto stat = log::EvalStat{.mode = log::EvalMode::Symmetrize,
                              .time = time,
                              .memory = log::bytes(pre, result)};
    log::eval(stat, node0(arg0(std::forward<Args>(args)...))->label());
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Calls @code evaluate followed by the anti-symmetrization function on
///        the bra indices and the ket indices. The bra and ket indices are
///        inferred from the evaluation node(s).
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
ResultPtr evaluate_antisymm(Args&&... args) {
  ResultPtr pre = evaluate<EvalTrace>(std::forward<Args>(args)...);
  assert(pre);

  auto const& n0 = node0(arg0(std::forward<Args>(args)...));

  ResultPtr result;
  auto time = timed_eval_inplace(
      [&]() { result = pre->antisymmetrize(n0->as_tensor().bra_rank()); });

  // logging
  if constexpr (trace(EvalTrace)) {
    auto stat = log::EvalStat{.mode = log::EvalMode::Antisymmetrize,
                              .time = time,
                              .memory = log::bytes(pre, result)};
    log::eval(stat, n0->label());
  }
  return result;
}
}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
