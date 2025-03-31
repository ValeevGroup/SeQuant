#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
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

namespace {

using Seconds = std::chrono::duration<double>;

///
/// Invokes @c fun that returns void on the arguments @c args and returns the
/// time duration as @c std::chrono::duration<double>.
template <
    typename F, typename... Args,
    std::enable_if_t<std::is_invocable_r_v<void, F, Args...>, bool> = true>
auto timed_eval_inplace(F&& fun, Args&&... args) {
  using Clock = std::chrono::high_resolution_clock;
  auto tstart = Clock::now();
  std::forward<F>(fun)(std::forward<Args>(args)...);
  auto tend = Clock::now();
  return Seconds(tend - tstart);
}

template <typename... Args>
void log_tag(std::string_view tag, Args const&... args) {
  auto& l = Logger::instance();
  if constexpr (sizeof...(Args))
    if (l.eval.level > 0)
      write_log(l, std::format("[{}]", tag), std::format(" | {}", args)...,
                '\n');
}

auto log_eval = [](auto const&... args) { log_tag("EVAL", args...); };
auto log_cache = [](auto const&... args) { log_tag("CACHE", args...); };
auto log_term = [](auto const&... args) { log_tag("TERM", args...); };

void log_cache_access(size_t key, CacheManager const& cm) {
  auto const cur_l = cm.life(key);
  auto const max_l = cm.max_life(key);
  bool const release = cur_l == 0;
  bool const store = cur_l + 1 == max_l;
  log_cache(release ? "RELEASE"
            : store ? "STORE"
                    : "ACCESS",
            key,                                 //
            std::format("{}/{}", cur_l, max_l),  //
            cm.alive_count(),                    //
            std::format("{}B", cm.size_in_bytes()));
}

[[maybe_unused]] std::string perm_groups_string(
    container::svector<std::array<size_t, 3>> const& perm_groups) {
  std::string result;
  for (auto const& g : perm_groups)
    result += "(" + std::to_string(g[0]) + "," + std::to_string(g[1]) + "," +
              std::to_string(g[2]) + ") ";
  result.pop_back();  // remove last space
  return result;
}

template <typename... Args>
concept last_type_is_cache_manager =
    std::same_as<CacheManager, std::remove_cvref_t<std::tuple_element_t<
                                   sizeof...(Args) - 1, std::tuple<Args...>>>>;

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

    auto mult_by_phase = [phase = node->canon_phase()](ResultPtr res) {
      return phase == 1 ? res : res->mult_by_phase(phase);
    };

    auto const h = hash::value(*node);
    if (auto ptr = cache.access(h); ptr) {
      if constexpr (trace(EvalTrace)) log_cache_access(h, cache);

      return mult_by_phase(ptr);
    } else if (cache.exists(h)) {
      auto ptr = cache.store(
          h, mult_by_phase(
                 evaluate<EvalTrace, CacheCheck::Unchecked>(node, le, cache)));
      if constexpr (trace(EvalTrace)) log_cache_access(h, cache);

      return mult_by_phase(ptr);
    } else {
      // do nothing
    }
  }

  ResultPtr result;
  ResultPtr left;
  ResultPtr right;

  Seconds seconds;

  if (node.leaf()) {
    seconds = timed_eval_inplace([&]() { result = le(node); });
  } else {
    left = evaluate<EvalTrace>(node.left(), le, cache);
    right = evaluate<EvalTrace>(node.right(), le, cache);
    assert(left);
    assert(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};
    if (node->op_type() == EvalOp::Sum) {
      timed_eval_inplace([&]() { result = left->sum(*right, ann); });
    } else {
      assert(node->op_type() == EvalOp::Product);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();
      timed_eval_inplace([&]() {
        result = left->prod(*right, ann,
                            de_nest ? TA::DeNest::True : TA::DeNest::False);
      });
    }
  }

  assert(result);

  // logging
  if constexpr (trace(EvalTrace)) {
    struct {
      std::string type, annot;
      size_t bytes;
    } log;
    if (node.leaf()) {
      log.type = node->is_constant()   ? "CONSTANT"
                 : node->is_variable() ? "VARIABLE"
                                       : "TENSOR";
      log.annot = node->label();
      log.bytes = result->size_in_bytes();
    } else {
      log.type = node->is_product() ? "PRODUCT" : node->is_sum() ? "SUM" : "ID";
      log.annot = node->is_primary()
                      ? node->label()
                      : std::format("{} {} {} -> {}",                  //
                                    node.left()->label(),              //
                                    (node->is_product() ? "*" : "+"),  //
                                    node.right()->label(), node->label());
      log.bytes = left->size_in_bytes()     //
                  + right->size_in_bytes()  //
                  + result->size_in_bytes();
    }

    log_eval(log.type,                       //
             std::format("{}", seconds),     //
             std::format("{}B", log.bytes),  //
             log.annot);
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
  std::string xpr;
  if constexpr (trace(EvalTrace)) {
    xpr = to_string(deparse(to_expr(node)));
    log_term("BEGIN", xpr);
  }

  struct {
    ResultPtr pre, post;
  } result;

  result.pre = evaluate<EvalTrace>(node, le, cache);

  auto seconds = timed_eval_inplace([&]() {
    result.post =
        result.pre->permute(std::array<std::any, 2>{node->annot(), layout});
  });

  assert(result.post);

  // logging
  if constexpr (trace(EvalTrace)) {
    auto bytes = result.pre->size_in_bytes() + result.post->size_in_bytes();
    log_eval("PERMUTE",                   //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node->label());

    log_term("END", xpr);
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
    auto seconds = timed_eval_inplace([&]() { result->add_inplace(*pre); });

    // logging
    if constexpr (trace(EvalTrace)) {
      auto bytes = result->size_in_bytes() + pre->size_in_bytes();
      log_eval("ADD_INPLACE",               //
               std::format("{}", seconds),  //
               std::format("{}B", bytes),   //
               n->label());
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
/// \brief Calls @code evaluate followd by the particle-symmetrization function.
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
  auto seconds = timed_eval_inplace([&]() { result = pre->symmetrize(); });

  // logging
  if constexpr (trace(EvalTrace)) {
    auto&& arg0 =
        std::get<0>(std::forward_as_tuple(std::forward<Args>(args)...));
    std::string node_label;
    if constexpr (meta::can_evaluate_range<decltype(arg0)>)
      node_label = ranges::front(arg0)->label();
    else
      node_label = arg0->label();

    size_t bytes = pre->size_in_bytes() + result->size_in_bytes();
    log_eval("SYMMETRIZE",                //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node_label);
  }

  return result;
}

///
/// \tparam EvalTrace If Trace::On, trace is written to the logger's stream.
///                   Default is to follow Trace::Default, which is itself
///                   equal to Trace::On or Trace::Off.
/// \brief Calls @code evaluate followd by the anti-symmetrization function on
///        the bra indices and the ket indices. The bra and ket indices are
///        inferred from the evaluation node(s).
/// \return Evaluated result as ResultPtr.
///
template <Trace EvalTrace = Trace::Default, typename... Args>
ResultPtr evaluate_antisymm(Args&&... args) {
  size_t bra_rank;
  std::string node_label;  // for logging
  auto&& arg0 = std::get<0>(std::forward_as_tuple(std::forward<Args>(args)...));
  if constexpr (meta::can_evaluate_range<decltype(arg0)>) {
    assert(!ranges::empty(arg0));
    bra_rank = ranges::front(arg0)->as_tensor().bra_rank();
    node_label = ranges::front(arg0)->label();
  } else {
    bra_rank = arg0->as_tensor().bra_rank();
    node_label = arg0->label();
  }

  ResultPtr pre = evaluate<EvalTrace>(std::forward<Args>(args)...);
  assert(pre);

  ResultPtr result;
  auto seconds =
      timed_eval_inplace([&]() { result = pre->antisymmetrize(bra_rank); });

  // logging
  if constexpr (trace(EvalTrace)) {
    size_t bytes = pre->size_in_bytes() + result->size_in_bytes();
    log_eval("ANTISYMMETRIZE",            //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node_label);
  }
  return result;
}
}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
