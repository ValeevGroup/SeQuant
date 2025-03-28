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
#include <SeQuant/domain/eval/eval_result.hpp>

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
void log_eval(Args const&... args) noexcept {
  auto& l = Logger::instance();
  if constexpr (sizeof...(Args))
    if (l.eval.level > 0)
      write_log(l, "[EVAL]", std::format(" | {}", args)..., '\n');
}

template <typename... Args>
void log_term(Args const&... args) noexcept {
  auto& l = Logger::instance();
  if constexpr (sizeof...(Args))
    if (l.eval.level > 0) write_log(l, "[TERM]", args..., '\n');
}

[[maybe_unused]] void log_cache_access(size_t key, CacheManager const& cm) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) {
    assert(cm.exists(key));
    auto max_l = cm.max_life(key);
    auto cur_l = cm.life(key);
    write_log(l,                                    //
              "[CACHE] Accessed key: ", key, ". ",  //
              cur_l, "/", max_l, " lives remain.\n");
    if (cur_l == 0) {
      write_log(l,  //
                "[CACHE] Released key: ", key, ".\n");
    }
  }
}

[[maybe_unused]] void log_cache_store(size_t key, CacheManager const& cm) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) {
    assert(cm.exists(key));
    write_log(l,  //
              "[CACHE] Stored key: ", key, ".\n");
    // because storing automatically implies immediately accessing it
    log_cache_access(key, cm);
  }
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

template <CacheCheck Cache = CacheCheck::Checked, meta::can_evaluate Node,
          typename F>
  requires meta::leaf_node_evaluator<Node, F>
ERPtr evaluate(Node const& node,  //
               F const& le,       //
               CacheManager& cache) {
  if constexpr (Cache == CacheCheck::Checked) {  // return from cache if found

    auto mult_by_phase = [phase = node->canon_phase()](ERPtr res) {
      return phase == 1 ? res : res->mult_by_phase(phase);
    };

    auto const h = hash::value(*node);
    if (auto ptr = cache.access(h); ptr) {
      log_cache_access(h, cache);
      return mult_by_phase(ptr);
    } else if (cache.exists(h)) {
      auto ptr = cache.store(
          h, mult_by_phase(evaluate<CacheCheck::Unchecked>(node, le, cache)));
      log_cache_store(h, cache);
      return mult_by_phase(ptr);
    } else {
      // do nothing
    }
  }

  ERPtr result;

  Seconds seconds;
  size_t bytes;

  if (node.leaf()) {
    seconds = timed_eval_inplace([&]() { result = le(node); });
  } else {
    ERPtr const left = evaluate(node.left(), le, cache);
    ERPtr const right = evaluate(node.right(), le, cache);

    assert(left);
    assert(right);

    bytes += left->size_in_bytes();
    bytes += right->size_in_bytes();

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};
    if (node->op_type() == EvalOp::Sum) {
      timed_eval_inplace([&]() { result = left->sum(*right, ann); });
    } else {
      assert(node->op_type() == EvalOp::Prod);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();
      timed_eval_inplace([&]() {
        result = left->prod(*right, ann,
                            de_nest ? TA::DeNest::True : TA::DeNest::False);
      });
    }
  }

  assert(result);

  bytes += result->size_in_bytes();

  {  // logging
    struct {
      std::string type, annot;
    } log;
    if (node.leaf()) {
      log.type = node->is_constant()   ? "CONSTANT"
                 : node->is_variable() ? "VARIABLE"
                                       : "TENSOR";
      log.annot = node->label();
    } else {
      log.type = node->is_prod() ? "PROD" : node->is_sum() ? "SUM" : "ID";
      log.annot = node->is_id()
                      ? node->label()
                      : std::format("{} {} {} -> {}",               //
                                    node.left()->label(),           //
                                    (node->is_prod() ? "*" : "+"),  //
                                    node.right()->label(), node->label());
    }

    log_eval(log.type,                    //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             log.annot);
  }

  return result;
}

template <meta::can_evaluate Node, typename F>
  requires meta::leaf_node_evaluator<Node, F>  //
ERPtr evaluate(Node const& node,               //
               auto const& layout,             //
               F const& le,                    //
               CacheManager& cache) {
  log_term("[TERM]", " ", to_string(deparse(to_expr(node))), '\n');
  struct {
    ERPtr pre, post;
  } result;
  result.pre = evaluate(node, le, cache);

  auto seconds = timed_eval_inplace([&]() {
    result.post =
        result.pre->permute(std::array<std::any, 2>{node->annot(), layout});
  });
  assert(result.post);

  {  // logging
    auto bytes = result.pre->size_in_bytes() + result.post->size_in_bytes();
    log_eval("PERMUTE",                   //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node->label());
  }

  return result.post;
}

template <meta::can_evaluate_range Nodes, typename F>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ERPtr evaluate(Nodes const& nodes,  //
               auto const& layout,  //
               F const& le, CacheManager& cache) {
  ERPtr result;

  for (auto&& n : nodes) {
    if (!result) {
      result = evaluate(n, layout, le, cache);
      continue;
    }

    size_t bytes;
    Seconds seconds;
    ERPtr pre = evaluate(n, layout, le, cache);
    seconds = timed_eval_inplace([&]() { result->add_inplace(*pre); });
    bytes = result->size_in_bytes() + pre->size_in_bytes();

    // logging
    {
      log_eval("ADD_INPLACE",               //
               std::format("{}", seconds),  //
               std::format("{}B", bytes),   //
               n->label());
    }
  }

  return result;
}

template <typename... Args>
  requires(!last_type_is_cache_manager<Args...>)
ERPtr evaluate(Args&&... args) {
  auto cache = CacheManager::empty();
  return evaluate(std::forward<Args>(args)..., cache);
}

template <typename... Args>
[[deprecated]] ERPtr evaluate_symm(Args&&... args) {
  ERPtr pre = evaluate(std::forward<Args>(args)...);
  assert(pre);
  ERPtr result;
  auto seconds = timed_eval_inplace([&]() { result = pre->symmetrize(); });
  size_t bytes = pre->size_in_bytes() + result->size_in_bytes();

  // logging
  {
    auto&& arg0 =
        std::get<0>(std::forward_as_tuple(std::forward<Args>(args)...));
    std::string node_label;
    if constexpr (meta::can_evaluate_range<decltype(arg0)>)
      node_label = ranges::front(arg0)->label();
    else
      node_label = arg0->label();

    log_eval("SYMMETRIZE",                //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node_label);
  }

  return result;
}

template <typename... Args>
[[deprecated]] ERPtr evaluate_antisymm(Args&&... args) {
  size_t bra_rank;
  std::string node_label;  // for logging
  auto&& arg0 = std::get<0>(std::forward_as_tuple(std::forward<Args>(args)...));
  if constexpr (meta::can_evaluate_range<decltype(arg0)>) {
    assert(!ranges::empty(arg0));
    bra_rank = ranges::front(arg0)->as_tensor().bra_rank();
  } else {
    bra_rank = arg0->as_tensor().bra_rank();
  }

  ERPtr pre = evaluate(std::forward<Args>(args)...);
  assert(pre);

  ERPtr result;
  auto seconds =
      timed_eval_inplace([&]() { result = pre->antisymmetrize(bra_rank); });
  size_t bytes = pre->size_in_bytes() + result->size_in_bytes();

  // logging
  {
    log_eval("ANTISYMMETRIZE",            //
             std::format("{}", seconds),  //
             std::format("{}B", bytes),   //
             node_label);
  }
  return result;
}
}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
