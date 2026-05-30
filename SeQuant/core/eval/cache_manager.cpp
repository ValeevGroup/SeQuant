//
// Created by Bimal Gaudel on 9/26/23.
//

#include <SeQuant/core/eval/cache_manager.hpp>

#include <SeQuant/core/eval/eval_node.hpp>

#include <range/v3/view/transform.hpp>

namespace sequant {

namespace {

template <meta::eval_node Node, typename N, bool F>
void max_cache(Node const& node,        //
               CacheManager<N, F>& cm,  //
               AsyCost& curr,           //
               AsyCost& max) {
  if (auto ptr = cm.access(node); ptr) {
    if (cm.life(node) == 0) {
      curr -= Memory{}(node);
    }
    return;
  }
  if (!node.leaf()) {
    max_cache(node.left(), cm, curr, max);
    max_cache(node.right(), cm, curr, max);
    if (cm.exists(node)) {
      curr += Memory{}(node);
      max = std::max(curr, max);
      // simulate cache store
      auto s = cm.store(node, nullptr);
    }
  }
}

}  // namespace

AsyCost peak_cache(Sum const& expr, std::optional<size_t> min_repeats) {
  // Materialize into a vector so that nodes have stable addresses for
  // pointer-based scanning in cache_manager(). per-summand binarize for
  // cache-cost analysis only; the positional head doesn't escape.
#if defined(__clang__) || defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
  auto const nodes_vec = expr | ranges::views::transform([](const auto& expr) {
                           return binarize<EvalExpr>(expr);
                         }) |
                         ranges::to_vector;
#if defined(__clang__) || defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  auto cm = min_repeats ? cache_manager(nodes_vec, min_repeats.value())
                        : cache_manager(nodes_vec);
  auto max = AsyCost::zero();
  auto curr = AsyCost::zero();
  for (auto const& n : nodes_vec) max_cache(n, cm, curr, max);
  return max;
}

}  // namespace sequant
