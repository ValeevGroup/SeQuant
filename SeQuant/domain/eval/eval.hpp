#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include "eval_result.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>

#include <btas/btas.h>
#include <tiledarray.h>

#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <any>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace sequant {
#if __cplusplus < 202002L
template <class T>
struct remove_cvref {
  typedef std::remove_cv_t<::std::remove_reference_t<T>> type;
};

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;
#else
template <typename T>
using remove_cvref = std::remove_cvref<T>;

template <typename T>
using remove_cvref_t = std::remove_cvref_t<T>;
#endif
}  // namespace sequant

namespace sequant::eval {

template <typename T, typename = void>
constexpr bool IsIterable{};

template <typename T>
constexpr bool IsIterable<
    T, std::void_t<
           decltype(std::begin(std::declval<std::remove_reference_t<T>>())),
           decltype(std::end(std::declval<std::remove_reference_t<T>>()))>> =
    true;

template <typename I, typename = std::enable_if_t<IsIterable<I>>>
using IteredT =
    std::remove_reference_t<decltype(*std::begin(std::declval<I>()))>;

template <typename, typename = void>
constexpr bool HasAnnotMethod{};

template <typename T>
constexpr bool HasAnnotMethod<
    T, std::void_t<decltype(std::declval<remove_cvref_t<T>>().annot())>> = true;

template <typename, typename = void>
constexpr bool IsEvaluable{};

template <typename T>
constexpr bool IsEvaluable<
    FullBinaryNode<T>,
    std::enable_if_t<std::is_convertible_v<T, EvalExpr> && HasAnnotMethod<T>>> =
    true;

template <typename T>
constexpr bool IsEvaluable<
    const FullBinaryNode<T>,
    std::enable_if_t<std::is_convertible_v<T, EvalExpr> && HasAnnotMethod<T>>> =
    true;

template <typename, typename = void>
constexpr bool IsIterableOfEvaluableNodes{};

template <typename Iterable>
constexpr bool IsIterableOfEvaluableNodes<
    Iterable, std::enable_if_t<IsEvaluable<IteredT<Iterable>>>> = true;

template <typename, typename, typename = void>
constexpr bool IsLeafEvaluator{};

template <typename NodeT>
constexpr bool IsLeafEvaluator<NodeT, CacheManager<ERPtr>, void>{};

template <typename NodeT, typename Le>
constexpr bool IsLeafEvaluator<
    NodeT, Le,
    std::enable_if_t<
        IsEvaluable<NodeT> &&
        std::is_same_v<
            ERPtr, std::remove_reference_t<std::invoke_result_t<Le, NodeT>>>>> =
    true;

template <typename NodesI,
          typename Pred = std::function<bool(IteredT<NodesI> const&)>,
          typename = std::enable_if_t<
              IsIterableOfEvaluableNodes<NodesI> &&
              std::is_invocable_r_v<bool, Pred, IteredT<NodesI> const&>>>
CacheManager<ERPtr> cache_manager(
    NodesI const& nodes, Pred&& pred = [](auto&&) { return true; },
    size_t min_repeats = 2) noexcept {
  auto imed_counts = container::map<size_t, size_t>{};

  // counts number of times each internal node appears in
  // all of @c nodes trees
  auto imed_visitor = [&imed_counts, &pred](auto&& n) {
    if (!std::invoke(std::forward<Pred>(pred), n)) return;

    auto&& end = imed_counts.end();
    auto&& h = n->hash_value();
    if (auto&& found = imed_counts.find(h); found != end)
      ++found->second;
    else
      imed_counts.emplace(h, 1);
  };  // imed_visitor

  // visit imeds
  ranges::for_each(nodes, [&imed_visitor](auto&& tree) {
    tree.visit_internal(imed_visitor);
  });

  // remove less repeating imeds
  auto less_repeating = [min_repeats](auto&& pair) {
    return pair.second < min_repeats;
  };
  ranges::actions::remove_if(imed_counts, less_repeating);

  return CacheManager<ERPtr>{imed_counts};
}

class EvalExprTA final : public EvalExpr {
 public:
  ///
  /// annotation for TiledArray
  ///
  [[nodiscard]] std::string const& annot() const;

  explicit EvalExprTA(Tensor const&);

  explicit EvalExprTA(Constant const&);

  EvalExprTA(EvalExprTA const&, EvalExprTA const&, EvalOp);

 private:
  std::string annot_;
};

class EvalExprBTAS final : public EvalExpr {
 public:
  using annot_t = container::svector<long>;

  ///
  /// annotation for BTAS tensor
  ///
  [[nodiscard]] annot_t const& annot() const noexcept;

  explicit EvalExprBTAS(Tensor const&) noexcept;

  explicit EvalExprBTAS(Constant const&) noexcept;

  EvalExprBTAS(EvalExprBTAS const&, EvalExprBTAS const&, EvalOp) noexcept;

 private:
  annot_t annot_;
};

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le&&);

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le&&, CacheManager<ERPtr>&);

template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_core(NodeT const& node, Le&& le, Args&&... args) {
  if (node.leaf()) {
    return std::invoke(std::forward<Le>(le), node);
  } else {
    ERPtr const left = evaluate_crust(node.left(), std::forward<Le>(le),
                                      std::forward<Args>(args)...);
    ERPtr const right = evaluate_crust(node.right(), std::forward<Le>(le),
                                       std::forward<Args>(args)...);

    assert(left);
    assert(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};

    if (node->op_type() == EvalOp::Sum) {
      return left->sum(*right, ann);
    } else {
      assert(node->op_type() == EvalOp::Prod);
      return left->prod(*right, ann);
    }
  }
}

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool>>
ERPtr evaluate_crust(NodeT const& node, Le&& le) {
  return evaluate_core(node, std::forward<Le>(le));
}

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool>>
ERPtr evaluate_crust(NodeT const& node, Le&& le, CacheManager<ERPtr>& cache) {
  auto const h = hash::value(*node);
  if (auto ptr = cache.access(h); ptr) {
    return *ptr;
  } else if (cache.exists(h)) {
    return *cache.store(h, evaluate_core(node, std::forward<Le>(le), cache));
  } else {
    return evaluate_core(node, std::forward<Le>(le), cache);
  }
}

template <typename NodeT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node,    //
              Annot const& layout,  //
              Le&& le, Args&&... args) {
  return evaluate_crust(node, std::forward<Le>(le), std::forward<Args>(args)...)
      ->permute(std::array<std::any, 2>{node->annot(), layout});
}

template <typename NodesT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>,
                           bool> = true>
auto evaluate(NodesT const& nodes,  //
              Annot const& layout,  //
              Le&& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter, layout, std::forward<Le>(le),
                         std::forward<Args>(args)...);
  for (++iter; iter != end; ++iter) {
    result->add_inplace(*evaluate(*iter, layout, std::forward<Le>(le),
                                  std::forward<Args>(args)...));
  }
  return result;
}

template <typename NodeT, typename Annot, typename Le, typename... Args>
auto evaluate_symm(NodeT const& node, Annot const& layout,
                   container::svector<std::array<size_t, 3>> const& perm_groups,
                   Le&& le, Args&&... args) {
  if (perm_groups.empty()) {
    // asked for symmetrization without specifying particle symmetric index
    // ranges

    ExprPtr expr_ptr{};
    if constexpr (IsIterableOfEvaluableNodes<NodeT>) {
      expr_ptr = (*std::begin(node))->expr();
    } else {
      expr_ptr = node->expr();
    }
    assert(expr_ptr->is<Tensor>());
    auto const& t = expr_ptr->as<Tensor>();
    assert(t.bra_rank() == t.ket_rank());

    size_t const half_rank = t.bra_rank();
    return evaluate(node, layout, std::forward<Le>(le),
                    std::forward<Args>(args)...)
        ->symmetrize({{0, half_rank, half_rank}});
  }

  return evaluate(node, layout, std::forward<Le>(le),
                  std::forward<Args>(args)...)
      ->symmetrize(perm_groups);
}

template <typename NodeT, typename Annot, typename Le,
          typename... Args>
auto evaluate_antisymm(
    NodeT const& node,                                             //
    Annot const& layout,                                           //
    container::svector<std::array<size_t, 2>> const& perm_groups,  //
    Le&& le,                                                       //
    Args&&... args) {
  if (perm_groups.empty()) {
    // Asked for antisymmetrization without specifying particle antisymmetric
    // index ranges.
    // Assume both bra indices and ket indices are antisymmetric in
    // the particle exchange.
    ExprPtr expr_ptr{};
    if constexpr (IsIterableOfEvaluableNodes<NodeT>) {
      expr_ptr = (*std::begin(node))->expr();
    } else {
      expr_ptr = node->expr();
    }
    assert(expr_ptr->is<Tensor>());
    auto const& t = expr_ptr->as<Tensor>();

    size_t const b = t.bra_rank();
    size_t const k = t.ket_rank();
    return evaluate(node, layout, std::forward<Le>(le),
                    std::forward<Args>(args)...)
        ->antisymmetrize({{0, b}, {b, k}});
  }
  return evaluate(node, layout, std::forward<Le>(le),
                  std::forward<Args>(args)...)
      ->antisymmetrize(perm_groups);
}

template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node, Le&& le, Args&&... args) {
  return evaluate_crust(node, std::forward<Le>(le),
                        std::forward<Args>(args)...);
}

template <typename NodesT, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>, bool> = true>
auto evaluate(NodesT const& nodes, Le&& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter,                 //
                         std::forward<Le>(le),  //
                         std::forward<Args>(args)...);
  for (++iter; iter != end; ++iter) {
    result->add_inplace(*evaluate(*iter,                 //
                                  std::forward<Le>(le),  //
                                  std::forward<Args>(args)...));
  }
  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
