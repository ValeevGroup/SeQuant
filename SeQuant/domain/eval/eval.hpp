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

namespace {

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

}  // namespace

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

///
/// \brief Make a cache manager from an iterable of evaluable nodes.
///
/// \param nodes An iterable of evaluable nodes.
///
/// \param pred A predicate to filter nodes. By default all nodes are treated to
///             take part in contributing to the number of repeats.
///
/// \param min_repeats Minimum number of repeats for a node to be cached. By
///                    default anything repeated twice or more is cached.
///
/// \return A cache manager.
///
/// \see CacheManager
///
template <typename NodesI,
          typename Pred = std::function<bool(IteredT<NodesI> const&)>,
          typename = std::enable_if_t<
              IsIterableOfEvaluableNodes<NodesI> &&
              std::is_invocable_r_v<bool, Pred, IteredT<NodesI> const&>>>
CacheManager<ERPtr> cache_manager(
    NodesI const& nodes, Pred const& pred = [](auto&&) { return true; },
    size_t min_repeats = 2) noexcept {
  auto imed_counts = container::map<size_t, size_t>{};

  // counts number of times each internal node appears in
  // all of @c nodes trees
  auto imed_visitor = [&imed_counts, &pred](auto&& n) {
    if (!pred(n)) return;

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

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
///        that it can be used to evaluate using TiledArray.
///
class EvalExprTA final : public EvalExpr {
 public:
  ///
  /// \return String annotation for TA::DistArray.
  ///
  [[nodiscard]] std::string const& annot() const;

  ///
  /// \brief Construct an EvalExprTA from a Tensor.
  ///
  /// \see EvalExpr(Tensor const&).
  ///
  explicit EvalExprTA(Tensor const&);

  ///
  /// \brief Construct an EvalExprTA from a Constant.
  ///
  /// \see EvalExpr(Constant const&).
  ///
  explicit EvalExprTA(Constant const&);

  ///
  /// \brief Construct an EvalExprTA from two EvalExprTA and an EvalOp.
  /// \see EvalExpr(EvalExpr const&, EvalExpr const&, EvalOp).
  ///
  EvalExprTA(EvalExprTA const&, EvalExprTA const&, EvalOp);

 private:
  std::string annot_;
};  // class EvalExprTA

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
///        that it can be used to evaluate using BTAS.
///
class EvalExprBTAS final : public EvalExpr {
 public:
  using annot_t = container::svector<long>;

  ///
  /// \return Annotation (container::svector<long>) for BTAS::Tensor.
  ///
  [[nodiscard]] annot_t const& annot() const noexcept;

  ///
  /// \brief Construct an EvalExprBTAS from a Tensor.
  ///
  /// \see EvalExpr(Tensor const&).
  ///
  explicit EvalExprBTAS(Tensor const&) noexcept;

  ///
  /// \brief Construct an EvalExprBTAS from a Constant.
  ///
  /// \see EvalExpr(Constant const&).
  ///
  explicit EvalExprBTAS(Constant const&) noexcept;

  ///
  /// \brief Construct an EvalExprBTAS from two EvalExprBTAS and an EvalOp.
  ///
  /// \see EvalExpr(EvalExpr const&, EvalExpr const&, EvalOp).
  ///
  EvalExprBTAS(EvalExprBTAS const&, EvalExprBTAS const&, EvalOp) noexcept;

 private:
  annot_t annot_;
};  // EvalExprBTAS

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le const&);

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le const&, CacheManager<ERPtr>&);

template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_core(NodeT const& node, Le const& le, Args&&... args) {
  if (node.leaf()) {
    return le(node);
  } else {
    ERPtr const left =
        evaluate_crust(node.left(), le, std::forward<Args>(args)...);
    ERPtr const right =
        evaluate_crust(node.right(), le, std::forward<Args>(args)...);

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
ERPtr evaluate_crust(NodeT const& node, Le const& le) {
  return evaluate_core(node, le);
}

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool>>
ERPtr evaluate_crust(NodeT const& node, Le const& le,
                     CacheManager<ERPtr>& cache) {
  auto const h = hash::value(*node);
  if (auto ptr = cache.access(h); ptr) {
    return *ptr;
  } else if (cache.exists(h)) {
    return *cache.store(h, evaluate_core(node, le, cache));
  } else {
    return evaluate_core(node, le, cache);
  }
}

///
/// \param node An EvalNode to be evaluated.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting EvalResult.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node, Le&& le, Args&&... args) {
  return evaluate_crust(node, le, std::forward<Args>(args)...);
}

///
/// \param nodes An iterable of EvalNode objects that will be evaluated turn by
///              turn and summed up.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting EvalResult.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodesT, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>, bool> = true>
auto evaluate(NodesT const& nodes, Le const& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter, le, std::forward<Args>(args)...);

  for (++iter; iter != end; ++iter) {
    result->add_inplace(*evaluate(*iter,  //
                                  le,     //
                                  std::forward<Args>(args)...));
  }
  return result;
}

///
/// \param node An EvalNode to be evaluated into a tensor.
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node,    //
              Annot const& layout,  //
              Le const& le, Args&&... args) {
  return evaluate_crust(node, le, std::forward<Args>(args)...)
      ->permute(std::array<std::any, 2>{node->annot(), layout});
}

///
/// \param nodes An iterable of EvalNode objects that will be evaluated turn by
///              turn and summed up into a tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodesT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>,
                           bool> = true>
auto evaluate(NodesT const& nodes,  //
              Annot const& layout,  //
              Le const& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter, layout, le, std::forward<Args>(args)...);

  for (++iter; iter != end; ++iter) {
    result->add_inplace(
        *evaluate(*iter, layout, le, std::forward<Args>(args)...));
  }
  return result;
}

///
/// \param node An EvalNode or an iterable of such nodes to be evaluated into a
///             tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param perm_groups A vector of 3-element arrays of size_t. Each array
///                    represents a group of indices that are particle
///                    symmetric. The first two elements of the array are the
///                    indices of the bra and ket of the resulting tensor,
///                    respectively, and the third element is the number of
///                    symmetric indices in the group.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le, typename... Args>
auto evaluate_symm(NodeT const& node, Annot const& layout,
                   container::svector<std::array<size_t, 3>> const& perm_groups,
                   Le const& le, Args&&... args) {
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
    return evaluate(node, layout, le, std::forward<Args>(args)...)
        ->symmetrize({{0, half_rank, half_rank}});
  }

  return evaluate(node, layout, le, std::forward<Args>(args)...)
      ->symmetrize(perm_groups);
}

///
/// \param node An EvalNode or an iterable of such nodes to be evaluated into a
///             tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param perm_groups A vector of 3-element arrays of size_t. Each array
///                    represents a group of indices that are particle
///                    anti-symmetric. The first two elements of the array are
///                    the indices of the bra and ket of the resulting tensor,
///                    respectively, and the third element is the number of
///                    symmetric indices in the group.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le,
          typename... Args>
auto evaluate_antisymm(
    NodeT const& node,                                             //
    Annot const& layout,                                           //
    container::svector<std::array<size_t, 3>> const& perm_groups,  //
    Le const& le,                                                  //
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
    assert(b == t.ket_rank());
    return evaluate(node, layout, le, std::forward<Args>(args)...)
        ->antisymmetrize({{0, b, b}});
  }
  return evaluate(node, layout, le, std::forward<Args>(args)...)
      ->antisymmetrize(perm_groups);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
