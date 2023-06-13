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

namespace {

template <typename N, typename = void>
struct IsFbNode_ : std::false_type {};

template <typename T>
struct IsFbNode_<FullBinaryNode<T>> : std::true_type {};

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

template <typename N>
constexpr bool IsFbNode = IsFbNode_<N>::value;

template <typename, typename = void>
constexpr bool IsEvaluable{};

template <typename T>
constexpr bool IsEvaluable<
    FullBinaryNode<T>, std::enable_if_t<std::is_convertible_v<T, EvalExpr>>> =
    true;

template <typename T>
constexpr bool
    IsEvaluable<const FullBinaryNode<T>,
                std::enable_if_t<std::is_convertible_v<T, EvalExpr>>> = true;

template <typename NodeT, typename E,
          typename = std::enable_if_t<IsEvaluable<NodeT>>>
using EvalResultT = remove_cvref_t<std::invoke_result_t<E, NodeT const&>>;

template <typename Iterable, typename = void>
constexpr bool IsIterableOfEvaluableNodes{};

template <typename Iterable>
constexpr bool IsIterableOfEvaluableNodes<
    Iterable, std::enable_if_t<IsEvaluable<IteredT<Iterable>>>> = true;

template <typename NodeT, typename Le, typename Cm,
          typename = std::enable_if_t<std::is_convertible_v<
              typename std::remove_reference_t<Cm>::cached_type, ERPtr>>>
constexpr bool IsCacheManager = true;

template <typename F,
          typename = std::enable_if_t<std::is_invocable_v<F, perm_type const&>>>
void symmetric_permutation(size_t half_rank, F&& symmetrizer) noexcept {
  using ranges::views::concat;
  // using ranges::views::iota;
  using ranges::views::transform;

  // this vector contains indices from 0...rank/2 where rank is that
  // of the tensor being symmetrized
  //
  // Caveat:
  // creating perm_vec this way is not allowed in gcc-11,
  // sth to do with container::svector (boost::svector)
  // clang-format off
  // auto perm_vec = iota(size_t{0}, half_rank) | ranges::to<perm_type>;
  // clang-format on
  //
  auto perm_vec = perm_type(half_rank);
  for (auto i = 0; i < half_rank; ++i) perm_vec[i] = i;

  auto add_half_rank = [half_rank](auto x) { return x + half_rank; };
  do {
    auto const total_perm =
        concat(perm_vec, perm_vec | transform(add_half_rank)) |
        ranges::to<perm_type>;
    std::forward<F>(symmetrizer)(total_perm);
  } while (std::next_permutation(std::begin(perm_vec), std::end(perm_vec)));
}

template <typename F, typename = std::enable_if_t<
                          std::is_invocable_v<F, double, perm_type const&>>>
void antisymmetric_permutation(size_t half_rank, F&& antisymmetrizer) noexcept {
  using phase_type = double;

  using ranges::views::concat;
  //  using ranges::views::iota;
  using ranges::views::join;

  int bra_parity = 0;
  //
  // Caveat:
  // creating bra_perm_vec this way is not allowed in gcc-11,
  // sth to do with container::svector (boost::svector)
  // clang-format off
  // auto bra_perm_vec = iota(size_t{0}, half_rank) | ranges::to<perm_type>;
  // clang-format on
  //
  auto bra_perm_vec = perm_type(half_rank);
  for (auto i = 0; i < half_rank; ++i) bra_perm_vec[i] = i;

  do {
    int ket_parity = 0;
    // same problem as with bra_perm_vec
    // clang-format off
    // auto ket_perm_vec = iota(half_rank, 2 * half_rank) | ranges::to<perm_type>;
    // clang-format on
    auto ket_perm_vec = perm_type(half_rank);
    for (auto i = 0; i < half_rank; ++i) ket_perm_vec[i] = i + half_rank;
    do {
      phase_type const phase_factor =
          (bra_parity + ket_parity) % 2 == 0 ? 1 : -1;
      auto const perm_vec =
          concat(bra_perm_vec, ket_perm_vec) | ranges::to<perm_type>;
      std::forward<F>(antisymmetrizer)(phase_factor, perm_vec);
    } while (next_permutation_parity(ket_parity,
                                     std::begin(ket_perm_vec),  //
                                     std::end(ket_perm_vec)));
  } while (next_permutation_parity(bra_parity,
                                   std::begin(bra_perm_vec),  //
                                   std::end(bra_perm_vec)));
}

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

///
/// Given an iterable of Index objects, generate a string annotation
/// that can be used for TiledArray tensor expressions.
/// Tensor-of-tensors also supported.
template <typename Indices>
std::string braket_to_annot(Indices const& indices) {
  using ranges::find;
  using ranges::views::filter;
  using ranges::views::intersperse;
  using ranges::views::join;

  // make a comma-separated string out of an iterable of strings
  auto add_commas = [](auto const& strs) -> std::string {
    return strs | intersperse(",") | join | ranges::to<std::string>;
  };

  container::svector<std::string> idxs{}, pidxs{};
  for (auto&& idx : indices) {
    idxs.emplace_back(sequant::to_string(idx.label()));
    for (auto&& pidx : idx.proto_indices())
      pidxs.emplace_back(sequant::to_string(pidx.label()));
  }

  if (pidxs.empty()) {
    // not a tensor-of-tensor type expression
    return add_commas(idxs);
  } else {
    ranges::stable_sort(pidxs);
    ranges::actions::unique(pidxs);
    auto not_in_pidx = [&pidxs](auto&& l) {
      return find(pidxs, l) == pidxs.end();
    };
    return add_commas(pidxs) + ";" +
           add_commas(idxs | filter(not_in_pidx) | ranges::to<decltype(idxs)>);
  }
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

template <typename NodeT, typename Le>
ERPtr evaluate_crust(NodeT const&, Le&&);

template <typename NodeT, typename Le, typename Cm>
ERPtr evaluate_crust(NodeT const&, Le&&, Cm&&);

template <typename NodeT, typename Le, typename... Args>
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

template <typename NodeT, typename Le>
ERPtr evaluate_crust(NodeT const& node, Le&& le) {
  return evaluate_core(node, std::forward<Le>(le));
}

template <typename NodeT, typename Le, typename Cm>
ERPtr evaluate_crust(NodeT const& node, Le&& le, Cm&& cm) {
  auto&& cache = std::forward<Cm>(cm);
  auto const h = hash::value(*node);
  if (auto ptr = cache.access(h); ptr) {
    return *ptr;
  } else if (cache.exists(h)) {
    return *cache.store(
        h, evaluate_core(node, std::forward<Le>(le), std::forward<Cm>(cm)));
  } else {
    return evaluate_core(node, std::forward<Le>(le), std::forward<Cm>(cm));
  }
}

template <typename NodeT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsEvaluable<NodeT>, bool> = true>
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
  iter = std::next(iter);
  for (; iter != end; ++iter) {
    result->add_inplace(*evaluate(*iter, layout, std::forward<Le>(le),
                                  std::forward<Args>(args)...));
  }
  return result;
}

template <typename NodeT, typename Annot, typename Le, typename... Args>
auto evaluate_symm(NodeT const& node,    //
                   Annot const& layout,  //
                   Le&& le, Args&&... args) {
  return evaluate(node, layout, std::forward<Le>(le),
                  std::forward<Args>(args)...)
      ->symmetrize();
}

template <typename NodeT, typename Annot, typename Le,
          typename... Args>
auto evaluate_antisymm(NodeT const& node,    //
                       Annot const& layout,  //
                       Le&& le, Args&&... args) {
  return evaluate(node, layout, std::forward<Le>(le),
                  std::forward<Args>(args)...)
      ->antisymmetrize();
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
