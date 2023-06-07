#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include "eval_result.hpp"

#include <SeQuant/core/algorithm.hpp>
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

template <typename T>
struct IsString_ : public std::disjunction<
                       std::is_same<char*, typename std::decay_t<T>>,
                       std::is_same<const char*, typename std::decay_t<T>>,
                       std::is_same<std::string, typename std::decay_t<T>>> {};

template <typename T, typename = void>
struct IsAnnot_ : std::false_type {};

template <typename T>
struct IsAnnot_<T, std::enable_if_t<IsString_<T>::value>> : std::true_type {};

template <typename T>
struct IsAnnot_<std::initializer_list<T>,
                std::enable_if_t<std::is_integral_v<T>>> : std::true_type {};

template <typename T>
struct IsAnnot_<container::svector<T>, std::enable_if_t<std::is_integral_v<T>>>
    : std::true_type {};

template <typename T, typename = void>
struct HasAnnotMethod_ : std::false_type {};

template <typename T>
struct HasAnnotMethod_<
    T, std::enable_if_t<std::is_same_v<
           std::string, remove_cvref_t<decltype(std::declval<T>().annot())>>>>
    : std::true_type {};

}  // namespace

using perm_type = container::svector<size_t>;

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

template <typename StrT>
constexpr bool IsString = IsString_<StrT>::value;

template <typename T>
constexpr bool IsAnnot = IsAnnot_<T>::value;

template <typename Iterable, typename = void>
constexpr bool IsIterableOfEvaluableNodes{};

template <typename Iterable>
constexpr bool IsIterableOfEvaluableNodes<
    Iterable, std::enable_if_t<IsEvaluable<IteredT<Iterable>>>> = true;

template <typename NodeT, typename Le, typename Cm,
          typename = std::enable_if_t<std::is_convertible_v<
              typename std::remove_reference_t<Cm>::cached_type,
              EvalResultT<NodeT, Le>>>>
constexpr bool IsCacheManager = true;

template <typename C,
          typename = std::enable_if_t<std::is_same_v<C, sequant::Constant>>>
void assert_imaginary_zero(C const& c) {
#ifndef NDEBUG
  assert(c.value().imag() == 0 && "complex scalar unsupported");
#endif
}

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

template <typename DataT, typename NodesI,
          typename Pred = std::function<bool(IteredT<NodesI> const&)>,
          typename = std::enable_if_t<
              IsIterableOfEvaluableNodes<NodesI> &&
              std::is_invocable_r_v<bool, Pred, IteredT<NodesI> const&>>>
CacheManager<DataT> cache_manager(
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

  return CacheManager<DataT>{imed_counts};
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

template <typename T>
constexpr bool HasAnnotMethod =
    HasAnnotMethod_<std::remove_reference_t<T>>::value;

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

template <typename N, typename = std::enable_if_t<IsEvaluable<N>>>
std::string annot(N const& node) {
  if constexpr (HasAnnotMethod<typename N::value_type>)
    return node->annot();
  else if (node->expr()->template is<Constant>())
    return "";
  assert(node->expr()->template is<Tensor>());
  return braket_to_annot(node->as_tensor().const_braket());
}

template <typename RngOfOrdinals>
std::string ords_to_annot(RngOfOrdinals const& ords) {
  using ranges::views::intersperse;
  using ranges::views::join;
  using ranges::views::transform;
  auto to_str = [](auto x) { return std::to_string(x); };
  return ords | transform(to_str) | intersperse(std::string{","}) | join |
         ranges::to<std::string>;
}

namespace {

template <typename... Args>
auto permute(TA::DistArray<Args...> const& pre,  //
             std::string const& pre_annot,       //
             std::string const& post_annot) noexcept {
  TA::DistArray<Args...> post;
  post(post_annot) = pre(pre_annot);
  decltype(post)::wait_for_lazy_cleanup(post.world());
  return post;
}

template <typename... Args>
auto symmetrize(TA::DistArray<Args...> const& arr) {
  auto result = TA::DistArray<Args...>{arr.world(), arr.trange()};
  result.fill(0);
  size_t rank = arr.trange().rank();
  auto const lannot = ords_to_annot(ranges::views::iota(size_t{0}, rank));

  auto symmetrizer = [&result, &lannot, &arr](auto const& permutation) {
    auto const rannot = ords_to_annot(permutation);
    result(lannot) += arr(rannot);
  };

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] symmetrizing rank-" << rank << " tensor" << std::endl;
#endif

  symmetric_permutation(rank / 2, symmetrizer);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
auto antisymmetrize(TA::DistArray<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.trange().rank();
  auto const lannot = ords_to_annot(iota(size_t{0}, rank));

  auto result = TA::DistArray<Args...>{arr.world(), arr.trange()};
  result.fill(0);

  auto antisymmetrizer = [&result, &lannot, &arr](auto phase,
                                                  auto const& permutation) {
    auto const rannot = ords_to_annot(permutation);
    result(lannot) += phase * arr(rannot);
  };

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] antisymmetrizing rank-" << rank << " tensor"
            << std::endl;
#endif

  antisymmetric_permutation(rank / 2, antisymmetrizer);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename... Args>
void add_to(TA::DistArray<Args...>& lhs, TA::DistArray<Args...> const& rhs) {
  assert(lhs.trange() == rhs.trange());

  auto const annot = TA::detail::dummy_annotation(lhs.trange().rank());

  lhs(annot) += rhs(annot);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(lhs.world());
}

}  // namespace

///
/// \tparam NodeT Node type. eg. FullBinaryNode<EvalExpr>
/// \tparam Le    Leaf evaluator type.
/// \param node   Node to evaluate
/// \param le     Leaf evaluator: returns result for leaf nodes.
/// \return
template <typename NodeT, typename Le,
          typename = std::enable_if_t<IsEvaluable<NodeT>>>
ERPtr evaluate_(NodeT const& node, Le&& le) {
  using numeric_type = typename EvalResultT<NodeT, Le>::numeric_type;
  using constant_type = EvalConstant<numeric_type>;
  if (node.leaf()) {
    if (node->result_type() == ResultType::Constant) {
      auto n = constant_type{node->as_constant().value().real()}.value();
      return eval_result<constant_type>(n);
    } else {
      assert(node->result_type() == ResultType::Tensor);
      return eval_result<EvalTensorTA<EvalResultT<NodeT, Le>>>(
          std::forward<Le>(le)(node));
    }
  } else {
    auto const left = evaluate_(node.left(), std::forward<Le>(le));
    auto const right = evaluate_(node.right(), std::forward<Le>(le));

    std::array<std::any, 3> ann{annot(node.left()),   //
                                annot(node.right()),  //
                                annot(node)};

    if (node->op_type() == EvalOp::Sum) {
      return left->sum(*right, ann);
    } else {
      assert(node->op_type() == EvalOp::Prod);
      return left->prod(*right, ann);
    }
  }
}

template <typename NodeT, typename Annot, typename Le,
          std::enable_if_t<IsAnnot<Annot> && IsEvaluable<NodeT>, bool> = true>
auto evaluate(NodeT const& node,    //
              Annot const& layout,  //
              Le&& le) {
  auto const pre = evaluate_(node, std::forward<Le>(le))
                       ->template get<EvalResultT<NodeT, Le>>();
  return permute(pre, annot(node), layout);
}

template <typename NodesT, typename Annot, typename Le,
          std::enable_if_t<IsAnnot<Annot> && IsIterableOfEvaluableNodes<NodesT>,
                           bool> = true>
auto evaluate(NodesT const& nodes,  //
              Annot const& layout,  //
              Le&& le) {
  auto beg = std::cbegin(nodes);
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter, layout, std::forward<Le>(le));

  iter = std::next(iter);
  for (; iter != end; ++iter) {
    add_to(result, evaluate(*iter, layout, std::forward<Le>(le)));
  }
  return result;
}

template <typename NodeT, typename Annot, typename Le>
auto evaluate_symm(NodeT const& node,    //
                   Annot const& layout,  //
                   Le&& le) {
  auto const pre = evaluate(node, layout, std::forward<Le>(le));
  return symmetrize(pre);
}

template <typename NodeT, typename Annot, typename Le>
auto evaluate_antisymm(NodeT const& node,    //
                       Annot const& layout,  //
                       Le&& le) {
  auto const pre = evaluate(node, layout, std::forward<Le>(le));
  return antisymmetrize(pre);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
