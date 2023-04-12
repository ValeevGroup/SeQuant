#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>

#include <btas/btas.h>
#include <tiledarray.h>

#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <iostream>
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

#ifndef NDEBUG
template <typename C,
          typename = std::enable_if_t<std::is_same_v<C, sequant::Constant>>>
void assert_imaginary_zero(C const& c) {
  assert(c.value().imag() == 0 && "complex scalar unsupported");
}
#endif

template <typename BidirIt,
          typename Comp = std::less<decltype(*std::declval<BidirIt>())>>
bool next_permutation_parity(int& init_parity, BidirIt first, BidirIt last,
                             Comp&& comp = {}) {
  BidirIt i = last;
  if (first == last || first == --i) return false;

  for (;;) {
    BidirIt ii = i;
    --i;
    if (std::invoke(std::forward<Comp>(comp), *i, *ii)) {
      BidirIt j = last;

      while (!std::invoke(std::forward<Comp>(comp), *i, *(--j))) {
        // do nothing here
      }

      int parity = init_parity + 1 /* for the single iter_swap */
                   + std::distance(ii, last) / 2;
      init_parity = parity % 2;
      std::iter_swap(i, j);
      std::reverse(ii, last);
      return true;
    }
    if (i == first) {
      std::reverse(first, last);
      return false;
    }
  }
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

//==============================================================================
//                     TA::DistArray kernel definitions
//==============================================================================
class EvalExprTA final : public EvalExpr {
 public:
  ///
  /// annotation for TiledArray
  ///
  [[nodiscard]] std::string const& annot() const;

  ///
  /// Whether this object represents a tensor-of-tensor kind expression
  ///
  [[nodiscard]] bool tot() const;

  explicit EvalExprTA(Tensor const&);

  EvalExprTA(EvalExprTA const&, EvalExprTA const&, EvalOp);

 private:
  std::string annot_;

  bool tot_;
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
  else
    return braket_to_annot(node->tensor().const_braket());
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

namespace kernel {
template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
auto sum(NodeT const& n, TA::DistArray<Args...> const& lhs,
         TA::DistArray<Args...> const& rhs) {
  auto lscal = n.left()->scalar().value().real();
  auto rscal = n.right()->scalar().value().real();
  TA::DistArray<Args...> result;
  auto const& lannot = annot(n.left());
  auto const& rannot = annot(n.right());
  auto const& pannot = annot(n);

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] " << n->label() << "(" << pannot << ") = " << lscal
            << "*" << n.left()->label() << "(" << lannot << ") + " << rscal
            << "*" << n.right()->label() << "(" << rannot << ")" << std::endl;
#endif

  result(pannot) = lscal * lhs(lannot) + rscal * rhs(rannot);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename NodeT, typename... Args>
void add_to(TA::DistArray<Args...>& lhs, TA::DistArray<Args...> const& rhs,
            NodeT const& lnode, NodeT const& rnode) {
  using ranges::views::intersperse;
  using ranges::views::iota;
  using ranges::views::join;
  using ranges::views::transform;

  assert(lhs.trange() == rhs.trange());

  // Caveat:
  // This ->
  //  auto const annot = iota(size_t{0}, lhs.trange().rank()) |
  //               transform([](auto x) { return std::to_string(x); }) |
  //               intersperse(",") | join | ranges::to<std::string>;
  // Or, this? ->
  auto const annot = TA::detail::dummy_annotation(lhs.trange().rank());

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] " << lnode->label() << " += " << rnode->label()
            << std::endl;
#endif

  lhs(annot) += rhs(annot);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(lhs.world());
}

template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
auto prod(NodeT const& n, TA::DistArray<Args...> const& lhs,
          TA::DistArray<Args...> const& rhs) {
  TA::DistArray<Args...> result;
  auto const& lannot = annot(n.left());
  auto const& rannot = annot(n.right());
  auto const& pannot = annot(n);

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] " << n->label() << "(" << pannot
            << ") = " << n.left()->label() << "(" << lannot << ") * "
            << n.right()->label() << "(" << rannot << ")" << std::endl;
#endif

  result(pannot) = lhs(lannot) * rhs(rannot);

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
}

template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
auto scale(NodeT const& n, std::string const& annot,
           TA::DistArray<Args...> const& arr) {
  TA::DistArray<Args...> result;
  auto scalar = n->scalar().value().real();

#ifdef SEQUANT_EVAL_TRACE
  std::cout << "[EVAL] " << n->label() << "(" << annot << ") = " << scalar
            << "*" << n->label() << "(" << eval::annot(n) << ")" << std::endl;
#endif

  result(annot) = scalar * arr(eval::annot(n));

  TA::DistArray<Args...>::wait_for_lazy_cleanup(result.world());

  return result;
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

}  // namespace kernel

//==============================================================================

//=============================================================================
//                    btas::Tensor kernel definitions
//=============================================================================
///
/// \param bk iterable of Index objects.
/// \return vector of long-type hash values
///         of the labels of indices in \c bk
///
template <typename Iterable>
auto index_hash(Iterable const& bk) {
  return ranges::views::transform(
             bk,
             [](auto const& idx) {
               //
               // WARNING!
               // The BTAS expects index types to be long by default.
               // There is no straight-forward way to turn the default.
               // Hence here we explicitly cast the size_t values to long
               // Which is a potentailly narrowing conversion leading to
               // integral overflow. Hence, the values in the returned
               // container are mixed negative and positive integers (long type)
               //
               return static_cast<long>(hash::value(idx.label()));
             }) |
         ranges::to<container::svector<long>>;
}  // index_hash

namespace kernel {

template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
auto sum(NodeT const& n, btas::Tensor<Args...> const& lhs,
         btas::Tensor<Args...> const& rhs) {
  ///
  auto const post_annot = index_hash(n->tensor().const_braket());
  auto permute_and_scale = [&post_annot](auto const& btensor,
                                         auto const& child_seqt, auto scal) {
    auto pre_annot = index_hash(child_seqt.const_braket());
    btas::Tensor<Args...> result;
    btas::permute(btensor, pre_annot, result, post_annot);
    btas::scal(scal, result);
    return result;
  };
  ///
  auto lscal = n.left()->scalar().value().real();
  auto rscal = n.right()->scalar().value().real();

  auto sum = permute_and_scale(lhs, n.left()->tensor(), lscal);
  sum += permute_and_scale(rhs, n.right()->tensor(), rscal);

  return sum;
}

template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
void add_to(btas::Tensor<Args...>& lhs, btas::Tensor<Args...> const& rhs,
            NodeT const&, NodeT const&) {
  lhs += rhs;
}

template <typename NodeT, typename = std::enable_if_t<IsEvaluable<NodeT>>,
          typename... Args>
auto prod(NodeT const& n, btas::Tensor<Args...> const& lhs,
          btas::Tensor<Args...> const& rhs) {
  auto lannot = index_hash(n.left()->tensor().const_braket());
  auto rannot = index_hash(n.right()->tensor().const_braket());
  auto this_annot = index_hash(n->tensor().const_braket());

  btas::Tensor<Args...> prod;

  btas::contract(   //
      1.0,          //
      lhs, lannot,  //
      rhs, rannot,  //
      0.0,          //
      prod, this_annot);
  return prod;
}

template <typename NodeT, typename AnnotT,
          typename = std::enable_if_t<IsEvaluable<NodeT>>, typename... Args>
auto scale(NodeT const& n, AnnotT const& annot,
           btas::Tensor<Args...> const& arr) {
  btas::Tensor<Args...> result;
  auto const pre_annot = index_hash(n->tensor().const_braket());
  btas::permute(arr, pre_annot, result, annot);
  btas::scal(n->scalar().value().real(), result);
  return result;
}

template <typename... Args>
auto symmetrize(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // clang-format off
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_type>;
  // clang-format on
  auto const lannot = [rank]() {
    auto p = perm_type(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto symmetrizer = [&result, &lannot, &arr](auto const& permutation) {
    auto const& rannot = permutation;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    result += temp;
  };

  symmetric_permutation(rank / 2, symmetrizer);

  return result;
}

template <typename... Args>
auto antisymmetrize(btas::Tensor<Args...> const& arr) {
  using ranges::views::iota;

  size_t const rank = arr.rank();
  // Caveat:
  // auto const lannot = iota(size_t{0}, rank) | ranges::to<perm_type>;
  //
  auto const lannot = [rank]() {
    auto p = perm_type(rank);
    for (auto i = 0; i < rank; ++i) p[i] = i;
    return p;
  }();

  auto result = btas::Tensor<Args...>{arr.range()};
  result.fill(0);

  auto antisymmetrizer = [&result, &lannot, &arr](auto phase,
                                                  auto const& permutation) {
    auto const& rannot = permutation;
    btas::Tensor<Args...> temp;
    btas::permute(arr, lannot, temp, rannot);
    btas::scal(phase, temp);
    result += temp;
  };

  antisymmetric_permutation(rank / 2, antisymmetrizer);

  return result;
}

}  // namespace kernel

//
//=============================================================================

// =============================================================================
//                           Evaluation backend
// =============================================================================
template <typename NodeT, typename Le>
EvalResultT<NodeT, Le> evaluate_core(NodeT const& n, Le&& lev);

template <typename NodeT, typename Le, typename Cm,
          typename = std::enable_if_t<IsCacheManager<NodeT, Le, Cm>>>
EvalResultT<NodeT, Le> evaluate_core(NodeT const& n, Le&& lev, Cm&& cm);

template <typename NodeT, typename Le, typename... Args>
EvalResultT<NodeT, Le> evaluate_impl(NodeT const& n, Le&& lev, Args&&... args) {
#ifndef NDEBUG
  assert_imaginary_zero(n->scalar());
#endif
  if (n.leaf()) return std::invoke(std::forward<Le>(lev), n);
  EvalResultT<NodeT, Le> lres = evaluate_core(n.left(), std::forward<Le>(lev),
                                              std::forward<Args>(args)...);
  EvalResultT<NodeT, Le> rres = evaluate_core(n.right(), std::forward<Le>(lev),
                                              std::forward<Args>(args)...);
  if (n->op() == EvalOp::Sum) {
    return kernel::sum(n, lres, rres);
  }
  assert(n->op() == EvalOp::Prod);
  return kernel::prod(n, lres, rres);
}

template <typename NodeT, typename Le>
EvalResultT<NodeT, Le> evaluate_core(NodeT const& n, Le&& lev) {
  return evaluate_impl(n, std::forward<Le>(lev));
}

template <typename NodeT, typename Le, typename Cm, typename>
EvalResultT<NodeT, Le> evaluate_core(NodeT const& n, Le&& lev, Cm&& cm) {
  auto h = hash::value(*n);
  if (auto ptr = std::forward<Cm>(cm).access(h); ptr) {
#ifdef SEQUANT_EVAL_TRACE
    auto const max_c = std::forward<Cm>(cm).max_life(h);
    auto const curr_c = std::forward<Cm>(cm).life(h);
    std::cout << "[EVAL] [ACCESSED] intermediate for " << n->label() << "("
              << annot(n) << ") using key: " << hash::value(*n) << " ["
              << curr_c << "/" << max_c << "] accesses remain" << std::endl;
    if (curr_c == 0)
      std::cout << "[EVAL] [RELEASED] intermediate for " << n->label() << "("
                << annot(n) << ") using key: " << hash::value(*n) << std::endl;
#endif

    return *ptr;
  }

  if (std::forward<Cm>(cm).exists(h)) {
#ifdef SEQUANT_EVAL_TRACE
    auto const max_c = std::forward<Cm>(cm).max_life(h);
    // curr_c will be reduced by one right when it is stored
    auto const curr_c = std::forward<Cm>(cm).life(h) - 1;
    std::cout << "[EVAL] [STORED] [ACCESSED] intermediate for " << n->label()
              << "(" << annot(n) << ") using key: " << hash::value(*n) << " ["
              << curr_c << "/" << max_c << "] accesses remain" << std::endl;
#endif

    return *std::forward<Cm>(cm).store(
        h, std::move(
               evaluate_impl(n, std::forward<Le>(lev), std::forward<Cm>(cm))));
  }

#ifndef NDEBUG
  if (auto const curr_c = std::forward<Cm>(cm).life(h); curr_c == 0) {
    auto const max_c = std::forward<Cm>(cm).max_life(h);
    std::cout << "[EVAL] [MISSED] intermediate for " << n->label() << "("
              << annot(n) << ") using key: " << hash::value(*n) << " ["
              << curr_c << "/" << max_c << "] accesses remain" << std::endl;
  }

  assert(!std::forward<Cm>(cm).zombie(h) &&
         "Cached data outlives lifetime. Cache manager is in invalid state!");
#endif

  return evaluate_impl(n, std::forward<Le>(lev), std::forward<Cm>(cm));
}
// =============================================================================

// =============================================================================
//                         Evaluation frontend
// =============================================================================
template <typename NodeT, typename AnnotT, typename Le, typename... Args,
          std::enable_if_t<IsAnnot<AnnotT>, bool> = true>
EvalResultT<NodeT, Le> evaluate(NodeT const& n, AnnotT const& annot, Le&& leval,
                                Args&&... args) {
#ifndef NDEBUG
  assert_imaginary_zero(n->scalar());
#endif
  return kernel::scale(
      n, annot,
      evaluate_core(n, std::forward<Le>(leval), std::forward<Args>(args)...));
}

template <
    typename NodesT, typename AnnotT, typename Le, typename... Args,
    std::enable_if_t<IsAnnot<AnnotT> && IsIterableOfEvaluableNodes<NodesT>,
                     bool> = true>
auto evaluate(NodesT const& nodes, AnnotT const& annot, Le&& leval,
              Args&&... args) {
  auto beg = std::cbegin(nodes);
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  auto result = evaluate(*iter, annot, std::forward<Le>(leval),
                         std::forward<Args>(args)...);

  iter = std::next(iter);
  for (; iter != end; ++iter) {
    kernel::add_to(result,
                   evaluate(*iter, annot, std::forward<Le>(leval),
                            std::forward<Args>(args)...),
                   *beg, *iter);
  }
  return result;
}

template <typename NodeT, typename AnnotT, typename Le, typename... Args>
auto evaluate_symm(NodeT const& nodes, AnnotT const& annot, Le&& le,
                   Args&&... args) {
  using ranges::views::concat;
  // using ranges::views::iota;
  using ranges::views::transform;
#ifndef NDEBUG
  auto even_rank_assert = [](Tensor const& t) {
    assert((t.bra_rank() + t.ket_rank()) % 2 == 0 &&
           "Odd ranked tensor symmetrization not supported");
  };
  if constexpr (IsIterableOfEvaluableNodes<NodeT>)
    even_rank_assert((*std::begin(nodes))->tensor());
  else {
    static_assert(IsEvaluable<NodeT>);
    even_rank_assert(nodes->tensor());
  }
#endif
  return kernel::symmetrize(evaluate(nodes, annot, std::forward<Le>(le),
                                     std::forward<Args>(args)...));
}

template <typename NodeT, typename StrT, typename Le, typename... Args>
auto evaluate_antisymm(NodeT const& nodes, StrT const& annot, Le&& le,
                       Args&&... args) {
#ifndef NDEBUG
  auto even_rank_assert = [](Tensor const& t) {
    assert((t.bra_rank() + t.ket_rank()) % 2 == 0 &&
           "Odd ranked tensor anti-symmetrization not supported");
  };
  if constexpr (IsIterableOfEvaluableNodes<NodeT>)
    even_rank_assert((*std::begin(nodes))->tensor());
  else {
    static_assert(IsEvaluable<NodeT>);
    even_rank_assert(nodes->tensor());
  }
#endif
  return kernel::antisymmetrize(evaluate(nodes, annot, std::forward<Le>(le),
                                         std::forward<Args>(args)...));
}

// =============================================================================

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
