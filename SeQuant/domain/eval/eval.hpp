#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>
#include <type_traits>

namespace sequant::eval {

/// Represents permutation of a sequence.
/// Equivalent to the range [0, N) where N is the length of the sequence to be
/// permuted.
using perm_type = container::svector<size_t>;

/// Even or odd permuation. Negative if odd, positive if even.
using phase_type = int;

///
/// Keeps a permutation and its phase together
///
struct PermWithPhase {
  phase_type phase;
  perm_type const& perm;
};

namespace detail {

template <typename F>
void permute_ords(perm_type& ords, F&& callback, size_t beg = 0,
                  size_t swaps = 0) {
  static_assert(std::is_invocable_v<F, PermWithPhase const&>,
                "F(perm_with_phase) not possible");

  if (beg + 1 == ords.size())
    callback(PermWithPhase{swaps % 2 == 0 ? 1 : -1, ords});

  for (auto ii = beg; ii < ords.size(); ++ii) {
    std::swap(ords[beg], ords[ii]);
    permute_ords(ords, std::forward<F>(callback), beg + 1,
                 ii == beg ? swaps : swaps + 1);
    std::swap(ords[ii], ords[beg]);
  }
}  // permute_ords

template <typename Maplike>
void count_imeds(EvalNode const& node, Maplike& cmap) {
  if (node.leaf()) return;

  if (auto&& found = cmap.find(node->hash()); found != cmap.end()) {
    ++found->second;
    return;
  } else {
    cmap.emplace(node->hash(), 1);
  }
  count_imeds<Maplike>(node.left(), cmap);
  count_imeds<Maplike>(node.right(), cmap);
}  // count_imeds

}  // namespace detail

///
/// \tparam F Type of the callable.
/// \param rank Rank of the tensor to be particle-symmetrized.
///             Needed for permutation vector length determination.
/// \param callback A function that takes @c perm_type object and does
///                 something. For example call TiledArray to perform actual
///                 symmetrization of a DistArray.
template <typename F>
void symmetrize_tensor(size_t rank, F&& callback) {
  static_assert(std::is_invocable_v<F, perm_type>);

  using ranges::views::iota;

  auto add_half_rank = [n = rank / 2](auto x) { return x + n; };

  auto perm_vec = iota(size_t{0}, rank / 2) | ranges::to<perm_type>;
  do {
    auto const rannot =
        ranges::views::concat(
            perm_vec, perm_vec | ranges::views::transform(add_half_rank)) |
        ranges::to<perm_type>;

    std::invoke(callback, rannot);

  } while (std::next_permutation(perm_vec.begin(), perm_vec.end()));
}

///
/// \tparam F Type of the callable.
/// \param rank Rank of the tensor to be particle-antisymmetrized.
///             Needed for permutation vector length determination.
/// \param callback A function that takes @c PermWithPhase object and does
///                 something. For example call TiledArray to perform actual
///                 antisymmetrization of a DistArray.
template <typename F>
void antisymmetrize_tensor(size_t rank, F&& callback) {
  static_assert(std::is_invocable_v<F, PermWithPhase const&>);

  auto asymm_impl = [](size_t rank, F&& callback) -> void {
    assert(rank % 2 == 0 &&
           "odd ranked tensor antisymmetrization not supported");
    auto const init_perm =
        ranges::views::iota(size_t{0}, rank / 2) | ranges::to<perm_type>;
    auto bra_perm = init_perm;
    detail::permute_ords(bra_perm, [&init_perm, &callback,
                                    n = rank / 2](auto const& bp) {
      auto ket_perm = init_perm |
                      ranges::views::transform([n](auto x) { return n + x; }) |
                      ranges::to<perm_type>;

      detail::permute_ords(ket_perm, [&callback, &bp](auto const& kp) {
        auto annot =
            ranges::views::concat(bp.perm, kp.perm) | ranges::to<perm_type>;
        auto total_perm = PermWithPhase{bp.phase * kp.phase, annot};
        std::invoke(callback, total_perm);
      });  // permute ket
    });    // permute bra
  };       // asymm_impl

  asymm_impl(rank, std::forward<F>(callback));
}

///
/// \tparam T type of the objects to be cached associated with an
///                  intermediate or a leaf node.
///
/// \param node An @c EvalNode object.
/// \param persistent_leaves Whether the cache manager should store the data
///                          associated with the leaf nodes indefinitely.
/// \return A cache manager object. @see CacheManager.
template <typename T>
CacheManager<T> make_cache_manager(EvalNode const& node,
                                   bool persistent_leaves) {
  container::map<size_t, size_t> hash_to_counts{};
  detail::count_imeds(node, hash_to_counts);

  auto less_repeating = [](auto const& pair) { return pair.second < 2; };
  ranges::actions::remove_if(hash_to_counts, less_repeating);

  if (!persistent_leaves) return CacheManager<T>{hash_to_counts, {}};

  container::svector<size_t> leaf_hashes{};
  node.visit_leaf([&leaf_hashes](auto const& node) {
    leaf_hashes.emplace_back(node->hash());
  });

  return CacheManager<T>{hash_to_counts, leaf_hashes};
}

///
/// Make a @c CacheManager object that caches the data associated with an
/// intermediate for a certain number of accesses. After that many accesses, the
/// associated resource is freed. Also, optionally store the data associated
/// with the leaf nodes for indefinitely many accesses. This is a specific
/// use-case function intended to be used in Coupled-Cluster calculations.
/// Because the amplitude tensors (the leaf nodes whose tensor has the label
/// 't') are updated every iteration, their data should never be reused -- and
/// thus stored. This fact is taken into consideration by this function.
/// \tparam  T type of the objects to be cached associated with an intermediate
///            or a leaf node.
///
/// \param nodes An iterable of @c EvalNode objects.
/// \param persistent_leaves Whether the cache manager should store the data
///                          associated with the leaf nodes indefinitely.
/// \return A cache manager object. @see CacheManager.
///
template <typename T, typename Iterable>
CacheManager<T> make_cache_manager(Iterable const& nodes,
                                   bool persistent_leaves) {
  container::map<size_t, size_t> hash_to_counts{};

  for (auto const& n : nodes) detail::count_imeds(n, hash_to_counts);

  auto less_repeating = [](auto const& pair) { return pair.second < 2; };
  ranges::actions::remove_if(hash_to_counts, less_repeating);

  if (!persistent_leaves) return CacheManager<T>{hash_to_counts, {}};

  container::set<size_t> leaf_hashes{};
  for (auto const& n : nodes) {
    n.visit_leaf([&leaf_hashes](auto const& node) {
      if (node->tensor().label() != L"t") leaf_hashes.insert(node->hash());
    });
  }

  return CacheManager<T>{hash_to_counts, leaf_hashes};
}

///
/// Make @c CacheManager object that persistently stores the data associated
/// with the leaf nodes of an  @c EvalNode object. This is a specific use-case
/// function intended to be used in Coupled-Cluster calculations. Because the
/// amplitude tensors (the leaf nodes whose tensor has the label 't') are
/// updated every iteration, their data should never be reused -- and thus
/// stored. This fact is taken into consideration by this function.
///
/// \tparam T type of the objects to be cached associated with the leaf nodes.
///
/// \param nodes An iterable of @c EvalNode objects.
/// \return A cache manager object. @see CacheManager.
///
template <typename T, typename Iterable>
CacheManager<T> make_cache_manager_leaves_only(Iterable const& nodes) {
  container::set<EvalExpr::hash_t> leaf_hashes{};
  for (auto const& n : nodes) {
    n.visit_leaf([&leaf_hashes](auto const& node) {
      if (node->tensor().label() != L"t") leaf_hashes.insert(node->hash());
    });
  }
  return CacheManager<T>{{}, leaf_hashes};
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
