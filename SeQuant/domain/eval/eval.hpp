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
/// Make a @c CacheManager object that caches the data associated with an
/// intermediate for a certain number of accesses. After that many accesses, the
/// associated resource is freed. Aditionally, store the data associated with
/// the leaf nodes for indefinitely many accesses if @c persistent_leaves in on.
///
/// \tparam T type of the objects to be cached associated with an
///                  intermediate or a leaf node.
///
/// \param nodes Range of FullBinaryNode objects. Elements of @c nodes should
///              have a method named <code> hash() </code> that returns size_t.
///
/// \param proj Only cache data for nodes that <code> proj(node) </code>
///             evaluates to true. Useful to avoid caching trivially computable
///             intermediates or the leaves that should not be cached at all.
///             For example the 't' amplitude tensors in coupled-cluster
///             expressions. @note @c proj is called on both internal and leaf
///             nodes.
/// \param persistent_leaves Whether the cache manager should store the data
///                          associated with the leaf nodes indefinitely.
/// \param min_repeats Only cache intermediates that repeats at least this many
///                    times.
///
/// \todo Restrictions on template params.
///
template <typename T, typename NodesRng, typename P>
CacheManager<T> cache_manager(NodesRng const& nodes, P&& proj,
                              bool persistent_leaves, size_t min_repeats = 2) {
  using key_t = typename CacheManager<T>::key_t;
  using count_t = typename CacheManager<T>::count_t;
  using map_t = container::map<key_t, count_t>;

  auto imed_counts = map_t{};

  // counts number of times each internal node appears in
  // all of @c nodes trees
  auto imed_visitor = [&proj, &imed_counts](auto&& n) {
    if (!std::invoke(std::forward<P>(proj), n)) return;

    auto&& end = imed_counts.end();
    auto&& h = n->hash();
    if (auto&& found = imed_counts.find(h); found != end)
      ++found->second;
    else
      imed_counts.emplace(h, 1);
  };

  ranges::for_each(nodes, [&imed_visitor](auto&& tree) {
    tree.visit_internal(imed_visitor);
  });

  // remove less repeating imeds
  auto less_repeating = [min_repeats](auto&& pair) {
    return pair.second < min_repeats;
  };
  ranges::actions::remove_if(imed_counts, less_repeating);

  auto leafs = container::set<key_t>{};
  auto leaf_visitor = [&proj, &leafs](auto&& n) {
    if (std::invoke(std::forward<P>(proj), n)) {
      leafs.emplace(n->hash());
    }
  };

  if (persistent_leaves)
    ranges::for_each(
        nodes, [&leaf_visitor](auto&& tree) { tree.visit_leaf(leaf_visitor); });

  return {imed_counts, leafs};
}

struct UncacheAmplitudeTensors {
  template <typename N>
  [[nodiscard]] bool operator()(N&& node) const noexcept {
    return node->tensor().label() != L"t";
  }
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_EVAL_HPP
