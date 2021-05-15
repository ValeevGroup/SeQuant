#ifndef SEQUANT_DOMAIN_EVAL_HPP
#define SEQUANT_DOMAIN_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/cache_manager.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>
#include <type_traits>

namespace sequant::eval {
using eval_node = utils::binary_node<utils::eval_expr>;

using perm_type = container::svector<size_t>;

using phase_type = int;

struct perm_with_phase {
  phase_type phase;
  perm_type const& perm;
};

namespace detail {

template <typename F>
void permute_ords(perm_type& ords, F&& callback, size_t beg = 0,
                  size_t swaps = 0) {
  static_assert(std::is_invocable_v<F, perm_with_phase const&>,
                "F(perm_with_phase) not possible");

  if (beg + 1 == ords.size())
    callback(perm_with_phase{swaps % 2 == 0 ? 1 : -1, ords});

  for (auto ii = beg; ii < ords.size(); ++ii) {
    std::swap(ords[beg], ords[ii]);
    permute_ords(ords, std::forward<F>(callback), beg + 1,
                 ii == beg ? swaps : swaps + 1);
    std::swap(ords[ii], ords[beg]);
  }
}  // permute_ords

template <typename Maplike>
void count_imeds(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Maplike& cmap) {
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

template <typename F>
void antisymmetrize_tensor(size_t rank, F&& callback) {
  static_assert(std::is_invocable_v<F, perm_with_phase const&>);

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
        auto total_perm = perm_with_phase{bp.phase * kp.phase, annot};
        std::invoke(callback, total_perm);
      });  // permute ket
    });    // permute bra
  };       // asymm_impl

  asymm_impl(rank, std::forward<F>(callback));
}

template <typename Tensor_t>
utils::cache_manager<Tensor_t> make_cache_man(
    utils::binary_node<utils::eval_expr> const& node, bool persistent_leaves) {
  container::map<size_t, size_t> hash_to_counts{};
  detail::count_imeds(node, hash_to_counts);

  auto less_repeating = [](auto const& pair) { return pair.second < 2; };
  ranges::actions::remove_if(hash_to_counts, less_repeating);

  if (!persistent_leaves) return utils::cache_manager<Tensor_t>{hash_to_counts};

  container::svector<size_t> leaf_hashes{};
  node.visit_leaf([&leaf_hashes](auto const& node) {
    leaf_hashes.emplace_back(node->hash());
  });

  return utils::cache_manager<Tensor_t>{hash_to_counts, leaf_hashes};
}

template <typename Tensor_t, typename Iterable>
utils::cache_manager<Tensor_t> make_cache_man(Iterable const& nodes,
                                              bool persistent_leaves) {
  container::map<size_t, size_t> hash_to_counts{};

  for (auto const& n : nodes) detail::count_imeds(n, hash_to_counts);

  auto less_repeating = [](auto const& pair) { return pair.second < 2; };
  ranges::actions::remove_if(hash_to_counts, less_repeating);

  if (!persistent_leaves) return utils::cache_manager<Tensor_t>{hash_to_counts};

  container::set<size_t> leaf_hashes{};
  for (auto const& n : nodes) {
    n.visit_leaf([&leaf_hashes](auto const& node) {
      if (node->tensor().label() != L"t") leaf_hashes.insert(node->hash());
    });
  }

  return utils::cache_manager<Tensor_t>{hash_to_counts, leaf_hashes};
}

}  // namespace sequant::eval

#endif  // SEQUANT_DOMAIN_EVAL_HPP
