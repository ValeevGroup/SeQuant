#ifndef SEQUANT_DOMAIN_EVAL_HPP
#define SEQUANT_DOMAIN_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>
#include <type_traits>

namespace sequant::eval {
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

template <typename Data_t>
class cache_manager {
 public:
  /**
   * Key type to access data.
   */
  using key_t = size_t;

  /**
   * Count type to allow max access of a stored datum.
   */
  using count_t = size_t;

  enum struct Lifetime { Persistent, Decaying };

 private:
  struct data_cache {
    Lifetime lifetime{Lifetime::Persistent};
    count_t max_access = 0;
    std::unique_ptr<Data_t> ptr{nullptr};
  };  //

  container::map<key_t, data_cache> cache{};

 public:
  cache_manager() = default;

  void clear() { cache.clear(); }

  void add_key_persistent(key_t k) { cache.try_emplace(k, data_cache{}); }

  void add_key_persistent(container::svector<key_t> keys) {
    for (auto k : keys) add_key_persistent(k);
  }

  void add_key_decaying(key_t k, count_t c = 1) {
    if (c == 0) return;
    if (auto&& found = cache.find(k); found != cache.end())
      ++found->second.max_access;
    else
      cache.emplace(k, data_cache{Lifetime::Decaying, c, nullptr});
  }

  void add_key_decaying(container::map<key_t, count_t> const& counts) {
    for (auto&& [k, c] : counts) add_key_decaying(k, c);
  }

  container::map<key_t, count_t> key_counts_decaying() const {
    return cache | ranges::views::filter([](auto const& pair) {
             return pair.second.lifetime == Lifetime::Decaying;
           }) |
           ranges::views::transform([](auto const& pair) {
             return std::pair{pair.first, pair.second.max_access};
           }) |
           ranges::to<container::map<key_t, count_t>>;
  }

  container::svector<key_t> keys_persistent() const {
    return cache | ranges::views::filter([](auto const& pair) {
             return pair.second.lifetime == Lifetime::Persistent;
           }) |
           ranges::views::transform(
               [](auto const& pair) { return pair.first; }) |
           ranges::to<container::svector<key_t>>;
  }

  void clear_persistent() {
    auto decaying = key_counts_decaying();
    clear();
    add_key_decaying(decaying);
  }

  void clear_decaying() {
    auto persistent = keys_persistent();
    clear();
    add_key_persistent(persistent);
  }

  void keep_more_repeating(count_t min_count = 2) {
    auto decaying = key_counts_decaying();
    ranges::actions::remove_if(decaying, [min_count](auto const& pair) {
      return pair.second <= min_count;
    });
    clear_decaying();
    add_key_decaying(decaying);
  }

  Data_t store(key_t h, Data_t val) {
    auto&& found = cache.find(h);
    if (found == cache.end()) return std::move(val);  // stores nothing

    if (found->second.lifetime == Lifetime::Decaying &&
        found->second.max_access == 1) {
      cache.erase(found);
      return std::move(val);
    }

    assert(found->second.ptr == nullptr &&
           "re-write attempted on existing data");

    found->second.ptr = std::make_unique<Data_t>(std::move(val));

    if (found->second.lifetime == Lifetime::Decaying)
      --found->second.max_access;

    return *found->second.ptr;
  }

  std::optional<Data_t> access(key_t h) {
    auto&& found = cache.find(h);
    if (found == cache.end() || found->second.ptr == nullptr)
      return std::nullopt;

    if (found->second.lifetime == Lifetime::Persistent)
      return *found->second.ptr;

    --found->second.max_access;

    if (found->second.max_access == 0) {
      auto temp = *found->second.ptr;
      cache.erase(found);
      return temp;
    } else {
      return *found->second.ptr;
    }
  }
};  // cache manager

}  // namespace sequant::eval

#endif  // SEQUANT_DOMAIN_EVAL_HPP
