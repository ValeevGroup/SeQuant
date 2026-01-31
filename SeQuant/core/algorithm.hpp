//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_ALGORITHM_HPP
#define SEQUANT_ALGORITHM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/meta.hpp>

#include <boost/dynamic_bitset.hpp>
#include <range/v3/view.hpp>

#include <algorithm>
#include <functional>
#include <iterator>
#include <tuple>
#include <type_traits>

namespace sequant {

template <typename Callable, typename... Args>
using suitable_call_operator =
    decltype(std::declval<Callable>()(std::declval<Args>()...));

/// @brief bubble sort that uses swap exclusively
/// @tparam ForwardIter a forward iterator
/// @tparam Sentinel a sentinel type
/// @tparam Compare a less-than relationship that specifies a weak order
/// @note does not assume `Compare` to be strict, unlike `std::sort`
template <typename ForwardIter, typename Sentinel,
          typename Compare = std::less<>>
void bubble_sort(ForwardIter begin, Sentinel end, Compare comp = {}) {
  if (begin == end) return;  // no-op for empty range
  bool swapped;
  do {
    swapped = false;
    auto i = begin;
    auto inext = i;

    using deref_type = decltype(*(std::declval<ForwardIter>()));
    constexpr const bool comp_works_for_range_type =
        meta::is_detected_v<suitable_call_operator, Compare, deref_type,
                            deref_type>;

    for (++inext; inext != end; ++i, ++inext) {
      if constexpr (comp_works_for_range_type) {
        const auto& val0 = *inext;
        const auto& val1 = *i;
        if (comp(val0, val1)) {
          // current assumption: whenever iter_swap from below does not fall
          // back to std::iter_swap, we are handling zipped ranges where the
          // tuple sizes is two (even) -> thus using a non-std swap
          // implementation won't mess with the information of whether or not an
          // even amount of swaps has occurred.
          using ranges::iter_swap;
          iter_swap(i, inext);
          swapped = true;
        }
      } else {
        const auto& val0 = *inext;
        const auto& val1 = *i;
        static_assert(std::tuple_size_v<std::decay_t<decltype(val0)>> == 2,
                      "need to generalize comparer to handle tuples");
        using lhs_type = decltype(std::get<0>(val0));
        using rhs_type = decltype(std::get<1>(val0));
        constexpr const bool comp_works_for_tuple_entries =
            meta::is_detected_v<suitable_call_operator, Compare, lhs_type,
                                rhs_type>;
        static_assert(comp_works_for_tuple_entries,
                      "Provided comparator not suitable for entries in "
                      "tuple-like objects (in zipped range?)");

        auto composite_compare = [&comp](auto&& c0, auto&& c1) {
          if (comp(std::get<0>(c0), std::get<0>(c1))) {  // c0[0] < c1[0]
            return true;
          } else if (!comp(std::get<0>(c1),
                           std::get<0>(c0))) {  // c0[0] == c1[0]
            return comp(std::get<1>(c0), std::get<1>(c1));
          } else {  // c0[0] > c1[0]
            return false;
          }
        };
        auto composite_swap = [](auto& c0, auto& c1) {
          using std::swap;
          swap(std::get<0>(c0), std::get<0>(c1));
          swap(std::get<1>(c0), std::get<1>(c1));
        };
        if (composite_compare(val0, val1)) {
          composite_swap(val1, val0);
          swapped = true;
        }
      }
    }
  } while (swapped);
}

///
/// Just like std::next_permutation but capture the parity (even/odd-ness)
/// of the next permutation in the parameter @c parity.
///
/// Set the initial value of @c parity to an odd number only if the initial
/// sequence [first, last) is an odd permutation of the lexicographically
/// first sequence constructible from [first, last).
/// @see
/// https://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order
///
template <typename BidirIt,
          typename Comp = std::less<decltype(*std::declval<BidirIt>())>>
bool next_permutation_parity(int& parity, BidirIt first, BidirIt last,
                             Comp const& comp = {}) {
  BidirIt i = last;
  if (first == last || first == --i) return false;

  for (;;) {
    BidirIt ii = i;
    --i;
    if (comp(*i, *ii)) {
      BidirIt j = last;

      while (!comp(*i, *(--j))) {
        // do nothing here
      }

      int p = parity + 1 /* for the single iter_swap */
              + std::distance(ii, last) / 2;
      parity = p % 2;
      std::iter_swap(i, j);
      std::reverse(ii, last);
      return true;
    }
    if (i == first) {
      std::reverse(first, last);
      int p = parity + std::distance(first, last) / 2;
      parity = p % 2;
      return false;
    }
  }
}

///
/// All elements in the vector belong to the integral range [-1,N)
/// where N is the length of the [Expr] (ie. the iterable of expressions)
///   * only applicable for binary evaluations
///   * the integer -1 can appear in certain places: it implies the binary
///     operation between the last two expressions
///   * eg.
///         * {0,1,-1,2,-1} => ( (e[0], e[1]), e[2])
///         * {0,1,-1,2,3,-1,-1} => ((e[0], e[1]), (e[2],e[3]))
///
using EvalSequence = container::svector<int>;

///
/// Given a range [a_0, a_1, ..., a_n] returns a range of sliced views of
/// lengths 1 through n: [[a_0], [a_0, a_1], ..., [a_0, a_1, ..., a_n]].
///
template <typename Rng>
auto inits(Rng const& rng) {
  using ranges::views::slice;
  using ranges::views::transform;
  return rng | transform([n = 0, &rng](auto&&) mutable {
           return slice(rng, 0, ++n);
         });
}

namespace bits {

///
/// @brief Generates all ordered bipartitions of the set bits in the mask.
///
/// For a given mask representing a set S, this function returns a view of pairs
/// (A, B) such that A ∪ B = S, A ∩ B = {}, and A, B are subsets of S. The
/// iteration order results in B descending and A ascending.
///
/// @tparam T The unsigned integer type of the mask.
/// @param mask The bitmask representing the set.
/// @return A view of pairs of disjoint submasks.
///
template <std::unsigned_integral T>
constexpr auto bipartitions_ordered(T mask) {
  // number of possible bipartitions
  size_t const nparts = (1 << std::popcount(mask));
  return std::views::iota(size_t{0}, nparts) |
         std::views::transform([mask, p = mask](auto) mutable {
           auto p_ = p;
           p = (p - 1) & mask;
           return std::pair<T, T>{p_ ^ mask, p_};
         });
}

///
/// @brief Generates unordered bipartitions of the mask.
///
/// This function returns a subset of the ordered bipartitions such that for
/// every pair {A, B}, the pair {B, A} is excluded (unless A == B, which is
/// handled naturally). Specifically, it takes the first half of the ordered
/// bipartitions.
///
/// @param mask The bitmask representing the set.
/// @return A view of unordered bipartitions.
///
constexpr auto bipartitions_unordered(std::unsigned_integral auto mask) {
  auto bps = bipartitions_ordered(mask);
  size_t const nparts = std::ranges::distance(bps) / 2;
  return bps | std::views::take(nparts);
}

///
/// @brief Alias for bipartitions_unordered.
///
/// @param mask The bitmask representing the set.
/// @return A view of unordered bipartitions.
///
constexpr auto bipartitions(std::unsigned_integral auto mask) {
  return bipartitions_unordered(mask);
}

///
/// @brief Generates all subsets of the mask in ascending numerical order.
///
/// @param mask The bitmask representing the set.
/// @return A view of submasks sorted ascending.
///
constexpr auto subsets_ascending(std::unsigned_integral auto mask) {
  return bipartitions_ordered(mask) | std::views::elements<0>;
}

///
/// @brief Generates all subsets of the mask in descending numerical order.
///
/// @param mask The bitmask representing the set.
/// @return A view of submasks sorted descending.
///
constexpr auto subsets_descending(std::unsigned_integral auto mask) {
  return bipartitions_ordered(mask) | std::views::elements<1>;
}

///
/// @brief Generates all subsets of the mask in ascending numerical order.
///
/// @param mask The bitmask representing the set.
/// @return A view of submasks.
///
constexpr auto subsets(std::unsigned_integral auto mask) {
  return subsets_ascending(mask);
}

///
/// @brief Generates indices of set bits in the mask.
/// @param mask The bitmask representing the set.
/// @return A view of indices of set bits.
///
constexpr auto on_bits_index(std::unsigned_integral auto mask) {
  size_t const nbits = std::popcount(mask);
  return std::views::iota(size_t{0}, nbits) |
         std::views::transform([m = mask](auto) mutable {
           size_t i = std::countr_zero(m);
           m &= (m - 1);
           return i;
         });
}

///
/// @brief Creates a view adaptor that selects elements from the input range
/// based on indices.
/// This function returns a closure object that, when applied to a range of
/// indices, maps those indices to elements in the source range `rng`.
/// @param rng The source range from which elements are selected.
/// @return A view adaptor that transforms indices into elements from `rng`.
///
constexpr auto sieve(std::ranges::viewable_range auto&& rng) {
  auto elems = std::views::all(std::forward<decltype(rng)>(rng));
  return std::views::transform([elems](auto i) {
    return *std::ranges::next(std::ranges::begin(elems), i);
  });
}

}  // namespace bits

}  // namespace sequant

#endif  // SEQUANT_ALGORITHM_HPP
