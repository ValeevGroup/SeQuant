//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_ALGORITHM_HPP
#define SEQUANT_ALGORITHM_HPP

#include <SeQuant/core/meta.hpp>

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
template <typename ForwardIter, typename Sentinel, typename Compare>
void bubble_sort(ForwardIter begin, Sentinel end, Compare comp) {
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

}  // namespace sequant

#endif  // SEQUANT_ALGORITHM_HPP
