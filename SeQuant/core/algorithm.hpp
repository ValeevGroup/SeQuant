//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_ALGORITHM_HPP
#define SEQUANT_ALGORITHM_HPP

#include <type_traits>
#include <tuple>
#include <functional>
#include <iterator>
#include <algorithm>

namespace sequant {

/// @brief bubble sort that uses swap exclusively
template <typename ForwardIter, typename Sentinel, typename Compare>
void bubble_sort(ForwardIter begin, Sentinel end, Compare comp) {
  if (begin == end) return;  // no-op for empty range
  bool swapped;
  do {
    swapped = false;
    auto i = begin;
    auto inext = i;
    // iterators either dereference to a reference or to a composite of references
    constexpr const bool iter_deref_to_ref = std::is_reference_v<decltype(*(std::declval<ForwardIter>()))>;
    for (++inext; inext != end; ++i, ++inext) {
      if constexpr (iter_deref_to_ref) {
        auto& val0 = *inext;
        auto& val1 = *i;
        if (comp(val0, val1)) {
          using std::swap;
          swap(val1, val0);
          swapped = true;
        }
      }
      else {
        auto val0 = *inext;
        auto val1 = *i;
        static_assert(std::tuple_size_v<decltype(val0)> == 2,
                      "need to generalize comparer to handle tuples");
        auto composite_compare = [](auto&& c0, auto&& c1) {
          if (std::get<0>(c0) < std::get<0>(c1)) {  // c0[0] < c1[1]
            return true;
          } else if (!(std::get<0>(c1) <
                       std::get<0>(c0))) {  // c0[0] == c1[0]
            return std::get<1>(c0) < std::get<1>(c1);
          }
          else {  // c0[0] > c1[0]
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
/// @see https://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order
///
template <typename BidirIt,
          typename Comp = std::less<decltype(*std::declval<BidirIt>())>>
bool next_permutation_parity(int& parity, BidirIt first, BidirIt last,
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
