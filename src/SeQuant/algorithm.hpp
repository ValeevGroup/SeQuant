//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT_ALGORITHM_HPP
#define SEQUANT_ALGORITHM_HPP

namespace sequant {

/// @brief bubble sort that uses swap exclusively
template <typename ForwardIter, typename Compare>
void bubble_sort(ForwardIter begin, ForwardIter end, Compare comp) {
  if (begin == end) return;  // no-op for empty range
  bool swapped;
  do {
    swapped = false;
    for (auto i = begin, inext = std::next(begin); inext != end; ++i, ++inext) {
      auto& val0 = *inext;
      auto& val1 = *i;
      if (comp(val0, val1)) {
        using std::swap;
        swap(val1, val0);
        swapped = true;
      }
    }
  } while (swapped);
}

}  // namespace sequant

#endif  // SEQUANT_ALGORITHM_HPP
