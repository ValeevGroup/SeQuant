//
// Created by Eduard Valeyev on 2019-02-19.
//

#ifndef SEQUANT2_ALGORITHM_HPP
#define SEQUANT2_ALGORITHM_HPP

namespace sequant2 {

/// @brief bubble sort that uses swap exclusively
template <typename ForwardIter, typename Compare>
void bubble_sort(ForwardIter begin, ForwardIter end, Compare comp) {
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

}  // namespace sequant2

#endif  // SEQUANT2_ALGORITHM_HPP
