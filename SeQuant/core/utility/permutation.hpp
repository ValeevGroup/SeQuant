#ifndef SEQUANT_PERMUTATION_HPP
#define SEQUANT_PERMUTATION_HPP

#include <SeQuant/core/index.hpp>

#include <range/v3/algorithm.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <type_traits>
#include <utility>

namespace sequant {

/// @brief Returns the number of cycles

/// Counts the number of cycles of a permutation represented in a 2-line form
/// by stacking \p v0 and \p v1 on top of each other.
/// @tparam Seq0 (reference to) a container type
/// @tparam Seq1 (reference to) a container type
/// @param v0 first sequence; if passed as an rvalue reference, it is moved from
/// @param[in] v1 second sequence
/// @pre \p v0 is a permutation of \p v1
/// @return the number of cycles
template <typename Seq0, typename Seq1>
std::size_t count_cycles(Seq0&& v0, const Seq1& v1) {
  std::remove_reference_t<Seq0> v(std::forward<Seq0>(v0));
  using T = std::decay_t<decltype(v[0])>;
  assert(ranges::is_permutation(v, v1));
  // This function can't deal with duplicate entries in v0 or v1
  assert(std::set(std::begin(v0), std::end(v0)).size() == v0.size());
  assert(std::set(std::begin(v1), std::end(v1)).size() == v1.size());

  auto make_null = []() -> T {
    if constexpr (std::is_arithmetic_v<T>) {
      return -1;
    } else if constexpr (std::is_same_v<T, Index>) {
      return L"p_50";
    } else  // unreachable
      abort();
  };

  const auto null = make_null();
  assert(ranges::contains(v, null) == false);
  assert(ranges::contains(v1, null) == false);

  std::size_t n_cycles = 0;
  for (auto it = v.begin(); it != v.end(); ++it) {
    if (*it != null) {
      n_cycles++;

      auto idx = std::distance(v.begin(), it);
      assert(idx >= 0);

      auto it0 = it;

      auto it1 = std::find(v1.begin(), v1.end(), *it0);
      assert(it1 != v1.end());

      auto idx1 = std::distance(v1.begin(), it1);
      assert(idx1 >= 0);

      do {
        it0 = std::find(v.begin(), v.end(), v[idx1]);
        assert(it0 != v.end());

        it1 = std::find(v1.begin(), v1.end(), *it0);
        assert(it1 != v1.end());

        idx1 = std::distance(v1.begin(), it1);
        assert(idx1 >= 0);

        *it0 = null;
      } while (idx1 != idx);
    }
  }
  return n_cycles;
};

}  // namespace sequant

#endif
