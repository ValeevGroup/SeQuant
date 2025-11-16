#ifndef SEQUANT_PERMUTATION_HPP
#define SEQUANT_PERMUTATION_HPP

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm.hpp>

#include <algorithm>
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
  SEQUANT_ASSERT(ranges::is_permutation(v, v1));
  // This function can't deal with duplicate entries in v0 or v1
  SEQUANT_ASSERT(std::set(std::begin(v0), std::end(v0)).size() == v0.size());
  SEQUANT_ASSERT(std::set(std::begin(v1), std::end(v1)).size() == v1.size());

  auto make_null = []() -> T {
    if constexpr (std::is_arithmetic_v<T>) {
      return -1;
    } else if constexpr (std::is_same_v<T, Index>) {
      return L"p_50";
    }

    SEQUANT_UNREACHABLE;
  };

  const auto null = make_null();
  SEQUANT_ASSERT(ranges::contains(v, null) == false);
  SEQUANT_ASSERT(ranges::contains(v1, null) == false);

  std::size_t n_cycles = 0;
  for (auto it = v.begin(); it != v.end(); ++it) {
    if (*it != null) {
      n_cycles++;

      auto idx = std::distance(v.begin(), it);
      SEQUANT_ASSERT(idx >= 0);

      auto it0 = it;

      auto it1 = std::find(v1.begin(), v1.end(), *it0);
      SEQUANT_ASSERT(it1 != v1.end());

      auto idx1 = std::distance(v1.begin(), it1);
      SEQUANT_ASSERT(idx1 >= 0);

      do {
        it0 = std::find(v.begin(), v.end(), v[idx1]);
        SEQUANT_ASSERT(it0 != v.end());

        it1 = std::find(v1.begin(), v1.end(), *it0);
        SEQUANT_ASSERT(it1 != v1.end());

        idx1 = std::distance(v1.begin(), it1);
        SEQUANT_ASSERT(idx1 >= 0);

        *it0 = null;
      } while (idx1 != idx);
    }
  }
  return n_cycles;
};

/// computes parity of a permutation of 0 ... N-1
///
/// @param p permutation
/// @param overwrite if true, will overwrite @p p
template <std::integral T>
int permutation_parity(std::span<T> p, bool overwrite = false) {
  // https://stackoverflow.com/a/20703469
  // compute cycles, mutating elements of p to mark used elements
  const std::size_t N = p.size();
  int parity = 1;
  // search for next element to start cycle with
  for (std::size_t k = 0; k != N; ++k) {
    if (p[k] >= N) continue;
    std::size_t i = k;
    std::size_t cycle_length = 1;
    do {
      i = p[i];
      p[i] += N;
      ++cycle_length;
    } while (p[i] < N);
    if (cycle_length % 2 == 0) parity *= -1;
  }

  if (overwrite) {
    ranges::for_each(p, [N](auto& e) { e -= N; });
  }

  return parity;
}

}  // namespace sequant

#endif
