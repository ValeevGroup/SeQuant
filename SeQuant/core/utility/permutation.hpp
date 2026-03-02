#ifndef SEQUANT_PERMUTATION_HPP
#define SEQUANT_PERMUTATION_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <vector>

namespace sequant {

/// @brief Returns the number of cycles

/// Counts the number of cycles of a permutation represented in a 2-line form
/// by stacking \p v0 and \p v1 on top of each other.
/// @tparam Seq0 a container type
/// @tparam Seq1 a container type
/// @param[in] v0 first sequence
/// @param[in] v1 second sequence
/// @pre \p v0 is a permutation of \p v1
/// @return the number of cycles
template <typename Seq0, typename Seq1>
std::size_t count_cycles(const Seq0& v0, const Seq1& v1) {
  SEQUANT_ASSERT(ranges::is_permutation(v0, v1));
  // This function can't deal with duplicate entries in v0 or v1
  SEQUANT_ASSERT(std::set(std::begin(v0), std::end(v0)).size() == v0.size());
  SEQUANT_ASSERT(std::set(std::begin(v1), std::end(v1)).size() == v1.size());

  const auto n = v0.size();
  std::vector<bool> visited(n, false);
  std::size_t n_cycles = 0;

  for (std::size_t i = 0; i < n; ++i) {
    if (!visited[i]) {
      ++n_cycles;
      auto j = i;
      do {
        visited[j] = true;
        auto it = std::find(std::begin(v1), std::end(v1), v0[j]);
        SEQUANT_ASSERT(it != std::end(v1));
        j = static_cast<std::size_t>(std::distance(std::begin(v1), it));
      } while (j != i);
    }
  }

  return n_cycles;
}

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
