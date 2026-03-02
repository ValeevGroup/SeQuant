#ifndef SEQUANT_PERMUTATION_HPP
#define SEQUANT_PERMUTATION_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <ranges>
#include <set>

namespace sequant {

/// @brief Returns the number of cycles

/// Counts the number of cycles of a permutation represented in a 2-line form
/// by stacking \p v0 and \p v1 on top of each other.
/// @tparam Seq0 (reference to) a container type
/// @tparam Seq1 (reference to) a container type
/// @param v0 first range
/// @param v1 second range
/// @pre \p v0 is a permutation of \p v1
/// @return the number of cycles
template <typename Seq0, typename Seq1>
std::size_t count_cycles(Seq0&& v0, Seq1&& v1) {
  using std::ranges::begin;
  using std::ranges::end;
  using std::ranges::size;
  SEQUANT_ASSERT(std::ranges::is_permutation(v0, v1));
  // This function can't deal with duplicate entries in v0 or v1
  SEQUANT_ASSERT(std::set(begin(v0), end(v0)).size() == size(v0));
  SEQUANT_ASSERT(std::set(begin(v1), end(v1)).size() == size(v1));

  container::set<decltype(begin(v0))> visited;
  visited.reserve(size(v0));

  std::size_t n_cycles = 0;
  for (auto it = begin(v0); it != end(v0); ++it) {
    if (visited.find(it) != end(visited)) {
      continue;
    }

    n_cycles++;

    auto idx = std::distance(begin(v0), it);
    SEQUANT_ASSERT(idx >= 0);

    auto it0 = it;

    auto it1 = std::find(begin(v1), end(v1), *it0);
    SEQUANT_ASSERT(it1 != end(v1));

    auto idx1 = std::distance(begin(v1), it1);
    SEQUANT_ASSERT(idx1 >= 0);

    do {
      auto val_it = begin(v0);
      std::advance(val_it, idx1);
      SEQUANT_ASSERT(val_it != end(v0));

      it0 = std::find(begin(v0), end(v0), *val_it);
      SEQUANT_ASSERT(it0 != end(v0));

      it1 = std::find(begin(v1), end(v1), *it0);
      SEQUANT_ASSERT(it1 != end(v1));

      idx1 = std::distance(begin(v1), it1);
      SEQUANT_ASSERT(idx1 >= 0);

      SEQUANT_ASSERT(visited.find(it0) == end(visited));
      visited.insert(it0);
    } while (idx1 != idx);
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
    std::ranges::for_each(p, [N](auto& e) { e -= N; });
  }

  return parity;
}

}  // namespace sequant

#endif
