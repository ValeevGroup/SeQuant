#ifndef SEQUANT_PERMUTATION_HPP
#define SEQUANT_PERMUTATION_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find.hpp>
#include <range/v3/algorithm/for_each.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <ranges>

namespace sequant {

/// @brief Returns the number of (spin) loops in a closed-shell trace

/// Interprets \p v0 and \p v1 as the column-paired index slots of a contracted
/// tensor network: slot @c i is one particle slot of one tensor, holding @c
/// v0[i] in one row and @c v1[i] in the other (e.g. ket and bra). Two slots
/// must carry the same spin if they are (a) the two rows of the same column (a
/// *column edge* @c v0[i]—v1[i]) or (b) the two occurrences of the same
/// (contracted) index value (a *contraction edge*). Every internal slot then
/// has degree two, so this graph is a disjoint union of cycles and the number
/// of connected components is the number of independent spin loops (each loop
/// contributes a factor of two in the closed-shell trace).
///
/// Unlike a 2-line-permutation reading, this is agnostic to which row a value
/// sits in: a contraction may be covariant (one occurrence in \p v0, one in
/// \p v1) or, for bra-ket-symmetric (Hermitian/real) tensors that have been
/// reoriented, both occurrences in the same row (bra-bra / ket-ket). The
/// component count is invariant under such bra<->ket swaps and reduces to the
/// permutation's cycle count whenever \p v0 is a permutation of \p v1.
/// @tparam Seq0 (reference to) a container type
/// @tparam Seq1 (reference to) a container type
/// @param v0 first row of index slots
/// @param v1 second row of index slots (\p v1[i] is column-paired with \p
/// v0[i])
/// @pre \p v0 and \p v1 have equal size and every value occurs in exactly two
///      slots across \p v0 and \p v1 combined
/// @return the number of connected components (spin loops)
template <typename Seq0, typename Seq1>
std::size_t count_cycles(Seq0&& v0, Seq1&& v1) {
  using std::ranges::begin;
  using std::ranges::end;
  using std::ranges::size;
  const std::size_t n0 = size(v0);
  const std::size_t n1 = size(v1);
  const std::size_t n = n0 + n1;

  // precondition (debug builds only): every value occurs in exactly two slots
  // across v0 and v1 combined. This is what makes the column+contraction graph
  // 2-regular on its internal slots, so it decomposes into disjoint cycles and
  // the component count is the loop count. It generalizes the former "v0 is a
  // permutation of v1" contract (which also implied exactly two occurrences,
  // one per row) to allow both occurrences in the same row (bra-bra / ket-ket
  // edges from reoriented bra-ket-symmetric tensors). Without this check
  // malformed input (a value appearing once, or 3+ times) is silently accepted
  // and returns a meaningless count.
  if constexpr (assert_enabled()) {
    container::map<std::ranges::range_value_t<Seq0>, std::size_t> counts;
    for (auto&& x : v0) ++counts[x];
    for (auto&& x : v1) ++counts[x];
    SEQUANT_ASSERT(
        ranges::all_of(counts, [](auto const& kv) { return kv.second == 2; }));
  }

  // slot ids: row-0 slot i -> i ; row-1 slot i -> n0 + i. Size the union-find
  // by n0 + n1 (not 2*n0) so it is safe even if the two rows differ in length.
  // union-find with path halving
  container::svector<std::size_t> parent(n);
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  auto find = [&parent](std::size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  auto unite = [&](std::size_t a, std::size_t b) { parent[find(a)] = find(b); };

  // column edges: the two rows of each column carry the same spin
  const std::size_t ncols = std::min(n0, n1);
  for (std::size_t i = 0; i < ncols; ++i) unite(i, n0 + i);

  // contraction edges: the slots sharing an index value carry the same spin
  container::map<std::ranges::range_value_t<Seq0>, std::size_t> first_slot;
  auto add_slot = [&first_slot, &unite](auto const& idx, std::size_t slot) {
    auto [it, inserted] = first_slot.try_emplace(idx, slot);
    if (!inserted)
      unite(it->second, slot);  // contraction with first occurrence
  };
  {
    std::size_t i = 0;
    for (auto it = begin(v0); it != end(v0); ++it, ++i) add_slot(*it, i);
    i = 0;
    for (auto it = begin(v1); it != end(v1); ++it, ++i) add_slot(*it, n0 + i);
  }

  // number of connected components = number of spin loops
  std::size_t n_cycles = 0;
  for (std::size_t i = 0; i < n; ++i)
    if (find(i) == i) ++n_cycles;
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
