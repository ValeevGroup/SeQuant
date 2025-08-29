#include <SeQuant/core/export/reordering_context.hpp>

#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cassert>
#include <ranges>

#include <range/v3/algorithm/sort.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/zip.hpp>

namespace sequant {

template <typename Range>
  requires std::ranges::range<Range> &&
           std::is_same_v<std::ranges::range_value_t<Range>, Index>
IndexSpace max_space_size_idx(Range &&rng) {
  IndexSpace max;

  for (const Index &idx : rng) {
    if (idx.space().approximate_size() > max.approximate_size()) {
      max = idx.space();
    }
  }

  return max;
}

template <typename Range>
  requires std::ranges::range<Range> &&
           std::is_same_v<std::ranges::range_value_t<Range>, Index>
IndexSpace relevant_space(Range &&rng, MemoryLayout layout) {
  if (std::ranges::empty(rng)) {
    return {};
  }

  switch (layout) {
    case MemoryLayout::RowMajor:
      // Get last element
      return (--std::ranges::end(rng))->space();
    case MemoryLayout::ColumnMajor:
      // Get first element
      return std::ranges::begin(rng)->space();
    case MemoryLayout::Unspecified:
      return {};
  }

  SEQUANT_UNREACHABLE;
}

template <typename Comp>
struct PairComparator {
  const Comp &cmp;
  bool prioritize_first;

  PairComparator(const Comp &cmp, bool prioritize_first)
      : cmp(cmp), prioritize_first(prioritize_first) {}

  bool operator()(const auto &lhs, const auto &rhs) const {
    const auto &[lhs_first, lhs_second] = lhs;
    const auto &[rhs_first, rhs_second] = rhs;

    if (prioritize_first) {
      if (cmp(lhs_first, rhs_first)) {
        return true;
      } else if (!cmp(rhs_first, lhs_first)) {
        // lhs_first and rhs_first considered equal -> check second
        if (cmp(lhs_second, rhs_second)) {
          return true;
        } else if (cmp(rhs_second, lhs_second)) {
          return false;
        }
      }

      // The pairings are equal
      return false;
    }

    if (cmp(lhs_second, rhs_second)) {
      return true;
    } else if (!cmp(rhs_second, lhs_second)) {
      // lhs_second and rhs_second considered equal -> check second
      if (cmp(lhs_first, rhs_first)) {
        return true;
      } else if (cmp(rhs_first, lhs_first)) {
        return false;
      }
    }

    // The pairings are equal
    return false;
  }
};

bool ReorderingContext::rewrite(Tensor &tensor) const {
  using std::ranges::begin;
  using std::ranges::end;

  auto comparator = [this](const Index &lhs, const Index &rhs) {
    return this->is_less(lhs.space(), rhs.space());
  };

  const IndexSpace max_bra_space = max_space_size_idx(tensor.bra());
  const IndexSpace max_ket_space = max_space_size_idx(tensor.ket());
  const bool bra_over_ket =
      max_bra_space.approximate_size() >= max_ket_space.approximate_size();

  const PairComparator pair_comp(comparator, bra_over_ket);

  // Note: In theory we could apply the same logic to antisymmetric index
  // groups but reordering those might incur a sign change, which we can't
  // represent by simply changing the Tensor object
  bool sort_bra = tensor.symmetry() == Symmetry::Symm &&
                  !std::ranges::is_sorted(tensor.bra(), comparator);
  bool sort_ket = tensor.symmetry() == Symmetry::Symm &&
                  !std::ranges::is_sorted(tensor.ket(), comparator);

  // Note: Sorting bra and ket individually is more powerful than doing
  // particle-wise sorting. Hence, if we do the former, we don't have to do the
  // latter
  bool sort_particles =
      !sort_bra && !sort_ket && !tensor.bra().empty() &&
      tensor.bra().size() == tensor.ket().size() &&
      (tensor.column_symmetry() == ColumnSymmetry::Symm ||
       m_column_permutability) &&
      !std::ranges::is_sorted(ranges::views::zip(tensor.bra(), tensor.ket()),
                              pair_comp);

  IndexSpace relevant_bra = sort_bra || (sort_particles && bra_over_ket)
                                ? max_space_size_idx(tensor.bra())
                                : relevant_space(tensor.bra(), m_layout);
  IndexSpace relevant_ket = sort_ket || (sort_particles && !bra_over_ket)
                                ? max_space_size_idx(tensor.ket())
                                : relevant_space(tensor.ket(), m_layout);

  // Same as for antisymmetry: We can't deal with conjugation at the tensor
  // level
  bool swap_braket = !tensor.bra().empty() && !tensor.ket().empty() &&
                     (tensor.braket_symmetry() == BraKetSymmetry::Symm ||
                      (m_real_orbitals && tensor.braket_symmetry() ==
                                              BraKetSymmetry::Conjugate)) &&
                     needs_swap(relevant_bra, relevant_ket);
  bool prioritize_aux = !tensor.aux().empty() &&
                        needs_swap((swap_braket ? relevant_ket : relevant_bra),
                                   relevant_space(tensor.aux(), m_layout));

  if (!sort_bra && !sort_ket && !swap_braket && !prioritize_aux &&
      !sort_particles) {
    return false;
  }

  container::svector<Index> indices;
  indices.reserve(tensor.num_indices());

  if (prioritize_aux) {
    indices.insert(end(indices), begin(tensor.aux()), end(tensor.aux()));
  }

  std::size_t first_begin_idx;
  std::size_t second_begin_idx;
  if (swap_braket) {
    first_begin_idx = indices.size();
    indices.insert(end(indices), begin(tensor.ket()), end(tensor.ket()));
    second_begin_idx = indices.size();
    indices.insert(end(indices), begin(tensor.bra()), end(tensor.bra()));
  } else {
    first_begin_idx = indices.size();
    indices.insert(end(indices), begin(tensor.bra()), end(tensor.bra()));
    second_begin_idx = indices.size();
    indices.insert(end(indices), begin(tensor.ket()), end(tensor.ket()));
  }

  auto bra_range = [&]() {
    if (swap_braket) {
      return ranges::make_subrange(indices.begin() + second_begin_idx,
                                   indices.end());
    }

    return ranges::make_subrange(indices.begin() + first_begin_idx,
                                 indices.begin() + second_begin_idx);
  }();

  auto ket_range = [&]() {
    if (swap_braket) {
      return ranges::make_subrange(indices.begin() + first_begin_idx,
                                   indices.begin() + second_begin_idx);
    }

    return ranges::make_subrange(indices.begin() + second_begin_idx,
                                 indices.end());
  }();

  if (sort_bra) {
    std::sort(bra_range.begin(), bra_range.end(), comparator);
  }
  if (sort_ket) {
    std::sort(ket_range.begin(), ket_range.end(), comparator);
  }

  if (sort_particles) {
    assert(!bra_range.empty());
    assert(bra_range.size() == ket_range.size());

    ranges::sort(ranges::views::zip(bra_range, ket_range), pair_comp);
  }

  if (!prioritize_aux) {
    indices.insert(end(indices), begin(tensor.aux()), end(tensor.aux()));
  }

  tensor = Tensor(tensor.label(), bra(), ket(), aux(std::move(indices)),
                  Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                  ColumnSymmetry::Nonsymm);

  return true;
}

bool ReorderingContext::is_less(const IndexSpace &lhs,
                                const IndexSpace &rhs) const {
  switch (m_layout) {
    case MemoryLayout::RowMajor:
      // Fastest index is right-most in [abcd]
      return lhs.approximate_size() < rhs.approximate_size();
    case MemoryLayout::ColumnMajor:
      // Fastest index is left-most in [abcd]
      return lhs.approximate_size() > rhs.approximate_size();
    case MemoryLayout::Unspecified:
      return false;
  }

  SEQUANT_UNREACHABLE;
}

bool ReorderingContext::is_ordered(const IndexSpace &lhs,
                                   const IndexSpace &rhs) const {
  return is_less(lhs, rhs) || lhs.approximate_size() == rhs.approximate_size();
}

bool ReorderingContext::needs_swap(const IndexSpace &lhs,
                                   const IndexSpace &rhs) const {
  return lhs != rhs && !is_ordered(lhs, rhs);
}

}  // namespace sequant
