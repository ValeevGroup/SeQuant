#include <SeQuant/core/export/reordering_context.hpp>

#include <algorithm>
#include <cassert>
#include <ranges>

namespace sequant {

struct IndexSpaceSizeComparer {
  IndexSpaceSizeComparer(MemoryLayout layout) : layout(layout) {}

  bool operator()(const Index &lhs, const Index &rhs) const {
    return (*this)(lhs.space(), rhs.space());
  }

  // If index spaces are sorted based on this comparison, the largest index
  // space will end up in the slot with greatest cache-locality
  bool operator()(const IndexSpace &lhs, const IndexSpace &rhs) const {
    switch (layout) {
      case MemoryLayout::RowMajor:
        // Fastest index is right-most in [abcd]
        return lhs.approximate_size() < rhs.approximate_size();
      case MemoryLayout::ColumnMajor:
        // Fastest index is left-most in [abcd]
        return lhs.approximate_size() > rhs.approximate_size();
      case MemoryLayout::Unspecified:
        return false;
    }

    assert(false);
    return false;
  }

  MemoryLayout layout;
};

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

  assert(false);
  return {};
}

bool ReorderingContext::rewrite(Tensor &tensor) const {
  using std::ranges::begin;
  using std::ranges::end;

  const IndexSpaceSizeComparer cmp(m_layout);
  const auto needs_swap = [&](const IndexSpace &lhs, const IndexSpace &rhs) {
    return lhs != rhs && !cmp(lhs, rhs);
  };

  // Note: In theory we could apply the same logic to antisymmetric index
  // groups but reordering those might incur a sign change, which we can't
  // represent by simply changing the Tensor object
  bool sort_bra = tensor.symmetry() == Symmetry::symm &&
                  !std::ranges::is_sorted(tensor.bra(), cmp);
  bool sort_ket = tensor.symmetry() == Symmetry::symm &&
                  !std::ranges::is_sorted(tensor.ket(), cmp);

  IndexSpace relevant_bra = sort_bra ? max_space_size_idx(tensor.bra())
                                     : relevant_space(tensor.bra(), m_layout);
  IndexSpace relevant_ket = sort_ket ? max_space_size_idx(tensor.ket())
                                     : relevant_space(tensor.ket(), m_layout);

  // Same as for antisymmetry: We can't deal with conjugation at the tensor
  // level
  bool swap_braket = !tensor.bra().empty() && !tensor.ket().empty() &&
                     (tensor.braket_symmetry() == BraKetSymmetry::symm ||
                      (m_real_orbitals && tensor.braket_symmetry() ==
                                              BraKetSymmetry::conjugate)) &&
                     needs_swap(relevant_bra, relevant_ket);
  bool prioritize_aux = !tensor.aux().empty() &&
                        needs_swap((swap_braket ? relevant_ket : relevant_bra),
                                   relevant_space(tensor.aux(), m_layout));

  if (!sort_bra && !sort_ket && !swap_braket && !prioritize_aux) {
    return false;
  }

  container::svector<Index> indices;
  indices.reserve(tensor.const_indices().size());

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

  if (sort_bra) {
    if (swap_braket) {
      std::sort(indices.begin() + second_begin_idx, indices.end(), cmp);
    } else {
      std::sort(indices.begin() + first_begin_idx,
                indices.begin() + second_begin_idx, cmp);
    }
  }
  if (sort_ket) {
    if (swap_braket) {
      std::sort(indices.begin() + first_begin_idx,
                indices.begin() + second_begin_idx, cmp);
    } else {
      std::sort(indices.begin() + second_begin_idx, indices.end(), cmp);
    }
  }

  if (!prioritize_aux) {
    indices.insert(end(indices), begin(tensor.aux()), end(tensor.aux()));
  }

  tensor = Tensor(tensor.label(), bra(), ket(), aux(std::move(indices)));

  return true;
}

}  // namespace sequant
