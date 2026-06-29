#include <SeQuant/core/eval/slot_symmetry.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>

#include <algorithm>
#include <cstddef>

namespace sequant {

SlotSymmetry from_leaf_tensor(Tensor const& t) {
  SlotSymmetry ss;

  const std::size_t bra_rank = t.bra_rank();
  const std::size_t ket_rank = t.ket_rank();

  // Column group: invariance under permuting matched (bra[c], ket[c]) columns.
  // Only the paired columns over min(bra_rank, ket_rank) form columns; unpaired
  // bra/ket slots and aux are never part of a ColumnGroup (spec 1.1 R2).
  if (t.column_symmetry() == ColumnSymmetry::Symm) {
    const std::size_t ncols = std::min(bra_rank, ket_rank);
    if (ncols >= 2) {
      SlotSymmetry::ColumnGroup cg;
      cg.sign = 1;
      cg.cols.reserve(ncols);
      for (std::size_t c = 0; c < ncols; ++c) cg.cols.push_back(c);
      ss.column_groups.push_back(std::move(cg));
    }
  }

  // Within-bundle (bra-only / ket-only) permutational (anti)symmetry. SeQuant's
  // Symmetry attribute applies jointly to the bra and the ket bundles, so emit
  // a group for each bundle of rank >= 2.
  if (t.symmetry() == Symmetry::Symm || t.symmetry() == Symmetry::Antisymm) {
    const std::int8_t sign = (t.symmetry() == Symmetry::Antisymm) ? -1 : 1;

    if (bra_rank >= 2) {
      SlotSymmetry::SlotGroup bg;
      bg.sign = sign;
      bg.slots.reserve(bra_rank);
      for (std::size_t s = 0; s < bra_rank; ++s) bg.slots.push_back(s);
      ss.bra_groups.push_back(std::move(bg));
    }
    if (ket_rank >= 2) {
      SlotSymmetry::SlotGroup kg;
      kg.sign = sign;
      kg.slots.reserve(ket_rank);
      for (std::size_t s = 0; s < ket_rank; ++s) kg.slots.push_back(s);
      ss.ket_groups.push_back(std::move(kg));
    }
  }

  return ss;
}

}  // namespace sequant
