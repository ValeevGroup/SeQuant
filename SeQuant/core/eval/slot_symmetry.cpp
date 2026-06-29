#include <SeQuant/core/eval/slot_symmetry.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>

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

namespace {

/// Where an index sits in an operand tensor: which bundle and column.
struct SlotLoc {
  enum class Bundle { Bra, Ket } bundle;
  std::size_t column;   ///< position within the bra (or ket) bundle
  bool column_grouped;  ///< true iff this column lies in an operand ColumnGroup
};

/// Map index-label -> SlotLoc for the matched columns of an operand tensor.
/// Only bra/ket slots over min(bra_rank, ket_rank) are columns; aux and
/// unpaired bra/ket slots are excluded.
std::unordered_map<std::wstring, SlotLoc> column_locations(
    Tensor const& t, SlotSymmetry const& sym) {
  auto in_column_group = [&sym](std::size_t col) {
    for (auto const& cg : sym.column_groups)
      if (std::find(cg.cols.begin(), cg.cols.end(), col) != cg.cols.end())
        return true;
    return false;
  };

  std::unordered_map<std::wstring, SlotLoc> locs;
  const std::size_t ncols = std::min(t.bra_rank(), t.ket_rank());
  auto const& bra = t.bra();
  auto const& ket = t.ket();
  for (std::size_t c = 0; c < ncols; ++c) {
    const bool grouped = in_column_group(c);
    if (bra[c].nonnull())
      locs.emplace(std::wstring{bra[c].label()},
                   SlotLoc{SlotLoc::Bundle::Bra, c, grouped});
    if (ket[c].nonnull())
      locs.emplace(std::wstring{ket[c].label()},
                   SlotLoc{SlotLoc::Bundle::Ket, c, grouped});
  }
  return locs;
}

}  // namespace

SlotSymmetry deduce_slot_symmetry(EvalExpr const& left, EvalExpr const& right,
                                  Tensor const& result) {
  SlotSymmetry ss;

  // Only tensor*tensor contributes column inheritance here.
  if (!left.is_tensor() || !right.is_tensor()) return ss;

  Tensor const& lt = left.as_tensor();
  Tensor const& rt = right.as_tensor();
  auto lloc = column_locations(lt, left.slot_symmetry());
  auto rloc = column_locations(rt, right.slot_symmetry());

  // Trace a result index to the operand (0 = left, 1 = right) and slot that
  // supplies it. Externals appear in exactly one operand slot.
  auto trace = [&](Index const& idx) -> std::optional<std::pair<int, SlotLoc>> {
    const std::wstring key{idx.label()};
    if (auto it = lloc.find(key); it != lloc.end()) return {{0, it->second}};
    if (auto it = rloc.find(key); it != rloc.end()) return {{1, it->second}};
    return std::nullopt;
  };

  const std::size_t ncols = std::min(result.bra_rank(), result.ket_rank());
  if (ncols < 2) return ss;

  auto const& rbra = result.bra();
  auto const& rket = result.ket();

  // Identify the single bra-supplier and ket-supplier operands, and require
  // every result column's bra/ket index to trace into that operand's
  // column-grouped slots.
  std::optional<int> bra_supplier, ket_supplier;
  bool all_columns_inherit = true;
  for (std::size_t c = 0; c < ncols && all_columns_inherit; ++c) {
    if (!rbra[c].nonnull() || !rket[c].nonnull()) {
      all_columns_inherit = false;
      break;
    }
    auto b = trace(rbra[c]);
    auto k = trace(rket[c]);
    if (!b || !k || !b->second.column_grouped || !k->second.column_grouped) {
      all_columns_inherit = false;
      break;
    }
    if (!bra_supplier)
      bra_supplier = b->first;
    else if (*bra_supplier != b->first)
      all_columns_inherit = false;
    if (!ket_supplier)
      ket_supplier = k->first;
    else if (*ket_supplier != k->first)
      all_columns_inherit = false;
  }

  if (all_columns_inherit && bra_supplier && ket_supplier) {
    // The contraction between the two supplying operands is symmetric: each
    // supplier carries a full column group over its supplying columns, so the
    // contracted indices (the other bundle of each supplier) sit in matched
    // grouped columns -- the PPL / giant matched-pair-swap pattern. (When
    // bra_supplier == ket_supplier the whole column comes from one operand's
    // ColumnGroup, the leaf-passthrough case.)
    SlotSymmetry::ColumnGroup cg;
    cg.sign = 1;
    cg.cols.reserve(ncols);
    for (std::size_t c = 0; c < ncols; ++c) cg.cols.push_back(c);
    ss.column_groups.push_back(std::move(cg));
  }

  return ss;
}

}  // namespace sequant
