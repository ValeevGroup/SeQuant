#include <SeQuant/core/eval/slot_symmetry.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <tuple>
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
  std::size_t column;  ///< position within the bra (or ket) bundle
  /// Index of the operand ColumnGroup this column belongs to, or nullopt if
  /// the column is not part of any ColumnGroup.
  std::optional<std::size_t> column_group_idx;
};

/// Map index full_label -> SlotLoc for the matched columns of an operand
/// tensor. Only bra/ket slots over min(bra_rank, ket_rank) are columns; aux
/// and unpaired bra/ket slots are excluded.
///
/// @note Keys use Index::full_label() (not label()) so that proto-indexed
///       indices with the same base label but distinct proto-indices are
///       stored as separate entries (I2 fix).
std::unordered_map<std::wstring, SlotLoc> column_locations(
    Tensor const& t, SlotSymmetry const& sym) {
  // Return the index of the ColumnGroup in sym that contains col, or nullopt.
  auto which_column_group =
      [&sym](std::size_t col) -> std::optional<std::size_t> {
    for (std::size_t gi = 0; gi < sym.column_groups.size(); ++gi)
      if (std::find(sym.column_groups[gi].cols.begin(),
                    sym.column_groups[gi].cols.end(),
                    col) != sym.column_groups[gi].cols.end())
        return gi;
    return std::nullopt;
  };

  std::unordered_map<std::wstring, SlotLoc> locs;
  const std::size_t ncols = std::min(t.bra_rank(), t.ket_rank());
  auto const& bra = t.bra();
  auto const& ket = t.ket();
  for (std::size_t c = 0; c < ncols; ++c) {
    const auto grp = which_column_group(c);
    if (bra[c].nonnull())
      locs.emplace(std::wstring{bra[c].full_label()},
                   SlotLoc{SlotLoc::Bundle::Bra, c, grp});
    if (ket[c].nonnull())
      locs.emplace(std::wstring{ket[c].full_label()},
                   SlotLoc{SlotLoc::Bundle::Ket, c, grp});
  }
  return locs;
}

}  // namespace

SlotSymmetry adjoint(SlotSymmetry const& s) {
  SlotSymmetry result;
  result.column_groups = s.column_groups;
  result.bra_groups = s.ket_groups;
  result.ket_groups = s.bra_groups;
  return result;
}

SlotSymmetry intersect(SlotSymmetry const& a, SlotSymmetry const& b) {
  SlotSymmetry result;

  // Column groups: keep a group from a iff b contains a group with the same
  // sign and the same set of column positions (order-insensitive).
  auto col_match = [](SlotSymmetry::ColumnGroup const& ga,
                      SlotSymmetry::ColumnGroup const& gb) {
    if (ga.sign != gb.sign) return false;
    auto ac = ga.cols;
    auto bc = gb.cols;
    std::sort(ac.begin(), ac.end());
    std::sort(bc.begin(), bc.end());
    return ac == bc;
  };
  for (auto const& ga : a.column_groups) {
    for (auto const& gb : b.column_groups) {
      if (col_match(ga, gb)) {
        result.column_groups.push_back(ga);
        break;
      }
    }
  }

  // Bra groups: keep a group from a iff b contains a group with the same sign
  // and the same set of slot positions (order-insensitive).
  auto slot_match = [](SlotSymmetry::SlotGroup const& ga,
                       SlotSymmetry::SlotGroup const& gb) {
    if (ga.sign != gb.sign) return false;
    auto as = ga.slots;
    auto bs = gb.slots;
    std::sort(as.begin(), as.end());
    std::sort(bs.begin(), bs.end());
    return as == bs;
  };
  for (auto const& ga : a.bra_groups) {
    for (auto const& gb : b.bra_groups) {
      if (slot_match(ga, gb)) {
        result.bra_groups.push_back(ga);
        break;
      }
    }
  }

  // Ket groups.
  for (auto const& ga : a.ket_groups) {
    for (auto const& gb : b.ket_groups) {
      if (slot_match(ga, gb)) {
        result.ket_groups.push_back(ga);
        break;
      }
    }
  }

  return result;
}

SlotSymmetry deduce_slot_symmetry(EvalExpr const& left, EvalExpr const& right,
                                  Tensor const& result) {
  SlotSymmetry ss;

  // Only tensor*tensor contributes column inheritance here.
  if (!left.is_tensor() || !right.is_tensor()) return ss;

  Tensor const& lt = left.as_tensor();
  Tensor const& rt = right.as_tensor();

  // ---- Soundness guards ----
  // The flat index->slot trace below assumes each external occupies exactly one
  // slot and that the two factors are distinguishable. Two cases violate that;
  // both bail to an empty descriptor (a false negative is safe, a false
  // positive would corrupt a consumer). See design note
  // doc/dev/specs/2026-07-01-symmetry-deduction-via-tn-canonicalization.md.
  //
  // Guard 1: proto-indexed externals. An index that also appears as a proto-
  // index of another slot (e.g. i1 in F{a1<i1,i2>; i1}) couples slots the trace
  // treats as independent, which can yield a WRONG symmetry. If any index of
  // either operand or the result carries proto-indices, decline. (This declines
  // on CSV/PNO intermediates until the graph-canonicalization deducer lands.)
  auto any_proto = [](Tensor const& t) {
    auto has = [](auto const& bundle) {
      for (auto const& idx : bundle)
        if (idx.has_proto_indices()) return true;
      return false;
    };
    return has(t.bra()) || has(t.ket()) || has(t.aux());
  };
  if (any_proto(lt) || any_proto(rt) || any_proto(result)) return ss;

  // NB: repeated identical factors (e.g. A{a1;i1} A{a2;i2}) carry an emergent
  // exchange symmetry the flat rules cannot see, but that is a SAFE false
  // negative -- the deduced group is always a subset of the true symmetry (the
  // 4-tuple supplier key prevents wrong cross-copy merges), so no guard is
  // warranted. Recovering the emergent symmetry needs the
  // graph-canonicalization deducer (Phase 0.5), not a guard.

  auto lloc = column_locations(lt, left.slot_symmetry());
  auto rloc = column_locations(rt, right.slot_symmetry());

  // Trace a result index to the operand (0 = left, 1 = right) and slot that
  // supplies it. Externals appear in exactly one operand slot.
  // Keys use full_label() to avoid aliasing proto-indexed indices that share
  // a base label (I2 fix).
  auto trace = [&](Index const& idx) -> std::optional<std::pair<int, SlotLoc>> {
    const std::wstring key{idx.full_label()};
    if (auto it = lloc.find(key); it != lloc.end()) return {{0, it->second}};
    if (auto it = rloc.find(key); it != rloc.end()) return {{1, it->second}};
    return std::nullopt;
  };

  const std::size_t ncols = std::min(result.bra_rank(), result.ket_rank());
  auto const& rbra = result.bra();
  auto const& rket = result.ket();

  // ---- Column-group inheritance (PPL / giant / n-column / maximal-subset)
  // ----
  if (ncols >= 2) {
    // A result column c inherits iff both rbra[c] and rket[c] are nonnull and
    // both trace to column-grouped operand slots with known group identities.
    // Cluster inheriting columns by their 4-tuple
    //   (bra_supplier, bra_group_idx, ket_supplier, ket_group_idx)
    // so that columns from distinct source groups in the same operand are
    // never merged into one oversized ColumnGroup (C1 fix). Emit one
    // ColumnGroup per cluster of size >= 2, sign +1.
    using ClusterKey = std::tuple<int, std::size_t, int, std::size_t>;
    std::map<ClusterKey, container::svector<std::size_t>> clusters;
    for (std::size_t c = 0; c < ncols; ++c) {
      if (!rbra[c].nonnull() || !rket[c].nonnull()) continue;
      auto b = trace(rbra[c]);
      auto k = trace(rket[c]);
      if (!b || !k || !b->second.column_group_idx ||
          !k->second.column_group_idx)
        continue;
      clusters[{b->first, *b->second.column_group_idx, k->first,
                *k->second.column_group_idx}]
          .push_back(c);
    }
    for (auto& [key, cols] : clusters) {
      if (cols.size() < 2) continue;
      SlotSymmetry::ColumnGroup cg;
      cg.sign = 1;
      cg.cols = std::move(cols);
      ss.column_groups.push_back(std::move(cg));
    }
  }

  // ---- Bra-only / ket-only group inheritance ----
  // Build full_label -> result-bundle-position maps for fast membership checks.
  // Using full_label() (not label()) avoids aliasing proto-indexed indices
  // that share a base label but differ in proto-indices (I2 fix).
  const std::size_t rbra_rank = result.bra_rank();
  const std::size_t rket_rank = result.ket_rank();

  std::unordered_map<std::wstring, std::size_t> rbra_pos, rket_pos;
  for (std::size_t p = 0; p < rbra_rank; ++p)
    if (rbra[p].nonnull())
      rbra_pos.emplace(std::wstring{rbra[p].full_label()}, p);
  for (std::size_t p = 0; p < rket_rank; ++p)
    if (rket[p].nonnull())
      rket_pos.emplace(std::wstring{rket[p].full_label()}, p);

  // Try to inherit one operand slot-group into one result bundle. All member
  // indices of the operand group must appear in the result bundle
  // (whole-group-survives guard). Emits a result SlotGroup over the result
  // positions, with the operand sign carried verbatim.
  auto try_inherit =
      [&](SlotSymmetry::SlotGroup const& og, auto const& ot_bundle,
          std::size_t ot_bundle_rank,
          std::unordered_map<std::wstring, std::size_t> const& res_pos,
          std::size_t res_rank,
          container::svector<SlotSymmetry::SlotGroup>& res_groups) {
        if (res_rank < 2) return;
        container::svector<std::size_t> result_positions;
        result_positions.reserve(og.slots.size());
        for (std::size_t s : og.slots) {
          if (s >= ot_bundle_rank || !ot_bundle[s].nonnull()) return;
          auto it = res_pos.find(std::wstring{ot_bundle[s].full_label()});
          if (it == res_pos.end()) return;
          result_positions.push_back(it->second);
        }
        if (result_positions.empty()) return;
        SlotSymmetry::SlotGroup rg;
        rg.sign = og.sign;
        rg.slots = std::move(result_positions);
        res_groups.push_back(std::move(rg));
      };

  // Check each operand's bra_groups and ket_groups for whole-group survival
  // into the result bra or ket bundle.
  auto inherit_from_operand = [&](Tensor const& ot, SlotSymmetry const& oss) {
    for (auto const& og : oss.bra_groups) {
      try_inherit(og, ot.bra(), ot.bra_rank(), rbra_pos, rbra_rank,
                  ss.bra_groups);
      try_inherit(og, ot.bra(), ot.bra_rank(), rket_pos, rket_rank,
                  ss.ket_groups);
    }
    for (auto const& og : oss.ket_groups) {
      try_inherit(og, ot.ket(), ot.ket_rank(), rbra_pos, rbra_rank,
                  ss.bra_groups);
      try_inherit(og, ot.ket(), ot.ket_rank(), rket_pos, rket_rank,
                  ss.ket_groups);
    }
  };
  inherit_from_operand(lt, left.slot_symmetry());
  inherit_from_operand(rt, right.slot_symmetry());

  return ss;
}

}  // namespace sequant
