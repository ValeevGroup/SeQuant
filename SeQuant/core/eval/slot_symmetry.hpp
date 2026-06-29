#ifndef SEQUANT_EVAL_SLOT_SYMMETRY_HPP
#define SEQUANT_EVAL_SLOT_SYMMETRY_HPP

#include <SeQuant/core/container.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace sequant {

///
/// \brief Out-of-band descriptor of the permutational symmetry of an
///        intermediate result tensor produced by an EvalExpr node.
///
/// \details This struct records the exploitable index-permutation symmetry of
///          an EvalExpr's result tensor AFTER deduction (a later task).  It is
///          carried out-of-band on EvalExpr (i.e. NOT embedded in the result
///          Tensor's symmetry tag) so that it does NOT influence hashing,
///          canonicalization, bliss-graph construction, or export.  Default-
///          constructed instances represent "no exploitable symmetry" and are
///          considered empty.
///
struct SlotSymmetry {
  ///
  /// \brief A group of matched (bra[c], ket[c]) column pairs that may be
  ///        permuted together (possibly with a sign).
  ///
  struct ColumnGroup {
    /// Zero-based column indices whose (bra, ket) pairs may be permuted.
    container::svector<std::size_t> cols;
    /// +1 for symmetric, -1 for antisymmetric permutation within this group.
    std::int8_t sign{1};
  };

  ///
  /// \brief A group of bra (or ket) slot indices that may be permuted
  ///        (possibly with a sign).
  ///
  struct SlotGroup {
    /// Zero-based slot indices (within the bra or ket) that may be permuted.
    container::svector<std::size_t> slots;
    /// +1 for symmetric, -1 for antisymmetric permutation within this group.
    std::int8_t sign{1};
  };

  /// Permutation symmetry over matched (bra[c], ket[c]) column pairs.
  container::svector<ColumnGroup> column_groups;

  /// Permutation symmetry within the bra slots only.
  container::svector<SlotGroup> bra_groups;

  /// Permutation symmetry within the ket slots only.
  container::svector<SlotGroup> ket_groups;

  ///
  /// \return true if this descriptor records no exploitable symmetry
  ///         (all group containers are empty).
  ///
  [[nodiscard]] bool empty() const noexcept {
    return column_groups.empty() && bra_groups.empty() && ket_groups.empty();
  }

  ///
  /// \brief Order-insensitive equality: two SlotSymmetry objects are equal if
  ///        their group containers have the same elements regardless of order.
  ///
  friend bool operator==(SlotSymmetry const& lhs,
                         SlotSymmetry const& rhs) noexcept {
    auto col_eq = [](ColumnGroup const& a, ColumnGroup const& b) {
      return a.sign == b.sign && a.cols == b.cols;
    };
    auto slot_eq = [](SlotGroup const& a, SlotGroup const& b) {
      return a.sign == b.sign && a.slots == b.slots;
    };

    if (lhs.column_groups.size() != rhs.column_groups.size()) return false;
    if (lhs.bra_groups.size() != rhs.bra_groups.size()) return false;
    if (lhs.ket_groups.size() != rhs.ket_groups.size()) return false;

    // Check column_groups: same multiset of ColumnGroup elements.
    {
      auto l = lhs.column_groups;
      auto r = rhs.column_groups;
      auto cmp = [](ColumnGroup const& a, ColumnGroup const& b) {
        if (a.sign != b.sign) return a.sign < b.sign;
        return a.cols < b.cols;
      };
      std::sort(l.begin(), l.end(), cmp);
      std::sort(r.begin(), r.end(), cmp);
      for (std::size_t i = 0; i < l.size(); ++i)
        if (!col_eq(l[i], r[i])) return false;
    }

    // Check bra_groups.
    {
      auto l = lhs.bra_groups;
      auto r = rhs.bra_groups;
      auto cmp = [](SlotGroup const& a, SlotGroup const& b) {
        if (a.sign != b.sign) return a.sign < b.sign;
        return a.slots < b.slots;
      };
      std::sort(l.begin(), l.end(), cmp);
      std::sort(r.begin(), r.end(), cmp);
      for (std::size_t i = 0; i < l.size(); ++i)
        if (!slot_eq(l[i], r[i])) return false;
    }

    // Check ket_groups.
    {
      auto l = lhs.ket_groups;
      auto r = rhs.ket_groups;
      auto cmp = [](SlotGroup const& a, SlotGroup const& b) {
        if (a.sign != b.sign) return a.sign < b.sign;
        return a.slots < b.slots;
      };
      std::sort(l.begin(), l.end(), cmp);
      std::sort(r.begin(), r.end(), cmp);
      for (std::size_t i = 0; i < l.size(); ++i)
        if (!slot_eq(l[i], r[i])) return false;
    }

    return true;
  }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_SLOT_SYMMETRY_HPP
