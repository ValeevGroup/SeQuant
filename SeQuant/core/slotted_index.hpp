#ifndef SEQUANT_SLOTTED_INDEX_H
#define SEQUANT_SLOTTED_INDEX_H

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/index.hpp>

namespace sequant {

/// Combination of an Index along with information about the Slot it occupies
///
/// @note Usage of this class only makes sense, if the index in question does
/// not appear in multiple slots
///
/// @sa Index
/// @sa Slot
class SlottedIndex {
 public:
  SlottedIndex(Index idx, Slot slot)
      : idx_(std::move(idx)), slot_(std::move(slot)) {}

  Index &index() { return idx_; }

  const Index &index() const { return idx_; }

  Slot slot() const { return slot_; }

  friend bool operator==(const SlottedIndex &lhs, const SlottedIndex &rhs) {
    return lhs.index() == rhs.index() && lhs.slot() == rhs.slot();
  }

  friend bool operator!=(const SlottedIndex &lhs, const SlottedIndex &rhs) {
    return !(lhs == rhs);
  }

 private:
  Index idx_;
  Slot slot_;
};

}  // namespace sequant

#endif
