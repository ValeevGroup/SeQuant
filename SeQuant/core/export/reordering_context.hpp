#ifndef SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP

#include <SeQuant/core/export/context.hpp>

namespace sequant {

class Tensor;

enum class MemoryLayout {
  RowMajor,
  ColumnMajor,
  Unspecified,
};

struct assume_real_orbitals {};

/// An export context class that implements functionality for reordering indices
/// in tensors for better cache locality during access, assuming the full tensor
/// will be held in memory and all elements will be iterated over.
///
/// In other words, this class will use existing tensor symmetries to move the
/// slot belonging to the largest index space into the slot with highest
/// cache-locality (depending on the chosen memory layout).
class ReorderingContext : public ExportContext {
 public:
  explicit ReorderingContext(MemoryLayout layout)
      : m_layout(layout), m_real_orbitals(false) {}
  ReorderingContext(MemoryLayout layout, assume_real_orbitals tag)
      : m_layout(layout), m_real_orbitals(true) {}

  MemoryLayout memory_layout() const { return m_layout; }
  void set_memory_layout(MemoryLayout layout) { m_layout = layout; }

  bool assumes_real_orbitals() const { return m_real_orbitals; }
  void set_assume_real_orbitals(bool real) { m_real_orbitals = real; }

  bool rewrite(Tensor &tensor) const override;

 private:
  MemoryLayout m_layout;
  bool m_real_orbitals;

 protected:
  bool is_ordered(const IndexSpace &lhs, const IndexSpace &rhs) const;
  bool needs_swap(const IndexSpace &lhs, const IndexSpace &rhs) const;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP
