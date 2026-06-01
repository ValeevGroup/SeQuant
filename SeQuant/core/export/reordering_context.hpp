#ifndef SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP
#define SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP

#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/memory_layout.hpp>

namespace sequant {

class Tensor;

/// An export context class that implements functionality for reordering indices
/// in tensors for better cache locality during access, assuming the full tensor
/// will be held in memory and all elements will be iterated over.
///
/// In other words, this class will use existing tensor symmetries to move the
/// slot belonging to the largest index space into the slot with highest
/// cache-locality (depending on the chosen memory layout).
class ReorderingContext : public ExportContext {
 public:
  explicit ReorderingContext(MemoryLayout layout) : m_layout(layout) {}

  MemoryLayout memory_layout() const { return m_layout; }
  void set_memory_layout(MemoryLayout layout) { m_layout = layout; }

  bool rewrite(Tensor &tensor) const override;

 private:
  MemoryLayout m_layout;

 protected:
  bool is_less(const IndexSpace &lhs, const IndexSpace &rhs) const;
  bool is_ordered(const IndexSpace &lhs, const IndexSpace &rhs) const;
  bool needs_swap(const IndexSpace &lhs, const IndexSpace &rhs) const;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_REORDERING_CONTEXT_HPP
