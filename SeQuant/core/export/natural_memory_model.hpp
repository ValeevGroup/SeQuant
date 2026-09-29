#ifndef SEQUANT_CORE_EXPORT_NATURAL_MEMORY_MODEL_HPP
#define SEQUANT_CORE_EXPORT_NATURAL_MEMORY_MODEL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/export/memory_model.hpp>

namespace sequant {

/// For every non-Leaf cell of \p table, the MemoryModel its own live span
/// (from production to last use) naturally satisfies on its own, given that
/// \p table.cells is already in execution order (see validate_cell_table's
/// doc comment on that invariant, which this function also relies on):
/// MemoryModel::Stack iff, at the point its span closes (its last use),
/// nothing else that opened after it is still open -- i.e. it would sit at
/// the top of a single, real stack at that moment. Every other cell
/// naturally satisfies only MemoryModel::RandomAccess, which imposes no
/// nesting requirement at all, so that is the result whenever the stricter
/// Stack test fails. Leaf cells are exempt: per CellTable's own model they
/// are "fetched on demand from outside" the table rather than tracked with a
/// refcounted residency, so a read of one is always independently loadable
/// and never forces anything else open.
///
/// A caller driving MemoryModel::Stack (ScheduleWalkingGenerationVisitor,
/// schedule_export.hpp) uses this to find every cell whose natural model
/// falls short of Stack, and force-persists exactly those -- storing each
/// right after its own production, then reloading it independently at every
/// later read -- so the exported code's alloc/dealloc sequence stays
/// properly stack-nested despite what the cell's own span would otherwise
/// require.
///
/// \note Assumes an unbatched schedule (every cell's CellScope is the root
/// scope, i.e. \c table.cells never carries loop-nested forms). A batched
/// schedule additionally needs loop-boundary information (from
/// OrderedSchedule) to correctly force-persist a value whose sole reader
/// sits inside a loop that reads it once per batch -- not handled here yet.
[[nodiscard]] container::vector<MemoryModel> natural_memory_model(
    const eval::CellTable &table);

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_NATURAL_MEMORY_MODEL_HPP
