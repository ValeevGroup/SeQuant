#ifndef SEQUANT_CORE_EXPORT_VALUE_RESOLVER_HPP
#define SEQUANT_CORE_EXPORT_VALUE_RESOLVER_HPP

#include <SeQuant/core/eval/eval_expr.hpp>

#include <cstddef>

namespace sequant {

/// Resolves a schedule value id (CellTable/OrderedSchedule's \c value_id) to
/// the EvalExpr that computes it.
///
/// This exists so ScheduleWalkingGenerationVisitor (schedule_export.hpp)
/// does not need to know how a concrete pipeline bridges a schedule back to
/// the eval forest it was built from. In the real eval pipeline that bridge
/// is a RichSchedule plus a ValueNodeMap over the original forest (see
/// SeQuant/core/eval/value_node_map.hpp) -- assembling that bridge is a
/// pipeline-wiring concern (migration step 5), not something the
/// schedule-walking visitor itself should own; a ValueResolver
/// implementation is where that assembly happens.
class ValueResolver {
 public:
  virtual ~ValueResolver() = default;

  /// @returns The EvalExpr that computes the value with the given id.
  [[nodiscard]] virtual const EvalExpr &node_of(std::size_t value_id) const = 0;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_VALUE_RESOLVER_HPP
