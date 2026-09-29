#ifndef SEQUANT_CORE_EXPORT_SCHEDULE_EXPORT_HPP
#define SEQUANT_CORE_EXPORT_SCHEDULE_EXPORT_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/memory_model.hpp>
#include <SeQuant/core/export/natural_memory_model.hpp>
#include <SeQuant/core/export/value_resolver.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cstddef>
#include <type_traits>
#include <variant>

// This translation unit (and everything else under SeQuant/core/export/) may
// depend on eval's pure scheduling data (cell_table.hpp, ordered_schedule.hpp,
// ...) but must never depend on its runtime execution machinery
// (cell_registry.hpp, ordered_executor.hpp), which has no meaning for a
// one-shot compile. A CI check (.github/workflows/formatting_check.yml) and a
// matching pre-commit hook enforce this by rejecting those two includes
// anywhere under this directory.

namespace sequant {

/// Drives an (unmodified) Generator<Context> from eval's schedule IR
/// (OrderedSchedule/CellTable) instead of from export's own ExportNode tree
/// walk (see GenerationVisitor in export.hpp). This is the schedule-driven
/// alternative to that tree walk, so that future improvements to eval's
/// scheduling (contraction order, batching, sharing detection) are reused by
/// export instead of being independently re-implemented there.
///
/// Scope of this first cut (documented gaps, to be closed before this can be
/// wired to a real Generator pipeline in migration step 5):
///  - EvalOp::Adjoint is not supported; encountering one throws.
///  - Sum's accumulation-folding optimization (export's ComputeSelection,
///    see compute_selection.hpp) has no equivalent here yet: every Sum-typed
///    cell is computed in full from all of its operands in one compute()
///    call, which is correct but, unlike the legacy path, never elides a
///    redundant intermediate via an implicit "+=" fold.
///  - Scalar-factor folding (export's
///    PreprocessVisitor::prune_scalar_factor, export.hpp) has no equivalent
///    here yet: a Variable- or Power-of-Variable-valued operand is
///    created/loaded/computed/unloaded like any other operand, rather than
///    being folded into its consumer's compute() call.
///  - Only an unbatched schedule is accepted (every step of
///    OrderedSchedule::root is a BuildStep, never a nested ScopeBlock); a
///    batched schedule throws. Batching support is deliberately deferred
///    (see the migration plan's "must not preclude" checklist).
///  - A leaf operand is (re)loaded independently at every read, with no
///    reuse across reads of the same leaf within one export -- unlike the
///    legacy path's GenerationVisitor::m_tensorUses/m_variableUses, which
///    de-duplicates repeated leaf occurrences within one tree. Always
///    correct, merely not as I/O-efficient; a documented minor gap.
///  - Every compute() call currently corresponds to exactly one CellTable
///    Build/Assemble cell (one EvalExpr node, i.e. one former binary
///    contraction). Operands are deliberately gathered generically from
///    however many Reads a cell has, rather than assuming exactly two, so
///    that a future preprocessing pass could fuse a run of single-consumer
///    cells into one combined compute() call spanning a larger subtree
///    (export of an arbitrary (sub)tree expression, not just a binary
///    contraction) without requiring a rewrite of this operand-gathering
///    logic -- see the corresponding item in the migration plan.
template <typename Context>
class ScheduleWalkingGenerationVisitor {
  static_assert(
      std::is_base_of_v<ExportContext, Context>,
      "Generator context class must inherit from sequant::ExportContext");

 public:
  /// The memory model is not a parameter here: it is a fixed property of
  /// \p generator's target (see Generator::memory_model(), generator.hpp),
  /// queried from it rather than chosen by whoever constructs this visitor.
  ScheduleWalkingGenerationVisitor(const eval::CellTable &table,
                                   const eval::OrderedSchedule &schedule,
                                   const ValueResolver &resolver,
                                   Generator<Context> &generator, Context &ctx)
      : m_table(table),
        m_schedule(schedule),
        m_resolver(resolver),
        m_generator(generator),
        m_ctx(ctx),
        m_model(generator.memory_model()),
        m_remaining_life(table.cells.size()) {
    // Only MemoryModel::Stack ever needs to know a cell's natural model (to
    // find the ones that fall short and must be force-persisted); computing
    // it is pointless work under any other model. Written as a
    // switch-without-default (rather than `if (m_model == Stack)`) so that
    // adding a third MemoryModel value later forces a compile warning here
    // instead of silently reusing the "nothing to do" branch for it.
    switch (m_model) {
      case MemoryModel::RandomAccess:
        break;
      case MemoryModel::Stack:
        m_natural = natural_memory_model(table);
        break;
    }

    for (eval::CellId id = 0; id < table.cells.size(); ++id)
      m_remaining_life.at(id) = table.cells.at(id).life;
  }

  /// Walks the schedule once, driving the generator.
  void run() {
    require_unbatched();

    for (eval::CellId id = 0; id < m_table.cells.size(); ++id) {
      const eval::TableCell &cell = m_table.cells.at(id);
      if (cell.production.kind == eval::ProductionKind::Leaf) continue;
      produce(id);
    }
  }

 private:
  void require_unbatched() const {
    for (const eval::Step &step : m_schedule.root.steps) {
      if (!std::holds_alternative<eval::BuildStep>(step.value)) {
        throw Exception(
            "ScheduleWalkingGenerationVisitor: batched schedules (nested "
            "ScopeBlocks) are not yet supported");
      }
    }
  }

  const EvalExpr &node(eval::CellId cell) const {
    return m_resolver.node_of(m_table.cells.at(cell).value_id);
  }

  bool is_leaf(eval::CellId cell) const {
    return m_table.cells.at(cell).production.kind == eval::ProductionKind::Leaf;
  }

  /// Under MemoryModel::Stack, whether \p cell's natural model
  /// (natural_memory_model.hpp) falls short of Stack -- i.e. it must be
  /// force-persisted: reloaded independently at every read rather than kept
  /// resident. Checking "falls short of Stack" rather than "equals
  /// RandomAccess" (and, below, dispatching on the natural model via a
  /// switch with no default) means a future third MemoryModel value is
  /// handled conservatively (or flagged by a compiler warning for an
  /// unhandled case) rather than silently treated as already
  /// stack-compatible. Always false under MemoryModel::RandomAccess (the
  /// natural model is never even computed there, see the constructor).
  bool requires_persist(eval::CellId cell) const {
    if (m_model != MemoryModel::Stack) return false;
    switch (m_natural.at(cell)) {
      case MemoryModel::Stack:
        return false;
      case MemoryModel::RandomAccess:
        return true;
    }
    SEQUANT_UNREACHABLE;
  }

  /// Whether \p cell must be (re)loaded at every read site rather than being
  /// already resident from its own production: true for a genuine leaf, and,
  /// under MemoryModel::Stack, for a force-persisted non-leaf cell too.
  bool reloaded_at_every_read(eval::CellId cell) const {
    return is_leaf(cell) || requires_persist(cell);
  }

  /// Ensures the operand at \p cell is loaded, immediately before use. Used
  /// both for genuine leaves and, under MemoryModel::Stack, for a
  /// force-persisted non-leaf cell's reload at each of its reads.
  void load_operand(eval::CellId cell) {
    const EvalExpr &n = node(cell);
    if (n.is_constant()) return;  // Constants need no load/create.
    if (n.is_tensor()) {
      const Tensor &t = n.as_tensor();
      m_generator.load(
          t,
          static_cast<bool>(m_ctx.zeroStrategy(t) & ZeroStrategy::ZeroOnLoad),
          m_ctx);
    } else if (n.is_variable()) {
      const Variable &v = n.as_variable();
      m_generator.load(
          v,
          static_cast<bool>(m_ctx.zeroStrategy(v) & ZeroStrategy::ZeroOnLoad),
          m_ctx);
    } else {
      throw Exception(
          "ScheduleWalkingGenerationVisitor: unsupported leaf expression "
          "type");
    }
  }

  /// Closes a reload-at-every-read operand immediately after use. A genuine
  /// leaf is always merely unloaded (never destroyed: it remains reloadable
  /// from outside the schedule forever, and its life is not tracked). A
  /// force-persisted non-leaf cell's life is tracked exactly like a
  /// resident cell's: unloaded on a non-draining read, destroyed on the one
  /// that drains it (its persisted copy, made right after its own
  /// production, is what makes every read but the last one safely
  /// reloadable).
  void close_reloaded_operand(eval::CellId cell) {
    const EvalExpr &n = node(cell);
    if (n.is_constant()) return;

    bool destroy_now = false;
    if (!is_leaf(cell)) {
      SEQUANT_ASSERT(m_remaining_life.at(cell) > 0);
      --m_remaining_life.at(cell);
      destroy_now = m_remaining_life.at(cell) == 0;
    }

    if (n.result_type() == ResultType::Tensor) {
      const Tensor &t = n.as_tensor();
      destroy_now ? m_generator.destroy(t, m_ctx)
                  : m_generator.unload(t, m_ctx);
    } else {
      const Variable &v = n.as_variable();
      destroy_now ? m_generator.destroy(v, m_ctx)
                  : m_generator.unload(v, m_ctx);
    }
  }

  /// Creates (or loads, per LoadStrategy) the result placeholder for the
  /// cell about to be computed.
  void create_result(eval::CellId cell) {
    const EvalExpr &n = node(cell);
    if (n.result_type() == ResultType::Tensor) {
      const Tensor &t = n.as_tensor();
      const ZeroStrategy zero = m_ctx.zeroStrategy(t);
      if (m_ctx.loadStrategy(t) == LoadStrategy::Create) {
        m_generator.create(
            t, static_cast<bool>(zero & ZeroStrategy::ZeroOnCreate), m_ctx);
      } else {
        m_generator.load(t, static_cast<bool>(zero & ZeroStrategy::ZeroOnLoad),
                         m_ctx);
      }
    } else {
      const Variable &v = n.as_variable();
      const ZeroStrategy zero = m_ctx.zeroStrategy(v);
      if (m_ctx.loadStrategy(v) == LoadStrategy::Create) {
        m_generator.create(
            v, static_cast<bool>(zero & ZeroStrategy::ZeroOnCreate), m_ctx);
      } else {
        m_generator.load(v, static_cast<bool>(zero & ZeroStrategy::ZeroOnLoad),
                         m_ctx);
      }
    }
  }

  /// Decrements a resident (not reload-at-every-read) operand's remaining
  /// life after it has been read. A cell with further reads ahead simply
  /// stays resident (no generator call at all -- it was created once, at
  /// its own production, and every read up to and including this one just
  /// references that same, still-live result); only the read that drains
  /// its life to zero destroys it.
  void close_operand(eval::CellId source) {
    SEQUANT_ASSERT(m_remaining_life.at(source) > 0);
    --m_remaining_life.at(source);
    if (m_remaining_life.at(source) > 0) return;

    const EvalExpr &n = node(source);
    if (n.result_type() == ResultType::Tensor)
      m_generator.destroy(n.as_tensor(), m_ctx);
    else
      m_generator.destroy(n.as_variable(), m_ctx);
  }

  void persist_root(eval::CellId cell) {
    const EvalExpr &n = node(cell);
    if (n.result_type() == ResultType::Tensor)
      m_generator.persist(n.as_tensor(), m_ctx);
    else
      m_generator.persist(n.as_variable(), m_ctx);
  }

  /// MemoryModel::Stack, cells whose natural model falls short of Stack
  /// only: immediately after this cell's own compute() call, store it and
  /// free the in-memory copy that compute() just filled -- every later read
  /// reloads it independently (see load_operand/close_reloaded_operand), so
  /// its live span is always immediately closed rather than held open
  /// across the schedule.
  void persist_and_unload_self(eval::CellId cell) {
    const EvalExpr &n = node(cell);
    if (n.result_type() == ResultType::Tensor) {
      const Tensor &t = n.as_tensor();
      m_generator.persist(t, m_ctx);
      m_generator.unload(t, m_ctx);
    } else {
      const Variable &v = n.as_variable();
      m_generator.persist(v, m_ctx);
      m_generator.unload(v, m_ctx);
    }
  }

  void produce(eval::CellId id) {
    const eval::TableCell &cell = m_table.cells.at(id);
    const EvalExpr &self = node(id);

    SEQUANT_ASSERT(self.op_type().has_value());
    if (*self.op_type() == EvalOp::Adjoint) {
      throw Exception(
          "ScheduleWalkingGenerationVisitor: EvalOp::Adjoint is not yet "
          "supported");
    }

    create_result(id);

    // Gather this cell's operands, in table order, (re)loading each operand
    // that isn't already resident just-in-time: a genuine leaf, or, under
    // MemoryModel::Stack, a force-persisted non-leaf cell. Any other
    // non-leaf operand is already resident, having been created at its own,
    // earlier production step.
    container::svector<eval::CellId> operand_cells;
    container::svector<ExprPtr> operand_exprs;
    for (const eval::Read &r : m_table.reads) {
      if (r.consumer != id) continue;
      if (reloaded_at_every_read(r.source)) load_operand(r.source);
      operand_cells.push_back(r.source);
      operand_exprs.push_back(node(r.source).expr());
    }
    SEQUANT_ASSERT(!operand_exprs.empty());

    ExprPtr computation =
        *self.op_type() == EvalOp::Product
            ? ex<Product>(operand_exprs.begin(), operand_exprs.end(),
                          Product::Flatten::No)
            : ex<Sum>(operand_exprs.begin(), operand_exprs.end());

    if (self.result_type() == ResultType::Tensor)
      m_generator.compute(*computation, self.as_tensor(), m_ctx);
    else
      m_generator.compute(*computation, self.as_variable(), m_ctx);

    for (eval::CellId operand : operand_cells) {
      if (reloaded_at_every_read(operand))
        close_reloaded_operand(operand);
      else
        close_operand(operand);
    }

    if (cell.life == 0)
      persist_root(id);
    else if (requires_persist(id))
      persist_and_unload_self(id);
  }

  const eval::CellTable &m_table;
  const eval::OrderedSchedule &m_schedule;
  const ValueResolver &m_resolver;
  Generator<Context> &m_generator;
  Context &m_ctx;
  MemoryModel m_model;
  container::vector<std::size_t> m_remaining_life;
  /// Empty under MemoryModel::RandomAccess (never consulted there); one
  /// entry per cell under MemoryModel::Stack, from natural_memory_model().
  container::vector<MemoryModel> m_natural;
};

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_SCHEDULE_EXPORT_HPP
