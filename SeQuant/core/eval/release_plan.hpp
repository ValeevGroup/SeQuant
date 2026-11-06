#ifndef SEQUANT_EVAL_RELEASE_PLAN_HPP
#define SEQUANT_EVAL_RELEASE_PLAN_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <variant>

namespace sequant::eval {

/// When each cell of a table is let go, read off the schedule.
///
/// A cell lives for one epoch: from its production until the loop instance
/// that re-produces it advances (\c detail::epoch_depth). Within an epoch the
/// block of that loop -- the root block for a cell bound to no loop -- runs
/// its timeline once: its steps in order, then its outputs in order. Every
/// read of the cell happens at one position of that timeline, either directly
/// (a build of that block, or one of its outputs) or inside the nested block
/// that is the step at that position. The cell is released when the last
/// position holding a read of it is done. Nothing is counted: a step that is
/// skipped, or replays its reads once per batch of loops nested in it, is
/// done all the same.
struct ReleasePlan {
  /// Per block of the schedule (the root block included), per position of
  /// its timeline (steps, then outputs): the cells released once that
  /// position is done.
  std::unordered_map<ScopeBlock const*,
                     container::vector<container::svector<CellId>>>
      after;
  /// Per \c Read, by its position in \c CellTable::reads: whether it is the
  /// last use of its source, so that the reader takes the value over. Only a
  /// read made once per epoch of its source can be: one by a build or an
  /// output of the source's own epoch block.
  container::vector<char> last_read;
  /// Per cell: for an \c Assemble, whether its read of \c Production::source
  /// is the last use of that source.
  container::vector<char> last_assemble;
  /// Per cell: whether it is read from outside the table once the walk is
  /// over (\c build_release_plan's \p keep): never released, and never dead
  /// work, whatever becomes of its readers in the table.
  container::vector<char> kept;

  /// \return the cells released after position \p pos of \p block, or null
  [[nodiscard]] container::svector<CellId> const* released_after(
      ScopeBlock const& block, std::size_t pos) const {
    auto const it = after.find(&block);
    if (it == after.end() || pos >= it->second.size() ||
        it->second[pos].empty())
      return nullptr;
    return &it->second[pos];
  }
};

namespace detail {

/// The length of the prefix of \p c's scope path that ends at the loop
/// instance whose advance re-produces \p c: the deepest enclosing instance
/// \p c is bound to (\c bound_instances), or its innermost enclosing loop
/// when the cell is produced on every batch of it (not \c produce_if_absent).
/// 0 for a cell no loop re-produces (a leaf, or a cell bound to none of its
/// enclosing loops and produced on first visit only). This is the scope \c
/// CellRegistry::clear_bound_to empties the cell on.
[[nodiscard]] inline std::size_t epoch_depth(TableCell const& c) {
  if (c.production.kind == ProductionKind::Leaf) return 0;
  auto const& path = c.scope.path;
  if (!c.produce_if_absent) return path.size();
  std::size_t depth = 0;
  for (LoopKey const& k : bound_instances(c))
    for (std::size_t i = depth; i < path.size(); ++i)
      if (same_key(k, path[i].first)) depth = i + 1;
  return depth;
}

}  // namespace detail

/// Derives the \c ReleasePlan of \p table on the schedule tree \p root it was
/// built from. The cells in \p keep are never released (the schedule's
/// results, which the caller reads once the walk is over); neither are
/// persistent cells, nor cells nothing in the table reads.
/// \throw Exception when a cell's scope or a reader's place is not in the
///        tree, i.e. the table and the schedule disagree
[[nodiscard]] inline ReleasePlan build_release_plan(
    CellTable const& table, ScopeBlock const& root,
    container::svector<CellId> const& keep = {}) {
  std::size_t const n = table.cells.size();
  ReleasePlan plan;
  plan.last_read.assign(table.reads.size(), 0);
  plan.last_assemble.assign(n, 0);

  // Per block: where its nested blocks, its builds and its outputs sit in its
  // timeline.
  using LoopAt = std::tuple<std::size_t, int, int>;  // depth, slot, latitude
  struct BlockIndex {
    std::map<LoopAt, std::size_t> child;                  // -> step
    std::unordered_map<std::size_t, std::size_t> build;   // value -> step
    std::unordered_map<std::size_t, std::size_t> output;  // value -> position
  };
  std::unordered_map<ScopeBlock const*, BlockIndex> index;
  auto const index_of = [&](ScopeBlock const& b) -> BlockIndex const& {
    auto const it = index.find(&b);
    if (it != index.end()) return it->second;
    BlockIndex bi;
    for (std::size_t i = 0; i < b.steps.size(); ++i) {
      if (auto const* build = std::get_if<BuildStep>(&b.steps[i].value)) {
        bi.build.emplace(build->value_id, i);
      } else {
        auto const& c = std::get<ScopeBlock>(b.steps[i].value);
        bi.child.emplace(
            LoopAt{c.level.depth, c.level.loop_slot, c.latitude_ordinal}, i);
      }
    }
    for (std::size_t k = 0; k < b.outputs.size(); ++k)
      bi.output.emplace(b.outputs[k].first, b.steps.size() + k);
    return index.emplace(&b, std::move(bi)).first->second;
  };
  auto const loop_at = [](std::pair<LoopKey, int> const& e) {
    return LoopAt{e.first.depth, e.first.loop_slot, e.second};
  };
  auto const disagree = [&](CellId c, char const* what) {
    return Exception("build_release_plan: " + detail::cell_str(table, c) + " " +
                     what + " (cell table/schedule disagreement)");
  };

  // The reads of each cell: operand reads by position in table.reads, and the
  // Assemble cells sourced from it.
  container::vector<container::svector<std::size_t>> reads_of(n);
  for (std::size_t r = 0; r < table.reads.size(); ++r)
    if (table.reads[r].source < n) reads_of[table.reads[r].source].push_back(r);
  container::vector<container::svector<CellId>> assembles_of(n);
  for (CellId a = 0; a < n; ++a)
    if (table.cells[a].production.kind == ProductionKind::Assemble &&
        table.cells[a].production.source < n)
      assembles_of[table.cells[a].production.source].push_back(a);
  plan.kept.assign(n, 0);
  for (CellId c : keep)
    if (c < n) plan.kept[c] = 1;
  container::vector<char> const& kept = plan.kept;

  for (CellId s = 0; s < n; ++s) {
    TableCell const& cell = table.cells[s];
    if (cell.persistent || kept[s]) continue;
    if (reads_of[s].empty() && assembles_of[s].empty()) continue;

    // The block whose timeline one epoch of this cell spans.
    std::size_t const depth = detail::epoch_depth(cell);
    ScopeBlock const* block = &root;
    for (std::size_t i = 0; i < depth; ++i) {
      auto const& bi = index_of(*block);
      auto const it = bi.child.find(loop_at(cell.scope.path[i]));
      if (it == bi.child.end())
        throw disagree(s, "has a scope not in the tree");
      block = &std::get<ScopeBlock>(block->steps[it->second].value);
    }
    BlockIndex const& bi = index_of(*block);

    // The timeline position of a read made at scope `at`: by the build or
    // the output of value `vid` when `at` is the epoch block's own scope,
    // else somewhere inside the nested block on the way to `at`.
    auto const position_of =
        [&](container::svector<std::pair<LoopKey, int>> const& at,
            std::size_t vid, bool assemble,
            CellId reader) -> std::pair<std::size_t, bool> {
      if (at.size() < depth) throw disagree(reader, "reads outside the epoch");
      if (at.size() > depth) {
        auto const it = bi.child.find(loop_at(at[depth]));
        if (it == bi.child.end())
          throw disagree(reader, "reads from a scope not in the tree");
        return {it->second, false};
      }
      if (!assemble)
        if (auto const it = bi.build.find(vid); it != bi.build.end())
          return {it->second, true};
      auto const it = bi.output.find(vid);
      if (it == bi.output.end())
        throw disagree(reader, "is neither a build nor an output of its block");
      return {it->second, true};
    };

    // The last position holding a read, and the read made last within it --
    // an operand read (by position in table.reads) or an Assemble's read of
    // its source -- unless a read there is replayed once per batch of a
    // nested loop, in which case no reader takes the value.
    struct {
      bool any = false;
      std::size_t pos = 0;
      bool replayed = false;
      std::optional<std::size_t> operand_read = {};
      std::optional<CellId> assemble = {};
    } last;
    auto const note = [&](std::pair<std::size_t, bool> where, bool assemble,
                          std::size_t id) {
      auto const [pos, once] = where;
      if (last.any && pos < last.pos) return;
      if (!last.any || pos > last.pos) last = {.any = true, .pos = pos};
      if (!once)
        last.replayed = true;
      else if (assemble)
        last.assemble = id;
      else
        last.operand_read = id;  // table.reads is in leg order per consumer
    };
    for (std::size_t r : reads_of[s]) {
      CellId const consumer = table.reads[r].consumer;
      TableCell const& c = table.cells[consumer];
      note(position_of(c.scope.path, c.value_id, /*assemble=*/false, consumer),
           /*assemble=*/false, r);
    }
    for (CellId a : assembles_of[s]) {
      // An output reads its source inside the block that lists it: the loop
      // just below the Assemble's own scope along the source's path.
      TableCell const& ac = table.cells[a];
      if (cell.scope.path.size() <= ac.scope.path.size())
        throw disagree(a, "assembles a source not nested below it");
      container::svector<std::pair<LoopKey, int>> const at(
          cell.scope.path.begin(),
          cell.scope.path.begin() + ac.scope.path.size() + 1);
      note(position_of(at, ac.value_id, /*assemble=*/true, a),
           /*assemble=*/true, a);
    }

    auto& slots = plan.after[block];
    if (slots.empty())
      slots.resize(block->steps.size() + block->outputs.size());
    slots[last.pos].push_back(s);
    if (last.replayed) continue;
    // An output reads its source after its own operands.
    if (last.assemble)
      plan.last_assemble[*last.assemble] = 1;
    else if (last.operand_read)
      plan.last_read[*last.operand_read] = 1;
  }
  return plan;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_RELEASE_PLAN_HPP
