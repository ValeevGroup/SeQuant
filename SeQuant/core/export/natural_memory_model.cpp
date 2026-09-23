#include <SeQuant/core/export/natural_memory_model.hpp>

#include <algorithm>
#include <optional>

namespace sequant {

namespace {

/// For every cell, the position (CellId) of the event that exhausts its
/// life: the greatest consumer id among the Reads naming it as source and
/// the Assemble cells whose production.source names it. table.cells is in
/// execution order, so this is exactly the position where the cell is last
/// needed; nullopt means the cell is never consumed (a forest root).
container::vector<std::optional<eval::CellId>> last_consumer_positions(
    const eval::CellTable &table) {
  container::vector<std::optional<eval::CellId>> last(table.cells.size());
  auto update = [&](eval::CellId source, eval::CellId consumer) {
    if (!last.at(source) || consumer > *last.at(source))
      last.at(source) = consumer;
  };
  for (const eval::Read &r : table.reads) update(r.source, r.consumer);
  for (eval::CellId id = 0; id < table.cells.size(); ++id) {
    const eval::TableCell &c = table.cells.at(id);
    if (c.production.kind == eval::ProductionKind::Assemble)
      update(c.production.source, id);
  }
  return last;
}

}  // namespace

container::vector<MemoryModel> natural_memory_model(
    const eval::CellTable &table) {
  const auto n = table.cells.size();
  container::vector<MemoryModel> result(n, MemoryModel::Stack);
  const auto last_consumer = last_consumer_positions(table);

  // The simulated stack of cell ids whose production has been observed but
  // whose closing (draining) read has not yet been processed.
  container::vector<eval::CellId> stack;

  auto close = [&](eval::CellId source) {
    if (!stack.empty() && stack.back() == source) {
      stack.pop_back();
    } else {
      // Either source is buried under something still open, or it was
      // already flagged and removed earlier (e.g. read by more than one
      // consumer, both attempting to close it -- only the true draining
      // read should ever reach here per last_consumer_positions, so this
      // branch is the buried case). Either way it cannot sit at the top of
      // a single stack at this point.
      result.at(source) = MemoryModel::RandomAccess;
      stack.erase(std::remove(stack.begin(), stack.end(), source), stack.end());
    }
  };

  for (eval::CellId id = 0; id < n; ++id) {
    const eval::TableCell &cell = table.cells.at(id);
    if (cell.production.kind == eval::ProductionKind::Leaf) continue;

    // Resolve this cell's own operand reads first: an operand's span closes
    // (if this is its last use) before this cell's own value becomes live.
    for (const eval::Read &r : table.reads) {
      if (r.consumer != id) continue;
      if (table.cells.at(r.source).production.kind ==
          eval::ProductionKind::Leaf)
        continue;  // Leaves are always independently loadable; never tracked.
      if (last_consumer.at(r.source) && *last_consumer.at(r.source) == id)
        close(r.source);
    }
    if (cell.production.kind == eval::ProductionKind::Assemble &&
        last_consumer.at(cell.production.source) &&
        *last_consumer.at(cell.production.source) == id) {
      close(cell.production.source);
    }

    stack.push_back(id);
  }

  return result;
}

}  // namespace sequant
