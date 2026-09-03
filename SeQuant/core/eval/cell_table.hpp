#ifndef SEQUANT_EVAL_CELL_TABLE_HPP
#define SEQUANT_EVAL_CELL_TABLE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <utility>

namespace sequant::eval {

using LoopKey = sequant::LoopKey;

/// Explicit value cells (see the explicit-value-cells design document). A cell
/// is one FORM of a value resident at one scope; identity is by cell id,
/// carried position and loop instance -- never by canonical index label.
using CellId = std::size_t;

/// Enclosing loop instances outermost-first, each with its latitude (pass).
/// Empty = the root scope.
struct CellScope {
  container::svector<std::pair<LoopKey, int>> path;

  /// true iff this scope is the root or a strict/equal prefix of \p inner
  /// (same loop instances AND same passes along the prefix).
  [[nodiscard]] bool encloses(CellScope const& inner) const noexcept {
    if (path.size() > inner.path.size()) return false;
    for (std::size_t i = 0; i < path.size(); ++i)
      if (path[i] != inner.path[i]) return false;
    return true;
  }
  [[nodiscard]] friend bool operator==(CellScope const& a,
                                       CellScope const& b) noexcept {
    return a.path.size() == b.path.size() && a.encloses(b);
  }
};

enum class ProductionKind { Build, Assemble, Leaf };
enum class AssembleKind { Sum, Scatter };

struct Production {
  ProductionKind kind = ProductionKind::Build;
  AssembleKind assemble = AssembleKind::Sum;  // Assemble only
  CellId source = 0;                          // Assemble only: the partial
  /// Assemble/Scatter only: target position -> the loop instance whose
  /// batches are scattered into it.
  container::svector<std::pair<std::size_t, LoopKey>> scatter_map;
};

/// Represents one form of a value at one scope. Named TableCell to avoid
/// collision with the existing sequant::eval::Cell in peak_profile.hpp.
struct TableCell {
  std::size_t value_id = 0;
  /// carried position -> loop instance slicing it; empty = whole.
  container::svector<std::pair<std::size_t, LoopKey>> sliced;
  CellScope scope;
  Production production;
  bool produce_if_absent = false;
  bool persistent = false;
  std::size_t life = 0;
};

/// One per DAG edge into a consumer cell's production tree.
struct Read {
  CellId consumer = 0;
  std::size_t operand_value_id = 0;
  CellId source = 0;
  /// operand positions to slice to the current batch of that loop instance.
  container::svector<std::pair<std::size_t, LoopKey>> slice;
};

struct CellTable {
  container::vector<TableCell> cells;
  container::vector<Read> reads;
};

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_HPP
