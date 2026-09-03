#ifndef SEQUANT_EVAL_CELL_TABLE_HPP
#define SEQUANT_EVAL_CELL_TABLE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>

#include <cstddef>
#include <functional>
#include <stdexcept>
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

struct CellViolation {
  std::string rule;  // "visibility" | "form" | "chain" | "life" | "uniqueness"
  std::string what;
};

namespace detail {
[[nodiscard]] inline std::string cell_str(CellTable const& t, CellId id) {
  auto const& c = t.cells[id];
  std::string s = "cell#" + std::to_string(id) + "(value " +
                  std::to_string(c.value_id) + ", scope depth " +
                  std::to_string(c.scope.path.size()) + ", ";
  s += c.production.kind == ProductionKind::Build      ? "Build"
       : c.production.kind == ProductionKind::Assemble ? "Assemble"
                                                       : "Leaf";
  return s + ")";
}
[[nodiscard]] inline bool same_key(LoopKey const& a, LoopKey const& b) {
  return a.depth == b.depth && a.loop_slot == b.loop_slot;
}
}  // namespace detail

/// Static validation of a cell table (spec section 3). Production order is
/// the order of non-Leaf cells in \p table.cells (Task 3's builder emits
/// them in execution order). The block tree \p root is accepted for the
/// Task 3 integration (it is not consulted here beyond being passed through).
[[nodiscard]] inline container::vector<CellViolation> validate_cell_table(
    CellTable const& table, ScopeBlock const& /*root*/) {
  container::vector<CellViolation> out;
  auto const n = table.cells.size();
  for (Read const& r : table.reads)
    if (r.consumer >= n || r.source >= n) {
      out.push_back({"visibility", "read with out-of-range cell id"});
      return out;
    }

  // (1) visibility: a source is resident at the consumer's production iff
  // the source was produced earlier (or is a Leaf) and its scope encloses
  // the consumer's scope (a per-batch cell of an inner loop is dead outside
  // it), or it is persistent / produce_if_absent and at an enclosing scope.
  {
    container::vector<bool> produced(n, false);
    for (CellId id = 0; id < n; ++id)
      if (table.cells[id].production.kind == ProductionKind::Leaf)
        produced[id] = true;
    for (CellId id = 0; id < n; ++id) {
      TableCell const& c = table.cells[id];
      if (c.production.kind == ProductionKind::Leaf) continue;
      for (Read const& r : table.reads) {
        if (r.consumer != id) continue;
        TableCell const& s = table.cells[r.source];
        bool const ok = produced[r.source] && s.scope.encloses(c.scope);
        if (!ok)
          out.push_back({"visibility", detail::cell_str(table, r.source) +
                                           " not resident when " +
                                           detail::cell_str(table, id) +
                                           " is produced"});
      }
      if (c.production.kind == ProductionKind::Assemble) {
        CellId const src = c.production.source;
        if (src >= n || !produced[src])
          out.push_back({"visibility", detail::cell_str(table, id) +
                                           " assembles a cell not produced"});
      }
      produced[id] = true;
    }
  }

  // (2) form: per consumer and per loop group it is sliced by, every operand
  // bound to that group must be bound to the SAME instance; an operand whole
  // on a group that another operand is bound to is a mismatch.
  for (CellId id = 0; id < n; ++id) {
    TableCell const& c = table.cells[id];
    if (c.production.kind != ProductionKind::Build) continue;
    for (auto const& [cpos, L] : c.sliced) {
      (void)cpos;
      bool any_bound = false, any_whole = false;
      for (Read const& r : table.reads) {
        if (r.consumer != id) continue;
        TableCell const& s = table.cells[r.source];
        bool bound_here = false, bound_other = false;
        for (auto const& [p, k] : s.sliced)
          if (k.depth == L.depth)
            (k.loop_slot == L.loop_slot ? bound_here : bound_other) = true;
        for (auto const& [p, k] : r.slice)
          if (k.depth == L.depth)
            (k.loop_slot == L.loop_slot ? bound_here : bound_other) = true;
        if (bound_other)
          out.push_back({"form", detail::cell_str(table, r.source) +
                                     " bound to another instance of loop "
                                     "group " +
                                     std::to_string(L.depth) + " than " +
                                     detail::cell_str(table, id)});
        (bound_here ? any_bound : any_whole) = true;
      }
      if (any_bound && any_whole)
        out.push_back({"form", detail::cell_str(table, id) +
                                   ": one operand whole and another bound "
                                   "on loop group " +
                                   std::to_string(L.depth)});
    }
  }

  // (3) chain: Assemble source/scope/map consistency.
  for (CellId id = 0; id < n; ++id) {
    TableCell const& c = table.cells[id];
    if (c.production.kind != ProductionKind::Assemble) continue;
    CellId const src = c.production.source;
    if (src >= n) continue;  // reported by (1)
    TableCell const& s = table.cells[src];
    if (!(c.scope.encloses(s.scope) &&
          s.scope.path.size() > c.scope.path.size()))
      out.push_back({"chain", detail::cell_str(table, id) +
                                  " does not enclose its source scope"});
    if (c.production.assemble == AssembleKind::Scatter) {
      for (auto const& [pos, k] : c.production.scatter_map) {
        bool found = false;
        for (auto const& [sp, sk] : s.sliced)
          if (sp == pos && detail::same_key(sk, k)) found = true;
        if (!found)
          out.push_back({"chain", detail::cell_str(table, id) +
                                      " scatters position " +
                                      std::to_string(pos) +
                                      " from a source not sliced by that "
                                      "instance"});
      }
    } else {
      bool reduces_one = false;
      for (auto const& [sp, sk] : s.sliced) {
        bool kept = false;
        for (auto const& [cp, ck] : c.sliced)
          if (detail::same_key(sk, ck)) kept = true;
        if (!kept) reduces_one = true;
      }
      if (!reduces_one)
        out.push_back({"chain", detail::cell_str(table, id) +
                                    " sums a source sliced by no instance "
                                    "the sum removes"});
    }
  }

  // (4) life: reads (multiplicity carried in Read count) == life.
  for (CellId id = 0; id < n; ++id) {
    TableCell const& c = table.cells[id];
    if (c.production.kind == ProductionKind::Leaf) continue;
    std::size_t reads = 0;
    for (Read const& r : table.reads)
      if (r.source == id) ++reads;
    for (CellId o = 0; o < n; ++o)
      if (table.cells[o].production.kind == ProductionKind::Assemble &&
          table.cells[o].production.source == id)
        ++reads;
    if (reads != c.life)
      out.push_back({"life", detail::cell_str(table, id) + " life " +
                                 std::to_string(c.life) + " != reads " +
                                 std::to_string(reads)});
  }

  // (5) uniqueness: one Build per (value, scope).
  for (CellId a = 0; a < n; ++a)
    for (CellId b = a + 1; b < n; ++b) {
      TableCell const& x = table.cells[a];
      TableCell const& y = table.cells[b];
      if (x.production.kind == ProductionKind::Build &&
          y.production.kind == ProductionKind::Build &&
          x.value_id == y.value_id && x.scope == y.scope)
        out.push_back({"uniqueness", detail::cell_str(table, a) + " and " +
                                         detail::cell_str(table, b) +
                                         " build one value in one scope"});
    }
  return out;
}

inline void assert_valid_cell_table(CellTable const& table,
                                    ScopeBlock const& root) {
  auto const v = validate_cell_table(table, root);
  if (v.empty()) return;
  std::string msg =
      "cell table invalid (" + std::to_string(v.size()) + " violation(s)):";
  for (auto const& x : v) msg += "\n  [" + x.rule + "] " + x.what;
  throw std::runtime_error(msg);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_HPP
