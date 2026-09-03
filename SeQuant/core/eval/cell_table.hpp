#ifndef SEQUANT_EVAL_CELL_TABLE_HPP
#define SEQUANT_EVAL_CELL_TABLE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

/// The scope over which a produced cell \p s remains resident, once
/// produced: (1) a persistent cell is resident everywhere (the root scope);
/// (2) else, if \p s is bound (via \c s.sliced) to one instance of a loop
/// that also appears in \c s.scope.path, it is resident only through that
/// loop's current batch, i.e. through the prefix of \c s.scope.path ending
/// at the DEEPEST such bound loop (it dies when that loop's batch ends); (3)
/// else (a volatile whole cell, including a produce_if_absent cell whole on
/// its innermost loop) it is resident exactly at \c s.scope.
[[nodiscard]] inline CellScope residency_scope(TableCell const& s) {
  if (s.persistent) return CellScope{};
  std::size_t deepest = 0;
  bool bound = false;
  for (std::size_t i = 0; i < s.scope.path.size(); ++i) {
    LoopKey const& here = s.scope.path[i].first;
    for (auto const& [pos, k] : s.sliced) {
      (void)pos;
      if (same_key(k, here)) {
        bound = true;
        deepest = i;
      }
    }
  }
  if (!bound) return s.scope;
  CellScope r;
  r.path.assign(s.scope.path.begin(), s.scope.path.begin() + deepest + 1);
  return r;
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
  // the source was produced earlier (or is a Leaf) and the source's
  // RESIDENCY scope (detail::residency_scope; persistent -> root, bound to
  // an enclosing loop instance -> the prefix of its scope ending at that
  // loop, else its own scope) encloses the consumer's scope.
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
        bool const ok =
            produced[r.source] && detail::residency_scope(s).encloses(c.scope);
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
    for (auto const& entry : c.sliced) {
      LoopKey const& L = entry.second;
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
        if (bound_here)
          any_bound = true;
        else if (!bound_other)
          any_whole = true;
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

/// Inputs to \c build_cell_table: the finished ordered schedule plus the
/// per-value callbacks the table derivation needs but does not itself own
/// the source of (sliced-mode set, volatility, batch count of a loop
/// instance).
struct CellTableInputs {
  OrderedSchedule const* ordered = nullptr;
  RichSchedule const* rich = nullptr;
  LegalitySchedule const* legality = nullptr;
  SlicedModeAssignment const* sliced = nullptr;
  std::function<container::svector<Index>(std::size_t)> sliced_modes_of;
  std::function<bool(std::size_t)> volatile_of;
  std::function<std::size_t(LoopKey const&)> n_batches_of;
};

namespace detail {
struct CellBuildState {
  CellTable table;
  // per value: the cell id of its form visible at each scope depth during the
  // walk (Build inside a block, or the Assemble at that block's parent)
  std::unordered_map<std::size_t, container::svector<CellId>> forms_of;
};

/// The loop instance among \p path (with its axis spaces) that slices carried
/// position \p p of value \p vid: the production occurrence's fusion slot for
/// p, matched against an enclosing entry of the same index space.
[[nodiscard]] inline std::optional<LoopKey> slicing_instance(
    RichSchedule const& rich, std::size_t vid, std::size_t p,
    container::svector<std::pair<LoopKey, int>> const& path,
    container::svector<std::wstring> const& path_spaces) {
  ValueCell const& c = rich.cells[vid];
  if (c.occurrences.empty()) return std::nullopt;
  OccurrenceRec const& prod = c.occurrences.front();
  if (p >= prod.loop_slot.size() || prod.loop_slot[p] < 0) return std::nullopt;
  std::wstring const space{c.carried[p].space().base_key()};
  for (std::size_t i = 0; i < path.size(); ++i)
    if (path_spaces[i] == space && path[i].first.loop_slot == prod.loop_slot[p])
      return path[i].first;
  return std::nullopt;
}

/// The cell to assemble FROM for an escaping value \p ovid closing loop
/// instance \p inst: \p ovid's own cell inside the closing block if one was
/// registered there (a Build, or a nested escape already assembled at this
/// scope); otherwise \p ovid never got its own cell here (its production
/// fuses the reduction/scatter with an operand contraction, so the escaped
/// value's hash differs from any single operand's by exactly the closing
/// axis) and the true partial is whichever DIRECT operand's own registered
/// cell is itself sliced by \p inst.
[[nodiscard]] inline std::optional<CellId> resolve_escape_source(
    CellTableInputs const& in, CellBuildState const& st, std::size_t ovid,
    LoopKey const& inst) {
  if (auto const fit = st.forms_of.find(ovid);
      fit != st.forms_of.end() && !fit->second.empty())
    return fit->second.back();
  auto const oit = in.ordered->operand_vids.find(ovid);
  if (oit == in.ordered->operand_vids.end()) return std::nullopt;
  for (std::size_t op : oit->second) {
    auto const fit = st.forms_of.find(op);
    if (fit == st.forms_of.end() || fit->second.empty()) continue;
    CellId const cand = fit->second.back();
    for (auto const& [p, k] : st.table.cells[cand].sliced)
      if (detail::same_key(k, inst)) return cand;
  }
  return std::nullopt;
}

inline void emit_cells(CellTableInputs const& in, ScopeBlock const& block,
                       container::svector<std::pair<LoopKey, int>>& path,
                       container::svector<std::wstring>& path_spaces,
                       CellBuildState& st) {
  RichSchedule const& rich = *in.rich;
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value)) {
      std::size_t const vid = b->value_id;
      TableCell cell;
      cell.value_id = vid;
      cell.scope.path = path;
      cell.production.kind = ProductionKind::Build;
      auto const sm = in.sliced_modes_of(vid);
      ValueCell const& vc = rich.cells[vid];
      for (std::size_t p = 0; p < vc.carried.size(); ++p) {
        bool const in_sm =
            std::find(sm.begin(), sm.end(), vc.carried[p]) != sm.end();
        if (!in_sm) continue;
        if (auto k = slicing_instance(rich, vid, p, path, path_spaces))
          cell.sliced.push_back({p, *k});
      }
      bool sliced_on_any_enclosing = false, sliced_on_innermost = false;
      for (auto const& [p, k] : cell.sliced) {
        for (auto const& [pk, lat] : path)
          if (detail::same_key(pk, k)) sliced_on_any_enclosing = true;
        if (!path.empty() && detail::same_key(path.back().first, k))
          sliced_on_innermost = true;
      }
      cell.persistent = !in.volatile_of(vid) && !sliced_on_any_enclosing;
      cell.produce_if_absent = !path.empty() && !sliced_on_innermost;
      st.table.cells.push_back(std::move(cell));
      st.forms_of[vid].push_back(st.table.cells.size() - 1);
      continue;
    }
    auto const& child = std::get<ScopeBlock>(step.value);
    path.push_back({child.level.key(), child.latitude_ordinal});
    path_spaces.push_back(std::wstring{child.axis.space().base_key()});
    emit_cells(in, child, path, path_spaces, st);
    // the child's outputs assemble at THIS scope
    for (auto const& [ovid, okind] : child.outputs) {
      LoopKey const inst = child.level.key();
      auto const src_opt = resolve_escape_source(in, st, ovid, inst);
      if (!src_opt) continue;
      CellId const src = *src_opt;
      TableCell const& s = st.table.cells[src];
      TableCell a;
      a.value_id = ovid;
      a.production.kind = ProductionKind::Assemble;
      a.production.assemble = okind == OutputKind::AccumulateSum
                                  ? AssembleKind::Sum
                                  : AssembleKind::Scatter;
      a.production.source = src;
      for (auto const& [p, k] : s.sliced) {
        if (detail::same_key(k, inst)) {
          if (a.production.assemble == AssembleKind::Scatter)
            a.production.scatter_map.push_back({p, k});
        } else {
          a.sliced.push_back({p, k});
        }
      }
      path.pop_back();
      path_spaces.pop_back();
      a.scope.path = path;
      bool any_enclosing = false, innermost = false;
      for (auto const& [p, k] : a.sliced) {
        for (auto const& [pk, lat] : path)
          if (detail::same_key(pk, k)) any_enclosing = true;
        if (!path.empty() && detail::same_key(path.back().first, k))
          innermost = true;
      }
      a.persistent = !in.volatile_of(ovid) && !any_enclosing;
      a.produce_if_absent = !path.empty() && !innermost;
      st.table.cells.push_back(std::move(a));
      st.forms_of[ovid].push_back(st.table.cells.size() - 1);
      path.push_back({child.level.key(), child.latitude_ordinal});
      path_spaces.push_back(std::wstring{child.axis.space().base_key()});
    }
    path.pop_back();
    path_spaces.pop_back();
  }
}
}  // namespace detail

/// Derives the cell table's PRODUCTIONS from an ordered schedule already
/// built by \c build_ordered_schedule: one \c Build cell per \c BuildStep,
/// one \c Assemble cell per block output entry (at the block's parent
/// scope), and one \c Leaf cell for every leaf value that is an operand of
/// some computed value. Reads and lives are left empty/zero here; Task 4
/// fills them in.
[[nodiscard]] inline CellTable build_cell_table(CellTableInputs const& in) {
  detail::CellBuildState st;
  container::svector<std::pair<LoopKey, int>> path;
  container::svector<std::wstring> path_spaces;
  detail::emit_cells(in, in.ordered->root, path, path_spaces, st);
  // Leaf cells: every leaf value that is an operand of some value.
  auto const g = detail::ordered_schedule_dep_graph(*in.rich);
  std::unordered_set<std::size_t> leaf_done;
  for (auto const& [vid, ops] : g.depends_on)
    for (std::size_t op : ops)
      if (in.rich->cells[op].is_leaf && leaf_done.insert(op).second) {
        TableCell leaf;
        leaf.value_id = op;
        leaf.production.kind = ProductionKind::Leaf;
        leaf.persistent = !in.volatile_of(op);
        st.table.cells.push_back(std::move(leaf));
        st.forms_of[op].push_back(st.table.cells.size() - 1);
      }
  return std::move(st.table);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_HPP
