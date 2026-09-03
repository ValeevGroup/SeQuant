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
  /// loop instances this cell is a partial sum over (a reduced axis has no
  /// carried position, so it cannot appear in \c sliced); empty = complete
  /// over every reduced axis.
  container::svector<LoopKey> partial_over;
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

/// The loop instances \p c is bound to: its sliced positions' instances plus
/// the instances it is a partial sum over (\c partial_over) -- the single
/// place every "is this cell bound to instance k" test should read from.
[[nodiscard]] inline container::svector<LoopKey> bound_instances(
    TableCell const& c) {
  container::svector<LoopKey> out;
  for (auto const& [pos, k] : c.sliced) {
    (void)pos;
    out.push_back(k);
  }
  for (LoopKey const& k : c.partial_over) out.push_back(k);
  return out;
}

/// The scope over which a produced cell \p s remains resident, once
/// produced: (1) a persistent cell is resident everywhere (the root scope);
/// (2) else, if \p s is bound (\c bound_instances) to one instance of a loop
/// that also appears in \c s.scope.path, it is resident only through that
/// loop's current batch, i.e. through the prefix of \c s.scope.path ending
/// at the DEEPEST such bound loop (it dies when that loop's batch ends); (3)
/// else (a volatile whole cell, including a produce_if_absent cell whole on
/// its innermost loop) it is resident exactly at \c s.scope.
[[nodiscard]] inline CellScope residency_scope(TableCell const& s) {
  if (s.persistent) return CellScope{};
  auto const bi = bound_instances(s);
  std::size_t deepest = 0;
  bool bound = false;
  for (std::size_t i = 0; i < s.scope.path.size(); ++i) {
    LoopKey const& here = s.scope.path[i].first;
    for (LoopKey const& k : bi) {
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

/// Multiplicity of one read of a source cell by a consumer cell: one factor
/// of \c max(1, n_batches_of(key)) per loop instance on the consumer's scope
/// path that is NOT also on the source's scope path -- a consumer nested one
/// extra loop deeper than its source re-reads that source once per batch of
/// each such loop. The single place both the builder (life bookkeeping) and
/// the validator (rule 4) compute this, so they never drift apart. A null \p
/// n_batches_of (the validator's default) is treated as returning 1
/// everywhere.
[[nodiscard]] inline std::size_t read_multiplicity(
    CellScope const& source_scope, CellScope const& consumer_scope,
    std::function<std::size_t(LoopKey const&)> const& n_batches_of) {
  std::size_t mult = 1;
  for (auto const& [pk, lat] : consumer_scope.path) {
    (void)lat;
    bool in_source = false;
    for (auto const& [sk, slat] : source_scope.path) {
      (void)slat;
      if (same_key(sk, pk)) in_source = true;
    }
    if (!in_source)
      mult *= std::max<std::size_t>(1, n_batches_of ? n_batches_of(pk) : 1);
  }
  return mult;
}
}  // namespace detail

/// Static validation of a cell table (spec section 3). Production order is
/// the order of non-Leaf cells in \p table.cells (Task 3's builder emits
/// them in execution order). The block tree \p root is accepted for the
/// Task 3 integration (it is not consulted here beyond being passed through).
[[nodiscard]] inline container::vector<CellViolation> validate_cell_table(
    CellTable const& table, ScopeBlock const& /*root*/,
    std::function<std::size_t(LoopKey const&)> const& n_batches_of = {}) {
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
        for (LoopKey const& k : detail::bound_instances(s))
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
    if (s.value_id != c.value_id)
      out.push_back({"chain", detail::cell_str(table, id) +
                                  " assembles a source of a different "
                                  "value"});
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
    } else if (!s.scope.path.empty()) {
      LoopKey const closing = s.scope.path.back().first;
      bool found = false;
      for (LoopKey const& k : s.partial_over)
        if (detail::same_key(k, closing)) found = true;
      if (!found)
        out.push_back({"chain", detail::cell_str(table, id) +
                                    " sums a source whose partial_over "
                                    "lacks the closing instance"});
    }
  }

  // (4) life: reads (weighted by the consumer's extra enclosing batches via
  // detail::read_multiplicity, same rule the builder uses) plus one per
  // Assemble that consumes the cell == life.
  for (CellId id = 0; id < n; ++id) {
    TableCell const& c = table.cells[id];
    if (c.production.kind == ProductionKind::Leaf) continue;
    std::size_t reads = 0;
    for (Read const& r : table.reads)
      if (r.source == id)
        reads += detail::read_multiplicity(
            c.scope, table.cells[r.consumer].scope, n_batches_of);
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

inline void assert_valid_cell_table(
    CellTable const& table, ScopeBlock const& root,
    std::function<std::size_t(LoopKey const&)> const& n_batches_of = {}) {
  auto const v = validate_cell_table(table, root, n_batches_of);
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

/// Whether block \p b lists \p ovid as an \c AccumulateSum output (a
/// reduction escape realized at \p b's own loop instance).
[[nodiscard]] inline bool has_accumulate_sum_output(ScopeBlock const& b,
                                                    std::size_t ovid) {
  for (auto const& [vid, kind] : b.outputs)
    if (vid == ovid && kind == OutputKind::AccumulateSum) return true;
  return false;
}

/// Sets the \c persistent / \c produce_if_absent flags of \p cell from its
/// bound instances (\c bound_instances) against the current \p path: a cell
/// is persistent iff it carries no volatile leaf and is bound to none of its
/// enclosing loops; it is produced only on first visit (\c produce_if_absent)
/// iff its scope is nested inside a loop it is NOT bound to on that loop's
/// own instance.
inline void set_residency_flags(
    TableCell& cell, bool volatile_value,
    container::svector<std::pair<LoopKey, int>> const& path) {
  bool bound_on_any_enclosing = false, bound_on_innermost = false;
  for (LoopKey const& k : bound_instances(cell)) {
    for (auto const& [pk, lat] : path)
      if (same_key(pk, k)) bound_on_any_enclosing = true;
    if (!path.empty() && same_key(path.back().first, k))
      bound_on_innermost = true;
  }
  cell.persistent = !volatile_value && !bound_on_any_enclosing;
  cell.produce_if_absent = !path.empty() && !bound_on_innermost;
}

/// Builds and registers the \c Build cell for value \p vid at the current
/// scope (\p path / \p path_spaces / \p block_stack): \c sliced via \c
/// slicing_instance for each of its own carried positions; \c partial_over =
/// the \c level.key() of every block on \p block_stack (its enclosing
/// blocks, innermost included) that lists \c (vid, AccumulateSum) among its
/// own outputs -- so a value with a genuine \c BuildStep gets the same
/// treatment as one synthesized for an escape with no in-block form of its
/// own; residency flags via \c set_residency_flags. Returns the new cell's
/// id.
[[nodiscard]] inline CellId emit_build_cell(
    CellTableInputs const& in, std::size_t vid,
    container::svector<std::pair<LoopKey, int>> const& path,
    container::svector<std::wstring> const& path_spaces,
    container::svector<ScopeBlock const*> const& block_stack,
    CellBuildState& st) {
  RichSchedule const& rich = *in.rich;
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
  for (ScopeBlock const* b : block_stack)
    if (has_accumulate_sum_output(*b, vid))
      cell.partial_over.push_back(b->level.key());
  set_residency_flags(cell, in.volatile_of(vid), path);
  st.table.cells.push_back(std::move(cell));
  st.forms_of[vid].push_back(st.table.cells.size() - 1);
  return st.table.cells.size() - 1;
}

inline void emit_cells(CellTableInputs const& in, ScopeBlock const& block,
                       container::svector<std::pair<LoopKey, int>>& path,
                       container::svector<std::wstring>& path_spaces,
                       container::svector<ScopeBlock const*>& block_stack,
                       CellBuildState& st) {
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value)) {
      (void)emit_build_cell(in, b->value_id, path, path_spaces, block_stack,
                            st);
      continue;
    }
    auto const& child = std::get<ScopeBlock>(step.value);
    path.push_back({child.level.key(), child.latitude_ordinal});
    path_spaces.push_back(std::wstring{child.axis.space().base_key()});
    block_stack.push_back(&child);
    emit_cells(in, child, path, path_spaces, block_stack, st);
    // the child's outputs assemble at THIS scope
    CellScope const child_scope{path};  // path currently has child on top
    for (auto const& [ovid, okind] : child.outputs) {
      LoopKey const inst = child.level.key();
      // a form of ovid registered EXACTLY at the child's own scope (a
      // BuildStep there, an implicit Build synthesized below, or an
      // Assemble registered at the child's scope for a grandchild escape)
      // is ovid's own cell inside the child; a sibling block's forms live
      // at a different scope (a different loop_slot or latitude even at
      // the same depth), so full CellScope equality -- not just depth --
      // rules those out.
      CellId src = 0;
      bool found = false;
      if (auto const fit = st.forms_of.find(ovid); fit != st.forms_of.end())
        for (CellId id : fit->second)
          if (st.table.cells[id].scope == child_scope) {
            src = id;
            found = true;
            break;
          }
      if (!found) {
        // ovid never got its own cell inside the child: its production
        // fuses the reduction/scatter with an operand contraction (the
        // escaped value's hash differs from any single operand's by
        // exactly the closing axis), so the legacy schedule never emits a
        // plain BuildStep for it. The per-batch partial is nonetheless a
        // real form of ovid itself, realized once per batch at the child's
        // own scope; synthesize its Build cell here.
        src = emit_build_cell(in, ovid, path, path_spaces, block_stack, st);
      }
      TableCell const& s = st.table.cells[src];
      // Unreachable by construction: every entry under forms_of[ovid] is a
      // cell this builder itself registered under key ovid with
      // value_id == ovid (emit_build_cell's own cell, or an Assemble's
      // a.value_id = ovid below); kept as a loud guard against a future
      // emission bug rather than a silent cross-value assemble.
      if (s.value_id != ovid)
        throw std::logic_error(
            "cell table: value_id mismatch between escaping value " +
            std::to_string(ovid) + " and its registered form (value_id " +
            std::to_string(s.value_id) + ")");
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
      if (a.production.assemble == AssembleKind::Sum) {
        for (LoopKey const& k : s.partial_over)
          if (!detail::same_key(k, inst)) a.partial_over.push_back(k);
      } else {
        a.partial_over = s.partial_over;
      }
      path.pop_back();
      path_spaces.pop_back();
      block_stack.pop_back();
      a.scope.path = path;
      set_residency_flags(a, in.volatile_of(ovid), path);
      st.table.cells.push_back(std::move(a));
      st.forms_of[ovid].push_back(st.table.cells.size() - 1);
      path.push_back({child.level.key(), child.latitude_ordinal});
      path_spaces.push_back(std::wstring{child.axis.space().base_key()});
      block_stack.push_back(&child);
    }
    path.pop_back();
    path_spaces.pop_back();
    block_stack.pop_back();
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
  if (!in.ordered || !in.rich || !in.sliced_modes_of || !in.volatile_of)
    throw std::invalid_argument(
        "build_cell_table: ordered, rich, sliced_modes_of and volatile_of "
        "are all required");
  detail::CellBuildState st;
  container::svector<std::pair<LoopKey, int>> path;
  container::svector<std::wstring> path_spaces;
  container::svector<ScopeBlock const*> block_stack;
  detail::emit_cells(in, in.ordered->root, path, path_spaces, block_stack, st);
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

  // Reads: one per (consumer Build cell, operand value) edge.
  auto const key_of_lid = [&](LoopId lid) {
    return in.sliced->levels[lid].key();
  };
  for (CellId cid = 0; cid < st.table.cells.size(); ++cid) {
    TableCell const& c = st.table.cells[cid];
    if (c.production.kind != ProductionKind::Build) continue;
    auto const dit = g.depends_on.find(c.value_id);
    if (dit == g.depends_on.end()) continue;
    for (std::size_t op : dit->second) {
      Read r;
      r.consumer = cid;
      r.operand_value_id = op;
      auto const fit = st.forms_of.find(op);
      if (fit == st.forms_of.end() || fit->second.empty()) continue;
      // Source = among this operand's forms, the one whose RESIDENCY scope
      // (detail::residency_scope: root if persistent, else the prefix of
      // its scope ending at the deepest loop it is bound to, else its own
      // scope) encloses the consumer's scope, deepest such form wins; if
      // none is resident, fall back to the LAST form so the validator's
      // visibility rule surfaces the gap rather than the builder hiding it.
      std::optional<CellId> best;
      for (CellId f : fit->second) {
        TableCell const& fc = st.table.cells[f];
        if (!detail::residency_scope(fc).encloses(c.scope)) continue;
        if (!best ||
            fc.scope.path.size() > st.table.cells[*best].scope.path.size())
          best = f;
      }
      r.source = best ? *best : fit->second.back();
      TableCell const& s = st.table.cells[r.source];
      for (auto const& [w_vid, pos, lid, consumer_vid] : in.sliced->occ_facts) {
        if (w_vid != op || consumer_vid != c.value_id) continue;
        LoopKey const k = key_of_lid(lid);
        bool enclosing = false;
        for (auto const& [pk, lat] : c.scope.path)
          if (detail::same_key(pk, k)) enclosing = true;
        if (!enclosing) continue;
        bool already = false;
        for (auto const& [sp, sk] : s.sliced)
          if (sp == pos && detail::same_key(sk, k)) already = true;
        for (LoopKey const& pk : s.partial_over)
          if (detail::same_key(pk, k)) already = true;
        if (!already) r.slice.push_back({pos, k});
      }
      st.table.reads.push_back(std::move(r));
    }
  }
  // Lives: reads weighted by the consumer's extra enclosing batches
  // (detail::read_multiplicity), plus one per Assemble that consumes the
  // cell.
  for (TableCell& c : st.table.cells) c.life = 0;
  for (Read const& r : st.table.reads) {
    TableCell const& s = st.table.cells[r.source];
    TableCell const& c = st.table.cells[r.consumer];
    st.table.cells[r.source].life +=
        detail::read_multiplicity(s.scope, c.scope, in.n_batches_of);
  }
  for (TableCell const& c : st.table.cells)
    if (c.production.kind == ProductionKind::Assemble)
      st.table.cells[c.production.source].life += 1;
  return std::move(st.table);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_HPP
