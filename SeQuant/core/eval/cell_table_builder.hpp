#ifndef SEQUANT_EVAL_CELL_TABLE_BUILDER_HPP
#define SEQUANT_EVAL_CELL_TABLE_BUILDER_HPP

#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>

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

/// Inputs to \c build_cell_table: the finished ordered schedule plus the
/// per-value callbacks the table derivation needs but does not itself own
/// the source of (sliced-mode set, volatility, batch count of a loop
/// instance).
struct CellTableInputs {
  OrderedSchedule const* ordered = nullptr;
  RichSchedule const* rich = nullptr;
  /// the per-occurrence seam facts the reads loop consults; required.
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
  if (!in.ordered || !in.rich || !in.sliced || !in.sliced_modes_of ||
      !in.volatile_of)
    throw std::invalid_argument(
        "build_cell_table: ordered, rich, sliced, sliced_modes_of and "
        "volatile_of are all required");
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

  // Reads: one per (consumer Build cell, operand value) edge. The
  // occ_facts/occ_invariant scans below are linear per read (O(reads x
  // facts)); acceptable at this stage.
  auto const key_of_lid = [&](LoopId lid) {
    return in.sliced->levels[lid].key();
  };
  for (CellId cid = 0; cid < st.table.cells.size(); ++cid) {
    TableCell const& c = st.table.cells[cid];
    if (c.production.kind != ProductionKind::Build) continue;
    auto const dit = g.depends_on.find(c.value_id);
    if (dit == g.depends_on.end()) continue;
    auto const c_bound = detail::bound_instances(c);
    for (std::size_t op : dit->second) {
      auto const fit = st.forms_of.find(op);
      if (fit == st.forms_of.end() || fit->second.empty())
        throw std::logic_error("cell table: operand value " +
                               std::to_string(op) + " of consumer value " +
                               std::to_string(c.value_id) +
                               " has no registered form");
      Read r;
      r.consumer = cid;
      r.operand_value_id = op;
      // Source = among this operand's forms, the one whose RESIDENCY scope
      // (detail::residency_scope: the prefix of its scope ending at the
      // deepest loop it is bound to, else -- a whole form -- the root
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
        if (already) continue;
        bool dup = false;
        for (auto const& [ep, ek] : r.slice)
          if (ep == pos && detail::same_key(ek, k)) dup = true;
        if (!dup) r.slice.push_back({pos, k});
      }
      // Explicit invariant records: the seam found no slicing decision
      // binding this consumer's read of op on a loop the CONSUMER itself is
      // bound to, so the form rule accepts it as whole even though another
      // operand may be bound to that same loop (see Read::invariant_on).
      for (auto const& [w_vid, lid, consumer_vid] : in.sliced->occ_invariant) {
        if (w_vid != op || consumer_vid != c.value_id) continue;
        LoopKey const k = key_of_lid(lid);
        bool in_consumer_bound = false;
        for (LoopKey const& ck : c_bound)
          if (detail::same_key(ck, k)) in_consumer_bound = true;
        if (!in_consumer_bound) continue;
        bool dup = false;
        for (LoopKey const& ek : r.invariant_on)
          if (detail::same_key(ek, k)) dup = true;
        if (!dup) r.invariant_on.push_back(k);
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
        detail::read_multiplicity(s, c.scope, in.n_batches_of);
  }
  for (TableCell const& c : st.table.cells)
    if (c.production.kind == ProductionKind::Assemble)
      st.table.cells[c.production.source].life += 1;
  return std::move(st.table);
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_BUILDER_HPP
