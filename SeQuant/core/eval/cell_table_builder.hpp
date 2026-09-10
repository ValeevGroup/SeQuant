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
  /// OPTIONAL: the value's direct operand value ids WITH REPETITION -- one
  /// entry per LEG of its production tree, so a value contracted with itself
  /// appears twice. Both \c ordered_schedule_dep_graph's \c depends_on and
  /// the \c OrderedSchedule::operand_vids copied from it DE-DUPLICATE their
  /// operand lists, and a de-duplicated list loses exactly the read the
  /// runtime does perform (a home access per leg); a caller that can see the
  /// production tree (the forest node's two children, resolved back to value
  /// ids) supplies it here. When unset, the builder falls back to the
  /// de-duplicated \c depends_on and a self-contracting consumer gets one
  /// read instead of two.
  std::function<container::svector<std::size_t>(std::size_t)> operands_of;
};

namespace detail {
struct CellBuildState {
  CellTable table;
  // per value: the cell id of its form visible at each scope depth during the
  // walk (Build inside a block, or the Assemble at that block's parent)
  std::unordered_map<std::size_t, container::svector<CellId>> forms_of;
};

/// The loop instance among \p path (with its axis spaces) that slices carried
/// position \p p of value \p vid: an enclosing entry whose index space is the
/// position's own AND whose slot is the fusion slot ANY occurrence of the
/// value records for that position. Every occurrence is tested, not just the
/// production one, because that is the runtime's own test: a value CSE-shared
/// by several consumers is sliced on this loop if any of its occurrences was
/// bound to that slot, and the first occurrence is not privileged.
///
/// \note The match is on (space, slot), and so relies on fusion never giving
/// two levels of the SAME space one slot -- were that to happen, a position
/// could match the wrong level of its own space. Outermost enclosing entry
/// wins (\p path is outermost-first).
[[nodiscard]] inline std::optional<LoopKey> slicing_instance(
    RichSchedule const& rich, std::size_t vid, std::size_t p,
    container::svector<std::pair<LoopKey, int>> const& path,
    container::svector<std::wstring> const& path_spaces) {
  ValueCell const& c = rich.cells[vid];
  if (c.occurrences.empty() || p >= c.carried.size()) return std::nullopt;
  std::wstring const space{c.carried[p].space().base_key()};
  for (std::size_t i = 0; i < path.size(); ++i) {
    if (path_spaces[i] != space) continue;
    for (OccurrenceRec const& occ : c.occurrences)
      if (p < occ.loop_slot.size() && occ.loop_slot[p] >= 0 &&
          occ.loop_slot[p] == path[i].first.loop_slot)
        return path[i].first;
  }
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

/// Sets the \c produce_if_absent flag and the persistence CANDIDACY of \p
/// cell from its bound instances (\c bound_instances) against the current \p
/// path: a cell is a persistence candidate iff it carries no volatile leaf and
/// is bound to none of its enclosing loops; it is produced only on first visit
/// (\c produce_if_absent) iff its scope is nested inside a loop it is NOT
/// bound to on that loop's own instance.
///
/// \note Candidacy, not the final answer: persistence is the FRONTIER of the
/// invariant region, not all of it, and the remaining condition (at least one
/// volatile consumer, or no consumer at all) can only be decided once the
/// table's reads exist. \c apply_persistence_frontier demotes the rest at the
/// end of \c build_cell_table.
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
/// own; residency flags via \c set_residency_flags. A sliced mode whose
/// position matches no enclosing instance is recorded WHOLE and reported in
/// \c CellTable::unresolved. Returns the new cell's id.
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
  container::svector<std::size_t> unresolved_positions;
  ValueCell const& vc = rich.cells[vid];
  for (std::size_t p = 0; p < vc.carried.size(); ++p) {
    bool const in_sm =
        std::find(sm.begin(), sm.end(), vc.carried[p]) != sm.end();
    if (!in_sm) continue;
    if (auto k = slicing_instance(rich, vid, p, path, path_spaces))
      cell.sliced.push_back({p, *k});
    else
      unresolved_positions.push_back(p);
  }
  for (ScopeBlock const* b : block_stack)
    if (has_accumulate_sum_output(*b, vid))
      cell.partial_over.push_back(b->level.key());
  set_residency_flags(cell, in.volatile_of(vid), path);
  st.table.cells.push_back(std::move(cell));
  CellId const id = st.table.cells.size() - 1;
  for (std::size_t p : unresolved_positions)
    st.table.unresolved.push_back({id, p});
  st.forms_of[vid].push_back(id);
  return id;
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
    // The child's outputs assemble at THIS (the parent's) scope, while an
    // implicit Build synthesized for one of them belongs at the CHILD's. Leave
    // the child ONCE here, keeping a snapshot of the child-inclusive scope,
    // rather than pushing and popping all three stacks per output entry.
    CellScope const child_scope{path};  // path still has child on top
    auto const child_path = path;
    auto const child_spaces = path_spaces;
    auto const child_stack = block_stack;
    path.pop_back();
    path_spaces.pop_back();
    block_stack.pop_back();
    for (auto const& [ovid, okind] : child.outputs) {
      LoopKey const inst = child.level.key();
      // The child's IN-BLOCK form of ovid: the LAST form registered at a scope
      // the child's own scope ENCLOSES (the child's path is a prefix of it) --
      // the child's own scope qualifies, and so does any scope deeper inside
      // it. Enclosure, not equality: an escape chain may SKIP a level the
      // value is INVARIANT to, so the form that is live at the child's close
      // can be one the child's grandchild registered (an Assemble one or more
      // levels further in), and the runtime's home walk carries it out through
      // the skipped level. Taking the LAST such form takes the outermost
      // already-assembled one when several exist. A SIBLING block's forms live
      // at a scope the child does not enclose (a different loop_slot or
      // latitude even at the same depth), so they are still ruled out.
      CellId src = 0;
      bool found = false;
      if (auto const fit = st.forms_of.find(ovid); fit != st.forms_of.end())
        for (CellId id : fit->second)
          if (child_scope.encloses(st.table.cells[id].scope)) {
            src = id;
            found = true;
          }
      if (!found) {
        // ovid has NO form anywhere inside the child: its production
        // fuses the reduction/scatter with an operand contraction (the
        // escaped value's hash differs from any single operand's by
        // exactly the closing axis), so the legacy schedule never emits a
        // plain BuildStep for it. The per-batch partial is nonetheless a
        // real form of ovid itself, realized once per batch at the child's
        // own scope; synthesize its Build cell here.
        src = emit_build_cell(in, ovid, child_path, child_spaces, child_stack,
                              st);
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
      a.scope.path = path;
      set_residency_flags(a, in.volatile_of(ovid), path);
      st.table.cells.push_back(std::move(a));
      st.forms_of[ovid].push_back(st.table.cells.size() - 1);
    }
  }
}
/// Demotes every persistence CANDIDATE (\c set_residency_flags: no volatile
/// leaf, bound to no enclosing loop) that is not on the FRONTIER of the
/// invariant region, so that
///
///   \c persistent == carries no volatile leaf
///                 AND is bound to no enclosing loop instance
///                 AND (some CONSUMER of it holds a volatile value
///                      OR it has no consumer at all -- a forest root;
///                      validator rule 4 admits a zero-consumer cell only at
///                      the root scope).
///
/// The CONSUMERS of a cell are the consumer cells of every \c Read whose \c
/// source is it, plus every \c Assemble whose \c production.source is it --
/// the SAME edge set the runtime's cache-halt closure walks, by SOURCE CELL
/// and not by value, so the two cannot disagree about what reads a form.
///
/// Why the frontier and not the whole invariant region: a persistent cell is
/// held across evaluations and skipped by cache-halt on every later one. A
/// non-volatile cell all of whose consumers are non-volatile is therefore
/// never READ again -- each such consumer is itself either persistent-and-held
/// (skipped by cache-halt rule 1) or demoted here and, by induction from the
/// roots, skipped by rule 2 because ALL of ITS consumers are skipped. Keeping
/// such a cell resident across evaluations buys nothing and costs its bytes in
/// the persistent value store for the life of the cache handle. Only the cells
/// feeding volatile work -- and the forest's own results -- have a reader on a
/// later evaluation.
///
/// Leaves follow the same rule: an input tensor whose every consumer is
/// non-volatile is re-fetched on the (rare) evaluation that actually needs it
/// rather than held forever.
inline void apply_persistence_frontier(
    CellTable& table, std::function<bool(std::size_t)> const& volatile_of) {
  std::size_t const n = table.cells.size();
  container::vector<char> has_consumer(n, 0), has_volatile_consumer(n, 0);
  // volatile_of walks the value's whole production subtree; memoize it, since
  // one value is the consumer of many reads.
  std::unordered_map<std::size_t, bool> vol_of;
  auto const is_vol = [&](std::size_t vid) {
    auto const it = vol_of.find(vid);
    if (it != vol_of.end()) return it->second;
    bool const v = volatile_of(vid);
    vol_of.emplace(vid, v);
    return v;
  };
  auto const note = [&](CellId source, CellId consumer) {
    if (source >= n || consumer >= n) return;
    has_consumer[source] = 1;
    if (is_vol(table.cells[consumer].value_id))
      has_volatile_consumer[source] = 1;
  };
  for (Read const& r : table.reads) note(r.source, r.consumer);
  for (CellId c = 0; c < n; ++c)
    if (table.cells[c].production.kind == ProductionKind::Assemble)
      note(table.cells[c].production.source, c);
  for (CellId c = 0; c < n; ++c) {
    TableCell& cell = table.cells[c];
    if (cell.persistent && has_consumer[c] && !has_volatile_consumer[c])
      cell.persistent = false;
  }
}
}  // namespace detail

/// Derives the cell table's PRODUCTIONS from an ordered schedule already
/// built by \c build_ordered_schedule: one \c Build cell per \c BuildStep,
/// one \c Assemble cell per block output entry (at the block's parent
/// scope, sourced from the deepest form of that value registered inside the
/// block -- an escape chain may SKIP a level the value is invariant to), and
/// one \c Leaf cell for every leaf value that is an operand of
/// some computed value; then one \c Read per leg of every Build cell's
/// production tree, and every cell's \c life from those reads (weighted by
/// batch multiplicity) plus the Assembles that consume it. Sliced-mode
/// positions that match no enclosing loop instance are reported in \c
/// CellTable::unresolved.
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
  auto const g = detail::ordered_schedule_dep_graph(*in.rich);
  // The per-LEG operand list of a value (see CellTableInputs::operands_of),
  // falling back to the de-duplicated dependency graph when the caller cannot
  // supply one. Every operand question below -- which leaves are used, and
  // which reads exist -- goes through this one accessor, so the two never
  // disagree about what an operand is.
  auto const operands_of =
      [&](std::size_t vid) -> container::svector<std::size_t> {
    if (in.operands_of) return in.operands_of(vid);
    auto const it = g.depends_on.find(vid);
    if (it == g.depends_on.end()) return {};
    return it->second;
  };
  // Leaf cells: every leaf value that is an operand of some value. Emitted by
  // walking rich.cells in VALUE-ID order (the dependency map's iteration order
  // is unspecified), so cell ids are a deterministic function of the schedule.
  std::unordered_set<std::size_t> used_as_operand;
  for (ValueCell const& vc : in.rich->cells)
    for (std::size_t op : operands_of(vc.value_id)) used_as_operand.insert(op);
  for (ValueCell const& vc : in.rich->cells)
    if (vc.is_leaf && used_as_operand.count(vc.value_id)) {
      TableCell leaf;
      leaf.value_id = vc.value_id;
      leaf.production.kind = ProductionKind::Leaf;
      // Candidacy only, like every other cell's: apply_persistence_frontier
      // below demotes a leaf whose every consumer is non-volatile.
      leaf.persistent = !in.volatile_of(vc.value_id);
      st.table.cells.push_back(std::move(leaf));
      st.forms_of[vc.value_id].push_back(st.table.cells.size() - 1);
    }

  // Reads: one per (consumer Build cell, production-tree LEG). The
  // occ_facts/occ_invariant scans below are linear per read (O(reads x
  // facts)); acceptable at this stage.
  auto const key_of_lid = [&](LoopId lid) {
    return in.sliced->levels[lid].key();
  };
  for (CellId cid = 0; cid < st.table.cells.size(); ++cid) {
    TableCell const& c = st.table.cells[cid];
    if (c.production.kind != ProductionKind::Build) continue;
    auto const c_bound = detail::bound_instances(c);
    auto const ops = operands_of(c.value_id);
    std::size_t constexpr npos = static_cast<std::size_t>(-1);
    for (std::size_t oi = 0; oi < ops.size(); ++oi) {
      std::size_t const op = ops[oi];
      // The k-th read of value `op` by this consumer pairs with the k-th
      // distinct LEG the seam recorded for (op, consumer): a value on both
      // legs of one node has per-leg facts; a single-leg operand (or a seam
      // that recorded no leg) keeps every fact.
      std::size_t k = 0;
      for (std::size_t oj = 0; oj < oi; ++oj)
        if (ops[oj] == op) ++k;
      container::svector<std::size_t> legs;
      for (auto const& [w_vid, pos, lid, consumer_vid, f_leg] :
           in.sliced->occ_facts) {
        (void)pos;
        (void)lid;
        if (w_vid != op || consumer_vid != c.value_id || f_leg == npos)
          continue;
        if (std::find(legs.begin(), legs.end(), f_leg) == legs.end())
          legs.push_back(f_leg);
      }
      for (auto const& [w_vid, lid, consumer_vid, f_leg] :
           in.sliced->occ_invariant) {
        (void)lid;
        if (w_vid != op || consumer_vid != c.value_id || f_leg == npos)
          continue;
        if (std::find(legs.begin(), legs.end(), f_leg) == legs.end())
          legs.push_back(f_leg);
      }
      std::sort(legs.begin(), legs.end());
      std::optional<std::size_t> const my_leg =
          k < legs.size() ? std::optional<std::size_t>{legs[k]} : std::nullopt;
      auto const leg_matches = [&](std::size_t f_leg) {
        return !my_leg || f_leg == npos || f_leg == *my_leg;
      };
      auto const fit = st.forms_of.find(op);
      if (fit == st.forms_of.end() || fit->second.empty())
        throw std::logic_error("cell table: operand value " +
                               std::to_string(op) + " of consumer value " +
                               std::to_string(c.value_id) +
                               " has no registered form");
      Read r;
      r.consumer = cid;
      r.operand_value_id = op;
      // Source = the operand's form VISIBLE at the consumer's scope, by the
      // one shared rule (detail::deepest_visible_form, which the runtime's
      // CellRegistry::cell_of uses too): deepest resident form wins, ties
      // broken by emission order. If none is resident, fall back to the LAST
      // form so the validator's visibility rule surfaces the gap rather than
      // the builder hiding it.
      auto const best =
          detail::deepest_visible_form(st.table, fit->second, c.scope);
      r.source = best ? *best : fit->second.back();
      TableCell const& s = st.table.cells[r.source];
      for (auto const& [w_vid, pos, lid, consumer_vid, f_leg] :
           in.sliced->occ_facts) {
        if (w_vid != op || consumer_vid != c.value_id) continue;
        if (!leg_matches(f_leg)) continue;
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
      for (auto const& [w_vid, lid, consumer_vid, f_leg] :
           in.sliced->occ_invariant) {
        if (w_vid != op || consumer_vid != c.value_id) continue;
        if (!leg_matches(f_leg)) continue;
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
  // Persistence is the FRONTIER of the invariant region, and the frontier is
  // only knowable once the reads exist: demote every candidate whose every
  // consumer is itself non-volatile (see apply_persistence_frontier).
  detail::apply_persistence_frontier(st.table, in.volatile_of);
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
