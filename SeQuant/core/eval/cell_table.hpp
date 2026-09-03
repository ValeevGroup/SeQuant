#ifndef SEQUANT_EVAL_CELL_TABLE_HPP
#define SEQUANT_EVAL_CELL_TABLE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace sequant::eval {

using LoopKey = sequant::LoopKey;

/// The block tree of an ordered schedule (ordered_schedule.hpp). Only named
/// here, never dereferenced: \c validate_cell_table takes it by reference and
/// reserves it for the block-tree walk (see its \c root parameter), so the
/// types and the validator stay independent of the schedule machinery. The
/// builder (cell_table_builder.hpp) is where the schedule is actually read.
struct ScopeBlock;

/// Explicit value cells (see the explicit-value-cells design document). A cell
/// is one FORM of a value resident at one scope; identity is by cell id,
/// carried position and loop instance -- never by canonical index label.
using CellId = std::size_t;

/// Enclosing loop instances outermost-first, each with its latitude (pass).
/// Empty = the root scope.
struct CellScope {
  container::svector<std::pair<LoopKey, int>> path;

  /// true iff this scope is the root or a strict/equal prefix of \p inner
  /// (same loop instances AND same passes along the prefix). Compares full
  /// (LoopKey, latitude) entries -- LAYOUT -- by design: residency and
  /// multiplicity instead compare loop IDENTITY alone via \c same_key,
  /// ignoring latitude.
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
  /// loop instances the consumer is sliced by on which the schedule seam
  /// (\c SlicedModeAssignment::occ_invariant) recorded an explicit fact:
  /// an explicit record that no slicing decision bound this read on that
  /// instance; the form rule accepts it as whole, so this table check is
  /// necessary, not sufficient -- the dry-run range check remains the
  /// ground truth. Empty = no such record.
  container::svector<LoopKey> invariant_on;
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
/// produced -- decided by \c s.production.kind, exactly three cases:
/// (1) \c Assemble: the prefix of \c s.scope.path ending at the DEEPEST
/// instance in \c bound_instances(s) (it dies when that loop's batch ends);
/// the root scope (empty path) when \p s is bound to NONE of its enclosing
/// loops -- the runtime's close-store walk homes a whole block output that
/// far out, all the way to the chain root; (2) \c Build (plain and implicit
/// per-batch alike): always \c s.scope itself, bound or whole -- a step's
/// value is stored in its OWN block's cache and dies with that block, so
/// even a whole Build cell spans only that block's batches, never further
/// out; (3) \c Leaf: the root scope. \c persistent and \c produce_if_absent
/// are SEPARATE properties (survival across evaluations / first-visit
/// production, not where within one evaluation a cell lives) and are NOT
/// consulted here.
[[nodiscard]] inline CellScope residency_scope(TableCell const& s) {
  if (s.production.kind == ProductionKind::Leaf) return CellScope{};
  if (s.production.kind == ProductionKind::Build) return s.scope;
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
  if (!bound) return CellScope{};
  CellScope r;
  r.path.assign(s.scope.path.begin(), s.scope.path.begin() + deepest + 1);
  return r;
}

/// Multiplicity of one read of a source cell by a consumer cell: one factor
/// of \c max(1, n_batches_of(key)) per loop instance on the consumer's scope
/// path that is NOT on the source's RESIDENCY scope path (\c
/// residency_scope(source), not its raw \c scope -- an Assemble's residency
/// can extend past its own scope, and a read must be weighted against where
/// the source actually still lives, not merely where it was produced) -- a
/// consumer nested one extra (non-resident) loop deeper re-reads the source
/// once per batch of each such loop. Takes the source \p source as a whole
/// \c TableCell (never a bare scope) so this is the only place residency is
/// applied to multiplicity; the single place both the builder (life
/// bookkeeping) and the validator (rule 4) compute this, so they never
/// drift apart. A null \p n_batches_of (the validator's default) is treated
/// as returning 1 everywhere.
[[nodiscard]] inline std::size_t read_multiplicity(
    TableCell const& source, CellScope const& consumer_scope,
    std::function<std::size_t(LoopKey const&)> const& n_batches_of) {
  CellScope const source_residency = residency_scope(source);
  std::size_t mult = 1;
  for (auto const& [pk, lat] : consumer_scope.path) {
    (void)lat;
    bool in_source = false;
    for (auto const& [sk, slat] : source_residency.path) {
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
  // RESIDENCY scope (detail::residency_scope: bound to an enclosing loop
  // instance -> the prefix of its scope ending at that loop, else -- a
  // whole cell -- the root scope) encloses the consumer's scope.
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
  // bound to that group must be bound to the SAME instance; an UNDECIDED
  // whole operand (one with no explicit invariant record) on a group that
  // another operand is bound to is a mismatch. A read the seam recorded as
  // invariant on this instance (\c Read::invariant_on, from \c
  // SlicedModeAssignment::occ_invariant: an explicit record that no slicing
  // decision bound this read on that instance) is neither "bound" nor
  // "whole" for this rule -- the form rule accepts it as whole, so this
  // check is necessary, not sufficient (the dry-run range check remains the
  // ground truth), and it never contributes to the mismatch. A read bound
  // to another instance of the SAME group is still always its own
  // violation, invariant or not.
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
        bool invariant_here = false;
        for (LoopKey const& k : r.invariant_on)
          if (detail::same_key(k, L)) invariant_here = true;
        if (bound_here)
          any_bound = true;
        else if (!bound_other && !invariant_here)
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
        reads += detail::read_multiplicity(c, table.cells[r.consumer].scope,
                                           n_batches_of);
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

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_CELL_TABLE_HPP
