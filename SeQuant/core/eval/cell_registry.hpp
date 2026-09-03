#ifndef SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP
#define SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>  // BatchContext
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/result.hpp>

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace sequant::eval {

/// \c CacheManager<TreeNode, force_hash_collisions>::BatchContext does not
/// depend on either template parameter (see cache_manager.hpp), so this is
/// literally the same type under every instantiation -- naming it here lets
/// this header (and its tests) spell it without a CacheManager instance.
using BatchContextEntry = sequant::BatchContextEntry;
using BatchContext = container::svector<BatchContextEntry>;
using DagScopeLevel = sequant::DagScopeLevel;

/// Runtime side of the cell table: the current result of each cell and its
/// remaining life. Storage ownership stays with the scope caches in this
/// stage; the registry holds a second reference so that reads resolve by cell
/// id. Bound cells are cleared at the start of every batch of a loop instance
/// they are bound to (the per-batch scratch reset, expressed on cells).
class CellRegistry {
 public:
  explicit CellRegistry(CellTable const& table) : table_(&table) {
    slots_.resize(table.cells.size());
    for (CellId c = 0; c < table.cells.size(); ++c)
      slots_[c].life = table.cells[c].life;
    for (CellId c = 0; c < table.cells.size(); ++c) {
      auto const& cell = table.cells[c];
      switch (cell.production.kind) {
        case ProductionKind::Build:
          build_at_[cell.value_id].push_back(c);
          break;
        case ProductionKind::Assemble:
          assemble_at_[cell.value_id].push_back(c);
          break;
        case ProductionKind::Leaf:
          leaf_of_[cell.value_id] = c;
          break;
      }
    }
  }

  /// Production: overwrites the cell's current result and restores its life
  /// from the table (a new batch's/iteration's production of a cell whose
  /// prior life was drained starts fresh).
  void set(CellId c, ResultPtr v) {
    auto& s = slot(c);
    s.value = std::move(v);
    s.life = table_->cells[c].life;
  }

  /// Non-decrementing peek: the cell's current result, or null if unset.
  [[nodiscard]] ResultPtr peek(CellId c) const { return slot(c).value; }

  /// Decrementing read: throws if the cell has no current result, or (for a
  /// non-persistent cell) its life is already exhausted.
  [[nodiscard]] ResultPtr read(CellId c) {
    auto& s = slot(c);
    if (!s.value)
      throw std::runtime_error("CellRegistry::read: cell#" + std::to_string(c) +
                               " (value " +
                               std::to_string(table_->cells[c].value_id) +
                               ") has no current result");
    if (!table_->cells[c].persistent) {
      if (s.life == 0)
        throw std::runtime_error("CellRegistry::read: cell#" +
                                 std::to_string(c) + " read past its life " +
                                 std::to_string(table_->cells[c].life));
      --s.life;
    }
    return s.value;
  }

  /// Batch start: drops every cell bound to loop instance \p k (see
  /// detail::bound_instances) -- the per-batch scratch reset, expressed on
  /// cells instead of on a whole scope's storage.
  void clear_bound_to(LoopKey const& k) {
    for (CellId c = 0; c < slots_.size(); ++c)
      for (LoopKey const& b : detail::bound_instances(table_->cells[c]))
        if (detail::same_key(b, k)) {
          slots_[c].value.reset();
          break;
        }
  }

  [[nodiscard]] CellTable const& table() const { return *table_; }

  [[nodiscard]] std::optional<CellId> build_cell_at(
      std::size_t vid, CellScope const& scope) const {
    return find_at(build_at_, vid, scope);
  }
  [[nodiscard]] std::optional<CellId> assemble_cell_at(
      std::size_t vid, CellScope const& scope) const {
    return find_at(assemble_at_, vid, scope);
  }
  [[nodiscard]] std::optional<CellId> leaf_cell(std::size_t vid) const {
    auto it = leaf_of_.find(vid);
    if (it == leaf_of_.end()) return std::nullopt;
    return it->second;
  }

  /// The value's form VISIBLE at \p scope, by RESIDENCY (not exact scope
  /// equality) -- the table's own visibility contract (cell_table.hpp's
  /// validator rule 1 / \c detail::residency_scope's own doc comment): the
  /// DEEPEST-scoped (largest \c scope.path) Build or Assemble cell of \p vid
  /// whose \c detail::residency_scope ENCLOSES \p scope. A tie at equal
  /// depth prefers a Build cell (a Build and an Assemble of one value never
  /// actually coexist at one scope on a well-formed table, so this never
  /// actually arbitrates; kept for determinism). Leaf cells are not
  /// candidates here (see \c leaf_cell): a leaf has no scope in this sense.
  [[nodiscard]] std::optional<CellId> cell_of(std::size_t vid,
                                              CellScope const& scope) const {
    std::optional<CellId> best;
    auto const consider = [&](container::svector<CellId> const& ids) {
      for (CellId c : ids) {
        TableCell const& cell = table_->cells[c];
        if (!detail::residency_scope(cell).encloses(scope)) continue;
        if (!best) {
          best = c;
          continue;
        }
        TableCell const& b = table_->cells[*best];
        bool const deeper = cell.scope.path.size() > b.scope.path.size();
        bool const tie_prefers_build =
            cell.scope.path.size() == b.scope.path.size() &&
            cell.production.kind == ProductionKind::Build &&
            b.production.kind != ProductionKind::Build;
        if (deeper || tie_prefers_build) best = c;
      }
    };
    if (auto it = build_at_.find(vid); it != build_at_.end())
      consider(it->second);
    if (auto it = assemble_at_.find(vid); it != assemble_at_.end())
      consider(it->second);
    return best;
  }

 private:
  struct Slot {
    ResultPtr value;
    std::size_t life = 0;
  };
  Slot& slot(CellId c) {
    if (c >= slots_.size())
      throw std::out_of_range("CellRegistry: cell id out of range");
    return slots_[c];
  }
  Slot const& slot(CellId c) const {
    return const_cast<CellRegistry*>(this)->slot(c);
  }
  std::optional<CellId> find_at(
      std::unordered_map<std::size_t, container::svector<CellId>> const& m,
      std::size_t vid, CellScope const& scope) const {
    auto it = m.find(vid);
    if (it == m.end()) return std::nullopt;
    for (CellId c : it->second)
      if (table_->cells[c].scope == scope) return c;
    return std::nullopt;
  }

  CellTable const* table_;
  container::vector<Slot> slots_;
  std::unordered_map<std::size_t, container::svector<CellId>> build_at_,
      assemble_at_;
  std::unordered_map<std::size_t, CellId> leaf_of_;
};

/// Resolves one consumer cell's operand fetches to table reads. Installed on
/// a scratch cache for one consumer cell at a time (\c begin_consumer resets
/// the per-operand read cursors); \c fetch is consulted by \c evaluate_impl
/// ahead of every other probe.
class CellReadResolver {
 public:
  CellReadResolver(
      CellRegistry& reg,
      std::function<std::optional<std::size_t>(std::size_t)> vid_of_hash)
      : reg_(&reg), vid_of_hash_(std::move(vid_of_hash)) {}

  /// Resets the per-operand read cursors to every \c Read of \p consumer, in
  /// table order.
  void begin_consumer(CellId consumer) {
    consumer_ = consumer;
    cursor_.clear();
    auto const& t = reg_->table();
    for (std::size_t r = 0; r < t.reads.size(); ++r)
      if (t.reads[r].consumer == consumer)
        cursor_[t.reads[r].operand_value_id].push_back(r);
  }
  [[nodiscard]] CellId consumer() const { return consumer_; }

  /// \return nullopt when \p operand_node_hash is not a value of the table
  /// (a transient of this production tree, evaluated in place by the
  /// caller), OR when the matched Read's source is a LEAF cell with no
  /// current result yet (first touch: the caller must evaluate the leaf and
  /// call \c record_leaf; the cursor entry is left UNCONSUMED so the SAME
  /// Read is served -- and consumed -- by a later fetch once the leaf is
  /// recorded). Any OTHER matched Read whose source has no current result
  /// THROWS naming the consumer cell, the source cell and the value (spec
  /// section 4: "missing entry or non-resident source: throw with both
  /// ids") -- a well-formed table guarantees a Build cell's own source is
  /// always resident when read (registry lookups go by RESIDENCY, see \c
  /// CellRegistry::cell_of, and persistent cross-call values are seeded into
  /// the registry at entry, see \c run_ordered_schedule_pre_results), so
  /// this is a genuine table/tree or recording gap, never deferred.
  /// Otherwise the sliced source: each \c (pos, key) of the matched Read's
  /// \c slice is applied as \c slice_mode(pos, range.first, range.second)
  /// with \c range taken from the \p ctx entry whose \c level.key() equals
  /// \c key. Throws when the consumer has no remaining Read of that value at
  /// all (table/tree disagreement), or when a declared slice names a loop
  /// instance absent from \p ctx.
  [[nodiscard]] std::optional<ResultPtr> fetch(std::size_t operand_node_hash,
                                               BatchContext const& ctx) {
    auto const vid = vid_of_hash_(operand_node_hash);
    if (!vid) return std::nullopt;  // not a value of the table: a transient
    auto it = cursor_.find(*vid);
    if (it == cursor_.end() || it->second.empty())
      throw std::runtime_error(
          "CellReadResolver: consumer cell#" + std::to_string(consumer_) +
          " has no remaining read of value " + std::to_string(*vid));
    Read const& r = reg_->table().reads[it->second.front()];
    TableCell const& src_cell = reg_->table().cells[r.source];
    if (!reg_->peek(r.source)) {
      if (src_cell.production.kind == ProductionKind::Leaf)
        return std::nullopt;  // leaf first touch: cursor left UNCONSUMED
      throw std::runtime_error(
          "CellReadResolver: consumer cell#" + std::to_string(consumer_) +
          " source cell#" + std::to_string(r.source) + " (value " +
          std::to_string(*vid) + ") has no current result");
    }
    it->second.erase(it->second.begin());
    ResultPtr v = reg_->read(r.source);
    for (auto const& [pos, key] : r.slice) {
      std::optional<std::pair<std::size_t, std::size_t>> range;
      for (auto const& e : ctx)
        if (detail::same_key(e.level.key(), key)) range = e.range;
      if (!range)
        throw std::runtime_error(
            "CellReadResolver: read of cell#" + std::to_string(r.source) +
            " by cell#" + std::to_string(consumer_) + " binds position " +
            std::to_string(pos) +
            " to a loop instance not in the batch context");
      v = v->slice_mode(pos, range->first, range->second);
    }
    ++served_;
    return v;
  }

  /// Records a leaf's freshly evaluated result in the registry (called after
  /// the leaf evaluator runs on \p hash's first touch): a no-op if \p hash is
  /// not a value of the table or names no Leaf cell.
  void record_leaf(std::size_t hash, ResultPtr r) {
    if (auto vid = vid_of_hash_(hash))
      if (auto leaf = reg_->leaf_cell(*vid)) reg_->set(*leaf, std::move(r));
  }

  /// Diagnostic: the number of operand fetches this resolver has actually
  /// served (table-driven, sliced-per-Read) rather than deferring to the
  /// caller (a transient, a leaf's first touch, or a value this call's
  /// runtime cache-halt gate skipped -- see \c fetch's own doc comment).
  [[nodiscard]] std::size_t served() const noexcept { return served_; }

 private:
  CellRegistry* reg_;
  std::function<std::optional<std::size_t>(std::size_t)> vid_of_hash_;
  CellId consumer_ = 0;
  std::unordered_map<std::size_t, container::svector<std::size_t>> cursor_;
  std::size_t served_ = 0;
};

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP
