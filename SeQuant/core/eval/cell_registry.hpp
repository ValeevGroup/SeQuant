#ifndef SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP
#define SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>  // BatchContext
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/utility/macros.hpp>

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

/// Wiring \c CellRegistry needs from its owner (the ordered executor's
/// Stage-3 storage move) but must not construct itself: the cross-call \c
/// PersistentValueStore (Task 1's addition to \c CacheManager, one instance
/// per top-level evaluation call, outliving any one \c CellRegistry), the map
/// from a table value id to the canonical hash the store and the legacy scope
/// caches both key persistence on, and an optional metering hook. Every field
/// defaults to "off": a \c CellRegistry built with a default-constructed \c
/// CellRegistryHooks behaves exactly as it did before persistence/bytes
/// tracking existed (persistent cells are simply never seeded or published,
/// and \c live_bytes stays a bookkeeping-only counter nobody observes).
struct CellRegistryHooks {
  /// Cross-call persistence store; \c seed_persistent reads it, \c set
  /// publishes to it. Null disables both.
  PersistentValueStore* persistent = nullptr;
  /// A table value id's canonical hash, the key \c persistent is addressed
  /// by. Required (alongside \c persistent) for seeding/publishing; a null
  /// function disables both exactly as a null \c persistent does.
  std::function<std::size_t(std::size_t value_id)> hash_of;
  /// Metering: invoked with the new \c live_bytes total after every change.
  /// May be empty -- \c CellRegistry does not require a metering consumer.
  std::function<void(std::size_t live_bytes)> on_bytes_changed;
};

/// Runtime side of the cell table: the current result of each cell and its
/// remaining life. This is the OWNER of those results (Stage 3 of the
/// explicit-value-cells design): it tracks the live byte total across held
/// slots, enforces fill-once (a non-persistent cell produced twice without an
/// intervening clear is a duplicate producer -- a bug, not a legitimate
/// replay), and seeds/publishes persistent cells through a \c
/// PersistentValueStore so they survive across top-level evaluation calls
/// without depending on the legacy scope caches' own persistence. Bound cells
/// are cleared at the start of every batch of a loop instance they are bound
/// to (the per-batch scratch reset, expressed on cells).
///
/// ORIENTATION CONVENTION (the legacy scope cache's, kept verbatim): a cell
/// holds the value in its node's CANONICAL orientation, and every reader
/// applies that node's own \c canon_phase once. A producer therefore stores
/// \c apply_phase(node, result) -- \c evaluate_impl hands back the ORIENTED
/// result, and the phase is an involution -- exactly as \c
/// CacheManager::store did (\c cache.store(node, apply_phase(node, rb)),
/// readers \c apply_phase again, roots \c mult_by_phase). Storing the
/// oriented value here instead would double-apply the phase at every read
/// and flip the sign of every value whose canonicalization needed an odd
/// permutation.
class CellRegistry {
 public:
  explicit CellRegistry(CellTable const& table, CellRegistryHooks hooks = {})
      : table_(&table), hooks_(std::move(hooks)) {
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

  /// Seeds every persistent cell whose canonical hash is currently held by
  /// \c hooks_.persistent with the store's value, WITHOUT spending life (the
  /// value just arrived from a prior top-level call; nothing has read it
  /// yet this call). A no-op for a cell the store does not (yet) hold, and
  /// entirely a no-op when \c hooks_.persistent or \c hooks_.hash_of is
  /// unset. Intended to run once, right after construction, before the first
  /// read of this evaluation call (see \c run_ordered_schedule_pre_results).
  void seed_persistent() {
    if (!hooks_.persistent || !hooks_.hash_of) return;
    for (CellId c = 0; c < slots_.size(); ++c) {
      TableCell const& cell = table_->cells[c];
      if (!cell.persistent) continue;
      ResultPtr v = hooks_.persistent->get(hooks_.hash_of(cell.value_id));
      if (!v) continue;
      Slot& s = slots_[c];
      if (s.value)
        account(-static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
      s.value = std::move(v);
      s.life = cell.life;
      s.filled_since_clear = true;
      account(static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
    }
  }

  /// Production: overwrites the cell's current result and restores its life
  /// from the table (a new batch's/iteration's production of a cell whose
  /// prior life was drained starts fresh). FILL-ONCE: a NON-persistent cell
  /// already holding a value it was not read past nor cleared since (\c
  /// filled_since_clear) is a duplicate producer -- see \c
  /// eval::strict_fill_once (throws a named \c std::runtime_error under \c
  /// SEQUANT_UT_STRICT_FILL_ONCE; a plain \c SEQUANT_ASSERT otherwise,
  /// compiled out in Release). PERSISTENT cells are excluded from this check,
  /// exactly as \c CacheManager::entry::store excludes its own persistent
  /// entries: they legitimately re-store across batch replays and repeated
  /// top-level evaluation calls (e.g. successive CC iterations), often with
  /// no table-declared clear between productions. A persistent cell is
  /// published to \c hooks_.persistent (when set) on every production,
  /// mirroring what it now holds.
  void set(CellId c, ResultPtr v) {
    auto& s = slot(c);
    TableCell const& cell = table_->cells[c];
    if (!cell.persistent && s.filled_since_clear) {
      if (strict_fill_once())
        throw std::runtime_error(
            "CellRegistry::set: cell#" + std::to_string(c) + " (value " +
            std::to_string(cell.value_id) +
            ") filled twice without an intervening clear (duplicate "
            "producer) -- cache-fill-once violated");
      SEQUANT_ASSERT(!s.filled_since_clear &&
                     "CellRegistry::set: fill-once violated");
    }
    if (s.value)
      account(-static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
    s.value = std::move(v);
    s.life = cell.life;
    s.filled_since_clear = true;
    if (s.value) account(static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
    if (cell.persistent && hooks_.persistent && hooks_.hash_of)
      hooks_.persistent->put(hooks_.hash_of(cell.value_id), s.value);
  }

  /// Non-decrementing peek: the cell's current result, or null if unset.
  [[nodiscard]] ResultPtr peek(CellId c) const { return slot(c).value; }

  /// Decrementing read: throws if the cell has no current result, or (for a
  /// non-persistent cell) its life is already exhausted. The read that spends
  /// a non-persistent cell's LAST life also DROPS the registry's own
  /// reference (the cell has no reader left this evaluation, and holding on
  /// would both pin the memory and make the buffer look shared to the reader
  /// that just took it -- see \c CacheManager::entry::release, which the
  /// scope cache's own \c access() has always done at the same point). A
  /// later production of the cell restores both value and life via \c set.
  [[nodiscard]] ResultPtr read(CellId c) { return read(c, nullptr); }

  /// \overload As \c read(CellId), and reports through \p exhausted (when
  /// non-null) whether this read spent the cell's last life.
  ///
  /// Stage 3 re-derives in-place eligibility from this table's own \c life /
  /// \c persistent directly (see \c CellReadResolver::operand_drained,
  /// which queries \c drained() rather than this flag). The \p exhausted
  /// output tells a reader whether it took OWNERSHIP of the buffer: the
  /// ordered executor's Assemble step seeds its running sum with the first
  /// batch's partial when this read exhausted the source, and with a \c
  /// Result::clone of it otherwise (a source with life left, or a persistent
  /// one, is still going to be read again from that very buffer, so
  /// accumulating into it would corrupt every later read).
  [[nodiscard]] ResultPtr read(CellId c, bool* exhausted) {
    auto& s = slot(c);
    if (exhausted) *exhausted = false;
    if (!s.value)
      throw std::runtime_error("CellRegistry::read: cell#" + std::to_string(c) +
                               " (value " +
                               std::to_string(table_->cells[c].value_id) +
                               ") has no current result");
    if (table_->cells[c].persistent) return s.value;
    if (s.life == 0)
      throw std::runtime_error("CellRegistry::read: cell#" + std::to_string(c) +
                               " read past its life " +
                               std::to_string(table_->cells[c].life));
    --s.life;
    if (s.life != 0) return s.value;
    if (exhausted) *exhausted = true;
    // The value is fully consumed and about to be released, so a later set()
    // of this cell is a fresh re-production (a new batch's/iteration's), not
    // a duplicate -- clear the fill-once mark along with the byte accounting,
    // mirroring CacheManager::entry::access's stored_this_eval_ reset on the
    // same draining access.
    s.filled_since_clear = false;
    account(-static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
    return std::move(s.value);  // the slot is left null by the move
  }

  /// Batch start: drops every cell bound to loop instance \p k -- the
  /// per-batch scratch reset, expressed on cells instead of on a whole
  /// scope's storage. A cell is bound to \p k either explicitly (\c
  /// detail::bound_instances: a carried \c sliced position or a \c
  /// partial_over reduction on \p k) or IMPLICITLY, for a non-persistent
  /// cell whose own DECLARED home scope's deepest loop instance is \p k: the
  /// executor re-runs whatever this cell's tree position computes fresh
  /// every batch of that position's innermost enclosing loop (a step's own
  /// block for a Build; a block-close aggregation for an Assemble),
  /// independent of whether the value it holds happens to carry \p k as a
  /// mode or a reduced axis -- a "whole" cell unsliced and unsummed on its
  /// own home loop is still re-run every batch of it, only its declared
  /// life/persistence say how long the RESULT is then read for. A PERSISTENT
  /// cell is never cleared here (by definition it is bound to no loop
  /// instance -- see \c TableCell::persistent -- so neither check ever
  /// matches it; the explicit \c continue is a defensive redundant guard).
  /// A \c produce_if_absent cell is the ONE exception to the implicit rule:
  /// it is precisely a cell whose home scope sits inside a loop it is NOT
  /// bound to, and the table says it is produced on first visit and REUSED on
  /// every later batch of that loop (the executor's Build step skips its
  /// production while it is resident), so the implicit "its scope's innermost
  /// loop is \p k" clear would defeat the flag by dropping it at every batch
  /// boundary. Such a cell is cleared only on an instance it is genuinely
  /// bound to. Also resets the fill-once mark of every cell it clears: the
  /// boundary this crosses is exactly the one \c set's fill-once check must
  /// not treat as a duplicate producer across.
  void clear_bound_to(LoopKey const& k) {
    for (CellId c = 0; c < slots_.size(); ++c) {
      TableCell const& cell = table_->cells[c];
      if (cell.persistent) continue;
      bool bound = false;
      for (LoopKey const& b : detail::bound_instances(cell))
        if (detail::same_key(b, k)) {
          bound = true;
          break;
        }
      if (!bound && !cell.produce_if_absent && !cell.scope.path.empty() &&
          detail::same_key(cell.scope.path.back().first, k))
        bound = true;
      if (!bound) continue;
      Slot& s = slots_[c];
      if (s.value)
        account(-static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
      s.value.reset();
      s.filled_since_clear = false;
    }
  }

  /// Spends \p count of cell \p c's declared reads WITHOUT taking the value:
  /// the accounting half of a read the executor decided not to perform. A
  /// consumer the runtime cache-halt skipped (nothing left this evaluation
  /// reads its result, or a \c produce_if_absent cell that is resident and so
  /// is not re-produced) never calls \c read, so without this its source's
  /// life would never reach zero -- the source would stay resident to the end
  /// of the evaluation (pinning the memory) and keep looking shared to every
  /// later reader, which is exactly what disables in-place accumulation.
  ///
  /// Releases the value when the life reaches zero, exactly as the draining
  /// \c read does (bytes accounted, fill-once mark cleared so a later
  /// production of the cell is a fresh one). A no-op for a persistent cell,
  /// which never drains, and for \p count == 0. THROWS when \p count exceeds
  /// the remaining life: forgoing more reads than the table declared is an
  /// accounting bug in the caller, and silently saturating would release a
  /// value another consumer is still owed (see \c remaining_life, which the
  /// ordered executor consults for the one case where a skipped consumer's
  /// visits legitimately outnumber one production's budget).
  void forgo(CellId c, std::size_t count) {
    if (count == 0) return;
    auto& s = slot(c);
    if (table_->cells[c].persistent) return;
    if (count > s.life)
      throw std::runtime_error(
          "CellRegistry::forgo: cell#" + std::to_string(c) + " (value " +
          std::to_string(table_->cells[c].value_id) + ") forgoes " +
          std::to_string(count) + " reads but only " + std::to_string(s.life) +
          " of its declared life " + std::to_string(table_->cells[c].life) +
          " remain");
    s.life -= count;
    if (s.life != 0) return;
    s.filled_since_clear = false;
    if (s.value) {
      account(-static_cast<std::ptrdiff_t>(s.value->size_in_bytes()));
      s.value.reset();
    }
  }

  /// \return cell \p c's remaining declared life (0 for a spent cell; a
  /// persistent cell's life is never spent, so this is informational there).
  [[nodiscard]] std::size_t remaining_life(CellId c) const {
    return slot(c).life;
  }

  /// \return the sum of \c Result::size_in_bytes() over every slot currently
  /// holding a value (persistent and non-persistent alike).
  [[nodiscard]] std::size_t live_bytes() const { return live_bytes_; }

  /// Invokes \p f(CellId, ResultPtr const&) for every slot currently holding
  /// a value -- the registry-side source \c on_peak_liveset walks once
  /// storage moves fully onto the table.
  template <typename F>
  void for_each_live(F&& f) const {
    for (CellId c = 0; c < slots_.size(); ++c)
      if (slots_[c].value) f(c, slots_[c].value);
  }

  /// \return whether cell \p c is a spent non-persistent cell: no value held
  /// and no life left to spend. A persistent cell is never drained, so it is
  /// never reported as drained even when it happens to hold no value yet
  /// (e.g. before its first production this evaluation).
  [[nodiscard]] bool drained(CellId c) const {
    Slot const& s = slot(c);
    return !table_->cells[c].persistent && s.life == 0 && !s.value;
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
  /// equality) -- the table's own visibility contract, decided by the one
  /// shared rule \c detail::deepest_visible_form states (deepest resident
  /// form wins; at equal depth the earlier candidate does), which the table
  /// builder's read-source selection uses too. Build cells are offered
  /// FIRST, which is how a tie prefers a Build. Leaf cells are not
  /// candidates here (see \c leaf_cell): a leaf has no scope in this sense.
  [[nodiscard]] std::optional<CellId> cell_of(std::size_t vid,
                                              CellScope const& scope) const {
    container::svector<CellId> candidates;
    if (auto it = build_at_.find(vid); it != build_at_.end())
      candidates.assign(it->second.begin(), it->second.end());
    if (auto it = assemble_at_.find(vid); it != assemble_at_.end())
      candidates.insert(candidates.end(), it->second.begin(), it->second.end());
    return detail::deepest_visible_form(*table_, candidates, scope);
  }

 private:
  struct Slot {
    ResultPtr value;
    std::size_t life = 0;
    /// Fill-once tripwire: true once \c value has been \c set() since the
    /// last time this slot was emptied (by \c clear_bound_to, or by a \c
    /// read that spent the cell's last life). See \c set.
    bool filled_since_clear = false;
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

  /// Applies \p delta (signed) to \c live_bytes_ and, when \c
  /// hooks_.on_bytes_changed is set, reports the new total.
  void account(std::ptrdiff_t delta) {
    live_bytes_ = static_cast<std::size_t>(
        static_cast<std::ptrdiff_t>(live_bytes_) + delta);
    if (hooks_.on_bytes_changed) hooks_.on_bytes_changed(live_bytes_);
  }

  CellTable const* table_;
  CellRegistryHooks hooks_;
  container::vector<Slot> slots_;
  std::unordered_map<std::size_t, container::svector<CellId>> build_at_,
      assemble_at_;
  std::unordered_map<std::size_t, CellId> leaf_of_;
  std::size_t live_bytes_ = 0;
};

/// The value and exhaustion flag of one \c table_read.
struct TableRead {
  ResultPtr value;
  /// Whether this read spent \c source's LAST declared life (always \c false
  /// for a persistent cell, which never exhausts).
  bool exhausted = false;
};

/// The OWNERSHIP half of one table-driven read of \p source: spend one life
/// of \p source in \p reg and report whether that read spent the cell's
/// LAST life. \c CellReadResolver::fetch calls it for a consumer's operand
/// reads; an \c Assemble step reads its \c production.source through \c
/// CellRegistry::read directly, for the same accounting. A read served from
/// somewhere other than the registry would still owe the table that life:
/// skipping it leaves the source holding a fully consumed buffer, which pins
/// the memory and makes every later reader see the value as shared.
///
/// \c CellReadResolver::fetch no longer consults \c TableRead::exhausted
/// itself (Stage 3: \c CellReadResolver::operand_drained re-derives in-place
/// eligibility straight from \c CellRegistry::drained instead); an Assemble
/// step does, to decide whether the partial it just read is its own to
/// accumulate into (see \c CellRegistry::read's overload).
[[nodiscard]] inline TableRead table_read(CellRegistry& reg, CellId source) {
  TableRead r;
  r.value = reg.read(source, &r.exhausted);
  return r;
}

/// Resolves one consumer cell's operand fetches to table reads. Installed on
/// a scratch cache for one consumer cell at a time (\c begin_consumer resets
/// the per-operand read cursors); \c fetch is consulted by \c evaluate_impl
/// ahead of every other probe.
///
/// POSITIONAL MATCHING -- the invariant that makes a cursor per operand VALUE
/// sufficient, and the reason a \c Read carries no leg number:
/// 1. the table emits a consumer's \c Read entries in PRODUCTION-TREE LEG
///    order (cell_table_builder.hpp walks \c CellTableInputs::operands_of,
///    which is per-leg WITH repetition, left leg then right leg);
/// 2. the executor fetches a consumer's legs in that SAME order
///    (\c evaluate_impl requests the left operand, then the right);
/// 3. therefore the i-th surviving \c Read of one value in \c cursor_ is the
///    i-th leg of that value, and popping the front matches legs to reads
///    with no leg index anywhere. \c begin_consumer asserts (1) by checking
///    each per-value cursor is ordered by table position;
/// 4. every table value in a consumer's production tree is a DIRECT leg. A
///    non-leg intermediate node of that tree is a transient (not a value of
///    the table) and \c fetch reports it as such; a table value reached from
///    an intermediate node instead of a leg has no \c Read of its own and
///    \c fetch throws rather than borrowing another leg's read.
class CellReadResolver {
 public:
  CellReadResolver(
      CellRegistry& reg,
      std::function<std::optional<std::size_t>(std::size_t)> vid_of_hash)
      : reg_(&reg), vid_of_hash_(std::move(vid_of_hash)) {
    // Index the table's reads by consumer ONCE: begin_consumer runs before
    // every build step and every per-batch output evaluation, and scanning
    // the whole read list there made each step cost O(all reads).
    auto const& t = reg_->table();
    reads_of_.resize(t.cells.size());
    for (std::size_t r = 0; r < t.reads.size(); ++r)
      reads_of_[t.reads[r].consumer].push_back(r);
  }

  /// Resets the per-operand read cursors to every \c Read of \p consumer, in
  /// table order -- i.e. in the consumer's production-tree LEG order, which
  /// is what makes the front of a per-value cursor the next leg to fetch
  /// (see the positional-matching invariant on this class).
  void begin_consumer(CellId consumer) {
    consumer_ = consumer;
    cursor_.clear();
    auto const& t = reg_->table();
    if (consumer >= reads_of_.size()) return;  // a cell nothing reads for
    for (std::size_t r : reads_of_[consumer]) {
      auto& c = cursor_[t.reads[r].operand_value_id];
      // Leg order == table order (invariant 1): the index is built by
      // ascending read position, so each cursor must come out ordered.
      SEQUANT_ASSERT(
          (c.empty() || c.back() < r) &&
          "CellReadResolver: a consumer's reads of one value are out of "
          "production-tree leg order");
      c.push_back(r);
    }
  }
  [[nodiscard]] CellId consumer() const { return consumer_; }

  /// The registry this resolver reads from -- the storage the table owns.
  /// Read-only: a caller that wants to OBSERVE what a cell currently holds
  /// (a test probe, a diagnostic) goes through here; production is the
  /// executor's business.
  [[nodiscard]] CellRegistry const& registry() const noexcept { return *reg_; }

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
    // Stage 3: record WHICH cell this vid's read was most recently served
    // from -- operand_drained() below re-derives exhaustion from the
    // registry's own current state (life/value) at that cell instead of a
    // cached bool, so it can never drift from what the registry actually
    // holds.
    last_served_source_[*vid] = r.source;
    ResultPtr v = table_read(*reg_, r.source).value;
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

  /// Stage 3: whether \p operand_node_hash's operand is currently safe for
  /// in-place accumulation -- the registry-derived replacement for the
  /// legacy \c CacheManager::chain_holds_shared(f.left) check \c eval.hpp's
  /// in-place gate used before storage moved onto the table.
  ///
  /// True in exactly two cases:
  ///  - \p operand_node_hash is not a value of this table at all (a private
  ///    transient of the current production tree, e.g. an intermediate
  ///    running total of an accumulate-in-place \c Sum chain that was never
  ///    promoted to its own cell) -- no table cell could possibly be sharing
  ///    it, exactly as \c chain_holds_shared() reports "not held" for the
  ///    analogous untracked case; or
  ///  - it IS a table value, has been \c fetch()'d at least once, and its
  ///    MOST RECENT fetch's source cell (\c last_served_source_) is
  ///    currently \c CellRegistry::drained -- i.e. that read spent the
  ///    source's last declared life (never true for a persistent cell, which
  ///    never drains) and the registry has already let go of its own
  ///    reference (\c CellRegistry::read moves the value out on that read),
  ///    so nothing this evaluation will read it again.
  ///
  /// A table value never yet \c fetch()'d through this resolver (no entry in
  /// \c last_served_source_) is reported NOT drained -- the safe default. The
  /// registry may well be holding that value for other readers (it is a
  /// table cell, so some cell owns it), and this resolver has no evidence
  /// either way; answering "drained" would license an in-place mutation of a
  /// buffer this evaluation is going to read again. Out-of-place costs one
  /// allocation; the wrong answer corrupts a shared value.
  [[nodiscard]] bool operand_drained(std::size_t operand_node_hash) const {
    auto const vid = vid_of_hash_(operand_node_hash);
    if (!vid) return true;
    auto const it = last_served_source_.find(*vid);
    if (it == last_served_source_.end()) return false;
    return reg_->drained(it->second);
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
  /// consumer cell id -> its \c Read positions, ascending (built once in the
  /// constructor; see \c begin_consumer).
  container::vector<container::svector<std::size_t>> reads_of_;
  std::unordered_map<std::size_t, container::svector<std::size_t>> cursor_;
  std::size_t served_ = 0;
  /// value id -> the CellId its most recent \c fetch was served from (see
  /// \c operand_drained). Never reset per-consumer (begin_consumer leaves it
  /// alone): "most recent" is global across the whole resolver's lifetime,
  /// not scoped to one consumer's reads.
  std::unordered_map<std::size_t, CellId> last_served_source_;
};

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_CELL_REGISTRY_HPP
