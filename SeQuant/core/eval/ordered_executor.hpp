#ifndef SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
#define SEQUANT_EVAL_ORDERED_EXECUTOR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backend_array_ops.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/cell_table_builder.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/forest_combine.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/value_node_map.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace sequant::eval {

namespace detail {

///
/// \brief SP4 Task 4: per-LOOP-INSTANCE (not per-TYPE) batch count over the
/// whole \c ScopeBlock tree -- the \c CellTableInputs::n_batches_of the cell
/// table's life computation needs (\c detail::read_multiplicity,
/// cell_table.hpp): a source cell read by a consumer nested inside a loop it
/// is not resident on is re-read once per REAL batch of that loop (not once
/// per loop, the placeholder \c n_batches_of that always returns 1 --
/// self-consistent for Task 3's static-only validation, since the builder
/// and the validator used the SAME stub, but an undercount once Task 4
/// enforces \c life at runtime: a value read from OUTSIDE two nested batch
/// loops of \c n \& \c m real batches is read \c n*m times, not once). Keyed
/// by \c LoopKey::color() (the loop's stable depth+loop_slot identity, not
/// its canonical axis label, which a TYPE-keyed count would collapse across
/// distinct same-space instances) so it is exact per realized loop.
///
[[nodiscard]] inline std::function<std::size_t(LoopKey const&)>
ordered_n_batches_by_loop(
    OrderedSchedule const& ordered,
    std::function<std::size_t(Index const&)> const& target,
    BackendArrayOps const* aops) {
  auto by_loop =
      std::make_shared<std::unordered_map<std::size_t, std::size_t>>();
  auto const add = [&](auto&& self, ScopeBlock const& b) -> void {
    std::size_t const color = b.level.key().color();
    if (!by_loop->count(color)) {
      SEQUANT_ASSERT(aops &&
                     "ordered_n_batches_by_loop: batched schedule requires "
                     "backend array-ops (CacheManager::set_array_ops)");
      by_loop->emplace(color,
                       aops->axis_batches(b.axis, target(b.axis)).size());
    }
    for (Step const& s : b.steps)
      if (auto const* child = std::get_if<ScopeBlock>(&s.value))
        self(self, *child);
  };
  for (Step const& s : ordered.root.steps)
    if (auto const* child = std::get_if<ScopeBlock>(&s.value)) add(add, *child);
  return [by_loop](LoopKey const& k) -> std::size_t {
    auto const it = by_loop->find(k.color());
    return it == by_loop->end() ? std::size_t{1} : it->second;
  };
}

/// SP4 Task 4: the \c CellScope of the current point in the schedule walk --
/// enclosing loop instances outermost-first, one \c {level.key(),
/// level.latitude_ordinal} entry per already-opened \p ectx level, plus (the
/// two-argument overload) \p block's own \c {level.key(), latitude_ordinal}
/// entry appended last -- the ONE construction every table lookup this file
/// makes (\c build_cell_at / \c assemble_cell_at) goes through, so a scope a
/// lookup uses here and the scope \c cell_table_builder.hpp's \c emit_cells
/// recorded for the same point (an identical per-level construction) can
/// never drift apart. \c BatchContextEntry has no separate flat latitude
/// field (only \c level), so an enclosing entry reads
/// \c level.latitude_ordinal; \p block's own entry reads its flat
/// \c latitude_ordinal, mirroring \c emit_cells' \c child.latitude_ordinal
/// exactly (the two are the same value by construction, but this keeps the
/// two constructions textually identical rather than merely numerically
/// equal).
///
/// \note Qualified as \c eval::BatchContext (cell_registry.hpp), not the bare
/// \c BatchContext this namespace's OWN \c detail scope already binds to the
/// unrelated legacy alias in peak_profile.hpp (a
/// \c svector<pair<Index,pair<size_t,size_t>>>) -- unqualified lookup from
/// inside \c sequant::eval::detail finds that closer name first, and it is
/// NOT the same type as \c CacheManager<N, FHC>::BatchContext (the type \p
/// ectx actually has at every call site here).
[[nodiscard]] inline CellScope current_scope(eval::BatchContext const& ectx) {
  CellScope s;
  for (auto const& e : ectx)
    s.path.push_back({e.level.key(), e.level.latitude_ordinal});
  return s;
}
[[nodiscard]] inline CellScope current_scope(eval::BatchContext const& ectx,
                                             ScopeBlock const& block) {
  CellScope s = current_scope(ectx);
  s.path.push_back({block.level.key(), block.latitude_ordinal});
  return s;
}

/// Stage 3 cache-halt, computed ONCE per evaluation call over the table
/// (replacing the forest BFS over the legacy cache's alive entries): the set
/// of cells whose production this call may skip.
///
/// A cell is skipped when
///  1. it is PERSISTENT and the registry already HOLDS it -- it survived from
///     a previous evaluation of this schedule through the persistent value
///     store and was seeded back in (\c CellRegistry::seed_persistent), so
///     re-producing it would be pure waste (and, for an accumulating
///     Assemble, would corrupt the held value by summing it into itself); or
///  2. every CONSUMER of it is itself skipped -- nothing left this call will
///     read it. The consumers of a cell are the consumer cells of every
///     \c Read whose \c source is it, plus every \c Assemble whose \c
///     production.source is it. The edge is by SOURCE CELL, not by value:
///     a read names the exact form it consumes and the resolver serves that
///     form and no other, so an in-block partial whose only reader is the
///     Assemble that closes it is dead as soon as that Assemble is -- even
///     though the assembled form of the SAME value is still read elsewhere.
///
/// A cell with NO consumer at all is never skipped BY RULE 2 -- the closure
/// below only ever adds a cell all of whose consumers are skipped, and a cell
/// with no consumer has none to skip. Those cells are the schedule's own
/// results (validator rule 4 admits a zero-read cell only at the root scope),
/// and RULE 1 does reach them: a non-volatile forest root -- a constant term,
/// held by the persistent store from the previous evaluation -- is skipped
/// exactly as any other held persistent cell is, which is the whole point of
/// cache-halt. What the caller then receives for such a root is a private
/// COPY of the held value, never the stored buffer itself (see the root
/// results in \c run_ordered_schedule_pre_results).
///
/// Rule 2 is a fixpoint over the table's dependency edges: skipping a
/// resident persistent composite makes its own prerequisites dead, and so on
/// down. This is the table-side statement of what forest descent does by
/// halting its descent at a cache hit.
[[nodiscard]] inline container::vector<char> ordered_skip_closure(
    CellTable const& table, container::vector<char> seed) {
  std::size_t const n = table.cells.size();
  SEQUANT_ASSERT(seed.size() == n);
  container::vector<container::svector<CellId>> consumers(n);
  for (Read const& r : table.reads)
    if (r.source < n) consumers[r.source].push_back(r.consumer);
  for (CellId c = 0; c < n; ++c)
    if (table.cells[c].production.kind == ProductionKind::Assemble &&
        table.cells[c].production.source < n)
      consumers[table.cells[c].production.source].push_back(c);

  for (bool changed = true; changed;) {
    changed = false;
    for (CellId c = 0; c < n; ++c) {
      if (seed[c] || consumers[c].empty()) continue;
      bool all = true;
      for (CellId x : consumers[c])
        if (!seed[x]) {
          all = false;
          break;
        }
      if (all) {
        seed[c] = 1;
        changed = true;
      }
    }
  }
  return seed;
}

/// \overload The call-wide skip set: seeded by rule 1 (persistent and already
/// held after \c CellRegistry::seed_persistent) and closed under rule 2.
[[nodiscard]] inline container::vector<char> ordered_cache_halt_skip(
    CellTable const& table, CellRegistry const& registry) {
  container::vector<char> seed(table.cells.size(), 0);
  for (CellId c = 0; c < table.cells.size(); ++c)
    if (table.cells[c].persistent && registry.peek(c)) seed[c] = 1;
  return ordered_skip_closure(table, std::move(seed));
}

/// The sources one production of \c consumer reads, one entry per table \c
/// Read of it plus -- for an \c Assemble cell -- the one read it declares of
/// its \c production.source, together with what the count of elided reads
/// needs: the table and the per-loop batch counts. Built once per evaluation
/// next to the skip set; consumed only when a production is SKIPPED (see \c
/// ordered_forgo_reads).
struct ForgoPlan {
  CellTable const* table = nullptr;
  container::vector<container::svector<CellId>> sources;
  std::function<std::size_t(LoopKey const&)> n_batches_of;
};

[[nodiscard]] inline ForgoPlan ordered_forgo_plan(
    CellTable const& table,
    std::function<std::size_t(LoopKey const&)> n_batches_of) {
  ForgoPlan plan;
  plan.table = &table;
  plan.n_batches_of = std::move(n_batches_of);
  plan.sources.resize(table.cells.size());
  for (Read const& r : table.reads)
    if (r.consumer < plan.sources.size() && r.source < table.cells.size())
      plan.sources[r.consumer].push_back(r.source);
  for (CellId c = 0; c < table.cells.size(); ++c)
    if (table.cells[c].production.kind == ProductionKind::Assemble)
      plan.sources[c].push_back(table.cells[c].production.source);
  return plan;
}

/// How many reads of \p source one skip of \p consumer's production elides,
/// when the skip is taken at scope depth \p base_depth (the depth of the
/// scope the skip decision was made at): the product of the batch counts of
/// the loop instances on \p consumer's scope path FROM \p base_depth INWARD
/// that \p source is not resident on.
///
/// This is exactly \c detail::read_multiplicity restricted to a suffix of the
/// consumer's path -- with \p base_depth 0 the two agree literally -- and the
/// suffix is what makes it the count of ELIDED reads rather than of all
/// reads: a skip taken at \p base_depth recurs once per batch of every loop
/// OUTSIDE it, and each such recurrence is a fresh production epoch of a
/// source homed out there, whose life the table charged once per epoch. A
/// per-visit skip passes \p base_depth = the consumer's own scope depth and
/// so elides exactly one read per leg, which is what one visit performs.
[[nodiscard]] inline std::size_t ordered_elided_reads(ForgoPlan const& plan,
                                                      CellId consumer,
                                                      CellId source,
                                                      std::size_t base_depth) {
  CellScope const src_residency =
      detail::residency_scope(plan.table->cells[source]);
  auto const& path = plan.table->cells[consumer].scope.path;
  std::size_t m = 1;
  for (std::size_t i = base_depth; i < path.size(); ++i) {
    bool resident = false;
    for (auto const& [sk, slat] : src_residency.path) {
      (void)slat;
      if (detail::same_key(sk, path[i].first)) resident = true;
    }
    if (!resident)
      m *= std::max<std::size_t>(
          1, plan.n_batches_of ? plan.n_batches_of(path[i].first) : 1);
  }
  return m;
}

/// Spends the reads a SKIPPED production of \p c will not perform, so its
/// sources still reach the end of their declared lives and are released
/// there (a source nobody ever finishes reading stays resident to the end of
/// the evaluation and keeps looking shared, which disables in-place
/// accumulation for it). \p base_depth says how much the skip collapses: the
/// consumer's own scope depth for a skipped visit, the skipped block's parent
/// scope depth for a whole-block skip (see \c ordered_elided_reads).
///
/// The count is CLAMPED to the source's remaining life. A source that is
/// itself skipped is never produced, so its life is never restored, while its
/// skipped consumers keep being visited: their visits legitimately outnumber
/// the one production's budget the table charged. Clamping there is exact
/// (the budget is fully forgone, nothing is left to release); the registry's
/// own throw stays as the tripwire for a genuine over-forgo.
inline void ordered_forgo_reads(CellRegistry& registry, ForgoPlan const& plan,
                                CellId c, std::size_t base_depth) {
  if (c >= plan.sources.size()) return;
  for (CellId src : plan.sources[c]) {
    std::size_t const want = ordered_elided_reads(plan, c, src, base_depth);
    registry.forgo(src, std::min(want, registry.remaining_life(src)));
  }
}

inline std::size_t& ordered_last_block_skips_slot() {
  static std::size_t n = 0;
  return n;
}
/// Diagnostic: how many whole batch loops the most recent \c
/// run_ordered_schedule_pre_results call skipped outright -- a block every one
/// of whose productions was already resident, so it was not entered at all
/// (test-facing; not thread-safe).
[[nodiscard]] inline std::size_t ordered_last_block_skips() {
  return ordered_last_block_skips_slot();
}

/// Whether \p c may be SEEDED into a per-visit skip set (\c
/// run_ordered_contracted_block's own, at block entry).
///
/// Two conditions, and the second is the one that is easy to lose: the cell
/// must be \c produce_if_absent (its production is elided precisely while it
/// stays resident), and it must be bound to NO loop instance. A per-visit set
/// is computed once at a block's entry and consulted across every batch of
/// that block and inside every nested block, while \c
/// CellRegistry::clear_bound_to empties the cells bound to a loop instance
/// whenever that loop advances -- this block's own loop at each of its
/// batches, a nested loop at each of its. A bound cell can therefore lose its
/// value under a mark that still says "skip", and its production would be
/// elided while its consumers still need it. An unbound \c produce_if_absent
/// cell is reachable by no clear at all (the implicit "its scope's innermost
/// loop" rule exempts exactly these), so a mark on it stays true.
[[nodiscard]] inline bool ordered_visit_skip_seedable(TableCell const& c) {
  return c.produce_if_absent && detail::bound_instances(c).empty();
}

/// \overload A skipped VISIT of \p c: one read per leg.
inline void ordered_forgo_visit(CellRegistry& registry, ForgoPlan const& plan,
                                CellId c) {
  ordered_forgo_reads(registry, plan, c,
                      plan.table->cells[c].scope.path.size());
}

/// Whether the skip set \p skip covers EVERY production \p block realizes --
/// its own \c BuildStep cells at \p parent_scope + \p block, every nested
/// block's productions (recursively), and the \c Assemble cell of each of its
/// outputs at \p parent_scope. Such a block has nothing to do this visit: its
/// results are all resident already and its whole batch loop is skipped. A
/// production the table has no cell for is reported NOT skipped, so the block
/// still runs and the step's own lookup raises the table/schedule
/// disagreement with its full diagnostic.
///
/// \p skip is the VISIT's skip set (the call-wide one plus the
/// \c produce_if_absent cells the registry currently holds, closed under the
/// same consumer rule -- see \c run_ordered_contracted_block), so a block
/// whose only output is a loop-invariant Assemble that is already assembled
/// is skipped whole on the enclosing loop's later batches, rather than
/// re-running every step to feed an Assemble that will not be performed.
[[nodiscard]] inline bool ordered_block_fully_skipped(
    CellRegistry const& registry, container::vector<char> const& skip,
    ScopeBlock const& block, CellScope const& parent_scope) {
  CellScope inner = parent_scope;
  inner.path.push_back({block.level.key(), block.latitude_ordinal});
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value)) {
      auto const c = registry.build_cell_at(b->value_id, inner);
      if (!c || !skip[*c]) return false;
    } else if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
      if (!ordered_block_fully_skipped(registry, skip, *child, inner))
        return false;
    } else {
      return false;
    }
  }
  for (auto const& [ovid, okind] : block.outputs) {
    (void)okind;
    auto const a = registry.assemble_cell_at(ovid, parent_scope);
    if (!a || !skip[*a]) return false;
  }
  return true;
}

/// The accounting half of a whole-block skip (\c ordered_block_fully_skipped):
/// every production the block will not perform forgoes the reads it will not
/// make, with the collapsed per-epoch count (see \c ordered_forgo_reads).
/// Walks exactly the cells that walk enumerates.
inline void ordered_forgo_block(CellRegistry& registry, ForgoPlan const& plan,
                                ScopeBlock const& block,
                                CellScope const& parent_scope,
                                std::size_t base_depth) {
  CellScope inner = parent_scope;
  inner.path.push_back({block.level.key(), block.latitude_ordinal});
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value)) {
      if (auto const c = registry.build_cell_at(b->value_id, inner))
        ordered_forgo_reads(registry, plan, *c, base_depth);
    } else if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
      ordered_forgo_block(registry, plan, *child, inner, base_depth);
    }
  }
  bool built_here = false;
  for (auto const& [ovid, okind] : block.outputs) {
    (void)okind;
    auto const a = registry.assemble_cell_at(ovid, parent_scope);
    if (!a) continue;
    ordered_forgo_reads(registry, plan, *a, base_depth);
    // An Assemble whose per-batch source is an IMPLICIT build (a Build cell
    // at this block's own scope that no step of the block builds -- the
    // schedule fuses the reduction/scatter with an operand contraction, so it
    // emits no BuildStep) elides that production too when it is skipped: the
    // Assemble step is the only site that would have run it, so its own reads
    // are owed here and nowhere else.
    CellId const src = plan.table->cells[*a].production.source;
    TableCell const& sc = plan.table->cells[src];
    built_here = false;
    for (Step const& st : block.steps)
      if (auto const* b = std::get_if<BuildStep>(&st.value))
        if (b->value_id == ovid) built_here = true;
    if (!built_here && sc.production.kind == ProductionKind::Build &&
        sc.scope == inner)
      ordered_forgo_reads(registry, plan, src, base_depth);
  }
}

/// The current batch range of loop instance \p key, read off \p ctx (the
/// enclosing realized loops plus, inside a block's batch loop, that block's
/// own entry). Nullopt when the loop is not open here -- which every caller
/// turns into a throw naming what it was trying to bind.
[[nodiscard]] inline std::optional<std::pair<std::size_t, std::size_t>>
ordered_range_of(eval::BatchContext const& ctx, LoopKey const& key) {
  for (auto const& e : ctx)
    if (detail::same_key(e.level.key(), key)) return e.range;
  return std::nullopt;
}

///
/// \brief Realize one \c ScopeBlock's batch loop against \p parent_cache --
/// the ordered-schedule counterpart of \c scope_executor.hpp's \c
/// detail::walk_scope, executed ENTIRELY on the cell table (the
/// explicit-value-cells design, section 4).
///
/// \details \p block's own \c steps are a topologically ORDERED interleaving
/// of \c BuildStep's (each one produces the value's Build CELL at this
/// block's exact scope) and nested child \c ScopeBlock steps (realized
/// recursively, in full, once per batch of THIS loop). \p block's \c outputs
/// are its ASSEMBLE steps: each names the cell -- at the PARENT scope -- that
/// this block's batches assemble, either by summing the per-batch partials
/// (\c AccumulateSum: a reduced axis) or by scattering them into disjoint
/// slices of one destination (\c AccumulateScatter: a carried axis). Both
/// read the per-batch form from the table (\c Production::source), so a
/// block-close handoff spends the same declared life any other read does.
///
/// Per batch: every cell bound to this loop instance is dropped from the
/// registry (\c CellRegistry::clear_bound_to -- the per-batch reset,
/// expressed on cells), the batch context is extended by this block's axis
/// over the batch's element range, the steps run in schedule order, and each
/// output folds its per-batch partial into its running result. At close each
/// output's running result becomes its Assemble cell's value. There is no
/// home walk, no home slot, no reuse table and no store into a scope cache:
/// the registry owns every result, and the table says where each lives and
/// how long.
///
/// The per-block scratch is a BARE child cache: it holds no entries at all
/// (the table owns storage) and exists only to carry this block's batch
/// context and to let the hook lookups -- backend array-ops, the cell read
/// resolver, the peak monitor -- fall through to \p parent_cache.
///
/// \par A nest's pass blocks
/// A nest holding a forced-split axis realizes its loop as that nest's pass
/// blocks, one per pass (latitude = pass), as sibling \c ScopeBlock \c Step's
/// at the SAME nesting level rather than one nested inside the other, run in
/// schedule order. No special-casing is needed: sibling steps run
/// sequentially, the topological sort orders the pass blocks ascending by
/// pass, and a later pass's reads name an earlier pass's assembled cell
/// explicitly (scopes carry their latitude, so each pass's cells are
/// distinct).
///
/// \param skip The cache-halt skip set over cells (\c
///        ordered_cache_halt_skip), computed once per call; this block adds
///        the \c produce_if_absent cells currently held (its own visit's
///        seeds) and re-closes it.
/// \param forgo_plan The sources each consumer cell reads (\c
///        ordered_forgo_plan), so a SKIPPED production can still spend the
///        reads it will not perform and its sources reach the end of their
///        declared lives.
/// \param built The run-completeness ledger, marked at the exact site each
///        scheduled value is produced (or deliberately skipped).
///
template <Trace EvalTrace, typename node_t, typename F, typename N, bool FHC>
void run_ordered_contracted_block(
    ScopeBlock const& block,
    std::unordered_map<std::size_t, node_t> const& vmap,
    RichSchedule const& rich, OrderedSchedule const& ordered,
    F const& leaf_evaluator, CacheManager<N, FHC>& parent_cache,
    std::function<std::size_t(Index const&)> const& target,
    typename CacheManager<N, FHC>::BatchContext const& ectx,
    container::vector<char>& built,
    std::function<bool(node_t const&)> const& is_volatile,
    CellTable const* table, CellRegistry& registry, CellReadResolver& resolver,
    container::vector<char> const& skip, ForgoPlan const& forgo_plan) {
  using Cache = CacheManager<N, FHC>;
  using BatchContext = typename Cache::BatchContext;
  // Threaded for symmetry with the entry point and for the recursion below;
  // this block consults neither (the table already carries every
  // volatility-derived decision: persistence and lives are its own).
  (void)is_volatile;
  SEQUANT_ASSERT(table && &registry.table() == table &&
                 "run_ordered_contracted_block: the registry and the table "
                 "passed here must be the same table");

  // Backend array-ops (zero destination + axis chunking), sourced from the
  // cache chain (the backend -- mpqc's registries or a test's leaf source --
  // wires the root cache). A batched block cannot be realized without it.
  BackendArrayOps const* const aops = parent_cache.array_ops();
  SEQUANT_ASSERT(aops &&
                 "evaluate_ordered_schedule: batched eval requires backend "
                 "array-ops (CacheManager::set_array_ops)");

  // R4: loud guard on this block's batch-mode kind. The batch-loop primitive
  // below realizes Contracted and External blocks UNIFORMLY (their difference
  // is carried entirely by each Assemble cell's own kind), so both are
  // supported; any OTHER BatchModeType value is a schedule this executor
  // cannot interpret and must refuse loudly rather than silently mis-run.
  SEQUANT_ASSERT((block.kind == BatchModeType::Contracted ||
                  block.kind == BatchModeType::External) &&
                 "evaluate_ordered_schedule: unsupported ScopeBlock batch-mode "
                 "kind -- only Contracted/External are realizable");

  auto const resolve = [&](std::size_t vid) -> node_t const& {
    auto const hash = rich.cells[vid].hash;
    auto const it = vmap.find(hash);
    SEQUANT_ASSERT(it != vmap.end() &&
                   "evaluate_ordered_schedule: a loop-block value_id was not "
                   "found in the forest's value-node map");
    return it->second;
  };

  // A cell holds the CANONICAL orientation (see CellRegistry's own doc);
  // evaluate_impl returns the node's ORIENTED result, and the phase is an
  // involution, so a production converts by multiplying it back in -- exactly
  // what CacheManager::store did with apply_phase before storage moved onto
  // the table. Every reader (the resolver's fetch, and pre_results) applies
  // the node's phase once more and so sees the oriented value again.
  auto const canonical = [](node_t const& nd, ResultPtr r) -> ResultPtr {
    auto const ph = nd->canon_phase();
    // Null passes through untouched (never dereferenced): a missing result is
    // diagnosed where it is read, not here.
    return (!r || ph == 1) ? std::move(r) : r->mult_by_phase(ph);
  };

  CellScope const parent_scope = current_scope(ectx);
  CellScope const block_scope = current_scope(ectx, block);

  // THIS VISIT's skip set: the call-wide one plus the `produce_if_absent`
  // cells the registry currently holds -- such a cell is not re-produced
  // while it is resident, so it is a skipped consumer for the purpose of
  // deciding whether anything downstream of it still has to run -- closed
  // under the same by-source consumer rule.
  //
  // A seed must be a cell that CANNOT LOSE its value while this set is in
  // use -- and this set is in use for every batch of this block AND inside
  // every nested block, since it is handed down to them. That is exactly what
  // `ordered_visit_skip_seedable` decides (see its own doc comment). The
  // call-wide set is untouched by all this: it is decided from persistence
  // and consumer death, neither of which a per-batch clear can invalidate.
  container::vector<char> local_skip_storage;
  container::vector<char> const* skip_p = &skip;
  {
    bool seeded = false;
    for (CellId c = 0; c < table->cells.size(); ++c) {
      if (skip[c] || !ordered_visit_skip_seedable(table->cells[c])) continue;
      if (!registry.peek(c)) continue;
      if (!seeded) {
        local_skip_storage = skip;
        seeded = true;
      }
      local_skip_storage[c] = 1;
    }
    if (seeded) {
      local_skip_storage =
          ordered_skip_closure(*table, std::move(local_skip_storage));
      skip_p = &local_skip_storage;
    }
  }
  container::vector<char> const& vskip = *skip_p;

  // Cache-halt at BLOCK granularity: every production of this block (its own
  // steps', its descendants' and its outputs') is already resident, so the
  // whole batch loop is dead work this visit. Forgo the reads those
  // productions will not perform, and mark them accounted for so the
  // run-completeness ledger does not mistake the skip for a gap.
  if (ordered_block_fully_skipped(registry, vskip, block, parent_scope)) {
    ++ordered_last_block_skips_slot();
    ordered_forgo_block(registry, forgo_plan, block, parent_scope,
                        parent_scope.path.size());
    container::vector<std::size_t> ids;
    collect_production_ids(block, ids);
    for (std::size_t vid : ids)
      if (vid < built.size()) built[vid] = 1;
    if (std::getenv("SEQUANT_UT_BLOCK_DIAG"))
      std::cerr << "[BLOCK] axis=" << toUtf8(block.axis.full_label())
                << " depth=" << block.level.depth
                << " slot=" << block.level.loop_slot
                << " SKIPPED WHOLE (every production resident)" << std::endl;
    return;
  }

  // The value_ids this block BUILDS with a step of its own: an output whose
  // per-batch source cell is a Build at THIS scope that no step builds is an
  // IMPLICIT per-batch build (the schedule fuses the reduction/scatter with an
  // operand contraction, so it emits no BuildStep of its own) and is evaluated
  // by the Assemble step itself.
  std::unordered_set<std::size_t> built_here;
  for (Step const& s : block.steps)
    if (auto const* b = std::get_if<BuildStep>(&s.value))
      built_here.insert(b->value_id);

  // This block's ASSEMBLE steps: one cell per output entry, at the parent
  // scope, resolved once (a miss is a table/schedule disagreement, which must
  // be loud before any batch runs rather than half way through the loop).
  container::vector<CellId> out_cells(block.outputs.size(), 0);
  container::vector<char> out_skip(block.outputs.size(), 0);
  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    auto const vid = block.outputs[k].first;
    SEQUANT_ASSERT((block.outputs[k].second == OutputKind::AccumulateSum ||
                    block.outputs[k].second == OutputKind::AccumulateScatter) &&
                   "evaluate_ordered_schedule: unsupported escape OutputKind");
    auto const a = registry.assemble_cell_at(vid, parent_scope);
    if (!a)
      throw Exception(
          "evaluate_ordered_schedule: no Assemble cell for escaping value " +
          std::to_string(vid) +
          " at the parent scope (cell table/schedule disagreement), block "
          "depth " +
          std::to_string(block.level.depth) + " slot " +
          std::to_string(block.level.loop_slot));
    out_cells[k] = *a;
    // Skip this Assemble step when the cache-halt set covers its cell, and
    // ALSO when the cell is `produce_if_absent` and the registry already
    // holds it: such a cell is invariant to the loop its scope sits inside
    // (this block's parent loop), so the enclosing loop's later batches
    // re-enter this block and must REUSE the assembled value rather than
    // assemble it a second time (the Build step's rule below, applied to the
    // other production kind -- spec section 4 item 2).
    out_skip[k] = vskip[*a] || (table->cells[*a].produce_if_absent &&
                                registry.peek(*a) != nullptr);
  }

  // The per-block scratch: a BARE child cache. It registers nothing and
  // stores nothing (with a resolver wired, evaluate_impl neither probes nor
  // fills a scope cache); it carries this block's batch context and lets the
  // backend/resolver/monitor hooks fall through to the parent chain.
  auto bs_cache = Cache::empty();
  bs_cache.set_parent(&parent_cache);

  // Batch chunks over the loop axis, sourced per-space by the backend (no
  // carrier array is consulted).
  container::svector<std::pair<std::size_t, std::size_t>> const batches =
      aops->axis_batches(block.axis, target(block.axis));

  if (log::printing()) {
    BatchContext s = ectx;
    s.push_back({block.axis, block.level, {0, 0}, std::nullopt});
    log::log("BatchGroup", "Begin",
             std::format("{} steps over {} batches of {} {}",
                         block.steps.size(), batches.size(),
                         toUtf8(std::wstring(block.axis.space().base_key())),
                         log::scope_annot(s)));
  }
  {
    static bool const sched_dump = std::getenv("SEQUANT_SCHED_DUMP") != nullptr;
    if (sched_dump) {
      std::cerr << "ORDERED_RUN_BLOCK {\"kind\":\""
                << (block.kind == BatchModeType::External ? "external"
                                                          : "contracted")
                << "\",\"mode\":\"" << toUtf8(block.axis.full_label())
                << "\",\"depth\":" << block.level.depth
                << ",\"slot\":" << block.level.loop_slot
                << ",\"lat\":" << block.latitude_ordinal
                << ",\"blocks\":" << batches.size()
                << ",\"steps\":" << block.steps.size() << ",\"outs\":[";
      for (std::size_t k = 0; k != block.outputs.size(); ++k)
        std::cerr << (k ? "," : "") << "{\"h\":"
                  << (rich.cells[block.outputs[k].first].hash % 100000)
                  << ",\"cell\":" << out_cells[k]
                  << ",\"skip\":" << (int)out_skip[k] << "}";
      std::cerr << "]}\n";
    }
  }

  // Running results of this block's Assemble steps: a summed accumulator or a
  // scattered-into destination, one per output.
  container::vector<ResultPtr> acc(block.outputs.size());
  container::vector<ResultPtr> dest(block.outputs.size());

  for (auto const& [e_lo, e_hi] : batches) {
    if (e_lo == e_hi) continue;
    // The per-batch reset, expressed on cells: drop every cell bound to THIS
    // block's own loop instance (the batch just ended), so a stale prior-batch
    // cell is never read as this batch's.
    registry.clear_bound_to(block.level.key());
    BatchContext ctx = ectx;
    ctx.push_back({block.axis, block.level, {e_lo, e_hi}, std::nullopt});
    bs_cache.set_batch_context(ctx);

    for (Step const& step : block.steps) {
      if (auto const* build = std::get_if<BuildStep>(&step.value)) {
        // This BuildStep's own Build cell -- at this exact scope (this block,
        // entered from ectx) -- is the consumer every operand fetch inside the
        // evaluate_impl call below resolves against (CellReadResolver::fetch).
        // A miss means the table and the schedule tree disagree about what
        // this step builds.
        auto const build_cell =
            registry.build_cell_at(build->value_id, block_scope);
        if (!build_cell)
          throw Exception(
              "evaluate_ordered_schedule: no Build cell for value " +
              std::to_string(build->value_id) +
              " at this scope (cell table/schedule disagreement), block "
              "depth " +
              std::to_string(block.level.depth) + " slot " +
              std::to_string(block.level.loop_slot));
        built[build->value_id] = 1;
        // A loop-invariant cell homed inside a loop is produced on its FIRST
        // visit and reused by every later batch (the table's own flag; the
        // registry keeps it until an instance it IS bound to clears it).
        // RESIDENCY decides this FIRST, ahead of the skip set: a
        // `produce_if_absent` cell that a per-batch clear has emptied must be
        // re-produced, whatever an enclosing scope's per-visit set -- computed
        // before that clear -- still says about it.
        bool const pia_resident = table->cells[*build_cell].produce_if_absent &&
                                  registry.peek(*build_cell) != nullptr;
        // Cache-halt second: nothing reads this cell this visit. (A
        // `produce_if_absent` cell that is NOT resident still falls through to
        // here, so a cell the call-wide set has genuinely killed is not
        // rebuilt for nobody -- but a cell that is merely marked by some
        // enclosing scope's per-visit seed can no longer elide a production
        // the clear above made necessary, because such a cell is never seeded
        // any more: see the per-visit set at this block's entry.)
        if (pia_resident || vskip[*build_cell]) {
          // Resident and reused, or cache-halted. Either way the reads this
          // production will not perform are still owed to their sources.
          ordered_forgo_visit(registry, forgo_plan, *build_cell);
          continue;
        }
        // Read the env ONCE per translation unit, not once per step of every
        // batch of every block.
        static bool const block_diag_build =
            std::getenv("SEQUANT_UT_BLOCK_DIAG") != nullptr;
        if (block_diag_build)
          std::cerr << "[BLOCK] axis=" << toUtf8(block.axis.full_label())
                    << " batch=[" << e_lo << "," << e_hi
                    << ") BUILD vid=" << build->value_id
                    << " cell=" << *build_cell
                    << " hash=" << (rich.cells[build->value_id].hash % 100000u)
                    << std::endl;
        resolver.begin_consumer(*build_cell);
        node_t const& build_node = resolve(build->value_id);
        registry.set(
            *build_cell,
            canonical(build_node, evaluate_impl<EvalTrace>(
                                      build_node, leaf_evaluator, bs_cache)));
      } else if (auto const* child = std::get_if<ScopeBlock>(&step.value)) {
        run_ordered_contracted_block<EvalTrace>(
            *child, vmap, rich, ordered, leaf_evaluator, bs_cache, target, ctx,
            built, is_volatile, table, registry, resolver, vskip, forgo_plan);
      } else {
        // R4: the Step variant has exactly BuildStep/ScopeBlock alternatives;
        // a valueless-by-exception or future third alternative is a schedule
        // this executor cannot interpret.
        SEQUANT_ASSERT(false &&
                       "evaluate_ordered_schedule: unsupported Step variant");
      }
    }

    // ---- Assemble steps: fold this batch's partial into each output. ----
    for (std::size_t k = 0; k != block.outputs.size(); ++k) {
      if (out_skip[k]) {
        // The Assemble is not performed this batch, but the read it declares
        // of its per-batch source is still owed (the source may be produced
        // by a step of this very block, once per batch).
        ordered_forgo_visit(registry, forgo_plan, out_cells[k]);
        // ... and when that source is an IMPLICIT per-batch build -- a Build
        // cell at this block's scope that no step builds, whose only
        // production site IS this Assemble step -- the skip elides that
        // production too, so ITS reads are owed here and nowhere else.
        CellId const skipped_src = table->cells[out_cells[k]].production.source;
        TableCell const& skipped_src_cell = table->cells[skipped_src];
        if (skipped_src_cell.production.kind == ProductionKind::Build &&
            skipped_src_cell.scope == block_scope &&
            !built_here.count(block.outputs[k].first))
          ordered_forgo_visit(registry, forgo_plan, skipped_src);
        continue;
      }
      auto const vid = block.outputs[k].first;
      TableCell const& a = table->cells[out_cells[k]];
      CellId const src = a.production.source;
      TableCell const& s = table->cells[src];
      static bool const block_diag_assemble =
          std::getenv("SEQUANT_UT_BLOCK_DIAG") != nullptr;
      if (block_diag_assemble)
        std::cerr << "[BLOCK] axis=" << toUtf8(block.axis.full_label())
                  << " batch=[" << e_lo << "," << e_hi
                  << ") ASSEMBLE vid=" << vid << " cell=" << out_cells[k]
                  << " src=" << src << " kind="
                  << (a.production.assemble == AssembleKind::Sum ? "SUM"
                                                                 : "SCATTER")
                  << std::endl;
      // An IMPLICIT per-batch build: the source is this value's own Build cell
      // at this block's scope, but no step of this block builds it (the
      // schedule fuses the reduction/scatter with an operand contraction, so
      // it emits no BuildStep). Produce it here, as its own consumer's step
      // would have.
      if (s.production.kind == ProductionKind::Build &&
          s.scope == block_scope && !built_here.count(vid)) {
        SEQUANT_ASSERT(s.value_id == vid &&
                       "evaluate_ordered_schedule: an Assemble's per-batch "
                       "source names a different value than the Assemble");
        resolver.begin_consumer(src);
        node_t const& part_node = resolve(vid);
        registry.set(src, canonical(part_node,
                                    evaluate_impl<EvalTrace>(
                                        part_node, leaf_evaluator, bs_cache)));
      }
      // The read the Assemble DECLARES of its source (the table charged the
      // source +1 life for it): spending it here is what lets the source's
      // storage be released at its true last use.
      bool exhausted = false;
      ResultPtr part = registry.read(src, &exhausted);
      if (a.production.assemble == AssembleKind::Sum) {
        if (!acc[k]) {
          // The first batch's partial SEEDS the accumulator, which every later
          // batch then mutates in place. That is only safe when this read took
          // ownership: a source with life left, or a persistent one, is still
          // going to be read again from the very same buffer.
          acc[k] = exhausted ? std::move(part) : part->clone();
        } else {
          acc[k]->add_inplace(*part);
        }
      } else {
        if (!dest[k]) {
          // The destination is sized from the ASSEMBLE CELL'S OWN FORM: the
          // value's full index list, narrowed to the current batch of every
          // loop instance the cell itself is sliced by (its enclosing loops --
          // the instance being assembled is in the scatter map, not here). No
          // inference from the escaped axis's position on the node.
          dest[k] = aops->make_zeros(resolve(vid)->canon_indices());
          for (auto const& [pos, key] : a.sliced) {
            auto const range = ordered_range_of(ctx, key);
            if (!range)
              throw Exception("evaluate_ordered_schedule: Assemble cell#" +
                              std::to_string(out_cells[k]) + " (value " +
                              std::to_string(vid) + ") slices position " +
                              std::to_string(pos) +
                              " on a loop instance not in the batch context");
            dest[k] = dest[k]->slice_mode(pos, range->first, range->second);
          }
        }
        if (a.production.scatter_map.empty())
          throw Exception("evaluate_ordered_schedule: Assemble cell#" +
                          std::to_string(out_cells[k]) + " (value " +
                          std::to_string(vid) + ") scatters nothing");
        // Exactly ONE scattered position per Assemble. The loop below writes
        // the SAME per-batch partial at every position of the map, which is
        // only the right thing when there is one: a source sliced at two
        // positions by one loop instance is a joint sub-block, and writing it
        // twice (once per position, each time over the full extent of the
        // other) would be wrong. The builder emits one entry per instance and
        // an instance slices one position of a value here, so a second entry
        // means a schedule shape this executor does not implement.
        SEQUANT_ASSERT(a.production.scatter_map.size() == 1 &&
                       "evaluate_ordered_schedule: an Assemble scattering more "
                       "than one position would need a joint sub-block write "
                       "(unsupported)");
        for (auto const& [pos, key] : a.production.scatter_map) {
          auto const range = ordered_range_of(ctx, key);
          if (!range)
            throw Exception("evaluate_ordered_schedule: Assemble cell#" +
                            std::to_string(out_cells[k]) + " (value " +
                            std::to_string(vid) + ") scatters position " +
                            std::to_string(pos) +
                            " from a loop instance not in the batch context");
          dest[k]->write_into_slice(*part, pos, range->first, range->second);
        }
      }
    }
  }

  // ---- Block close: each Assemble step's running result IS its cell. ----
  for (std::size_t k = 0; k != block.outputs.size(); ++k) {
    auto const vid = block.outputs[k].first;
    built[vid] = 1;
    if (out_skip[k]) continue;  // cache-halt: the resident cell stands
    ResultPtr& out =
        table->cells[out_cells[k]].production.assemble == AssembleKind::Sum
            ? acc[k]
            : dest[k];
    SEQUANT_ASSERT(out &&
                   "evaluate_ordered_schedule: a loop block realized zero "
                   "batches for an Assemble step");
    registry.set(out_cells[k], std::move(out));
    if constexpr (::sequant::detail::trace(EvalTrace))
      log::cache(resolve(vid), parent_cache,
                 log::label(resolve(vid), parent_cache.batch_context()) +
                     " [assembled]");
  }
}

inline std::size_t& ordered_last_cell_table_size_slot() {
  static std::size_t n = 0;
  return n;
}
/// Diagnostic: number of cells in the table built by the most recent
/// \c run_ordered_schedule_pre_results call (test-facing; not thread-safe).
[[nodiscard]] inline std::size_t ordered_last_cell_table_size() {
  return ordered_last_cell_table_size_slot();
}

/// Diagnostic: what the cell registry was still holding when the most recent
/// \c run_ordered_schedule_pre_results call returned. \c live is the whole
/// live byte total; \c persistent is the part held by cells the table marks
/// persistent (they survive on purpose, into the next evaluation); \c roots
/// is the part held by the forest roots' own cells -- the exact CELLS the
/// results were taken from, not every cell of a root's value (the caller's \c
/// pre_results are private COPIES of them, and nobody in the table reads
/// them, so those cells keep holding their own until the registry dies with
/// the call). \c live beyond
/// those two is a non-persistent intermediate that never reached the end of its
/// declared life -- exactly what a missing \c CellRegistry::forgo at a
/// skipped production leaves behind.
struct OrderedRegistryResidency {
  std::size_t live = 0, persistent = 0, roots = 0;
};
inline OrderedRegistryResidency& ordered_last_registry_residency_slot() {
  static OrderedRegistryResidency r;
  return r;
}
/// \return the residency of the most recent \c
/// run_ordered_schedule_pre_results call (test-facing; not thread-safe).
[[nodiscard]] inline OrderedRegistryResidency
ordered_last_registry_residency() {
  return ordered_last_registry_residency_slot();
}

/// Assembles the \c CellTableInputs the cell table builder needs from what
/// \c run_ordered_schedule_pre_results already has in hand: the finished
/// schedule (\p ordered, \p rich), the seam it derives the per-occurrence
/// slicing facts from (\p sma), the value-id -> forest-node accessor every
/// homing site in this file already uses (\p resolve), and the node-level
/// volatility predicate (\p is_volatile, possibly empty). \c operands_of is
/// built here from the value's own production tree (its canonical node's two
/// children, resolved back to value ids via a hash -> value-id map built
/// once) so a value contracted with itself contributes one operand entry per
/// LEG, matching the runtime's own per-leg home accesses (see \c
/// CellTableInputs::operands_of's own note on why the de-duplicated
/// dependency graph cannot express that).
template <typename Resolve, typename IsVolatile>
CellTableInputs make_cell_table_inputs(OrderedSchedule const& ordered,
                                       RichSchedule const& rich,
                                       SlicedModeAssignment const& sma,
                                       Resolve const& resolve,
                                       IsVolatile const& is_volatile) {
  CellTableInputs in;
  in.ordered = &ordered;
  in.rich = &rich;
  in.sliced = &sma;
  in.sliced_modes_of = [&](std::size_t vid) {
    auto const& nd = resolve(vid);
    return container::svector<Index>(nd->sliced_modes().begin(),
                                     nd->sliced_modes().end());
  };
  in.volatile_of = [&](std::size_t vid) {
    return is_volatile && subtree_any(resolve(vid), is_volatile);
  };
  // Built ONCE here (not per call of the lambda below), then captured BY
  // VALUE: a reference into this function's own stack would dangle once it
  // returns, and both the caller's use of CellTableInputs (build_cell_table)
  // and CellTableInputs itself never outlive this call chain, so a value copy
  // is cheap and correct.
  std::unordered_map<std::size_t, std::size_t> vid_of_hash;
  for (auto const& vc : rich.cells) vid_of_hash.emplace(vc.hash, vc.value_id);
  in.operands_of = [&resolve, vid_of_hash](std::size_t vid) {
    container::svector<std::size_t> out;
    auto const& nd = resolve(vid);
    if (nd.leaf()) return out;
    for (auto const* child : {&nd.left(), &nd.right()})
      if (auto it = vid_of_hash.find((*child)->hash_value());
          it != vid_of_hash.end())
        out.push_back(it->second);
    return out;
  };
  in.n_batches_of = [](LoopKey const&) { return std::size_t{1}; };
  return in;
}

///
/// \brief Task 4 (multi-root single-DAG eval): the shared CORE every ordered
/// whole-forest entry point (\c evaluate_ordered_schedule's forest-wide SUM
/// and \c evaluate_ordered_multiroot's per-root MAP alike) delegates to --
/// walks \p ordered.root.steps exactly as \c evaluate_ordered_schedule always
/// has (Tasks 1-3: a root-level \c BuildStep built directly, a root-level \c
/// ScopeBlock realized via \c run_ordered_contracted_block) and returns each
/// forest root's own UNPERMUTED, already-built \c value_result, aligned
/// index-for-index with \p forest -- i.e. exactly the \c pre_results
/// \c combine_forest_roots expects, computed but NOT yet consumed by it. Pure
/// extraction of \c evaluate_ordered_schedule's original body (SP3 Tasks 1-4
/// through the original's own \c pre_results loop): no upstream logic
/// (schedule walk, cell table derivation, the run-completeness assert) is
/// touched -- see this file's own \note on why a
/// concatenated multi-root forest gets cross-root CSE for free from
/// \c compute_dag_boulevard's hash-keyed \c ValueCell bucketing (built
/// upstream of this function, in \p rich) with NO new dedup logic needed
/// here: a value shared across two \e independent root trees is just another
/// repeated hash, indistinguishable from a value shared across two summands
/// of one root's own forest, which this same walk already builds once.
///
/// \param forest Same requirement as \c evaluate_ordered_schedule's \p
///        forest: the roots whose results are computed by this call -- either
///        the summand terms of ONE equation (the pre-Task-4 caller) or
///        several INDEPENDENT equations' own root trees (Task 4's multi-root
///        caller); this function does not care which, since it produces one
///        unpermuted, unsummed result per element of \p forest either way.
/// \return Each element of \p forest's own already-built result, unpermuted,
///         same order and length as \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
[[nodiscard]] container::svector<ResultPtr> run_ordered_schedule_pre_results(
    Nodes const& forest, OrderedSchedule const& ordered,
    RichSchedule const& rich, F const& leaf_evaluator,
    CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // Task 5 (explicit value cells, stage 2): the per-(value,sliced-mode)->loop
  // facts (SlicedModeAssignment, value_id-keyed) feed the cell table builder
  // below (make_cell_table_inputs) -- the runtime no longer projects them onto
  // a hash-keyed seam for slice_to_use to consult; every ordered operand read
  // now goes through the cell table's CellReadResolver instead (already-sliced
  // cells, no runtime slice inference).
  SlicedModeAssignment const sliced_mode_assignment =
      compute_sliced_mode_assignment(ordered, rich);

  // hash -> node, resolving a BuildStep's value_id (via rich.cells[vid].hash)
  // to the forest node evaluate_impl builds.
  auto const vmap = build_value_node_map(forest);

  // value_id -> forest node: the same lookup run_ordered_contracted_block's
  // own `resolve` performs (see its definition above this function), built
  // again here since that one closes over the block function's own
  // parameters, not this function's locals.
  auto const resolve = [&](std::size_t vid) -> node_t const& {
    auto const hash = rich.cells[vid].hash;
    auto const it = vmap.find(hash);
    SEQUANT_ASSERT(it != vmap.end() &&
                   "evaluate_ordered_schedule: a value_id was not found in "
                   "the forest's value-node map");
    return it->second;
  };

  // Explicit value cells (SP4 Task 3): build and statically validate the cell
  // table once per call, consuming the sliced_mode_assignment seam just
  // built above -- before any block runs, so a schedule this table cannot
  // describe is refused up front rather than mis-executed. A fixture whose
  // table fails validation is a real finding on that fixture: the assertion
  // is not bypassed.
  //
  // SP4 Task 4: n_batches_of -- Task 3's own callback (make_cell_table_inputs)
  // stubbed this at a constant 1 (self-consistent for Task 3's static-only
  // validation: the builder and the validator both used that same stub, so
  // life came out internally consistent regardless of whether it matched the
  // true read count). Task 4's CellRegistry::read() is the first thing that
  // actually enforces life at runtime, which exposed the stub as a real
  // undercount (a value read from outside n*m nested real batches is read
  // n*m times, not the stub's 1) -- CellRegistry::read throwing "read past
  // its life" on a well-formed schedule. Fixed at its source
  // (ordered_n_batches_by_loop, exact per realized loop instance), passed to
  // BOTH the builder and the validator so their life computations (which
  // must use the identical n_batches_of, per read_multiplicity's own doc)
  // never diverge.
  CellTableInputs cell_table_inputs = make_cell_table_inputs(
      ordered, rich, sliced_mode_assignment, resolve, is_volatile);
  cell_table_inputs.n_batches_of =
      ordered_n_batches_by_loop(ordered, target, cache.array_ops());
  CellTable const cell_table = build_cell_table(cell_table_inputs);
  assert_valid_cell_table(cell_table, ordered.root,
                          cell_table_inputs.n_batches_of);
  if (!cell_table.unresolved.empty())
    throw Exception("evaluate_ordered_schedule: the cell table has " +
                    std::to_string(cell_table.unresolved.size()) +
                    " unresolved sliced positions (see CellTable::unresolved)");
  ordered_last_cell_table_size_slot() = cell_table.cells.size();
  ordered_last_block_skips_slot() = 0;

  // The runtime side of the table -- CellRegistry (the OWNER of every result,
  // with its remaining life) and CellReadResolver (table-driven operand reads
  // for evaluate_impl, which REPLACE the router/access_at probes there
  // entirely once wired), installed on the top-level cache for the duration
  // of this call. vid_of_hash resolves a fetched node's hash to its rich.cells
  // slot (its value_id); a hash absent from this map is not a value of the
  // table (a transient of some production tree), which CellReadResolver::fetch
  // reports as nullopt rather than mis-resolving.
  //
  // Stage 3: the registry OWNS every result, and it is wired to the cache
  // handle's cross-call PersistentValueStore -- it publishes a persistent
  // cell there on every production and seeds itself from it at entry, so a
  // value that survives between repeated evaluations of this schedule (e.g.
  // successive iterations of an iterative solver) comes back without any
  // scope cache having to hold it. The store is keyed by the canonical hash
  // every value id resolves to.
  CellRegistryHooks registry_hooks;
  registry_hooks.persistent = &cache.persistent_values();
  registry_hooks.hash_of = [&rich](std::size_t vid) -> std::size_t {
    return rich.cells[vid].hash;
  };
  CellRegistry registry(cell_table, std::move(registry_hooks));
  registry.seed_persistent();
  std::unordered_map<std::size_t, std::size_t> const cell_vid_of_hash = [&] {
    std::unordered_map<std::size_t, std::size_t> m;
    m.reserve(rich.cells.size());
    for (auto const& vc : rich.cells) m.emplace(vc.hash, vc.value_id);
    return m;
  }();
  CellReadResolver resolver(
      registry,
      [&cell_vid_of_hash](std::size_t h) -> std::optional<std::size_t> {
        auto const it = cell_vid_of_hash.find(h);
        if (it == cell_vid_of_hash.end()) return std::nullopt;
        return it->second;
      });
  struct CellReadResolverGuard {
    CacheManager<N, FHC>& c;
    CellReadResolver* prev;
    ~CellReadResolverGuard() { c.set_cell_read_resolver(prev); }
  } const cell_read_resolver_guard{cache, cache.cell_read_resolver()};
  cache.set_cell_read_resolver(&resolver);

  // Stage 3 (explicit value cells): the registry's bytes feed the peak
  // metering the same way the legacy scope caches' own cache_map_ entries
  // do -- CacheManager::chain_residency() folds in whatever
  // set_external_residency installs, and note_working_set()'s diagnostic
  // liveset capture folds in whatever set_external_liveset installs. Both
  // are looked up by resolve()'s own value_id -> hash convention (cell.
  // value_id indexes rich.cells, exactly as every other value_id in this
  // function does), so a live cell's reported hash matches the hash every
  // other liveset/coloring probe in this file already keys by.
  //
  // Double-counting guard (kept until the stage that removes the legacy
  // scope caches entirely): the table-driven path stores nothing in a scope
  // cache any more, but a caller's own cache handle may still hold a buffer
  // this registry also holds (e.g. a value the caller stored itself before
  // the run). Counting the registry's bytes unconditionally would double the
  // memory those values already contribute via current_residency()'s walk of
  // cache_map_, so a live cell's bytes are only added here when NO alive
  // legacy entry anywhere on the chain holds that same buffer by pointer
  // identity (cache.chain_holds,
  // read-only -- unlike chain_holds_shared it decays no lifetime and does
  // not care whether the legacy entry is shared).
  auto const cell_hash = [&](CellId c) -> std::size_t {
    return rich.cells[cell_table.cells[c].value_id].hash;
  };
  struct ExternalMeteringGuard {
    CacheManager<N, FHC>& c;
    std::function<std::size_t()> prev_residency;
    std::function<void(std::function<void(std::size_t, std::size_t)>)>
        prev_liveset;
    ~ExternalMeteringGuard() {
      c.set_external_residency(std::move(prev_residency));
      c.set_external_liveset(std::move(prev_liveset));
    }
  } const external_metering_guard{cache, cache.external_residency_hook(),
                                  cache.external_liveset()};
  cache.set_external_residency([&registry, &cache]() -> std::size_t {
    std::size_t bytes = 0;
    registry.for_each_live([&](CellId, ResultPtr const& v) {
      // chain_holds protects ONLY buffers a CALLER put in its own scope cache
      // before this run: under the cell-read resolver the executor never
      // stores into a scope cache itself, so on a run this executor drives
      // alone this test never fires and every live cell is counted here.
      if (!cache.chain_holds(v)) bytes += v->size_in_bytes();
    });
    return bytes;
  });
  cache.set_external_liveset(
      [&registry, &cache,
       cell_hash](std::function<void(std::size_t, std::size_t)> emit) {
        registry.for_each_live([&](CellId c, ResultPtr const& v) {
          if (!cache.chain_holds(v)) emit(cell_hash(c), v->size_in_bytes());
        });
      });

  // Stage 3 cache-halt: the skip set over CELLS, computed once here from the
  // table's own dependency edges against what the registry already holds
  // after the persistent seeding above (see ordered_cache_halt_skip). It
  // replaces the forest BFS over the legacy cache's alive entries: a
  // persistent cell that survived the previous evaluation is not re-produced,
  // and neither is any cell whose every consumer is itself skipped.
  container::vector<char> const skip = [&]() {
    PhaseTimer::Scope _pt("B.sched_setup");
    return ordered_cache_halt_skip(cell_table, registry);
  }();
  // The reads a skipped production will not perform, so their sources still
  // reach the end of their declared lives (see ordered_forgo_reads). Uses the
  // SAME per-loop batch counts the table's own life computation used.
  ForgoPlan const forgo_plan = [&]() {
    PhaseTimer::Scope _pt("B.sched_setup");
    return ordered_forgo_plan(cell_table, cell_table_inputs.n_batches_of);
  }();
  if (std::getenv("SEQUANT_UT_BLOCK_DIAG")) {
    std::size_t n_skip = 0;
    for (char const c : skip) n_skip += c ? 1 : 0;
    std::cerr << "[cell-registry] cache-halt skips " << n_skip << " of "
              << skip.size() << " cells" << std::endl;
  }

  // R3: executor-side run-completeness ledger. Set to 1 at the EXACT site
  // each scheduled value is produced -- a root-scope BuildStep below, a
  // block-local BuildStep, or a block's Assemble step at its close (see
  // run_ordered_contracted_block) -- or where the cache-halt skip set
  // deliberately leaves it to the resident cell it already has.
  container::vector<char> built(ordered.num_values, 0);

  // -------- Walk the root block's own steps, in order. --------
  // A root-scope BuildStep produces that value's Build cell at the EMPTY
  // scope; a child ScopeBlock (a realized batch loop, Contracted or External
  // alike) is run via detail::run_ordered_contracted_block, whose Assemble
  // steps produce the cells its outputs escape into. Every result lives in
  // the registry; nothing is homed in, or read back from, a scope cache.
  typename CacheManager<N, FHC>::BatchContext const root_ectx;

  for (Step const& step : ordered.root.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      std::size_t const vid = build->value_id;
      SEQUANT_ASSERT(vid < rich.cells.size());
      auto const it = vmap.find(rich.cells[vid].hash);
      SEQUANT_ASSERT(it != vmap.end() &&
                     "evaluate_ordered_schedule: BuildStep value not found "
                     "in the forest's value-node map");
      // A root-scope BuildStep is a Build cell at the EMPTY scope (the root,
      // never inside a batch loop). The resolver carries this cell as the
      // consumer for the operand reads inside evaluate_impl.
      auto const root_cell = registry.build_cell_at(vid, CellScope{});
      if (!root_cell)
        throw Exception(
            "evaluate_ordered_schedule: no root Build cell for value " +
            std::to_string(vid) + " (cell table/schedule disagreement)");
      built[vid] = 1;  // R3: produced, or deliberately skipped, either way
                       // accounted for.
      if (skip[*root_cell] || (cell_table.cells[*root_cell].produce_if_absent &&
                               registry.peek(*root_cell))) {
        // Cache-halt: nothing reads it this call. Its own reads are still
        // owed to their sources (see ordered_forgo_reads).
        ordered_forgo_visit(registry, forgo_plan, *root_cell);
        continue;
      }
      resolver.begin_consumer(*root_cell);
      // The registry holds the CANONICAL orientation (see CellRegistry's own
      // doc); evaluate_impl returns the oriented result and the phase is an
      // involution, so a production converts by multiplying it back in.
      {
        auto const ph = it->second->canon_phase();
        ResultPtr r =
            evaluate_impl<EvalTrace>(it->second, leaf_evaluator, cache);
        // A null result is recorded as null rather than dereferenced here, so
        // the diagnostic stays the "forest root was never produced" throw at
        // the combine below instead of a crash in the phase conversion.
        registry.set(*root_cell,
                     (!r || ph == 1) ? std::move(r) : r->mult_by_phase(ph));
      }
    } else if (auto const* block = std::get_if<ScopeBlock>(&step.value)) {
      run_ordered_contracted_block<EvalTrace>(
          *block, vmap, rich, ordered, leaf_evaluator, cache, target, root_ectx,
          built, is_volatile, &cell_table, registry, resolver, skip,
          forgo_plan);
    } else {
      // R4: the Step variant has exactly BuildStep/ScopeBlock alternatives; any
      // other state is a schedule this executor cannot interpret.
      SEQUANT_ASSERT(false &&
                     "evaluate_ordered_schedule: unsupported Step variant at "
                     "root scope");
    }
  }

  // R3: run-completeness. Every value_id the schedule PROMISES to produce --
  // every BuildStep (root-level or loop-local Transient) and every block escape
  // output, enumerated by collect_production_ids -- must have been actually
  // built by the walk above. Complements well_formed's STATIC single-producer
  // check, which checks no DUPLICATE production but NOT completeness (see
  // ordered_schedule.hpp well_formed's \note). A gap here means the executor
  // skipped a scheduled value: a silent under-execution to refuse, not ignore.
  {
    container::vector<std::size_t> production_ids;
    collect_production_ids(ordered.root, production_ids);
    for (std::size_t const vid : production_ids) {
      SEQUANT_ASSERT(vid < ordered.num_values &&
                     "evaluate_ordered_schedule: production value_id out of "
                     "range");
      SEQUANT_ASSERT(built[vid] &&
                     "evaluate_ordered_schedule: a scheduled value_id was "
                     "never produced during the run (incomplete execution)");
    }
  }

  // hash -> value_id, to resolve each forest root's own build above into the
  // per-root pre_results below.
  std::unordered_map<std::size_t, std::size_t> hash_to_vid;
  hash_to_vid.reserve(rich.cells.size());
  for (auto const& c : rich.cells) hash_to_vid.emplace(c.hash, c.value_id);

  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);

  container::svector<ResultPtr> pre_results(roots.size());
  // The CELLS the roots' results were taken from -- the residency diagnostic
  // below buckets by these, not by the roots' value ids: another cell of a
  // root's value (an in-block form of it, say) is an ordinary intermediate
  // and its retention must not be excused as "that is just the root".
  container::set<CellId> root_cells;
  for (std::size_t i = 0; i != roots.size(); ++i) {
    auto const vid_it = hash_to_vid.find(roots[i]->hash_value());
    SEQUANT_ASSERT(vid_it != hash_to_vid.end() &&
                   "evaluate_ordered_schedule: forest root not found in the "
                   "schedule's value map");
    std::size_t const vid = vid_it->second;
    // A forest root's ONLY reader is the combine below -- it has zero DAG
    // consumers (else it would not be a root), so its cell is read by nobody
    // in the table and its life is zero. PEEK it (non-decrementing): the cell
    // is the root's Build cell at the root scope, or, for a root produced
    // inside a batch loop and escaped out of it, the Assemble cell the
    // block's close produced at the same root scope.
    auto cell = registry.build_cell_at(vid, CellScope{});
    if (!cell) cell = registry.assemble_cell_at(vid, CellScope{});
    if (!cell)
      throw Exception(
          "evaluate_ordered_schedule: forest root value " +
          std::to_string(vid) +
          " has no root-scope cell (cell table/schedule disagreement)");
    root_cells.insert(*cell);
    ResultPtr ptr = registry.peek(*cell);
    if (!ptr)
      throw Exception("evaluate_ordered_schedule: forest root value " +
                      std::to_string(vid) + " (cell#" + std::to_string(*cell) +
                      ") holds no result at the combine read");
    // NEVER hand out a buffer the registry still holds. The read above is a
    // PEEK: it does not move the value out, so the cell -- and, for a
    // persistent cell, the PersistentValueStore it was published to, and any
    // caller cache entry holding the same buffer -- is still pointing at it.
    // The combine (forest_combine.hpp) accumulates the forest's roots IN
    // PLACE into the first one, so handing out the registry's own buffer
    // would mutate the stored value: a persistent root would then seed the
    // NEXT evaluation from root0 + root1, and two forest roots resolving to
    // one cell would double one of them. Hand out a private copy. (The
    // phase-shifting branch is not a substitute: a backend's \c
    // mult_by_phase may return a shallow handle onto the same tiles.) Only
    // the FIRST root's buffer is mutated (it becomes the accumulator); the
    // other roots are read-only addends, so they are handed out as-is: one
    // copy per evaluation, not one per root.
    if (i == 0) ptr = ptr->clone();
    // Orient the stored value to this root's phase, matching evaluate_impl's
    // own canonical->orientation return convention (apply_phase).
    auto const ph = roots[i]->canon_phase();
    pre_results[i] = (ph == 1) ? std::move(ptr) : ptr->mult_by_phase(ph);
    if (!pre_results[i])
      throw Exception(
          "evaluate_ordered_schedule: forest root was never produced");
  }

  // Diagnostic: what the registry is still holding now that the walk is over
  // (see OrderedRegistryResidency). Computed before the registry dies with
  // this call, and only from what it already tracks.
  {
    OrderedRegistryResidency res;
    registry.for_each_live([&](CellId c, ResultPtr const& v) {
      std::size_t const b = v->size_in_bytes();
      res.live += b;
      if (cell_table.cells[c].persistent)
        res.persistent += b;
      else if (root_cells.count(c))
        res.roots += b;
    });
    ordered_last_registry_residency_slot() = res;
    // Retention diagnostic: name every cell that is still holding a value it
    // should have finished with, with the reads that were supposed to spend
    // its life and whether each of their consumers was skipped -- which is
    // how a missing forgo at a skip site is localized (it found the one that
    // was missing at an Assemble's implicit per-batch source).
    if (std::getenv("SEQUANT_UT_RESIDENCY_DIAG"))
      registry.for_each_live([&](CellId c, ResultPtr const& v) {
        TableCell const& tc = cell_table.cells[c];
        if (tc.persistent || root_cells.count(c)) return;
        std::cerr << "[resid] cell#" << c << " value " << tc.value_id
                  << " kind=" << (int)tc.production.kind
                  << " scope_depth=" << tc.scope.path.size()
                  << " pia=" << (int)tc.produce_if_absent
                  << " life=" << registry.remaining_life(c)
                  << " declared=" << tc.life << " bytes=" << v->size_in_bytes()
                  << std::endl;
        for (Read const& r : cell_table.reads)
          if (r.source == c)
            std::cerr << "   [resid-read] consumer cell#" << r.consumer
                      << " (value " << cell_table.cells[r.consumer].value_id
                      << ", scope_depth "
                      << cell_table.cells[r.consumer].scope.path.size()
                      << ", kind "
                      << (int)cell_table.cells[r.consumer].production.kind
                      << ") skip=" << (int)skip[r.consumer] << " mult="
                      << detail::read_multiplicity(
                             tc, cell_table.cells[r.consumer].scope,
                             cell_table_inputs.n_batches_of)
                      << std::endl;
        for (CellId o = 0; o < cell_table.cells.size(); ++o)
          if (cell_table.cells[o].production.kind == ProductionKind::Assemble &&
              cell_table.cells[o].production.source == c)
            std::cerr << "   [resid-asm] assemble cell#" << o << " (value "
                      << cell_table.cells[o].value_id << ", scope_depth "
                      << cell_table.cells[o].scope.path.size()
                      << ", pia=" << (int)cell_table.cells[o].produce_if_absent
                      << ") skip=" << (int)skip[o] << std::endl;
      });
    if (std::getenv("SEQUANT_UT_BLOCK_DIAG"))
      std::cerr << "[cell-registry] residency live=" << res.live
                << " persistent=" << res.persistent << " roots=" << res.roots
                << std::endl;
  }

  // SP4 Task 4 diagnostic: how many operand reads the table-driven resolver
  // actually served this call (see CellReadResolver::served's own doc
  // comment for what does/doesn't count).
  if (std::getenv("SEQUANT_UT_BLOCK_DIAG"))
    std::cerr << "[cell-registry] resolver served " << resolver.served()
              << " reads" << std::endl;

  return pre_results;
}

}  // namespace detail

///
/// \brief SP3 Tasks 1-3 of the ordered-scope batched-eval design: the
/// ORDERED executor. Walks \c ordered.root.steps in sequence, building one
/// value per root-level \c BuildStep via \c evaluate_impl (which cache-probes
/// every operand at its own \c Stage::Enter, so a value already \c
/// cache.store'd by an earlier step short-circuits its own re-descent -- see
/// \c evaluate_impl's doc comment) and one value per root-level \c ScopeBlock
/// (Contracted or External alike) via \c detail::run_ordered_contracted_block
/// (Task 2: a realized batch loop, its own \c BuildStep's/nested child blocks
/// run per batch on a scratch cache, its \c AccumulateSum outputs summed
/// across batches; Task 3: its \c AccumulateScatter outputs written into a
/// disjoint slice of a pre-sized destination each batch -- both kinds stored
/// at the level the block itself sits in on close, including a nest's own
/// pass-ordered sibling blocks (one per pass, latitude = pass), run in
/// schedule order like any other set of sibling steps), then combines the
/// forest roots' results into the final \c ResultPtr exactly as \c
/// evaluate_whole_scope's shared root-combine loop does (permute-to-\p
/// layout + cross-root \c add_inplace, with the identical Term/Permute/
/// SumInplace trace bookkeeping).
///
/// \param forest The forest whose per-root results are evaluated and summed;
///        same requirement as \c evaluate_whole_scope's \p forest.
/// \param ordered The \c OrderedSchedule (\c build_ordered_schedule) for
///        \p forest, built from \p rich / a \c LegalitySchedule over it.
/// \param rich The \c RichSchedule that produced \p ordered (\c
///        compute_dag_boulevard) -- used to resolve a \c BuildStep's
///        \c value_id back to a forest node via its \c ValueCell::hash and
///        \c build_value_node_map, and to resolve each forest root's own
///        \c value_id for the final combine.
/// \param layout The layout each root's result is permuted to before being
///        summed; same meaning as \c evaluate_whole_scope's \p layout.
/// \param leaf_evaluator The leaf evaluator, as in \c sequant::evaluate.
/// \param cache The cache for common sub-expression elimination, as in
///        \c sequant::evaluate; a value repeated across \c BuildStep's
///        operand needs is deduped exactly as \c evaluate_impl's Checked
///        cache-probe already does for forest descent.
/// \param target Per-index batch partition size (elements): the source of
///        each realized loop block's batch partition (\c
///        detail::run_ordered_contracted_block's \c mode_batches call).
/// \param make_scope_guard Backend scope-guard factory; unused by Tasks 1-2
///        (no backend screening relaxation is threaded into the loop-block
///        walk yet), threaded for interface symmetry with later tasks.
/// \param is_volatile NODE-level volatility predicate (empty means
///        "never volatile"), the lift of \c BatchPolicy::is_volatile_leaf
///        exactly as \c make_evaluator's own \c is_volatile_node lift
///        (eval.hpp) computes it -- threaded down through \c
///        detail::run_ordered_contracted_block's recursion but not yet
///        CONSULTED here (a later task's job: classifying a home value
///        volatile-vs-persistent via \c subtree_any at the homing sites).
/// \return The summed, per-root-permuted result, as \c evaluate_whole_scope
///         (and, for an unbatched forest, forest descent itself) would
///         produce for the same \p forest.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate_ordered_schedule(
    Nodes const& forest, OrderedSchedule const& ordered,
    RichSchedule const& rich, auto const& layout, F const& leaf_evaluator,
    CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;

  // Task 4: the schedule walk itself (Tasks 1-3, plus the cell table
  // derivation and the run-completeness assert) is factored into
  // detail::run_ordered_schedule_pre_results, shared byte-for-byte with
  // evaluate_ordered_multiroot below -- see that function's own doc comment.
  // Nothing about the walk changes here; only what happens to its per-root
  // pre_results output differs (summed here, mapped there).
  container::svector<ResultPtr> pre_results =
      detail::run_ordered_schedule_pre_results<EvalTrace>(
          forest, ordered, rich, leaf_evaluator, cache, target,
          make_scope_guard, is_volatile);

  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);

  // -------- Shared combine: permute each root to layout and sum. --------
  // combine_forest_roots (forest_combine.hpp) is the SAME helper
  // evaluate_whole_scope's tail calls, so the two executors emit
  // byte-identical Term/Permute/SumInplace trace bookkeeping without a
  // hand-synced second copy.
  return combine_forest_roots<EvalTrace>(roots, pre_results, layout, cache);
}

///
/// \brief Task 4 of the multi-root single-DAG eval plan: the MULTI-ROOT
/// ordered entry point. Identical inputs to \c evaluate_ordered_schedule
/// (same \p roots/\p ordered/\p rich/\p leaf_evaluator/\p cache/\p target
/// contract -- \p ordered and \p rich must already have been built from the
/// SAME concatenated \p roots, e.g. via \c compute_dag_boulevard(roots, ...)
/// + \c build_ordered_schedule, exactly as any \c evaluate_ordered_schedule
/// caller already does for its own forest), EXCEPT \p layout widens to \p
/// layouts, one entry per root (Task 7 of the plan): unlike a single summed
/// forest, independent roots need not share a layout (e.g. distinct CC
/// residual annotations R1 `{a;i}` vs R2 `{ab;ij}`). Returns a MAP instead
/// of a SUM: each root's own (permuted to its own \p layouts entry, but NOT
/// cross-root-accumulated) result, aligned index-for-index with \p roots.
///
/// \details Runs the identical schedule walk \c evaluate_ordered_schedule
/// runs (via the same \c detail::run_ordered_schedule_pre_results this
/// function and \c evaluate_ordered_schedule both delegate to -- \c
/// the cell table derivation and the run-completeness assert
/// are UNCHANGED, exercised over the combined schedule exactly as they
/// already are for a multi-summand forest), so a value shared across two
/// \e independent roots is built exactly once for the identical reason a
/// value shared across two summands of one root's own forest already is:
/// \c compute_dag_boulevard buckets by structural hash into one \c ValueCell
/// regardless of which root(s) reference it, and \c well_formed's static
/// single-producer invariant (plus this walk's own run-completeness assert)
/// guarantees the schedule realizes exactly one \c BuildStep for it. The
/// ONLY change from \c evaluate_ordered_schedule is the tail: instead of one
/// forest-wide \c combine_forest_roots call (which sums every root into a
/// single \c ResultPtr), this calls \c combine_forest_roots once PER ROOT,
/// each with a singleton \c {roots[i]}/{pre_results[i]} pair -- reusing the
/// IDENTICAL per-root Term/Permute trace bookkeeping \c combine_forest_roots
/// already emits for each root of a summed forest, but never reaching its
/// cross-root \c add_inplace (a singleton call never has a second root to
/// sum against), so no roots are summed together. For a single-root \p roots
/// input this reduces to a singleton \c combine_forest_roots call byte-
/// identical to what \c evaluate_ordered_schedule would run for that same
/// one-root forest -- the regression anchor: \c evaluate_ordered_multiroot(
/// {r}, ...) == { evaluate_ordered_schedule({r}, ...) }.
///
/// \return One \c ResultPtr per element of \p roots, in \p roots' own order,
///         each permuted to its own \p layouts entry -- NOT summed together.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
[[nodiscard]] container::svector<ResultPtr> evaluate_ordered_multiroot(
    Nodes const& roots, OrderedSchedule const& ordered,
    RichSchedule const& rich, container::svector<std::string> const& layouts,
    F const& leaf_evaluator, CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> const& target,
    [[maybe_unused]] ScopeGuardFactory const& make_scope_guard = {},
    std::function<bool(std::ranges::range_value_t<Nodes> const&)> const&
        is_volatile = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;

  container::svector<ResultPtr> pre_results =
      detail::run_ordered_schedule_pre_results<EvalTrace>(
          roots, ordered, rich, leaf_evaluator, cache, target, make_scope_guard,
          is_volatile);

  container::svector<node_t> root_nodes;
  for (auto&& n : roots) root_nodes.push_back(n);
  SEQUANT_ASSERT(pre_results.size() == root_nodes.size());
  SEQUANT_ASSERT(layouts.size() == root_nodes.size());

  // Per-root combine: reuses combine_forest_roots' Term/Permute bookkeeping
  // one root at a time (a singleton call never reaches its cross-root
  // add_inplace), so each root is permuted to its OWN layouts[i] exactly as
  // it would be as part of a summed forest, but no two roots are ever
  // accumulated together -- the map, not the sum.
  container::svector<ResultPtr> results(root_nodes.size());
  for (std::size_t i = 0; i != root_nodes.size(); ++i) {
    container::svector<node_t> one_root{root_nodes[i]};
    container::svector<ResultPtr> one_pre{std::move(pre_results[i])};
    results[i] =
        combine_forest_roots<EvalTrace>(one_root, one_pre, layouts[i], cache);
  }
  return results;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_EXECUTOR_HPP
