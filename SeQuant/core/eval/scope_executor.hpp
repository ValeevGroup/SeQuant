#ifndef SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
#define SEQUANT_EVAL_SCOPE_EXECUTOR_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/forest_combine.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/member_axis.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/eval/value_node_map.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <any>
#include <array>
#include <cstddef>
#include <format>
#include <functional>
#include <initializer_list>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sequant::eval {

///
/// \brief The WEIGHTED in-block use count of a homed value -- the correct
/// NON-PERSISTENT cache life for a value stored at its home (built once per
/// home block by \c reset(), then drained by exactly its in-block reads).
/// Replaces the \c ensure_hoist_slot MAX-life hack: rather than living until
/// the home scope closes, the value frees as soon as its last in-block consumer
/// is done.
///
/// \details
/// \verbatim
///   life(V) = sum over consumers C of V
///               product over loops L on the path (home(V), scope(C)] of
///                 n_blocks(L)
/// \endverbatim
/// Each \c ValueCell::occurrences entry C is one consumer use-site; its \c ectx
/// holds the loops ENCLOSING that use. The loops strictly INNER to V's home --
/// the ectx modes whose TYPE (\c IndexSpace::base_key()) is not in \c
/// home_modes -- are exactly the path (home(V), scope(C)]: each fires \c
/// n_blocks times per home block, so an occurrence contributes their product
/// and the consumers sum. Every genuine occurrence's \c ectx contains V's home
/// loops (V is homed at the meet of its occurrences), so "ectx minus home" is
/// the inner path. A value all of whose occurrences sit AT its home (no inner
/// loop between) contributes 1 per occurrence -- the ordinary use count. An
/// empty \c home_modes (a whole-nest invariant homed at the root) charges every
/// enclosing loop of each consumer, so it is held across the entire nest, as
/// the old MAX life did, but freed by its true last read rather than never.
///
/// TYPE-keyed throughout (matching the narrow scope tree and \c
/// build_scope_schedule's value->node assignment), so an occurrence that binds
/// a loop mode under a different physical ordinal is still counted.
///
/// \note Used by the whole-scope executor's homing (below). The ordered
/// executor instead uses \c detail::ordered_home_reads (ordered_schedule.hpp),
/// which counts exact home reads from the ORDERED schedule's realized scopes
/// under the read-from-home discipline, rather than the boulevard occurrence
/// ectx this reads.
///
/// \param cell the value whose life is computed.
/// \param n_blocks the batch-block count of a loop mode (1 for a mode that is
///        not a realized loop -- an unbatched mode contributes a factor of 1).
///
[[nodiscard]] inline std::size_t weighted_use_count(
    ValueCell const& cell,
    std::function<std::size_t(Index const&)> const& n_blocks) {
  auto const in_home = [&cell](Index const& m) {
    auto const& bk = m.space().base_key();
    for (auto const& h : cell.home_modes)
      if (h.space().base_key() == bk) return true;
    return false;
  };
  std::size_t life = 0;
  for (auto const& occ : cell.occurrences) {
    std::size_t prod = 1;
    for (auto const& [mode, range] : occ.ectx)
      if (!in_home(mode)) prod *= n_blocks(mode);
    life += prod;
  }
  return life;
}

namespace detail {

// result_position_type / find_leaf_carrying_type / member_external_axis /
// member_contracted_axis -- the DAG-space -> node-physical-mode map -- were
// hoisted to member_axis.hpp so ordered_executor.hpp can apply the identical
// map (a circular include otherwise: scope_executor.hpp includes
// ordered_executor.hpp). Still sequant::eval::detail, so all uses below are
// unchanged.

///
/// \brief The recursive scope-tree walk (Task 4): realizes one batch loop per
/// \c ScopeNode, nesting into \c ScopeNode::children, and produces the result
/// of every \p member root assembled over \p node and all deeper loops, under
/// the enclosing batch context \p ectx held on \p parent_cache.
///
/// \details Three node shapes compose here, driven by \c ScopeNode::kind and
/// whether the node has a child:
///
/// - Contracted LEAF (a single aux loop, no nested loop): the Task-3 path,
///   BYTE-FOR-BYTE -- one shared \c make_batched_scratch over the whole member
///   group, per-block \c evaluate_impl with Enter-stage slice-on-use, and
///   cross-block \c add_inplace ACCUMULATE. Every K-carrying composite shared
///   across members is built once per block on the one scratch.
///
/// - Contracted with a child (aux OUTER, occ INNER): each block builds the
///   values HOMED at this node once (registered in the outer-level cache with
///   their \c weighted_use_count in-block life, so a value is REUSED across the
///   inner blocks that consume it and freed as soon as its last in-block use is
///   done, then rebuilt after \c reset() per outer block), then recurses into
///   the child loop for the members that carry the inner mode and
///   \c evaluate_impl's the members invariant to it, and ACCUMULATES each
///   block's per-member partial. The hoist cache IS the outer loop level (one
///   parent link == one loop), so slice-on-use crosses exactly the inner loops
///   when an inner body fetches an outer-homed value.
///
/// - External (occ scatter): realized PER MEMBER over the member's own physical
///   axis (\c member_external_axis), each on a solo \c make_batched_scratch,
///   the per-block partials SCATTERed by \c write_into_slice into a
///   \c pre_sized_zeros_over_mode destination -- the eval.hpp External
///   primitive (eval.hpp:1955-2045), driven from the schedule rather than
///   reimplemented. External blocks are disjoint slices, so they scatter (never
///   accumulate).
///
/// The backend scope guard is HELD across every loop here (contracted and
/// scatter alike), exactly as Task 3 and eval.hpp's External branch do.
///
/// \return per-member results, aligned with \p members.
///
template <Trace EvalTrace, meta::eval_node node_t, typename F, bool FHC,
          typename ScopeGuardFactory>
[[nodiscard]] container::svector<ResultPtr> walk_scope(
    ScopeNode const& node, container::svector<node_t const*> const& members,
    typename CacheManager<node_t, FHC>::BatchContext const& ectx,
    CacheManager<node_t, FHC>& parent_cache, F const& leaf_evaluator,
    std::function<std::size_t(Index const&)> const& target_batch_size,
    ScopeGuardFactory const& make_scope_guard,
    std::unordered_map<std::size_t, node_t> const& vmap,
    RichSchedule const& rich,
    std::unordered_set<std::size_t> const& root_hashes) {
  using Cache = CacheManager<node_t, FHC>;
  using BatchContext = typename Cache::BatchContext;
  using member_t = std::pair<node_t const*, Index>;

  container::svector<ResultPtr> out(members.size());
  ScopeNode const* const child =
      node.children.empty() ? nullptr : &node.children.front();
  std::wstring const base(node.mode.space().base_key());

  // Backend array-ops (zero destination + axis chunking), sourced from the
  // cache chain -- the SAME source the ordered (DAG) executor reads, so both
  // build identical arrays. nullptr (no backend wired, e.g. a raw unit test)
  // takes the legacy carrier path below.
  BackendArrayOps const* const aops = parent_cache.array_ops();

  // Recurse into the child loop for one member subset, forwarding all context.
  auto recurse = [&](ScopeNode const& c,
                     container::svector<node_t const*> const& subset,
                     BatchContext const& sub_ectx,
                     Cache& sub_parent) -> container::svector<ResultPtr> {
    return walk_scope<EvalTrace, node_t, F, FHC, ScopeGuardFactory>(
        c, subset, sub_ectx, sub_parent, leaf_evaluator, target_batch_size,
        make_scope_guard, vmap, rich, root_hashes);
  };

  if (node.kind == BatchModeType::External) {
    // -------- External SCATTER, realized PER MEMBER over its own axis.
    // --------
    for (std::size_t j = 0; j != members.size(); ++j) {
      node_t const* const m = members[j];
      auto const axopt = member_external_axis(*m, base);
      if (!axopt) {
        // This member does not batch the external mode: build it directly
        // (invariant to the scatter), recursing the child loop if any.
        if (child) {
          out[j] = recurse(*child, container::svector<node_t const*>{m}, ectx,
                           parent_cache)
                       .front();
        } else {
          parent_cache.set_batch_context(ectx);
          out[j] = evaluate_impl<EvalTrace>(*m, leaf_evaluator, parent_cache);
        }
        continue;
      }
      Index const axis = *axopt;
      auto const dm = index_position(*m, axis);  // exact: the member's own axis
      SEQUANT_ASSERT(dm &&
                     "external batch mode is not free on the member's result");
      // Batch chunks over the member's own axis, sourced per-space by the
      // backend (no carrier array is consulted).
      SEQUANT_ASSERT(aops &&
                     "batched whole-scope eval requires backend array-ops "
                     "(CacheManager::set_array_ops)");
      container::svector<std::pair<std::size_t, std::size_t>> const batches =
          aops->axis_batches(axis, target_batch_size(axis));

      // A solo scratch (an external mode is not a shared-final mode): dedups
      // repeats WITHIN the member subtree only. Reuses the same primitive the
      // Task-3 contracted path uses to seed above-homed values sliced-on-use.
      std::vector<member_t> solo{{m, axis}};
      auto bs = sequant::detail::make_batched_scratch(solo, parent_cache);
      for (auto const* s : bs.seeds)
        (void)bs.cache.store(*s, parent_cache.access(*s));
      bs.cache.set_parent(&parent_cache);

      // Loop-structure trace (mirrors eval.hpp's forest BatchGroup markers so
      // the two traces are directly comparable). Gated by the standard trace
      // level via log::printing() -- a first-class trace facility, not
      // env-gated.
      auto const scatter_scope = [&] {
        BatchContext s = ectx;
        s.push_back({axis, {0, 0}});
        return log::scope_annot(s);
      };
      if (log::printing()) {
        log::log("BatchGroup", "Begin",
                 std::format("1 member scattered over {} ext batches {}",
                             batches.size(), scatter_scope()));
        log::log("BatchMember",
                 toUtf8(io::serialization::to_string(to_expr(*m))));
        log::log("BatchAxes",
                 std::format("depth={} picked={}:ext nbatches={} {}",
                             ectx.size(), toUtf8(axis.full_label()),
                             batches.size(), scatter_scope()));
      }
      auto const scope_guard = make_scope_guard(batches.size());
      (void)scope_guard;

      ResultPtr dest;
      for (auto const& [e_lo, e_hi] : batches) {
        if (e_lo == e_hi) continue;
        bs.cache.reset();
        BatchContext ctx = ectx;
        ctx.push_back({axis, {e_lo, e_hi}});
        bs.cache.set_batch_context(ctx);
        ResultPtr part =
            child ? recurse(*child, container::svector<node_t const*>{m}, ctx,
                            bs.cache)
                        .front()
                  : evaluate_impl<EvalTrace>(*m, leaf_evaluator, bs.cache);
        if (!dest) dest = aops->make_zeros((*m)->canon_indices());
        dest->write_into_slice(*part, *dm, e_lo, e_hi);
      }
      SEQUANT_ASSERT(dest);
      if (log::printing()) log::log("BatchGroup", "End", scatter_scope());
      out[j] = std::move(dest);
    }
    return out;
  }

  // -------- Contracted ACCUMULATE over node.mode. --------
  Index const K = node.mode;

  // The schedule's canonical contracted mode K mapped to EACH member's own
  // physical axis (see member_contracted_axis): Index labels are meaningful
  // only within a tree, so a member is sliced on ITS OWN label, never the
  // schedule representative reused across members.
  container::svector<Index> axes;
  axes.reserve(members.size());
  for (auto const* m : members)
    axes.push_back(member_contracted_axis(*m, base).value_or(K));

  if (!child) {
    // Task-3 single aux loop: ONE shared scratch, whole-forest co-evaluation,
    // per-block accumulate. Members are guaranteed to carry K by the caller.
    std::vector<member_t> group;
    group.reserve(members.size());
    for (std::size_t m = 0; m != members.size(); ++m)
      group.emplace_back(members[m], axes[m]);

    auto bs = sequant::detail::make_batched_scratch(group, parent_cache);
    for (auto const* s : bs.seeds)
      (void)bs.cache.store(*s, parent_cache.access(*s));
    bs.cache.set_parent(&parent_cache);

    SEQUANT_ASSERT(aops &&
                   "batched whole-scope eval requires backend array-ops "
                   "(CacheManager::set_array_ops)");
    container::svector<std::pair<std::size_t, std::size_t>> const batches =
        aops->axis_batches(K, target_batch_size(K));

    // Loop-structure trace (mirrors eval.hpp's forest BatchGroup markers;
    // log::printing()-gated, first-class -- not env-gated).
    auto const con_scope = [&] {
      BatchContext s = ectx;
      s.push_back({K, {0, 0}});
      return log::scope_annot(s);
    };
    if (log::printing()) {
      log::log("BatchGroup", "Begin",
               std::format("{} members co-evaluated over {} aux batches {}",
                           members.size(), batches.size(), con_scope()));
      for (auto const* m : members)
        log::log("BatchMember",
                 toUtf8(io::serialization::to_string(to_expr(*m))));
      log::log("BatchAxes", std::format("depth={} picked={}:con nbatches={} {}",
                                        ectx.size(), toUtf8(K.full_label()),
                                        batches.size(), con_scope()));
    }
    container::svector<ResultPtr> acc(members.size());
    auto const scope_guard = make_scope_guard(batches.size());
    (void)scope_guard;
    for (auto const& [e_lo, e_hi] : batches) {
      if (e_lo == e_hi) continue;
      bs.cache.reset();
      for (std::size_t m = 0; m != members.size(); ++m) {
        // Slice each member on ITS OWN physical axis; the element range is the
        // same across members (one global aux tiling), only the label differs.
        BatchContext ctx = ectx;
        ctx.push_back({axes[m], {e_lo, e_hi}});
        bs.cache.set_batch_context(ctx);
        ResultPtr part =
            evaluate_impl<EvalTrace>(*members[m], leaf_evaluator, bs.cache);
        if (!acc[m])
          acc[m] = std::move(part);
        else
          acc[m]->add_inplace(*part);
      }
    }
    for (std::size_t m = 0; m != members.size(); ++m) {
      SEQUANT_ASSERT(acc[m]);
      out[m] = std::move(acc[m]);
    }
    if (log::printing()) log::log("BatchGroup", "End", con_scope());
    return out;
  }

  // Nested contracted loop (aux OUTER over an inner child loop). Members all
  // carry K (outer). Split by whether they carry the child (inner) mode.
  std::wstring const cbase(child->mode.space().base_key());
  auto carries_child = [&](node_t const& m) -> bool {
    if (child->kind == BatchModeType::External)
      return member_external_axis(m, cbase).has_value() ||
             result_position_type(m, cbase).has_value();
    return find_leaf_carrying_type(m, cbase).has_value();
  };
  container::svector<std::size_t> inner_pos, direct_pos;
  for (std::size_t m = 0; m != members.size(); ++m)
    (carries_child(*members[m]) ? inner_pos : direct_pos).push_back(m);

  // Group the inner members by their MAPPED physical outer axis: members that
  // share a physical contracted axis co-evaluate together in ONE child
  // recursion so the child level's make_batched_scratch preserves their
  // cross-member CSE (a sub-intermediate shared among them is built once, not
  // per member) -- the whole point of the grouped path, and what a per-member
  // (singleton) recursion would fragment. Members with distinct physical axes
  // form distinct groups (each sliced on its own label -- antipattern-1 fix).
  // Single-aux-label forests collapse to ONE group. Mirrors the Contracted-leaf
  // grouping above (make_batched_scratch over the members sharing the axis),
  // lifted to the nested level. Grouping is block-independent, so it is
  // computed once here.
  container::svector<Index> inner_group_axis;
  container::svector<container::svector<std::size_t>> inner_groups;
  for (auto p : inner_pos) {
    std::size_t g = 0;
    for (; g != inner_group_axis.size(); ++g)
      if (inner_group_axis[g] == axes[p]) break;
    if (g == inner_group_axis.size()) {
      inner_group_axis.push_back(axes[p]);
      inner_groups.emplace_back();
    }
    inner_groups[g].push_back(p);
  }

  // Block count per loop-mode TYPE in the subtree rooted at THIS node (this
  // loop and every deeper one), read off any member leaf carrying that type
  // under the one global per-type tiling. This is the \c n_blocks a homed
  // value's weighted in-block life is computed against (an inner consumer fires
  // once per block of each loop between the value's home and that consumer).
  std::unordered_map<std::wstring, std::size_t> nblocks_by_type;
  {
    auto add_scope = [&](auto&& self, ScopeNode const& sn) -> void {
      std::wstring const bk(sn.mode.space().base_key());
      if (!nblocks_by_type.count(bk)) {
        // Single source of truth with the actual batch loop (same
        // axis_batches).
        SEQUANT_ASSERT(aops &&
                       "batched whole-scope eval requires backend array-ops "
                       "(CacheManager::set_array_ops)");
        nblocks_by_type.emplace(
            bk, aops->axis_batches(sn.mode, target_batch_size(sn.mode)).size());
      }
      for (auto const& c : sn.children) self(self, c);
    };
    for (auto const& c : node.children) add_scope(add_scope, c);
  }
  std::function<std::size_t(Index const&)> const n_blocks =
      [&nblocks_by_type](Index const& m) -> std::size_t {
    auto const it = nblocks_by_type.find(std::wstring(m.space().base_key()));
    return it == nblocks_by_type.end() ? std::size_t{1} : it->second;
  };

  // Values HOMED at this node (shared composites carrying K but invariant to
  // the inner loop): built once per outer block and reused across every inner
  // block. Skip leaves (cheap, sliced on use) and forest roots (produced, not
  // pre-built). Each is registered in the outer cache with its WEIGHTED
  // in-block use count as its non-persistent life (see weighted_use_count):
  // reset() at each outer block rebuilds it, and it frees as soon as its last
  // in-block consumer is done -- rather than the old ensure_hoist_slot MAX life
  // that held it until the outer block closed. The +1 covers the homed build's
  // own store()-time access (the CacheManager stores by access, decaying once):
  // when a homed value is instead built as a byproduct of an earlier homed
  // value's subtree, no explicit build store happens and the +1 is a harmless
  // one-life overcount (the value still frees at the next reset()), so the life
  // is NEVER an undercount -- the invariant that keeps this a lifetime-only
  // change with byte-identical numerics and per-block build counts.
  using Hasher = sequant::TreeNodeHasher<node_t, FHC>;
  using Comp = sequant::TreeNodeEqualityComparator<node_t>;
  container::svector<node_t> homed;
  std::unordered_map<node_t, std::size_t, Hasher, Comp> homed_life;
  for (auto vid : node.homed_values) {
    auto const h = rich.cells[vid].hash;
    auto const it = vmap.find(h);
    if (it == vmap.end() || it->second.leaf() || root_hashes.count(h)) continue;
    homed.push_back(it->second);
    homed_life.emplace(it->second,
                       weighted_use_count(rich.cells[vid], n_blocks) + 1);
  }

  // The cache for THIS outer loop level, holding the per-block homed values
  // with their weighted lives (above). reset() rebuilds each entry's life and
  // drops its data for the next outer block; the ONE-parent-link-per-loop
  // invariant the Enter-stage slice-on-use relies on is preserved (each nested
  // loop wires exactly one parent link back to here).
  Cache outer{std::move(homed_life)};
  outer.set_parent(&parent_cache);

  SEQUANT_ASSERT(aops &&
                 "batched whole-scope eval requires backend array-ops "
                 "(CacheManager::set_array_ops)");
  container::svector<std::pair<std::size_t, std::size_t>> const batches =
      aops->axis_batches(K, target_batch_size(K));

  // Loop-structure trace (mirrors eval.hpp's forest BatchGroup markers;
  // log::printing()-gated, first-class -- not env-gated). Nested loops emit
  // their own Begin/End via the child recursion, so depth distinguishes levels.
  auto const nested_scope = [&] {
    BatchContext s = ectx;
    s.push_back({K, {0, 0}});
    return log::scope_annot(s);
  };
  if (log::printing()) {
    log::log("BatchGroup", "Begin",
             std::format(
                 "{} members ({} homed) over {} aux batches (nested) {}",
                 members.size(), homed.size(), batches.size(), nested_scope()));
    for (auto const* m : members)
      log::log("BatchMember",
               toUtf8(io::serialization::to_string(to_expr(*m))));
    log::log("BatchAxes",
             std::format("depth={} picked={}:con nbatches={} homed={} {}",
                         ectx.size(), toUtf8(K.full_label()), batches.size(),
                         homed.size(), nested_scope()));
  }
  container::svector<ResultPtr> acc(members.size());
  auto const scope_guard = make_scope_guard(batches.size());
  (void)scope_guard;
  for (auto const& [e_lo, e_hi] : batches) {
    if (e_lo == e_hi) continue;
    outer.reset();

    // Build the outer-homed values ONCE for this block, each sliced on ITS OWN
    // physical axis (a homed composite is a canonical node whose K label need
    // not equal the schedule representative).
    for (auto const& hv : homed) {
      if (outer.alive(hv)) continue;  // shared: already built this block
      BatchContext hctx = ectx;
      hctx.push_back(
          {member_contracted_axis(hv, base).value_or(K), {e_lo, e_hi}});
      outer.set_batch_context(hctx);
      (void)evaluate_impl<EvalTrace>(hv, leaf_evaluator, outer);
    }

    // Each member is looped on ITS OWN physical axis for the outer mode (same
    // element range across members; only the label differs). Direct members
    // (invariant to the inner loop) build on `outer`; inner members recurse the
    // child loop GROUPED by shared physical axis, so members sharing an axis
    // are co-evaluated in one child recursion (cross-member CSE preserved).
    container::svector<ResultPtr> block(members.size());
    for (auto p : direct_pos) {
      BatchContext mctx = ectx;
      mctx.push_back({axes[p], {e_lo, e_hi}});
      outer.set_batch_context(mctx);
      block[p] = evaluate_impl<EvalTrace>(*members[p], leaf_evaluator, outer);
    }
    for (std::size_t g = 0; g != inner_groups.size(); ++g) {
      BatchContext gctx = ectx;
      gctx.push_back({inner_group_axis[g], {e_lo, e_hi}});
      container::svector<node_t const*> gmembers;
      gmembers.reserve(inner_groups[g].size());
      for (auto p : inner_groups[g]) gmembers.push_back(members[p]);
      auto sub = recurse(*child, gmembers, gctx, outer);
      for (std::size_t k = 0; k != inner_groups[g].size(); ++k)
        block[inner_groups[g][k]] = std::move(sub[k]);
    }
    for (std::size_t m = 0; m != members.size(); ++m) {
      SEQUANT_ASSERT(block[m]);
      if (!acc[m])
        acc[m] = std::move(block[m]);
      else
        acc[m]->add_inplace(*block[m]);
    }
  }
  for (std::size_t m = 0; m != members.size(); ++m) {
    SEQUANT_ASSERT(acc[m]);
    out[m] = std::move(acc[m]);
  }
  if (log::printing()) log::log("BatchGroup", "End", nested_scope());
  return out;
}

}  // namespace detail

/// \brief Build a per-node trace-metadata provider from a \c RichSchedule.
///
/// \details Returns a callable, hash -> string, that annotates a value with its
/// SCHEDULE properties the running per-op annotation cannot see: the value's
/// remat HOME dag-scope and, for each of its use-sites, that use's dag-scope.
/// A dag-scope is the comma-joined \c IndexSpace::base_key() list of a mode
/// sequence (no trailing comma -- the same convention as \c log::scope_annot),
/// so the produced string is
///   `home={<home dag-scope>} uses=[{<use dag-scope>},{<use dag-scope>},...]`.
/// The map is built ONCE from \c rich.cells (keyed by \c ValueCell::hash, which
/// equals a node's \c hash_value()), then captured by the returned lambda; an
/// unknown hash yields the empty string. Intended to be installed on
/// \c Logger::instance().eval.node_meta so \c log::slice_home_annot appends it.
inline std::function<std::string(std::size_t)> make_node_meta(
    RichSchedule const& rich) {
  std::unordered_map<std::size_t, std::string> map;
  map.reserve(rich.cells.size());
  auto const dag_scope = [](auto const& modes) {
    std::string s;
    for (auto const& m : modes) {
      if (!s.empty()) s += ",";
      s += toUtf8(m.space().base_key());
    }
    return s;
  };
  for (auto const& cell : rich.cells) {
    std::string uses;
    for (auto const& occ : cell.occurrences) {
      if (!uses.empty()) uses += ",";
      std::string sc;
      for (auto const& e : occ.ectx) {
        if (!sc.empty()) sc += ",";
        sc += toUtf8(e.first.space().base_key());
      }
      uses += std::format("{{{}}}", sc);
    }
    map.emplace(cell.hash, std::format("home={{{}}} uses=[{}]",
                                       dag_scope(cell.home_modes), uses));
  }
  return [map = std::move(map)](std::size_t hash) -> std::string {
    auto const it = map.find(hash);
    return it != map.end() ? it->second : std::string{};
  };
}

///
/// \brief Task 2 of the whole-scope batched DAG execution design (see
/// `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md`):
/// the executor SKELETON -- a scope-tree walk that, for a forest with NO
/// batch loops (the scope tree is root-only), evaluates every forest root and
/// sums the results. Provably equivalent to the existing per-tree forest
/// descent (\c sequant::evaluate(Nodes const&, ...), eval.hpp) for that case.
/// Tasks 3-4 extend this to recurse into \c sched.root.children (one
/// accumulate/scatter block per realized batch loop); NONE of that batch-block
/// machinery lives here yet.
///
/// \details \c ScopeSchedule (Task 1) assigns every distinct forest VALUE --
/// not just forest roots, every hash-distinct node of the whole forest -- to
/// a \c ScopeNode by \c value_id, an index into the \c RichSchedule that built
/// it (see \c peak_profile.hpp's \c ValueCell::occurrences, which is the
/// value_id -> forest-node bridge in general). This function's signature,
/// though, is handed only the projected \c ScopeSchedule (no \c RichSchedule /
/// \c ValueCell in scope), so it cannot walk \c homed_values back to a forest
/// node in general here. For the UNBATCHED case that general bridge would
/// resolve to anyway is exactly "the forest's own top-level trees" (a root
/// occurrence, identified by \c OccurrenceRec::point == \c consumer_point, is
/// precisely a tree of the forest -- see \c scope_schedule.hpp's \c
/// detail::mode_is_external for the same identification), so this walk drives
/// directly off \p forest, per the Task-2 brief's ambiguity-resolution note.
/// \p sched is still consulted -- its shape gates dispatch (asserted below) --
/// so the entry point already has the scope-tree-walk SHAPE Tasks 3-4 extend
/// by recursing into \c sched.root.children; only the general value_id ->
/// node resolution is deferred to whichever later task plumbs \c
/// RichSchedule through this signature.
///
/// Reuses \c evaluate_impl (the single-node recursive engine, eval.hpp) for
/// every per-value build; no per-op evaluation is reimplemented here. The
/// per-root permute-to-\p layout and cross-root \c add_inplace accumulation
/// mirror \c sequant::evaluate(Node const&, layout, ...) and \c
/// sequant::evaluate(Nodes const&, layout, ...) respectively, inlined (rather
/// than delegated to those overloads) so that Tasks 3-4 have a natural seam
/// to inject batch-context management (slicing on entry, accumulate/scatter
/// on exit) around each per-value \c evaluate_impl call. The trace/hwmark
/// bookkeeping (\c log::term boundaries, the Permute \c EvalStat, the
/// SumInplace \c EvalStat) is reproduced line-for-line from those two
/// overloads too -- \c EvalTrace-gated, so it costs nothing when off, but
/// present so Tasks 3-4's batch-block trace/hwmark plumbing has a faithful
/// root-only baseline to extend rather than a silent gap to rediscover.
///
/// \param forest The forest whose per-root results are evaluated and summed;
///        same requirement as \c sequant::evaluate(Nodes const&, ...).
/// \param sched The scope-tree schedule (\c build_scope_schedule, Task 1)
///        for \p forest. Task 2 requires \c sched.root.children to be empty
///        (no realized batch loop): a non-empty scope tree needs the Task 3+
///        walk this function does not yet implement.
/// \param layout The layout each root's result is permuted to before being
///        summed; same meaning as \c sequant::evaluate(Node const&, layout,
///        ...)'s \p layout.
/// \param leaf_evaluator The leaf evaluator, as in \c sequant::evaluate.
/// \param cache The cache for common sub-expression elimination, as in \c
///        sequant::evaluate.
/// \param target_batch_size Per-index batch partition size (elements), the
///        source of the K partition for a realized batch loop; unused for the
///        root-only case. Same meaning as \c make_batched_custom_evaluator's
///        \p target_batch_size.
/// \param make_scope_guard Backend scope-guard factory, called with the batch
///        count and its RAII result HELD for the entire K-block loop -- exactly
///        as \c make_batched_custom_evaluator's contracted path does
///        (eval.hpp:2157). A screening backend uses it to relax block-sparse
///        screening scaled by the batch count so a contribution significant
///        over the FULL K range is not dropped within an individual batch (see
///        \c no_scope_guard). Defaults to \c make_no_scope_guard (a no-op for
///        the dense / DryRun backends), keeping every existing caller
///        byte-identical.
/// \return The summed, per-root-permuted result, as \c sequant::evaluate(
///         Nodes const&, layout, ...) would produce for the same \p forest.
///
/// \param rich The \c RichSchedule that produced \p sched (\c
///        compute_dag_boulevard). Used to resolve a \c ScopeNode's \c
///        homed_values (\c ValueCell::value_id's) back to forest nodes -- via
///        each cell's \c hash and \c build_value_node_map -- so the nested walk
///        (Task 4) can build a value HOMED at an outer loop once per outer
///        block (design integration point 3: drive off the type-keyed \c
///        homed_values, not a fragile exact-\c Index membership test).
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = ::sequant::make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate_whole_scope(
    Nodes const& forest, ScopeSchedule const& sched, RichSchedule const& rich,
    auto const& layout, F const& leaf_evaluator, CacheManager<N, FHC>& cache,
    std::function<std::size_t(Index const&)> target_batch_size = {},
    ScopeGuardFactory make_scope_guard = {}) {
  using node_t = std::ranges::range_value_t<Nodes>;
  static_assert(std::is_same_v<node_t, N>,
                "the forest's node type and the cache's node type must match");

  // Install a schedule-derived per-op trace annotation (remat home + use
  // scopes) for the duration of this replay, restored on ANY exit. Purely a
  // logging affordance -- gated on printing() so it is a no-op (byte-identical
  // trace, no map build) when tracing is off; the guard restores regardless.
  struct NodeMetaGuard {
    std::function<std::string(std::size_t)> prev;
    ~NodeMetaGuard() { Logger::instance().eval.node_meta = std::move(prev); }
  } node_meta_guard{std::move(Logger::instance().eval.node_meta)};
  if (log::printing()) Logger::instance().eval.node_meta = make_node_meta(rich);

  // The per-root, UNPERMUTED result, in forest order. The root-only case fills
  // each entry with a direct evaluate_impl; a batched scope tree fills the
  // loop-carrying roots from detail::walk_scope (accumulate/scatter per level)
  // and the invariant roots directly. The final permute-to-layout + cross-root
  // add_inplace combine (which reproduces sequant::evaluate(Nodes const&,
  // layout, ...)'s trace/hwmark bookkeeping line-for-line) is SHARED by both
  // cases, so every case emits identical Term/SumInplace records.
  container::svector<node_t> roots;
  for (auto&& n : forest) roots.push_back(n);
  container::svector<ResultPtr> pre_results(roots.size());

  if (sched.root.children.empty()) {
    // -------- Task 2 top-scope walk: every root built directly. --------
    for (std::size_t i = 0; i != roots.size(); ++i)
      pre_results[i] =
          evaluate_impl<EvalTrace>(roots[i], leaf_evaluator, cache);
  } else {
    // -------- Batched scope tree: recurse through sched.root.children.
    // --------
    //
    // Each child of the root is one realized batch loop; a Contracted loop
    // accumulates its per-block partials, an External loop scatters them into a
    // pre-sized result, and a Contracted loop with a child nests an inner loop
    // (aux OUTER over occ INNER) building its outer-homed values once per outer
    // block. All of that lives in detail::walk_scope; this driver only splits
    // the roots that participate in the loop nest from those invariant to it
    // and runs the shared combine below.
    SEQUANT_ASSERT(target_batch_size &&
                   "evaluate_whole_scope with a realized batch loop needs a "
                   "target_batch_size (the batch partition source).");

    // The value_id -> forest-node bridge (via ValueCell::hash) used to resolve
    // each ScopeNode's homed_values to nodes for the build-once hoist, plus the
    // forest-root hashes so a root is never pre-built as a homed value.
    auto const vmap = build_value_node_map(forest);
    std::unordered_set<std::size_t> root_hashes;
    for (auto const& r : roots) root_hashes.insert(r->hash_value());

    ScopeNode const& top = sched.root.children.front();
    std::wstring const top_base(top.mode.space().base_key());

    // Members = roots participating in the OUTERMOST loop (they carry its mode
    // TYPE -- structurally for a contracted mode, or free/annotated for an
    // external one). Every other root is invariant to the whole nest and built
    // directly, exactly like the top-scope path. Membership is TYPE-keyed (not
    // an exact-Index find_leaf_carrying) so a root that binds the mode type
    // under a different physical ordinal is not misclassified.
    container::svector<node_t const*> member_ptrs;
    container::svector<std::size_t> member_idx;
    for (std::size_t i = 0; i != roots.size(); ++i) {
      bool const member =
          top.kind == BatchModeType::Contracted
              ? detail::find_leaf_carrying_type(roots[i], top_base).has_value()
              : (detail::member_external_axis(roots[i], top_base).has_value() ||
                 detail::result_position_type(roots[i], top_base).has_value());
      if (member) {
        member_ptrs.push_back(&roots[i]);
        member_idx.push_back(i);
      } else {
        pre_results[i] =
            evaluate_impl<EvalTrace>(roots[i], leaf_evaluator, cache);
      }
    }

    if (!member_ptrs.empty()) {
      typename CacheManager<N, FHC>::BatchContext const empty_ctx;
      auto res =
          detail::walk_scope<EvalTrace, node_t, F, FHC, ScopeGuardFactory>(
              top, member_ptrs, empty_ctx, cache, leaf_evaluator,
              target_batch_size, make_scope_guard, vmap, rich, root_hashes);
      SEQUANT_ASSERT(res.size() == member_ptrs.size());
      for (std::size_t k = 0; k != member_ptrs.size(); ++k)
        pre_results[member_idx[k]] = std::move(res[k]);
    }
  }

  // -------- Shared combine: permute each root to layout and sum. --------
  // Factored into combine_forest_roots (forest_combine.hpp) so
  // evaluate_ordered_schedule (ordered_executor.hpp) can reuse the identical
  // bookkeeping rather than a hand-synced second copy; see that function's
  // doc comment.
  return combine_forest_roots<EvalTrace>(roots, pre_results, layout, cache);
}

}  // namespace sequant::eval

namespace sequant {

///
/// \brief Task 6 of the whole-scope batched DAG execution design (see
/// `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md`):
/// the COEXISTENCE driver entry. Routes a whole-forest evaluation to either
/// today's per-tree forest descent (\c sequant::evaluate(Nodes const&,
/// layout, leaf_evaluator, cache), eval.hpp:1105 -- called UNCHANGED) or the
/// whole-scope executor (\c eval::evaluate_whole_scope, Tasks 2-5), selected
/// by \p policy.scheduler.
///
/// \details \p policy.scheduler == BatchScheduler::forest_descent (the
/// default) is an unconditional, first-statement forward to the existing \c
/// sequant::evaluate(Nodes const&, layout, leaf_evaluator, cache) overload --
/// no schedule is built, nothing else runs on this path. Because this is an
/// ADDITIVE new overload (distinguished by the extra \p policy argument) and
/// the pre-existing overload is not modified, every caller of that overload
/// today stays byte-identical; a caller opts into this driver only by
/// supplying a \c BatchPolicy explicitly.
///
/// \p policy.scheduler == BatchScheduler::whole_scope builds the \c
/// eval::RichSchedule (\c eval::compute_dag_boulevard)
/// and the narrow \c eval::ScopeSchedule (\c eval::build_scope_schedule) from
/// \p forest's OWN placement -- the \c batched_here() annotations a prior
/// factorizer pass (e.g. \c optimize() with a batched objective, or a test
/// that stamps them directly) already recorded on \p forest's nodes -- then
/// drives \c eval::evaluate_whole_scope with \p policy.batch_target_size as
/// the batch-partition source and \p make_scope_guard forwarded verbatim. The
/// throwaway \c eval::dryrun::CostModel / \c SizeRegime built here are inert
/// for \c compute_dag_boulevard (its \c cm parameter is unused, kept only for
/// signature symmetry -- see its doc comment); they carry none of \p policy's
/// or the real backend's sizing information.
///
/// \param mode_order Ranks the scope-tree's canonical chain order (outermost
///        first; see \c eval::build_scope_schedule). Only consulted when the
///        flag is on. Empty (default) falls back to alphabetical order per
///        index type: per \c build_scope_schedule's own doc, chain order does
///        not affect NUMERIC correctness (only which axis nests outer vs.
///        inner -- a performance, not correctness, concern), so the default
///        is safe; a caller that cares about nesting order (e.g. aux outer /
///        occ inner, for peak) passes it explicitly.
/// \param make_scope_guard Backend scope-guard factory forwarded to \c
///        evaluate_whole_scope verbatim when the flag is on; unused (never
///        constructed or invoked) when it is off.
///
template <Trace EvalTrace = Trace::Default, meta::can_evaluate_range Nodes,
          typename F, typename N, bool FHC,
          typename ScopeGuardFactory = make_no_scope_guard>
  requires meta::leaf_node_evaluator<std::ranges::range_value_t<Nodes>, F>
ResultPtr evaluate(Nodes const& forest, BatchPolicy const& policy,
                   auto const& layout, F const& leaf_evaluator,
                   CacheManager<N, FHC>& cache,
                   std::initializer_list<std::wstring> mode_order = {},
                   ScopeGuardFactory make_scope_guard = {}) {
  // BatchPolicy docs: an empty batch_target_size means "no batching"; guard
  // the same way make_evaluator(policy, ...) does rather than let an empty
  // std::function throw std::bad_function_call out of compute_dag_boulevard /
  // evaluate_whole_scope / evaluate_ordered_schedule.
  auto const target_of = [&policy]() {
    return policy.batch_target_size
               ? policy.batch_target_size
               : std::function<std::size_t(Index const&)>(
                     [](Index const&) -> std::size_t { return 1; });
  };

  // SP3 gating switch: the ORDERED executor, driven by the SP2
  // OrderedSchedule IR (see BatchPolicy::scheduler's doc comment). Checked
  // FIRST so the two pre-existing dispatch arms below are reached,
  // byte-identically, only when the scheduler is NOT BatchScheduler::ordered.
  if (policy.scheduler == BatchScheduler::ordered) {
    std::function<std::size_t(Index const&)> const target = target_of();

    eval::dryrun::SizeRegime const regime;
    eval::dryrun::CostModel const cm{regime};
    eval::RichSchedule const rich =
        eval::compute_dag_boulevard(forest, cm, target);
    eval::LegalitySchedule const legality =
        eval::analyze_legality(rich, forest, policy);
    eval::OrderedSchedule const ordered =
        eval::build_ordered_schedule(rich, legality, policy, mode_order);

    // NODE-level lift of policy.is_volatile_leaf, exactly as make_evaluator's
    // own is_volatile_node lift (eval.hpp) computes it -- threaded into the
    // ordered executor for a later task's home-value volatile-vs-persistent
    // classification; not yet consulted there.
    using node_t = std::ranges::range_value_t<Nodes>;
    auto is_volatile_node =
        [p = policy.is_volatile_leaf](node_t const& n) -> bool {
      if (!n.leaf() || !n->is_tensor()) return false;
      return p && p(n->as_tensor());
    };

    return eval::evaluate_ordered_schedule<EvalTrace>(
        forest, ordered, rich, layout, leaf_evaluator, cache, target,
        make_scope_guard, is_volatile_node);
  }

  if (policy.scheduler != BatchScheduler::whole_scope)
    return evaluate<EvalTrace>(forest, layout, leaf_evaluator, cache);

  std::function<std::size_t(Index const&)> const target = target_of();

  eval::dryrun::SizeRegime const regime;
  eval::dryrun::CostModel const cm{regime};
  eval::RichSchedule const rich =
      eval::compute_dag_boulevard(forest, cm, target);
  eval::ScopeSchedule const sched =
      eval::build_scope_schedule(rich, mode_order);

  return eval::evaluate_whole_scope<EvalTrace>(forest, sched, rich, layout,
                                               leaf_evaluator, cache, target,
                                               make_scope_guard);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_SCOPE_EXECUTOR_HPP
