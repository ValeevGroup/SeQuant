#ifndef SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
#define SEQUANT_EVAL_ORDERED_SCHEDULE_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstdlib>
#include <iostream>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>

namespace sequant::eval {

///
/// \brief Task 2 of the ordered-scope batched-eval design (SP2): the \c
/// OrderedSchedule IR -- an ORDERED tree of loop blocks and build steps.
/// Purely data + a structural \c well_formed check; no sequencer/executor
/// here (that is a later task).
///
/// \details Unlike \c ScopeSchedule (\c scope_schedule.hpp), whose \c
/// ScopeNode groups all values homed at a node into one unordered \c
/// homed_values bag, \c OrderedSchedule threads builds and child loop blocks
/// through a single ORDERED sequence (\c ScopeBlock::steps): a value built
/// AFTER a child block in that sequence may read the child block's
/// accumulated output, so relative order among steps is load-bearing, not
/// incidental.
///

///
/// \brief What happens to a value at the close of its home \c ScopeBlock.
///
enum class OutputKind {
  Transient,      //!< the value's home IS this block; nothing carries out.
                  //!< NOTE (SP3 readers): a Transient value NEVER appears in
                  //!< any block's \c outputs list -- it is realized purely
                  //!< as a plain \c BuildStep, and an explicit Transient \c
                  //!< outputs entry would double-produce it (violating \c
                  //!< well_formed's single-producer check). Do NOT scan \c
                  //!< outputs for Transient; it is the ABSENCE of an escape
                  //!< output, not a recorded one.
  AccumulateSum,  //!< reduction: summed into an outer accumulator
  AccumulateScatter,  //!< loop-carried: scattered into a disjoint outer slice
};

///
/// \brief Build this value here: one contraction of already-cached operands.
///
struct BuildStep {
  std::size_t value_id;
};

struct ScopeBlock;  // fwd; see Step's doc comment (below ScopeBlock) for why
                    // both types are forward-declared here.
struct Step;        // fwd; see immediately below.

///
/// \brief One loop block: a batch loop over \c axis (sentinel/default on the
/// root block, which sits outside every loop), containing an ORDERED
/// sequence of build steps and nested child blocks, plus the set of values
/// that leave this block (and how) when it closes.
///
struct ScopeBlock {
  Index axis{};  //!< the loop axis; default (sentinel) on the root block.
  int latitude_ordinal =
      0;                  //!< layout: PROCON pass index (was: ordinal) --
                          //!< disambiguates recurring sibling blocks
                          //!< realizing the same axis (e.g. producer/consumer
                          //!< passes over the same axis TYPE at one level).
  DagScopeLevel level{};  //!< this block's DAG-scope nest position (mirrors
                          //!< \c axis/\c ordinal; see \c DagScopeLevel's doc
                          //!< comment). Default-valued (depth 0, empty
                          //!< space, ordinal 0) on the root block.
  BatchModeType kind =
      BatchModeType::Contracted;    //!< Contracted (accumulate on block exit)
                                    //!< or External (scatter on block exit);
                                    //!< meaningless on the root.
  container::vector<Step> steps{};  //!< ORDERED: build-or-child-block,
                                    //!< interleaved (see \c Step's doc
                                    //!< comment for why \c container::vector,
                                    //!< not \c svector).
  container::svector<std::pair<std::size_t, OutputKind>>
      outputs{};  //!< value_id -> how it leaves this block on close.
};

///
/// \brief One ORDERED step of a \c ScopeBlock: build a value at this scope,
/// or enter a nested child loop block.
///
/// \details The design brief's target shape is a bare alias, \c using Step =
/// std::variant<BuildStep, ScopeBlock>. That is not directly expressible:
/// \c ScopeBlock::steps must hold a sequence of \c Step (so a build
/// interleaves with child blocks in one ORDERED list, per the class doc
/// above), which means \c Step has to be at least forward-declared before \c
/// ScopeBlock; but \c std::variant requires every alternative type COMPLETE
/// at the point the variant specialization is instantiated (unlike \c
/// container::vector / \c std::vector, which the C++17 library tolerates
/// holding an incomplete element type up to first use -- the same relaxation
/// \c ScopeNode::children in \c scope_schedule.hpp relies on for ITS
/// self-reference). A bare \c using Step = std::variant<BuildStep,
/// ScopeBlock> declared before \c ScopeBlock therefore cannot compile
/// (\c ScopeBlock incomplete there), and a type alias cannot be
/// forward-declared separately from its definition (no "using Step;"
/// forward declaration exists in C++) to defer it to after \c ScopeBlock
/// either.
///
/// The fix: make \c Step a real (forward-declarable) class wrapping the
/// variant, not a bare alias -- \c ScopeBlock::steps holds \c
/// container::vector<Step> while \c Step is still only forward-declared
/// (legal, per the \c std::vector incomplete-type allowance above), and \c
/// Step's own definition (with the variant member) follows \c ScopeBlock,
/// where \c ScopeBlock is by then complete. This preserves the single
/// ORDERED sequence of interleaved build/child-block steps the design
/// requires, and preserves real \c std::variant semantics (\c
/// std::holds_alternative / \c std::get_if / \c std::visit all work on \c
/// Step::value) -- the only change from the brief's literal shape is the one
/// extra wrapping layer forced by the forward-declaration ordering.
///
struct Step {
  std::variant<BuildStep, ScopeBlock> value;

  Step(BuildStep b) : value(std::move(b)) {}   // NOLINT(*-explicit-*)
  Step(ScopeBlock b) : value(std::move(b)) {}  // NOLINT(*-explicit-*)
};

///
/// \brief The whole ordered schedule: the root block plus the total value
/// count (every \c BuildStep::value_id and every \c ScopeBlock::outputs
/// value_id is expected to be < \c num_values; see \c well_formed).
///
struct OrderedSchedule {
  ScopeBlock root{};
  std::size_t num_values = 0;
  /// Pillar 1 (slice-colored value identity): per value_id, the value's
  /// home-sliced modes paired with the DAG-scope DEPTH at which each is sliced
  /// -- recorded ONCE here, at home placement, in the value's OWN index-frame
  /// (keys are `rich.cells[vid].carried` indices, never a foreign frame or
  /// base_key). This is the authoritative mode->depth source the value-id
  /// coloring consumes; it is NOT re-derived from ectx. A value with no
  /// home-sliced modes has no entry (unsliced => empty coloring => value-id ==
  /// node-id). Distinct from Pillar-2 `occ_facts` (use-scope, per consumer).
  std::unordered_map<std::size_t, container::svector<std::pair<Index, int>>>
      home_mode_depth{};
  /// Pillar 1 / B-full: per value_id, the value_ids of its DIRECT operands --
  /// the value/occurrence DAG edges the value-driven ordered executor consumes
  /// to fetch each operand by its OWN home-colored key. Recorded here from
  /// `ordered_schedule_dep_graph(rich).depends_on`, whose edges come from every
  /// `OccurrenceRec`'s `consumer_point` (so split operands resolve to the
  /// specific consumed value, not an ambiguous node hash). A leaf value (no
  /// operands) has no entry.
  std::unordered_map<std::size_t, container::svector<std::size_t>>
      operand_vids{};
};

namespace detail {

///
/// \brief \c well_formed's recursive worker: checks \p block's own steps and
/// outputs, recurses into every child \c ScopeBlock step, and checks ordinal
/// uniqueness among \p block's own same-axis (\c IndexSpace::base_key())
/// sibling child blocks.
///
inline bool ordered_schedule_block_well_formed(ScopeBlock const& block,
                                               std::size_t num_values) {
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      if (build->value_id >= num_values) return false;
    } else {
      auto const& child = std::get<ScopeBlock>(step.value);
      if (!ordered_schedule_block_well_formed(child, num_values)) return false;
    }
  }

  // Ordinal uniqueness among same-axis sibling child blocks (direct children
  // of THIS block only; deeper levels are checked by the recursion above).
  for (std::size_t i = 0; i < block.steps.size(); ++i) {
    auto const* ci = std::get_if<ScopeBlock>(&block.steps[i].value);
    if (!ci) continue;
    for (std::size_t j = i + 1; j < block.steps.size(); ++j) {
      auto const* cj = std::get_if<ScopeBlock>(&block.steps[j].value);
      if (!cj) continue;
      // Two sibling blocks are the SAME realized loop only when their FULL
      // loop IDENTITY collides: (depth, loop_slot) AND the latitude (PROCON
      // pass). Keying on the axis SPACE (a fusion color, not identity) wrongly
      // rejected two DISTINCT same-space sibling loops the un-fuse legitimately
      // emits at different (depth, loop_slot) -- the w20 aux+occ case: two occ
      // (space "i") nests at (1,0) and (2,1), same latitude 0, are different
      // loops, not a duplicate. (See LoopKey::color / the same
      // space-vs-identity correction in the home-scope coloring.)
      if (ci->level.depth == cj->level.depth &&
          ci->level.loop_slot == cj->level.loop_slot &&
          ci->latitude_ordinal == cj->latitude_ordinal) {
        if (std::getenv("SEQUANT_DUMP_WF"))
          std::cerr << "[wf-fail] sibling-identity-collision space="
                    << toUtf8(std::wstring(ci->axis.space().base_key()))
                    << " depth=" << ci->level.depth
                    << " slot=" << ci->level.loop_slot
                    << " lat=" << ci->latitude_ordinal << std::endl;
        return false;
      }
    }
  }

  for (auto const& [value_id, kind] : block.outputs) {
    (void)kind;
    if (value_id >= num_values) return false;
  }

  return true;
}

///
/// \brief Append every value_id \p block PRODUCES -- as a \c BuildStep
/// (recursively, through every nested child block) or as a value_id in
/// \p block's own \c outputs -- to \p out.
///
/// \details Feeds the whole-schedule single-producer check in \c
/// well_formed: a \c BuildStep and a block \c outputs entry are both
/// "production sites" for a value_id in the SSA-like sense the schedule is
/// meant to hold, so both contribute to the same collected list.
///
inline void collect_production_ids(ScopeBlock const& block,
                                   container::vector<std::size_t>& out) {
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      out.push_back(build->value_id);
    } else {
      auto const& child = std::get<ScopeBlock>(step.value);
      collect_production_ids(child, out);
    }
  }
  for (auto const& [value_id, kind] : block.outputs) {
    (void)kind;
    out.push_back(value_id);
  }
}

///
/// \brief A single escape (output) site: a value_id and the root-to-block PATH
/// at which it escapes. Feeds \c well_formed's multi-level escape-chain check
/// (a value carried on an outer axis AND reduced on an inner one escapes at
/// BOTH -- see \c build_ordered_schedule's escape emission).
///
struct OutputSite {
  std::size_t value_id;
  container::svector<ScopeBlock const*>
      path;  //!< root .. this block, inclusive
};

///
/// \brief Collect every \c BuildStep value_id into \p builds and every output
/// escape site (value_id + its root-to-block path) into \p sites.
///
inline void collect_productions(ScopeBlock const& block,
                                container::svector<ScopeBlock const*>& path,
                                container::vector<std::size_t>& builds,
                                container::vector<OutputSite>& sites) {
  path.push_back(&block);
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value))
      builds.push_back(b->value_id);
    else
      collect_productions(std::get<ScopeBlock>(step.value), path, builds,
                          sites);
  }
  for (auto const& [value_id, kind] : block.outputs) {
    (void)kind;
    sites.push_back(OutputSite{value_id, path});
  }
  path.pop_back();
}

}  // namespace detail

///
/// \brief Structural sanity check on \p sched (no sequencer logic):
///   - every \c BuildStep::value_id is < \c sched.num_values;
///   - ordinals are unique among same-axis (\c IndexSpace::base_key())
///     sibling blocks within a parent;
///   - every \c ScopeBlock::outputs value_id is < \c sched.num_values;
///   - SINGLE-PRODUCER (SSA-like), with the multi-level escape chain allowed:
///     no value_id is built (\c BuildStep) more than once; no value_id is both
///     built AND escaped; and a value_id may escape (\c outputs) at MORE than
///     one block ONLY when those blocks lie on a single root-to-node nesting
///     path (distinct depths, each shallower one an ancestor of the deepest) --
///     the inner-sum / outer-scatter escape chain of \c build_ordered_schedule.
///     Escapes in unrelated (sibling) blocks, or two escapes at one depth, are
///     rejected as duplicate production.
///
/// \note This checks single-producer (no DUPLICATE production) but NOT
/// completeness (no value_id GAPS -- that every id in `[0, num_values)` is
/// produced somewhere). Completeness holds by construction of \c
/// build_ordered_schedule and is asserted in the Task-5 acceptance test, so an
/// SP3 reader must NOT assume \c well_formed implies every value_id is present.
///
[[nodiscard]] inline bool well_formed(OrderedSchedule const& sched) {
  if (!detail::ordered_schedule_block_well_formed(sched.root, sched.num_values))
    return false;

  container::svector<ScopeBlock const*> path;
  container::vector<std::size_t> builds;
  container::vector<detail::OutputSite> sites;
  detail::collect_productions(sched.root, path, builds, sites);

  // (a) no value_id built more than once.
  {
    container::vector<std::size_t> b = builds;
    std::sort(b.begin(), b.end());
    if (auto const it = std::adjacent_find(b.begin(), b.end()); it != b.end()) {
      if (std::getenv("SEQUANT_DUMP_WF"))
        std::cerr << "[wf-fail] double-build vid=" << *it << std::endl;
      return false;
    }
  }
  // (b) no value_id both built and escaped.
  {
    std::unordered_set<std::size_t> const build_set(builds.begin(),
                                                    builds.end());
    for (auto const& s : sites)
      if (build_set.count(s.value_id)) {
        if (std::getenv("SEQUANT_DUMP_WF"))
          std::cerr << "[wf-fail] built-and-escaped vid=" << s.value_id
                    << std::endl;
        return false;
      }
  }
  // (c) a value_id's escape sites (>1 => a multi-level chain) must lie on ONE
  // root-to-node nesting path: distinct depths, and every shorter path a
  // prefix of the deepest (so each is an ancestor of the next).
  {
    std::unordered_map<std::size_t,
                       container::svector<detail::OutputSite const*>>
        by_vid;
    for (auto const& s : sites) by_vid[s.value_id].push_back(&s);
    for (auto const& [vid, ss] : by_vid) {
      (void)vid;
      if (ss.size() == 1) continue;
      detail::OutputSite const* deepest = ss.front();
      for (auto const* s : ss)
        if (s->path.size() > deepest->path.size()) deepest = s;
      std::unordered_set<std::size_t> depths;
      bool const wfdbg = std::getenv("SEQUANT_DUMP_WF") != nullptr;
      for (auto const* s : ss) {
        if (!depths.insert(s->path.size()).second) {  // same depth
          if (wfdbg)
            std::cerr << "[wf-fail] escape-chain same-depth vid=" << vid
                      << " depth=" << s->path.size() << " nsites=" << ss.size()
                      << std::endl;
          return false;
        }
        if (s->path.size() > deepest->path.size()) return false;
        for (std::size_t k = 0; k < s->path.size(); ++k)
          if (s->path[k] != deepest->path[k]) {  // not an ancestor
            if (wfdbg)
              std::cerr << "[wf-fail] escape-chain not-ancestor vid=" << vid
                        << " depth=" << s->path.size()
                        << " deepest=" << deepest->path.size()
                        << " diverge_at=" << k << std::endl;
            return false;
          }
      }
    }
  }
  return true;
}

namespace detail {

///
/// \brief Per-axis-depth accumulator used while \c build_ordered_schedule
/// walks the canonical chain: which value_id's are plain \c BuildStep's at
/// this depth, and which are escape \c outputs (see the function's own doc
/// comment for what "escape" means).
///
struct OrderedScheduleDepthBucket {
  container::svector<std::size_t> build_ids;
  container::svector<std::pair<std::size_t, OutputKind>> outputs;
};

///
/// \brief Per-candidate-step metadata for \c ordered_schedule_topo_sort_steps:
/// which value_id's this step (a \c BuildStep or a nested child \c
/// ScopeBlock, already built) makes visible to ITS OWN siblings at this
/// SAME block level (\c produced), which value_id's its content directly
/// needs (\c requires_, unfiltered -- see \c ordered_schedule_topo_sort_steps
/// for how the irrelevant/external entries are dropped), and a tie-break key
/// for when the true dependency order leaves two ready steps unordered.
///
struct OrderedScheduleStepMeta {
  container::svector<std::size_t> produced;
  container::svector<std::size_t> requires_;
  std::size_t tie_key = 0;
};

///
/// \brief Topologically sort \p items (one already-built \c Step per entry,
/// paired index-for-index with \p meta) by the LOCAL dependency edges among
/// THIS block's own steps: step A must precede step B whenever B's \c
/// requires_ names a value_id that's in A's \c produced. Kahn's algorithm;
/// among simultaneously-ready steps, always picks the smallest \c tie_key
/// first, for a deterministic result when the true dependency order leaves
/// steps genuinely unordered relative to each other.
///
/// \details Replaces an earlier (unsound) scalar-key-only sort: a single
/// scalar per step can be MADE to sort a child block before every value that
/// reads its output (see \c build_ordered_schedule's own doc comment, part
/// 3), but it cannot ALSO guarantee a step sorts after every value its own
/// content reads as an input -- those are two independent constraints a
/// single total order can satisfy only when they happen to agree, which
/// water-20's own aux-only fixture never stresses (its \c {Κ} block's
/// content is leaf-only, needing no root-homed input). A real topological
/// sort over the ACTUAL per-step dependency edges satisfies both directions
/// by construction, superseding the scalar key -- the key survives only as
/// the tie-break \c build_ordered_schedule still needs for the (usual)
/// case of two steps with no dependency relation to each other at all.
///
/// \c SEQUANT_ASSERT's that every item is placed exactly once (a cycle in
/// this LOCAL edge set would be a bug -- these edges are a sub-relation of
/// the whole-forest DAG's edges, restricted to one block's own siblings, so
/// they inherit its acyclicity), then re-derives the local edges a SECOND
/// time against the FINAL order and \c SEQUANT_ASSERT's every one is
/// actually satisfied (a loud tripwire against any future violation of this
/// invariant, per the design review that requested it, rather than a silent
/// mis-order).
///
inline container::vector<Step> ordered_schedule_topo_sort_steps(
    container::vector<Step> items,
    container::vector<OrderedScheduleStepMeta> const& meta) {
  std::size_t const m = items.size();
  SEQUANT_ASSERT(meta.size() == m);

  // value_id -> which LOCAL item produces it. well_formed's whole-schedule
  // single-producer invariant guarantees at most one item at ANY level can
  // claim a given value_id; a value_id absent here is external to this
  // level (resolved at an ancestor level, not a local ordering constraint).
  std::unordered_map<std::size_t, std::size_t> produced_by;
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t vid : meta[i].produced) produced_by.emplace(vid, i);

  container::vector<container::svector<std::size_t>> prerequisites(m);
  container::vector<std::size_t> indegree(m, 0);
  container::vector<container::svector<std::size_t>> dependents(m);
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t vid : meta[i].requires_) {
      auto const it = produced_by.find(vid);
      if (it == produced_by.end() || it->second == i) continue;
      auto& preqs = prerequisites[i];
      if (std::find(preqs.begin(), preqs.end(), it->second) == preqs.end()) {
        preqs.push_back(it->second);
        dependents[it->second].push_back(i);
      }
    }
    indegree[i] = prerequisites[i].size();
  }

  container::svector<std::size_t> ready;
  for (std::size_t i = 0; i < m; ++i)
    if (indegree[i] == 0) ready.push_back(i);

  container::svector<std::size_t> order;
  order.reserve(m);
  while (!ready.empty()) {
    auto const best_it = std::min_element(
        ready.begin(), ready.end(), [&](std::size_t a, std::size_t b) {
          if (meta[a].tie_key != meta[b].tie_key)
            return meta[a].tie_key < meta[b].tie_key;
          return a < b;  // full determinism on an exact tie
        });
    std::size_t const cur = *best_it;
    ready.erase(best_it);
    order.push_back(cur);
    for (std::size_t dep : dependents[cur]) {
      SEQUANT_ASSERT(indegree[dep] > 0);
      if (--indegree[dep] == 0) ready.push_back(dep);
    }
  }
  // No cycle: see the function doc comment.
  SEQUANT_ASSERT(order.size() == m);

  // Post-sort validation (loud tripwire, see the function doc comment):
  // every local prerequisite must actually precede its dependent.
  container::vector<std::size_t> position(m);
  for (std::size_t pos = 0; pos < order.size(); ++pos)
    position[order[pos]] = pos;
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t p : prerequisites[i])
      SEQUANT_ASSERT(position[p] < position[i]);

  container::vector<Step> out_steps;
  out_steps.reserve(m);
  for (std::size_t idx : order) out_steps.push_back(std::move(items[idx]));
  return out_steps;
}

///
/// \brief The GLOBAL direct-dependency edges of a \c RichSchedule, recovered
/// from \c rich alone (no forest access): \c depends_on[p] lists every value_id
/// the value \c p directly READS, and \c consumers_of[c] lists every value_id
/// that directly reads \c c (the reverse). Same recovery \c
/// build_ordered_schedule uses inline (occurrence \c consumer_point ->
/// producing value_id via \c point_owner); factored here so the split-pass
/// partition and \c build_ordered_schedule agree edge-for-edge.
///
struct OrderedScheduleDepGraph {
  std::unordered_map<std::size_t, std::size_t> value_id_of;  //!< hash -> id
  std::unordered_map<std::size_t, container::svector<std::size_t>> depends_on;
  std::unordered_map<std::size_t, container::svector<std::size_t>> consumers_of;
};

inline OrderedScheduleDepGraph ordered_schedule_dep_graph(
    RichSchedule const& rich) {
  OrderedScheduleDepGraph g;
  g.value_id_of.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells)
    g.value_id_of.emplace(vc.hash, vc.value_id);

  std::unordered_map<std::size_t, std::size_t> point_owner;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences)
      point_owner[occ.point] = vc.value_id;

  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences) {
      if (occ.consumer_point == occ.point) continue;  // forest root
      auto const it = point_owner.find(occ.consumer_point);
      if (it == point_owner.end()) continue;  // defensive
      std::size_t const parent = it->second;
      auto& deps = g.depends_on[parent];
      if (std::find(deps.begin(), deps.end(), vc.value_id) == deps.end()) {
        deps.push_back(vc.value_id);
        g.consumers_of[vc.value_id].push_back(parent);
      }
    }
  return g;
}

///
/// \brief READ-FROM-HOME home use-count: the exact number of times a homed
/// value's HOME cache entry is accessed during ordered-scope eval, computed
/// from the ORDERED schedule's realized scopes (NOT the boulevard's occurrence
/// ectx that \c weighted_use_count reads).
///
/// \details Under the single read-from-home access discipline (see \c
/// make_batched_scratch: a batch-invariant home-resident operand is read from
/// the parent chain each batch, never seeded), a consumer \c W reads operand
/// \c V's home ONCE per batch of every loop on the path (home(V), scope(W)] --
/// exactly \c weighted_use_count's per-batch arithmetic, but each consumer's
/// enclosing loops come from where \c build_ordered_schedule actually PLACES
/// it, not from \c ValueCell::occurrences (the boulevard's natural nesting,
/// which the ordered schedule does not preserve -- the mismatch that made the
/// occurrence-fed count undercount). So:
/// \verbatim
///   home_reads(V) = 1 (the build store)
///     + Σ over DAG edges (W -> V), WITH multiplicity (V counted twice if it
///       is BOTH operands of W): Π n_blocks(L) over loops L in scope(W) whose
///       TYPE is not in scope(V) (the path strictly inner to V's home).
/// \endverbatim
/// TYPE-keyed (\c IndexSpace::base_key()) like \c weighted_use_count, and
/// multiplicity-aware from the DAG (a value used as both operands of one
/// contraction is two home reads), which the deduped \c consumers_of is not.
/// The \c +1 is the homing store's own decaying access (\c CacheManager stores
/// by access), matching the \c "+1" the seeded homing used.
///
/// \brief A value's home-scope key, per-cell (see \c ordered_home_reads). One
/// entry per enclosing loop of the home, as \c (depth, loop_slot,
/// latitude_ordinal) -- IDENTITY (depth, loop_slot) plus the PROCON LAYOUT
/// coordinate (latitude), so two cells of one value that differ ONLY by pass
/// are distinct. The executor builds this from the scope at each home site.
using HomeScopeKey = container::svector<std::tuple<std::size_t, int, int>>;

template <typename node_t>
[[nodiscard]] inline std::function<std::size_t(std::size_t,
                                               HomeScopeKey const&)>
ordered_home_reads(OrderedSchedule const& ordered, RichSchedule const& rich,
                   std::unordered_map<std::size_t, node_t> const& vmap,
                   std::function<std::size_t(Index const&)> const& n_blocks) {
  // Two scopes per value, both keyed off the ordered schedule's realized
  // nesting:
  //  - build_scope[vid]: the loops enclosing where vid is EVALUATED (reads its
  //    operands). This is what charges a consumer's reads of its operand.
  //  - home_scope[vid]: the loops enclosing where vid's cache entry LIVES.
  // They differ for an escape output: it is BUILT inside its block (per batch,
  // reading its operands there -- build_scope includes the block axis) but its
  // result HOMES one level OUT (stored at the parent on close -- home_scope
  // excludes the block axis). A plain BuildStep has build_scope == home_scope.
  // Each enclosing loop is kept as (axis, LoopKey) -- the axis for its block
  // count, the LoopKey (depth, loop_slot) for IDENTITY. Identity matters
  // because the un-fuse realizes the same SPACE as several DISTINCT loops in
  // separate nests; a space-only "is this loop above V's home" test would
  // wrongly treat a sibling nest's i-loop as V's own and drop its block factor
  // -> under-count -> premature eviction.
  using ScopeEntry = std::pair<Index, LoopKey>;
  auto build_scope = std::make_shared<
      std::unordered_map<std::size_t, container::svector<ScopeEntry>>>();
  auto home_scope = std::make_shared<
      std::unordered_map<std::size_t, container::svector<ScopeEntry>>>();
  auto const walk = [&build_scope, &home_scope](
                        auto&& self, ScopeBlock const& b,
                        container::svector<ScopeEntry> const& enc) -> void {
    for (Step const& s : b.steps) {
      if (auto const* build = std::get_if<BuildStep>(&s.value)) {
        (*build_scope)[build->value_id] = enc;
        (*home_scope)[build->value_id] = enc;
      } else if (auto const* child = std::get_if<ScopeBlock>(&s.value)) {
        container::svector<ScopeEntry> inner = enc;
        inner.push_back({child->axis, child->level.key()});
        self(self, *child, inner);
      }
    }
    // An escape output is BUILT inside b (build_scope == enc, reading its
    // operands once per batch) but HOMES one level OUT (home_scope == enc minus
    // b's own axis). Root (sentinel axis, enc=={}) has no outputs.
    for (auto const& [vid, kind] : b.outputs) {
      (void)kind;
      // build_scope = DEEPEST escape: a value carried on more than one nested
      // axis escapes at EACH, but the outer escapes only scatter an
      // already-built value -- they read NO operands; the contraction runs at
      // the innermost escape. Keep the deepest (longest enc); a shallower outer
      // escape must not clobber it, else `reads` loses that loop's block factor
      // and the value is evicted from its home before its last use. The walk
      // visits deeper escapes first (child recursion precedes this loop), so
      // guard on enc length.
      auto const bit = build_scope->find(vid);
      if (bit == build_scope->end() || enc.size() > bit->second.size())
        (*build_scope)[vid] = enc;
      // home_scope = the full result lives one level OUT of the OUTERMOST
      // (shallowest) escape; the shallowest own-outputs write is last, so
      // unconditional assignment is correct here.
      container::svector<ScopeEntry> out_scope = enc;
      if (!out_scope.empty()) out_scope.pop_back();
      (*home_scope)[vid] = out_scope;
    }
  };
  walk(walk, ordered.root, {});

  std::unordered_map<std::size_t, std::size_t> hash_to_vid;
  hash_to_vid.reserve(rich.cells.size());
  for (ValueCell const& c : rich.cells) hash_to_vid.emplace(c.hash, c.value_id);

  auto const scope_of =
      [](std::shared_ptr<std::unordered_map<
             std::size_t, container::svector<ScopeEntry>>> const& m,
         std::size_t vid) -> container::svector<ScopeEntry> {
    auto it = m->find(vid);
    return it == m->end() ? container::svector<ScopeEntry>{} : it->second;
  };
  // Is loop L (by IDENTITY, its LoopKey) one of the loops enclosing V's home?
  // Matched on (depth, loop_slot), NOT space: a sibling nest's same-space loop
  // is a DISTINCT loop, so V really is re-read across it.
  auto const loop_in = [](container::svector<ScopeEntry> const& s,
                          LoopKey const& k) {
    for (auto const& [ax, kk] : s)
      if (kk == k) return true;
    return false;
  };

  // Every value is cached at its scope (no min-use floor, member roots cached),
  // so evaluate_impl always stops at a consumer's DIRECT operands: V's home is
  // read exactly once per DIRECT parent build. Sum over V's direct DAG parents
  // W (with multiplicity -- V as BOTH operands of W is two reads), charging
  // each by W's build count over the path (home(V), build_scope(W)].
  auto reads = std::make_shared<std::unordered_map<std::size_t, std::size_t>>();
  for (ValueCell const& wc : rich.cells) {
    auto const wit = vmap.find(wc.hash);
    if (wit == vmap.end() || wit->second.leaf()) continue;
    container::svector<ScopeEntry> const scopeW =
        scope_of(build_scope, wc.value_id);
    auto const add_edge = [&](node_t const& child) {
      auto const cvit = hash_to_vid.find(child->hash_value());
      if (cvit == hash_to_vid.end()) return;
      std::size_t const cvid = cvit->second;
      container::svector<ScopeEntry> const homeC = scope_of(home_scope, cvid);
      std::size_t prod = 1;
      for (auto const& [Lax, Lkey] : scopeW)
        if (!loop_in(homeC, Lkey)) prod *= n_blocks(Lax);
      (*reads)[cvid] += prod;
      if (char const* dh = std::getenv("SEQUANT_DUMP_HOMEREADS");
          dh &&
          (child->hash_value() % 100000u) == std::strtoul(dh, nullptr, 10)) {
        auto const bsC = scope_of(build_scope, cvid);
        std::wcerr << L"[cons] consumed=" << (child->hash_value() % 100000u)
                   << L" W_hash=" << (wc.hash % 100000u) << L" prod=" << prod
                   << L" scopeW={";
        for (auto const& [ax, k] : scopeW)
          std::wcerr << ax.space().base_key() << L"#d" << k.depth << L"s"
                     << k.loop_slot << L" ";
        std::wcerr << L"} 36953.build={";
        for (auto const& [ax, k] : bsC)
          std::wcerr << ax.space().base_key() << L"#d" << k.depth << L"s"
                     << k.loop_slot << L" ";
        std::wcerr << L"}\n";
      }
    };
    add_edge(wit->second.left());
    add_edge(wit->second.right());
  }

  if (std::getenv("SEQUANT_DUMP_HOMEREADS")) {
    for (ValueCell const& c : rich.cells) {
      auto const it = reads->find(c.value_id);
      std::size_t const cnt = (it == reads->end() ? 0 : it->second) + 1;
      auto const hs = scope_of(home_scope, c.value_id);
      auto const bs = scope_of(build_scope, c.value_id);
      std::wcerr << L"[homereads] hash=" << (c.hash % 100000u) << L" reads="
                 << cnt << L" home={";
      for (auto const& [ax, k] : hs)
        std::wcerr << ax.space().base_key() << L"#d" << k.depth << L"s"
                   << k.loop_slot << L" ";
      std::wcerr << L"} build={";
      for (auto const& [ax, k] : bs)
        std::wcerr << ax.space().base_key() << L"#d" << k.depth << L"s"
                   << k.loop_slot << L" ";
      std::wcerr << L"}\n";
    }
  }

  // --- Per-cell INTERMEDIATE escape homes (multi-level escape) ---
  // A value carried on >1 nested axis escapes at EACH level, minting one CELL
  // per level: the OUTERMOST (shallowest home) is the full value -- its reads
  // are exactly `reads[vid]` above (external DAG consumers charged to the
  // collapsed home). Each INNER (partial) cell, homed one level out of a DEEPER
  // escape, is read by the NEXT-shallower escape once per that escape's own
  // build batch -- a STRUCTURAL read that is NOT a DAG-parent edge, so it never
  // appears in `reads` and the collapse (which keeps only the outermost home)
  // drops it. Record each inner cell's own home-scope (WITH latitude, so PROCON
  // passes stay distinct) and its structural read count, keyed by
  // (vid, home-scope-signature), so the executor can home each escape level
  // with its OWN cell's life instead of the collapsed full count. Single-home
  // values and the outermost cell are untouched (fall through to `reads[vid]`).
  using LatEntry = std::tuple<std::size_t, int, int>;  // depth, slot, latitude
  using LatScope = container::svector<LatEntry>;
  auto const sig_of = [](LatScope const& s) -> std::string {
    std::string r;
    for (auto const& [d, sl, lat] : s)
      r += std::to_string(d) + ":" + std::to_string(sl) + ":" +
           std::to_string(lat) + ";";
    return r;
  };
  // All escape home-scopes (latitude-aware) per vid, plus each level's build.
  struct LatCell {
    LatScope home;
    LatScope build;
  };
  std::unordered_map<std::size_t, container::svector<LatCell>> lat_homes;
  std::function<void(ScopeBlock const&, LatScope const&)> lwalk =
      [&](ScopeBlock const& b, LatScope const& enc) -> void {
    for (Step const& s : b.steps)
      if (auto const* child = std::get_if<ScopeBlock>(&s.value)) {
        LatScope inner = enc;
        inner.push_back({child->level.depth, child->level.loop_slot,
                         child->level.latitude_ordinal});
        lwalk(*child, inner);
      }
    for (auto const& [vid, kind] : b.outputs) {
      (void)kind;
      LatScope home = enc;
      if (!home.empty()) home.pop_back();
      lat_homes[vid].push_back({std::move(home), enc});
    }
  };
  lwalk(ordered.root, {});

  // Intermediate structural reads keyed by "vid|home_sig".
  auto inter = std::make_shared<std::unordered_map<std::string, std::size_t>>();
  for (auto& [vid, cells] : lat_homes) {
    if (cells.size() < 2) continue;  // single-level escape: no intermediate
    // Sort by home depth (shallowest = outermost = full, first).
    std::sort(cells.begin(), cells.end(),
              [](LatCell const& a, LatCell const& b) {
                return a.home.size() < b.home.size();
              });
    // Adjacent (shallower produces & reads deeper): the shallower escape's
    // build re-reads the deeper partial across the loops in shallower.build not
    // in deeper.home. Charge that to the deeper (inner) cell.
    for (std::size_t j = 1; j < cells.size(); ++j) {
      LatCell const& deeper = cells[j];
      LatCell const& shallower = cells[j - 1];
      std::size_t prod = 1;
      for (auto const& [d, sl, lat] : shallower.build) {
        (void)lat;
        bool in_home = false;
        for (auto const& [hd, hsl, hlat] : deeper.home)
          if (hd == d && hsl == sl) {
            in_home = true;
            break;
          }
        if (!in_home) {
          // block factor for identity (d, sl): find its axis via the collapsed
          // build_scope (any entry of matching identity carries a usable axis).
          Index ax;
          bool found = false;
          for (auto const& [wvid, sc] : *build_scope)
            for (auto const& [a, k] : sc)
              if (k.depth == d && k.loop_slot == sl) {
                ax = a;
                found = true;
                break;
              }
          prod *= found ? n_blocks(ax) : 1;
        }
      }
      (*inter)[std::to_string(vid) + "|" + sig_of(deeper.home)] += prod;
    }
  }

  if (std::getenv("SEQUANT_DUMP_HOMEREADS")) {
    for (ValueCell const& c : rich.cells) {
      auto const it = reads->find(c.value_id);
      std::size_t const cnt = (it == reads->end() ? 0 : it->second) + 1;
      auto const hs = scope_of(home_scope, c.value_id);
      std::wcerr << L"[homereads-cell] hash=" << (c.hash % 100000u)
                 << L" full_reads=" << cnt << L" ncells="
                 << (lat_homes.count(c.value_id) ? lat_homes[c.value_id].size()
                                                 : 0)
                 << L" home={";
      for (auto const& [ax, k] : hs)
        std::wcerr << ax.space().base_key() << L"#d" << k.depth << L"s"
                   << k.loop_slot << L" ";
      std::wcerr << L"}";
      if (auto const lit = lat_homes.find(c.value_id); lit != lat_homes.end())
        for (auto const& cell : lit->second) {
          auto const iit =
              inter->find(std::to_string(c.value_id) + "|" + sig_of(cell.home));
          std::string const s = sig_of(cell.home);
          std::wcerr << L" cell[" << std::wstring(s.begin(), s.end()) << L"]="
                     << (iit == inter->end() ? cnt : iit->second + 1);
        }
      std::wcerr << L"\n";
    }
  }

  return [reads, inter, sig_of](std::size_t vid,
                                HomeScopeKey const& home_key) -> std::size_t {
    // Per-cell: an intermediate escape level matches by (vid, home signature).
    std::string const ck = std::to_string(vid) + "|" + sig_of(home_key);
    if (auto const it = inter->find(ck); it != inter->end())
      return it->second + 1;
    // Outermost cell / single-home value: the collapsed external-consumer
    // count.
    auto const rit = reads->find(vid);
    return (rit == reads->end() ? std::size_t{0} : rit->second) + 1;
  };
}

///
/// \brief The producer/consumer pass partition of one forced-split axis TYPE
/// \p axis_key (SP2, Task 4): which values are LOOP-CARRIED on the axis (the
/// values that force the split), and which values belong to the CONSUMER pass
/// (a strict dependency-ancestor of some carried value -- it reads a carried
/// value's completed result, so it cannot run until the producer pass has
/// closed the axis loop and scattered that value to full). Every other value is
/// in the PRODUCER pass (it builds the carried values, or is unrelated).
///
struct ForcedSplitPasses {
  std::unordered_set<std::size_t> carried;        //!< LoopCarried-on-axis ids
  std::unordered_set<std::size_t> consumer_pass;  //!< strict ancestors thereof
};

inline ForcedSplitPasses forced_split_passes(std::wstring const& axis_key,
                                             LegalitySchedule const& legality,
                                             OrderedScheduleDepGraph const& g) {
  ForcedSplitPasses r;
  for (CellLegality const& cl : legality.cells) {
    bool const carried_here = std::any_of(
        cl.per_axis.begin(), cl.per_axis.end(), [&](AxisClass const& ac) {
          return ac.role == LoopRole::LoopCarried &&
                 ac.axis.space().base_key() == axis_key;
        });
    if (!carried_here) continue;
    auto const it = g.value_id_of.find(cl.hash);
    if (it != g.value_id_of.end()) r.carried.insert(it->second);
  }

  // consumer_pass = strict dependency-ancestors of the carried set: walk UP the
  // consumer edges from each carried value (a carried value that is itself an
  // ancestor of ANOTHER carried value is thereby included -- it consumes a
  // completed carried value -- while a carried value that is only a source of
  // others is not, and stays in the producer pass).
  container::svector<std::size_t> stack;
  auto const push = [&](std::size_t v) {
    if (r.consumer_pass.insert(v).second) stack.push_back(v);
  };
  for (std::size_t c : r.carried) {
    auto const it = g.consumers_of.find(c);
    if (it != g.consumers_of.end())
      for (std::size_t p : it->second) push(p);
  }
  while (!stack.empty()) {
    std::size_t const v = stack.back();
    stack.pop_back();
    auto const it = g.consumers_of.find(v);
    if (it != g.consumers_of.end())
      for (std::size_t p : it->second) push(p);
  }
  return r;
}

///
/// \brief The producer-side and consumer-side copies of a forked inner
/// sub-chain (see \c fork_subchain).
///
struct ForkedSubchain {
  container::vector<Step> producer;  //!< steps whose values are producer-side
  container::vector<Step> consumer;  //!< steps whose values are consumer-side
};

///
/// \brief Fork an already-built inner sub-chain (an ORDERED list of \c Step)
/// into a producer-side copy and a consumer-side copy, for a forced loop split
/// at a NON-innermost axis (\c build_ordered_schedule). \p
/// in_consumer(value_id) decides each value's side (true => the consumer pass).
///
/// \details A \c BuildStep goes wholly to one side by \p in_consumer of its
/// value. A nested \c ScopeBlock (an inner loop) is recursively forked; each
/// side that has surviving steps OR surviving escape \c outputs is rebuilt as a
/// per-side copy of the loop (same \c axis / \c ordinal / \c kind) carrying
/// only that side's steps and the escape \c outputs whose value lands on that
/// side; a side with neither is dropped (an empty loop is never emitted). NOTE:
/// a loop can be ALL-ESCAPE -- no \c BuildStep, only scatter \c outputs
/// contracted at the output step -- so "no surviving steps" does NOT imply "no
/// outputs to strand"; such a side must still be emitted for its outputs, else
/// a whole nested loop is silently dropped.
///
/// The relative order of the surviving steps is preserved. The input list is
/// already topologically valid (every \c ScopeBlock was built through \c
/// ordered_schedule_topo_sort_steps), and a subsequence of a valid order is
/// itself valid, so no re-sort is needed and no per-step meta is recomputed.
/// This makes the fork a pure structural transform of the \c Step tree.
///
inline ForkedSubchain fork_subchain(
    container::vector<Step> const& steps,
    std::function<bool(std::size_t)> const& in_consumer) {
  ForkedSubchain out;
  for (Step const& step : steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      (in_consumer(build->value_id) ? out.consumer : out.producer)
          .push_back(Step{BuildStep{build->value_id}});
      continue;
    }
    auto const& block = std::get<ScopeBlock>(step.value);
    ForkedSubchain sub = fork_subchain(block.steps, in_consumer);
    auto const make_side = [&](container::vector<Step>&& side_steps,
                               bool consumer_side) {
      // This side's escape outputs: those whose value lands on this side.
      container::svector<std::pair<std::size_t, OutputKind>> side_outputs;
      for (auto const& o : block.outputs)
        if (in_consumer(o.first) == consumer_side) side_outputs.push_back(o);
      // Emit the per-side loop when it has surviving STEPS *or* surviving
      // OUTPUTS. An ALL-ESCAPE loop -- one with no BuildStep, whose scatter
      // values are contracted at the output step itself -- is legitimate and
      // must NOT be dropped: doing so strands a whole nested loop (e.g. an
      // inner occ member loop or an aux loop) together with its escape outputs,
      // leaving the split axis as the only realized loop (the
      // is_range_set_congruent crash).
      if (side_steps.empty() && side_outputs.empty()) return;
      ScopeBlock fb;
      fb.axis = block.axis;
      fb.latitude_ordinal = block.latitude_ordinal;
      fb.level = block.level;
      fb.kind = block.kind;
      fb.steps = std::move(side_steps);
      fb.outputs = std::move(side_outputs);
      (consumer_side ? out.consumer : out.producer)
          .push_back(Step{std::move(fb)});
    };
    make_side(std::move(sub.producer), /*consumer_side=*/false);
    make_side(std::move(sub.consumer), /*consumer_side=*/true);
  }
  return out;
}

}  // namespace detail

///
/// \brief The SP2 \c DemotionSource (\c legality.hpp) for forced loop splits
/// (Task 4): every \c (hash, axis-base-key) that must be demoted \c LoopLocal
/// -> \c LoopCarried because a producer-side value is read across the split.
///
/// \details For each forced-split axis TYPE \c L (\c forced_split_types over
/// every cell) the split makes two ordered \c L-passes -- a PRODUCER pass that
/// builds the loop-carried values to full and a CONSUMER pass that reads them
/// (see \c build_ordered_schedule). A value \c V that is \c LoopLocal on \c L
/// (a single per-iteration copy, homed INSIDE the \c L loop) and HOMED on the
/// producer side (\c V itself is NOT in \c consumer_pass) is demoted the moment
/// ANY of its DIRECT consumers lands in the consumer pass: that consumer runs
/// in a later, disjoint \c L-loop and cannot see \c V's per-iteration
/// producer-side copy, so \c V must be materialized in full (\c LoopCarried ->
/// an \c AccumulateScatter escape output of the producer pass), lifting its
/// home floor out of \c L.
///
/// \par Why the trigger is SINGLE-SIDED (not "read from both passes")
/// The pass closure (\c forced_split_passes) is one-directional: if \c V were
/// in
/// \c consumer_pass, EVERY direct consumer of \c V would be forced into \c
/// consumer_pass too, so a consumer-homed value is never read from the producer
/// side and needs no demotion (hence the \c consumer_pass guard on \c V). The
/// only unsound case is therefore the asymmetric one -- \c V homed producer
/// whose SOLE consumer reads the carried leaf (that consumer is unconditionally
/// in \c consumer_pass). Requiring a producer-side consumer TOO would skip it,
/// leaving \c V an invisible producer-pass transient the consumer block
/// requires but no producer emits -- a silent mis-schedule. Demoting on ANY
/// consumer-pass reader closes that gap and is complete.
///
/// Wire this into \c analyze_legality as its \c demotion_source to grow the
/// monotone fixpoint (Task 4). Recomputed afresh each round on the CURRENT
/// schedule; a value already demoted is no longer \c LoopLocal, so it is not
/// re-reported, and the fixpoint converges.
///
[[nodiscard]] inline container::svector<std::pair<std::size_t, std::wstring>>
forced_split_demotions(RichSchedule const& rich,
                       LegalitySchedule const& legality) {
  auto const g = detail::ordered_schedule_dep_graph(rich);

  // forced-split axis TYPES, in first-discovery order across cells.
  container::svector<std::wstring> axes;
  for (CellLegality const& cl : legality.cells)
    for (Index const& ix : forced_split_types(cl)) {
      std::wstring key{ix.space().base_key()};
      if (std::find(axes.begin(), axes.end(), key) == axes.end())
        axes.push_back(std::move(key));
    }

  container::svector<std::pair<std::size_t, std::wstring>> out;
  for (std::wstring const& key : axes) {
    auto const passes = detail::forced_split_passes(key, legality, g);
    for (CellLegality const& cl : legality.cells) {
      bool const loop_local_here = std::any_of(
          cl.per_axis.begin(), cl.per_axis.end(), [&](AxisClass const& ac) {
            return ac.role == LoopRole::LoopLocal &&
                   ac.axis.space().base_key() == key;
          });
      if (!loop_local_here) continue;
      auto const vid_it = g.value_id_of.find(cl.hash);
      if (vid_it == g.value_id_of.end()) continue;
      std::size_t const vid = vid_it->second;
      // V must be HOMED on the producer side (not itself in consumer_pass): a
      // consumer-pass value's own consumers are, by the one-directional closure
      // in forced_split_passes, all in consumer_pass too, so it is never read
      // from the producer side and needs no demotion.
      if (passes.consumer_pass.count(vid)) continue;
      auto const cons_it = g.consumers_of.find(vid);
      if (cons_it == g.consumers_of.end()) continue;
      // SINGLE-SIDED trigger: demote as soon as ANY direct consumer of V lands
      // in the consumer pass. V is a producer-side LoopLocal transient, so
      // WITHOUT this demotion it would be an invisible per-iteration value the
      // consumer pass (a later, disjoint L-loop) cannot see -- a silent
      // mis-schedule. Requiring a producer-side consumer TOO (the old `&&`)
      // wrongly skipped the case where V's SOLE consumer reads the carried
      // value (and is thereby forced into consumer_pass): demoting V ->
      // LoopCarried materializes it as an AccumulateScatter escape output of
      // the producer pass, visible across the split.
      bool const read_by_consumer_pass = std::any_of(
          cons_it->second.begin(), cons_it->second.end(),
          [&](std::size_t c) { return passes.consumer_pass.count(c); });
      if (read_by_consumer_pass) out.push_back({cl.hash, key});
    }
  }
  return out;
}

///
/// \brief Task 3 of the ordered-scope batched-eval design (SP2): the
/// deterministic sequencer -- lowers SP1's \c LegalitySchedule (per-value
/// \c home_floor / \c per_axis roles) plus the \c RichSchedule (per-value
/// \c first_use / \c last_use over the forest's single post-order static-
/// point timeline) into an \c OrderedSchedule. Every batch axis TYPE realizes
/// ONE loop block, chained (not branched) exactly as \c build_scope_schedule's
/// single canonical chain -- EXCEPT the innermost forced-split axis, which
/// Task 4 realizes as two ordered producer/consumer passes (see step 2b and
/// \c forced_split_demotions).
///
/// \details Four-part algorithm, pure scheduling (no cost choice):
///
/// \par 1. The canonical chain
/// One loop block per distinct batch axis TYPE (\c IndexSpace::base_key())
/// appearing in ANY cell's \c CellLegality::per_axis (not just \c
/// home_floor -- a \c Reduction/\c LoopCarried-only axis, e.g. water-20's
/// \c Κ at the Κ-contraction RESULT cell, must still get a block to host its
/// escape output, even though NO cell is homed inside it in that role).
/// Ordered by \p mode_order (most-significant/outermost first), ties/
/// unlisted types alphabetical -- identical ranking to \c
/// build_scope_schedule's step 2.
///
/// \par 2. Per-value placement: home \c BuildStep vs. escape output
/// A value's \c CellLegality::per_axis has, by construction (see \c
/// CellLegality::per_axis's own doc comment), only \c LoopLocal, \c
/// Reduction, or \c LoopCarried entries (never the implicit \c
/// LoopInvariant). Two cases:
///   - EVERY \c per_axis entry is \c LoopLocal (this includes the empty
///     case: no batch-axis dependence at all, e.g. water-20's \c
///     I(i,i;a,a)) -- the value is a plain \c BuildStep. Its home block is
///     the depth whose accumulated (root-to-depth) TYPE SET equals \c
///     home_floor's TYPE SET (root if \c home_floor is empty), by the same
///     SET-equality rule \c build_scope_schedule's step 4 uses (an
///     unmatched/non-prefix \c home_floor falls back to root, mirroring
///     that function's existing behavior -- not "fixed" here). This value
///     gets NO \c outputs entry anywhere (see \c well_formed's single-
///     producer invariant): "Transient" (design point 4) is realized as
///     "produced by a \c BuildStep and nothing else", not as an explicit
///     \c OutputKind::Transient \c outputs record, since \c
///     well_formed::detail::collect_production_ids counts EVERY \c outputs
///     entry (regardless of \c OutputKind) as an independent production
///     site -- a \c Transient \c outputs entry alongside the \c BuildStep
///     would be flagged as double-production.
///   - AT LEAST ONE \c per_axis entry is \c Reduction or \c LoopCarried
///     ("escapes" that axis, per design point 4: \c Reduction ->
///     accumulate-summed out, \c LoopCarried -> accumulate-scattered out)
///     -- the value has NO \c BuildStep anywhere; instead it is recorded as
///     an \c outputs entry (kind \c AccumulateSum / \c AccumulateScatter)
///     of the DEEPEST such axis's block. "Deepest" is this function's
///     JUDGMENT CALL for the (untested-by-the-brief) case of a value
///     escaping MORE THAN ONE axis at its own node: the innermost escaped
///     loop is where the accumulation the value's own node performs
///     actually happens, so that block is its true production site: an
///     outer escape on a SHALLOWER axis would need this value already
///     complete before the outer loop can even close -- exactly consistent
///     with an outer accumulator reading an inner one. Ties (two \c
///     per_axis entries of the same axis TYPE, e.g. a same-space outer-
///     product's per-instance \c LoopCarried entries -- see \c
///     forced_split_axes's doc comment) resolve to the one axis TYPE they
///     both name, so there is no real tie to break.
///
/// \par 3. Topological order within a block -- a REAL topological sort
/// Each block's own \c steps interleave its \c BuildStep's (one per value
/// homed there) with, if the chain continues, ONE nested child \c
/// ScopeBlock \c Step for the next-deeper axis. These are ordered by
/// \c detail::ordered_schedule_topo_sort_steps against a per-step DEPENDENCY
/// GRAPH, not a scalar key alone (see below for why a scalar key cannot
/// suffice), reconstructed from \p rich alone -- no forest access needed:
///   - GLOBAL direct-dependency edges: for every \c OccurrenceRec of every
///     value, its \c consumer_point names the static point of its
///     structural PARENT node; resolving that point back to the value_id
///     whose OWN occurrence starts there (\c point_owner, built once up
///     front) recovers "this parent value directly reads that child value"
///     -- the exact same edges the forest itself encodes, without needing
///     the forest.
///   - Per LEVEL (one \c ScopeBlock's own \c steps list, including root),
///     each candidate step gets a \c detail::OrderedScheduleStepMeta:
///     - a \c BuildStep{v}'s \c produced = `{v}`; its \c requires_ = every
///       value_id \c v directly reads (raw, unfiltered -- irrelevant/
///       external entries are dropped by the topo-sort itself, since they
///       simply never match a LOCAL \c produced set).
///     - a nested child block's \c produced = that block's OWN top-level
///       \c outputs value_id's (what it makes visible to ITS OWN parent's
///       siblings -- its internal \c BuildStep's and any FURTHER-nested
///       child's content are never directly readable from outside it: by
///       construction, a value crossing a block boundary as an operand
///       must first have been resolved out of that axis, which is exactly
///       the escape/\c outputs case). Its \c requires_ is the FULL,
///       recursively bubbled external need of its WHOLE subtree (built
///       bottom-up alongside the block itself: \c requires_all(level) =
///       (this level's own direct needs UNION its child's already-bubbled
///       \c requires_all) MINUS \c produced_all(level), where
///       \c produced_all is everything ever produced ANYWHERE in the
///       subtree, recursively) -- so a need that is only satisfiable
///       several levels further out (e.g. a root-homed common factor
///       consumed by a value nested two axes deep) still surfaces at
///       whichever level can actually satisfy it.
///   - Ties (two ready steps with no dependency relation to each other) are
///     broken by \c tie_key ascending: a \c BuildStep's is its value's own
///     \c ValueCell::first_use; a child block's is the MIN \c first_use
///     over its own \c produced_all (deterministic, and -- though no longer
///     load-bearing for correctness, since the real edges now enforce both
///     directions -- still places a block as early as its own true
///     dependency slack allows).
///
/// A single SCALAR key alone cannot express both directions of this at
/// once: an earlier version of this function used the MIN-\c first_use
/// value itself as the sort key (not just a tie-break), which is provably
/// sound for "the block sorts before every true consumer" (see the MIN vs
/// MAX reasoning that was here, now superseded) but has NO corresponding
/// guarantee for "the block sorts after every true input it reads" -- a
/// value produced by a same-level sibling \c BuildStep (e.g. a root-homed
/// operand consumed by content nested inside a child block) could still
/// land, by raw point value, AFTER the block's MIN-derived key, silently
/// mis-ordering the schedule with no structural check to catch it. The real
/// topological sort above satisfies both directions by construction and is
/// checked twice (no-cycle placement count, then a second pass confirming
/// every edge survived the final order) -- see \c
/// ordered_schedule_topo_sort_steps's own doc comment.
///
/// \p policy is accepted for interface symmetry with the rest of the SP1/
/// SP2 pipeline (every stage from \c analyze_legality onward threads it)
/// and as a Task 4 hook (a future split threshold); this task's own logic
/// only consults \p rich and \p legality; the batchable-axis filtering
/// \p policy would otherwise provide is already baked into \c
/// CellLegality::per_axis by \c analyze_legality.
///
[[nodiscard]] inline OrderedSchedule build_ordered_schedule(
    RichSchedule const& rich, LegalitySchedule const& legality,
    [[maybe_unused]] BatchPolicy const& policy,
    std::initializer_list<std::wstring> mode_order = {}) {
  OrderedSchedule out;
  out.num_values = rich.cells.size();

  // hash -> value_id and the GLOBAL direct-dependency edges, recovered from
  // rich alone (see the function doc comment's part 3, and \c
  // ordered_schedule_dep_graph): for every occurrence of every value, its
  // consumer_point names its structural PARENT's own production point,
  // resolved back to the parent value_id.
  auto const g = detail::ordered_schedule_dep_graph(rich);
  auto const& value_id_of = g.value_id_of;
  static container::svector<std::size_t> const kNoDeps{};
  auto const requires_of =
      [&](std::size_t vid) -> container::svector<std::size_t> const& {
    auto const it = g.depends_on.find(vid);
    return it == g.depends_on.end() ? kNoDeps : it->second;
  };

  // The FUSION-assigned loop_slot (Task 2, compute_dag_boulevard) of a
  // legality cell's per_axis[pos]: which MEMBER of its same-space loop group
  // slices it -- an occurrence-invariant identity, established by producer->
  // consumer connectivity, NOT a within-cell position. Read off a
  // representative occurrence in the value's own frame (per_axis modes match
  // the occurrence's `carried` by label). Returns -1 if unavailable (a leaf, a
  // not-batched mode, or a divergent occurrence) -- callers fall back to slot
  // 0.
  std::unordered_map<std::size_t, ValueCell const*> hash_to_rich;
  hash_to_rich.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells) hash_to_rich.emplace(vc.hash, &vc);
  auto const fusion_slot = [&](CellLegality const& cl, std::size_t pos) -> int {
    auto const hit = hash_to_rich.find(cl.hash);
    if (hit == hash_to_rich.end() || hit->second->occurrences.empty())
      return -1;
    OccurrenceRec const& occ = hit->second->occurrences.front();
    Index const& m = cl.per_axis[pos].axis;
    auto const cit = std::find(occ.carried.begin(), occ.carried.end(), m);
    if (cit == occ.carried.end()) {
      // Not a carried mode: a Reduction mode is contracted at this value and
      // has no carried position. Its reduction loop shares loop identity with
      // the operand that home-slices it (carried->reduced propagation, stamped
      // by compute_dag_boulevard); read that slot so the reduction escape lands
      // in the SAME same-space nest as the operand it reduces, not the slot-0
      // default (a different nest -> the operand vanishes before the reduction
      // reaches it).
      for (auto const& [rm, rs] : occ.reduced_slot)
        if (rm == m) return rs;
      return -1;
    }
    std::size_t const p = static_cast<std::size_t>(cit - occ.carried.begin());
    return p < occ.loop_slot.size() ? occ.loop_slot[p] : -1;
  };

  // 1. The canonical chain: one representative Index per distinct axis TYPE
  // present in ANY cell's per_axis (LoopLocal, Reduction, OR LoopCarried --
  // NOT just home_floor; see the function doc comment's part 1).
  // Per-INSTANCE loop chain (2026-08-29 position-based de-collapse): for each
  // space, m_s = the MAX over cells of the number of same-space per_axis modes
  // in one cell; emit m_s consecutive depths, one per within-space SLOT
  // (loop_slot 0..m_s-1). A value carrying two same-space batched modes (a
  // doubles amplitude's two occ externals) thus gets TWO distinct loops instead
  // of one -- the collapse fix. `types[d]` is a representative Index of the
  // depth's space; `type_slot[d]` is its within-space slot (the DagScopeLevel
  // loop_slot). NOTE: same-space slots occupy distinct DEPTHS here (the
  // assembly is one loop per depth); depth carries both group and member
  // nesting for now.
  // Members present per space = the distinct FUSION loop_slots (Task 2) that
  // appear on that space across all cells. Each distinct slot becomes one
  // realized loop. (Was: the MAX same-space per_axis COUNT with local 0..m-1
  // slots -- e8bcee766's position-based numbering, which the atlas could not
  // match to a value's own frame; the fusion slot is that occurrence-invariant
  // identity.)
  std::map<std::wstring, std::set<int>> slots_of_space;
  std::map<std::wstring, Index> rep;  // space -> representative axis
  for (CellLegality const& cl : legality.cells)
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
      rep.emplace(bk, cl.per_axis[pos].axis);
      int const s = fusion_slot(cl, pos);
      slots_of_space[bk].insert(s >= 0 ? s : 0);
    }
  auto const rank_of = [&](std::wstring const& bk) -> std::size_t {
    std::size_t i = 0;
    for (auto const& key : mode_order) {
      if (bk == key) return i;
      ++i;
    }
    return static_cast<std::size_t>(-1);
  };
  container::svector<std::wstring> spaces;
  for (auto const& [bk, s] : slots_of_space) spaces.push_back(bk);
  std::sort(spaces.begin(), spaces.end(),
            [&](std::wstring const& a, std::wstring const& b) {
              auto const ra = rank_of(a), rb = rank_of(b);
              if (ra != rb) return ra < rb;
              return a < b;
            });
  container::svector<Index> types;
  container::svector<int> type_slot;  // the FUSION loop_slot of this loop
  for (auto const& bk : spaces)
    for (int s : slots_of_space.at(bk)) {  // ascending (std::set)
      types.push_back(rep.at(bk));
      type_slot.push_back(s);
    }

  // TEMP instrumentation (P1 Task 2 "before"): the realized loop chain is one
  // representative per SPACE (the collapse). Guarded by SEQUANT_DUMP_SCHEDULE.
  if (std::getenv("SEQUANT_DUMP_SCHEDULE")) {
    std::wcerr << L"[sched] per-instance loop chain: ";
    for (std::size_t d = 0; d < types.size(); ++d)
      std::wcerr << L"d" << d << L"=" << types[d].space().base_key() << L"#slot"
                 << type_slot[d] << L" ";
    std::wcerr << L"\n";
  }

  std::size_t const n = types.size();
  // Depth of the loop for a given (space, within-space slot).
  auto const depth_of_instance = [&](std::wstring const& bk,
                                     int slot) -> std::optional<std::size_t> {
    for (std::size_t d = 0; d < n; ++d)
      if (std::wstring(types[d].space().base_key()) == bk &&
          type_slot[d] == slot)
        return d;
    return std::nullopt;
  };

  // Co-occurrence clusters (un-fuse). Two loop members co-occur iff some value
  // is home-sliced on BOTH (carries both in one term, so its fusion loop_slots
  // -- Task 2, now home-based -- name both). Members that never co-occur live
  // in DISJOINT nests: e.g. two residual sub-DAGs that both batch occ i,j but
  // are connected only THROUGH a full symmetric intermediate (home meet empties
  // it, so it seeds no loop) -- they are genuinely separate loop groups.
  // Realizing them as one over-deep nested chain (slot0 superset ... superset
  // slotK) is the structure that deadlocks; each cluster must be a SEPARATE
  // sequential nest at root. Union-find over member depths.
  container::svector<std::size_t> mem_parent(n);
  for (std::size_t d = 0; d < n; ++d) mem_parent[d] = d;
  auto const mfind = [&](std::size_t x) {
    while (mem_parent[x] != x) {
      mem_parent[x] = mem_parent[mem_parent[x]];
      x = mem_parent[x];
    }
    return x;
  };
  for (CellLegality const& cl : legality.cells) {
    std::optional<std::size_t> anchor;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      int const fs = fusion_slot(cl, pos);
      if (fs < 0) continue;
      auto const dd = depth_of_instance(
          std::wstring{cl.per_axis[pos].axis.space().base_key()}, fs);
      if (!dd) continue;
      if (anchor)
        mem_parent[mfind(*dd)] = mfind(*anchor);
      else
        anchor = dd;
    }
  }
  container::svector<std::size_t> type_cluster(n);
  for (std::size_t d = 0; d < n; ++d) type_cluster[d] = mfind(d);

  // Assembly processing order: cluster by cluster (each cluster ordered by its
  // outermost = min depth, for determinism), and within a cluster innermost
  // (larger d) FIRST, so the loop wraps a cluster's members into one nest and
  // finalizes that nest at the cluster boundary.
  container::svector<std::size_t> order;
  {
    std::map<std::size_t, std::size_t> cluster_min;
    for (std::size_t d = 0; d < n; ++d) {
      auto const it = cluster_min.find(type_cluster[d]);
      if (it == cluster_min.end() || d < it->second)
        cluster_min[type_cluster[d]] = d;
    }
    container::svector<std::pair<std::size_t, std::size_t>> cl_by_min;
    for (auto const& [c, mn] : cluster_min) cl_by_min.push_back({mn, c});
    std::sort(cl_by_min.begin(), cl_by_min.end());
    for (auto const& [mn, c] : cl_by_min)
      for (std::size_t j = 0; j < n; ++j) {
        std::size_t const d = n - 1 - j;  // descending: innermost first
        if (type_cluster[d] == c) order.push_back(d);
      }
  }

  // 2. Per-value placement: home BuildStep (root-level bucket uses index n
  // as a sentinel "root" depth) vs. escape output at the deepest escaped
  // axis's depth.
  container::vector<detail::OrderedScheduleDepthBucket> buckets(n);
  container::svector<std::size_t> root_build_ids;

  for (CellLegality const& cl : legality.cells) {
    auto const vid_it = value_id_of.find(cl.hash);
    SEQUANT_ASSERT(vid_it != value_id_of.end());
    std::size_t const vid = vid_it->second;

    // A forest LEAF is an input fetched on demand by its consumers, never a
    // computed value: it must NOT be emitted as a BuildStep (doing so makes the
    // executor evaluate it as a standalone root -- one wasted leaf fetch per
    // iteration -- which the forest descent never does). Its consumers reach it
    // through the leaf evaluator exactly as in forest descent.
    if (rich.cells[vid].is_leaf) continue;

    // Emit an escape at EVERY non-local axis DEPTH (the multi-level escape
    // CHAIN, SP2 non-innermost split): a value that reduces an inner axis AND
    // is carried on an outer one escapes at BOTH -- AccumulateSum at the inner
    // (into the accumulator one level out) then AccumulateScatter at the outer
    // (that accumulator to full). The bottom-up assembly materializes them
    // inner -> outer. A value non-local on a SINGLE axis keeps exactly one
    // escape, unchanged. Same-depth axis-classes (e.g. two carried occ indices,
    // or a same-type reduce+carry pair) collapse to ONE escape at that depth,
    // with LoopCarried (AccumulateScatter) dominating Reduction: a carried axis
    // must materialize to full even if a same-type index reduces.
    container::svector<std::pair<std::size_t, OutputKind>> escapes;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      AxisClass const& ac = cl.per_axis[pos];
      if (ac.role == LoopRole::LoopLocal) continue;
      std::wstring const bk{ac.axis.space().base_key()};
      int const fs = fusion_slot(cl, pos);
      auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
      SEQUANT_ASSERT(d.has_value());  // types was built from this same union
      OutputKind const kind =
          (ac.role == LoopRole::Reduction)
              ? OutputKind::AccumulateSum
              : OutputKind::AccumulateScatter;  // LoopCarried
      // Each per-instance mode escapes to its OWN depth (distinct loop_slot),
      // so same-space modes no longer collapse to one escape; the dedup below
      // only merges a genuine reduce+carry pair that lands at ONE depth.
      auto it = std::find_if(escapes.begin(), escapes.end(),
                             [&](auto const& e) { return e.first == *d; });
      if (it == escapes.end())
        escapes.push_back({*d, kind});
      else if (kind == OutputKind::AccumulateScatter)
        it->second = OutputKind::AccumulateScatter;
    }

    // TEMP instrumentation (P1 Task 2 "before"): a value with >1 non-local mode
    // of ONE space (a doubles amplitude / PPL product carrying two occ
    // externals) has its distinct per-instance escapes COLLAPSED to fewer
    // escapes (one per depth == one per space). Dump those cells to expose the
    // collapse. Guarded by SEQUANT_DUMP_SCHEDULE.
    if (std::getenv("SEQUANT_DUMP_SCHEDULE")) {
      container::svector<Index> nonlocal;
      for (AxisClass const& ac : cl.per_axis)
        if (ac.role != LoopRole::LoopLocal) nonlocal.push_back(ac.axis);
      bool collapse = false;
      for (std::size_t a = 0; a < nonlocal.size() && !collapse; ++a)
        for (std::size_t b = a + 1; b < nonlocal.size(); ++b)
          if (nonlocal[a].space().base_key() ==
              nonlocal[b].space().base_key()) {
            collapse = true;
            break;
          }
      if (collapse) {
        std::wcerr << L"[sched-collapse] hash=" << cl.hash << L" nonlocal={";
        for (Index const& ix : nonlocal)
          std::wcerr << ix.full_label() << L":"
                     << (std::find_if(
                             cl.per_axis.begin(), cl.per_axis.end(),
                             [&](AxisClass const& ac) {
                               return ac.axis == ix;
                             })->role == LoopRole::Reduction
                             ? L"R"
                             : L"C")
                     << L" ";
        std::wcerr << L"} -> " << nonlocal.size() << L" modes COLLAPSE to "
                   << escapes.size() << L" escapes(depth:kind)={";
        for (auto const& [d, k] : escapes)
          std::wcerr << d << L":"
                     << (k == OutputKind::AccumulateScatter ? L"Scatter"
                                                            : L"Sum")
                     << L" ";
        std::wcerr << L"}\n";
      }
    }

    if (!escapes.empty()) {
      for (auto const& [d, kind] : escapes)
        buckets[d].outputs.push_back({vid, kind});
      continue;
    }

    // Plain BuildStep: home at the INNERMOST loop the value is LoopLocal on,
    // resolved PER-INSTANCE by fusion loop_slot -- NOT by shallowest same-space
    // count. A value local to occ slots 2,3 (its own fusion nest) homes in THAT
    // nest, not the FIRST {i,i} nest a space-multiset cover would pick; picking
    // the wrong same-space nest homes the value where its consumer's nest has
    // not opened (or has already closed), so an in-consumer-nest read misses
    // and the value vanishes. Resolve each LoopLocal mode's (space, fusion
    // slot) to its realized depth (exactly as the escape placement above does),
    // and home at the MAX such depth: within a co-occurrence cluster larger
    // depth nests inside smaller, so the innermost of the value's own home
    // slots is inside all of them. home_floor is the LoopLocal subset, but it
    // drops the pos->slot map, so walk per_axis directly for the slot.
    std::optional<std::size_t> target;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      if (cl.per_axis[pos].role != LoopRole::LoopLocal) continue;
      std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
      int const fs = fusion_slot(cl, pos);
      auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
      if (!d) continue;
      if (!target || *d > *target) target = *d;
    }

    if (target)
      buckets[*target].build_ids.push_back(vid);
    else
      root_build_ids.push_back(vid);
  }

  auto const union_into = [](container::svector<std::size_t>& acc,
                             container::svector<std::size_t> const& add) {
    for (std::size_t v : add)
      if (std::find(acc.begin(), acc.end(), v) == acc.end()) acc.push_back(v);
  };

  // 2b. FORCED LOOP SPLIT (Task 4). If the INNERMOST realized axis type
  // (types[n-1]) is forced to split -- some cell is LoopCarried on it -- and
  // its values divide into a non-empty PRODUCER pass and a non-empty CONSUMER
  // pass (a value is in the consumer pass iff it is a strict dependency-
  // ancestor of a loop-carried value, i.e. it reads that value's completed
  // result: see detail::forced_split_passes), the {L} loop is realized as TWO
  // ordered sibling blocks with distinct ordinals -- a producer pass (ordinal
  // 0, builds the loop-carried values to full via their AccumulateScatter
  // outputs) then a consumer pass (ordinal 1, the cross-iteration reads).
  //
  // Restricted to the INNERMOST axis (d == n-1) so there is no deeper loop
  // chain to fork across the two passes. A forced-split axis that is NOT the
  // innermost realized axis is left UNSPLIT here (a single block, exactly as
  // Task 3) -- a clearly-documented INERT hook deferred to a later task: the
  // general split forks the whole enclosed sub-chain, and no current fixture
  // (nor any consumer -- SP3's executor does not yet read this IR) exercises
  // it, so a sound conservative single block beats an unsound partial fork.
  // FORCED LOOP SPLIT detection (Task 4, generalized to any depth). Find the
  // OUTERMOST realized axis forced to split (some cell LoopCarried on it) whose
  // split is GENUINE -- both a pure producer (a carried value read only in a
  // later pass) and a consumer (a strict ancestor of a carried value). The
  // split is realized at that depth d*; when d* is not the innermost axis the
  // enclosed inner sub-chain is FORKED across the two passes (fork_subchain).
  // At most one forced axis is expected in practice (only an external axis is
  // LoopCarried; contracted axes are Reduction) -- a second, nested forced
  // split is a distinct feature and is rejected loudly rather than mis-emitted.
  std::optional<std::size_t> split_depth;
  std::optional<detail::ForcedSplitPasses> split_passes;
  {
    std::size_t forced_count = 0;
    container::svector<std::wstring> seen_split_spaces;
    for (std::size_t d = 0; d < n; ++d) {
      std::wstring const key{types[d].space().base_key()};
      // A space now spans several per-instance slot depths; the PROCON split is
      // per SPACE (base_key), so consider each space once, at its first (slot
      // 0) depth. (Per-slot PROCON of a multi-member group is deferred.)
      if (std::find(seen_split_spaces.begin(), seen_split_spaces.end(), key) !=
          seen_split_spaces.end())
        continue;
      seen_split_spaces.push_back(key);
      bool const forced =
          std::any_of(legality.cells.begin(), legality.cells.end(),
                      [&](CellLegality const& cl) {
                        for (Index const& ix : forced_split_types(cl))
                          if (ix.space().base_key() == key) return true;
                        return false;
                      });
      if (!forced) continue;
      auto passes = detail::forced_split_passes(key, legality, g);
      bool const has_cons = !passes.consumer_pass.empty();
      bool const has_prod = std::any_of(
          passes.carried.begin(), passes.carried.end(),
          [&](std::size_t v) { return !passes.consumer_pass.count(v); });
      if (!has_prod || !has_cons) continue;  // LoopCarried but no split needed
      ++forced_count;
      if (!split_depth) {  // outermost genuine split
        split_depth = d;
        split_passes = std::move(passes);
      }
    }
    SEQUANT_ASSERT(forced_count <= 1,
                   "build_ordered_schedule: multiple forced-split axes at "
                   "different depths (nested splits) not supported");
  }

  // Small helpers shared by the (usual) single-block path and the split path.
  auto const external_needs =
      [&](container::svector<std::size_t> const& produced)
      -> container::svector<std::size_t> {
    container::svector<std::size_t> req;
    for (std::size_t v : produced) union_into(req, requires_of(v));
    req.erase(std::remove_if(req.begin(), req.end(),
                             [&](std::size_t v) {
                               return std::find(produced.begin(),
                                                produced.end(),
                                                v) != produced.end();
                             }),
              req.end());
    return req;
  };
  auto const min_first_use =
      [&](container::svector<std::size_t> const& ids) -> std::size_t {
    std::size_t k = std::numeric_limits<std::size_t>::max();
    for (std::size_t v : ids) k = std::min(k, rich.cells[v].first_use);
    return k == std::numeric_limits<std::size_t>::max() ? std::size_t{0} : k;
  };
  auto const make_block =
      [&](Index const& axis, int latitude_ordinal, std::size_t depth,
          int loop_slot, container::svector<std::size_t> const& build_ids,
          container::svector<std::pair<std::size_t, OutputKind>> const& outputs,
          container::vector<Step>&& child_steps,
          container::vector<detail::OrderedScheduleStepMeta>&& child_metas)
      -> ScopeBlock {
    container::vector<Step> items;
    container::vector<detail::OrderedScheduleStepMeta> meta;
    for (std::size_t v : build_ids) {
      items.push_back(Step{BuildStep{v}});
      detail::OrderedScheduleStepMeta m;
      m.produced.push_back(v);
      m.requires_.assign(requires_of(v).begin(), requires_of(v).end());
      m.tie_key = rich.cells[v].first_use;
      meta.push_back(std::move(m));
    }
    for (std::size_t k = 0; k < child_steps.size(); ++k) {
      items.push_back(std::move(child_steps[k]));
      meta.push_back(std::move(child_metas[k]));
    }
    ScopeBlock block;
    block.axis = axis;
    block.latitude_ordinal = latitude_ordinal;
    block.level = DagScopeLevel{.depth = depth,
                                .space = std::wstring{axis.space().base_key()},
                                .loop_slot = loop_slot,
                                .latitude_ordinal = latitude_ordinal};
    block.kind = detail::mode_is_external(rich, axis)
                     ? BatchModeType::External
                     : BatchModeType::Contracted;
    block.steps =
        detail::ordered_schedule_topo_sort_steps(std::move(items), meta);
    block.outputs.assign(outputs.begin(), outputs.end());
    return block;
  };

  // 3. Assemble the chain bottom-up (innermost first). Thread up a LIST of
  // child block Steps to embed one level out -- normally ONE (the single block
  // for the deeper axis, exactly as Task 3), but TWO at the split depth (the
  // producer/consumer passes) -- each paired with the meta the outer topo-sort
  // needs (its OWN escape outputs as `produced`, its whole subtree's external
  // need as `requires_`, and a deterministic `tie_key`). child_produced_all /
  // child_requires_all carry the FULL recursive produced/external-need sets of
  // everything built at this depth (identical whether one block or two -- the
  // outer level sees the same production/need set either way) to grow the next
  // level's own sets, per the function doc comment's part 3.
  container::vector<Step> pending_steps;
  container::vector<detail::OrderedScheduleStepMeta> pending_metas;
  container::svector<std::size_t> child_produced_all;
  container::svector<std::size_t> child_requires_all;
  // Completed cluster nests, concatenated at root (each an independent nest).
  container::vector<Step> finished_steps;
  container::vector<detail::OrderedScheduleStepMeta> finished_metas;
  std::optional<std::size_t> prev_cluster;

  for (std::size_t oi = 0; oi < order.size(); ++oi) {
    std::size_t const d = order[oi];  // cluster-grouped, innermost-first
    // Cluster boundary: the previous cluster's nest is complete in `pending`;
    // move it to `finished` (a separate top-level nest) and start fresh.
    if (prev_cluster && type_cluster[d] != *prev_cluster) {
      for (std::size_t k = 0; k < pending_steps.size(); ++k) {
        finished_steps.push_back(std::move(pending_steps[k]));
        finished_metas.push_back(std::move(pending_metas[k]));
      }
      pending_steps.clear();
      pending_metas.clear();
      child_produced_all.clear();
      child_requires_all.clear();
    }
    prev_cluster = type_cluster[d];
    detail::OrderedScheduleDepthBucket const& bucket = buckets[d];

    container::svector<std::size_t> own_produced;
    for (std::size_t v : bucket.build_ids) own_produced.push_back(v);
    for (auto const& out_entry : bucket.outputs)
      own_produced.push_back(out_entry.first);

    container::svector<std::size_t> produced_all = own_produced;
    union_into(produced_all, child_produced_all);

    container::svector<std::size_t> requires_all;
    for (std::size_t v : own_produced) union_into(requires_all, requires_of(v));
    union_into(requires_all, child_requires_all);
    requires_all.erase(std::remove_if(requires_all.begin(), requires_all.end(),
                                      [&](std::size_t v) {
                                        return std::find(produced_all.begin(),
                                                         produced_all.end(),
                                                         v) !=
                                               produced_all.end();
                                      }),
                       requires_all.end());

    container::vector<Step> next_steps;
    container::vector<detail::OrderedScheduleStepMeta> next_metas;

    if (split_passes && d == *split_depth) {
      // Split at depth d*: partition this bucket's homed BuildSteps and escape
      // outputs into the producer (ordinal 0) and consumer (ordinal 1) passes
      // by consumer_pass membership, and FORK the enclosed inner sub-chain
      // (pending_steps) across the two passes. At the innermost axis pending is
      // empty and fork_subchain returns two empty inner lists, reducing to the
      // original childless two-block emission (byte-identical). The consumer
      // block reads the producer block's escaped (now-full) outputs, so its
      // `requires_` names them and the outer topo-sort orders producer first.
      auto const in_consumer = [&](std::size_t vid) {
        return split_passes->consumer_pass.count(vid) != 0;
      };

      detail::ForkedSubchain forked =
          detail::fork_subchain(pending_steps, in_consumer);

      // Recursive value_ids a forked child Step produces (builds + escape
      // outputs, through nested blocks), for its `requires_`/`tie_key` meta.
      std::function<void(ScopeBlock const&, container::svector<std::size_t>&)>
          collect_rec =
              [&](ScopeBlock const& blk, container::svector<std::size_t>& out) {
                for (Step const& s : blk.steps) {
                  if (auto const* b = std::get_if<BuildStep>(&s.value))
                    out.push_back(b->value_id);
                  else
                    collect_rec(std::get<ScopeBlock>(s.value), out);
                }
                for (auto const& o : blk.outputs) out.push_back(o.first);
              };
      auto const meta_for = [&](container::vector<Step> const& steps)
          -> container::vector<detail::OrderedScheduleStepMeta> {
        container::vector<detail::OrderedScheduleStepMeta> metas;
        for (Step const& s : steps) {
          detail::OrderedScheduleStepMeta m;
          if (auto const* b = std::get_if<BuildStep>(&s.value)) {
            m.produced.push_back(b->value_id);
            m.requires_.assign(requires_of(b->value_id).begin(),
                               requires_of(b->value_id).end());
            m.tie_key = rich.cells[b->value_id].first_use;
          } else {
            auto const& blk = std::get<ScopeBlock>(s.value);
            container::svector<std::size_t> rec;
            collect_rec(blk, rec);
            for (auto const& o : blk.outputs) m.produced.push_back(o.first);
            m.requires_ = external_needs(rec);
            m.tie_key = min_first_use(rec);
          }
          metas.push_back(std::move(m));
        }
        return metas;
      };

      auto const emit_pass =
          [&](int latitude_ordinal,
              container::svector<std::size_t> const& builds,
              container::svector<std::pair<std::size_t, OutputKind>> const&
                  outs,
              container::vector<Step>&& child_steps,
              container::vector<detail::OrderedScheduleStepMeta>&&
                  child_metas) {
            container::svector<std::size_t> pass_produced = builds;
            for (auto const& o : outs) pass_produced.push_back(o.first);
            for (Step const& s : child_steps) {
              if (auto const* b = std::get_if<BuildStep>(&s.value))
                pass_produced.push_back(b->value_id);
              else
                collect_rec(std::get<ScopeBlock>(s.value), pass_produced);
            }
            next_steps.push_back(Step{make_block(
                types[d], latitude_ordinal, d + 1, type_slot[d], builds, outs,
                std::move(child_steps), std::move(child_metas))});
            detail::OrderedScheduleStepMeta m;
            for (auto const& o : outs) m.produced.push_back(o.first);
            m.requires_ = external_needs(pass_produced);
            m.tie_key = min_first_use(pass_produced);
            next_metas.push_back(std::move(m));
          };

      container::svector<std::size_t> prod_builds, cons_builds;
      for (std::size_t v : bucket.build_ids)
        (in_consumer(v) ? cons_builds : prod_builds).push_back(v);
      container::svector<std::pair<std::size_t, OutputKind>> prod_outs,
          cons_outs;
      for (auto const& o : bucket.outputs)
        (in_consumer(o.first) ? cons_outs : prod_outs).push_back(o);

      auto prod_metas = meta_for(forked.producer);
      auto cons_metas = meta_for(forked.consumer);
      emit_pass(0, prod_builds, prod_outs, std::move(forked.producer),
                std::move(prod_metas));
      emit_pass(1, cons_builds, cons_outs, std::move(forked.consumer),
                std::move(cons_metas));
    } else {
      ScopeBlock block = make_block(
          types[d], 0, d + 1, type_slot[d], bucket.build_ids, bucket.outputs,
          std::move(pending_steps), std::move(pending_metas));
      next_steps.push_back(Step{std::move(block)});
      detail::OrderedScheduleStepMeta m;
      // RECURSIVE produced (not just this block's own outputs): a nest
      // advertises EVERYTHING it produces so the topo sort can order a sibling
      // nest that consumes a value produced by an INNER block of this one. With
      // the single chain this never mattered (one block per level, no
      // siblings); the un-fuse emits separate sibling nests at root, so an
      // under-reported `produced` orders a consumer nest before its producer ->
      // read-before-build.
      m.produced = produced_all;
      m.requires_ = requires_all;
      m.tie_key = min_first_use(produced_all);
      next_metas.push_back(std::move(m));
    }

    pending_steps = std::move(next_steps);
    pending_metas = std::move(next_metas);
    child_produced_all = std::move(produced_all);
    child_requires_all = std::move(requires_all);
  }
  // Finalize the LAST cluster's nest.
  for (std::size_t k = 0; k < pending_steps.size(); ++k) {
    finished_steps.push_back(std::move(pending_steps[k]));
    finished_metas.push_back(std::move(pending_metas[k]));
  }

  // Root assembly: root-level BuildStep's plus the SEPARATE cluster nests as
  // sibling top-level Steps (each an independent nest; the topo sort orders
  // them by dependency, and the executor runs sibling root steps sequentially).
  container::vector<Step> root_items;
  container::vector<detail::OrderedScheduleStepMeta> root_meta;
  for (std::size_t v : root_build_ids) {
    root_items.push_back(Step{BuildStep{v}});
    detail::OrderedScheduleStepMeta m;
    m.produced.push_back(v);
    m.requires_.assign(requires_of(v).begin(), requires_of(v).end());
    m.tie_key = rich.cells[v].first_use;
    root_meta.push_back(std::move(m));
  }
  for (std::size_t k = 0; k < finished_steps.size(); ++k) {
    root_items.push_back(std::move(finished_steps[k]));
    root_meta.push_back(std::move(finished_metas[k]));
  }
  out.root.steps = detail::ordered_schedule_topo_sort_steps(
      std::move(root_items), root_meta);

  // Pillar 1: record each value's home-sliced-mode -> DAG-scope depth, ONCE, in
  // the value's OWN frame. Walk the realized block tree accumulating the
  // enclosing levels (axis + depth); a BuildStep homes inside all of them, an
  // escape output homes one level OUT (mirrors ordered_home_reads' home_scope).
  // For each value pair its OWN `carried` (value frame) with those levels by
  // space + nest position -- the same frame-safe routing occ_facts uses -- and
  // keep only modes that are the value's `home_modes` (its home-sliced set).
  {
    auto const record =
        [&out, &rich](std::size_t vid,
                      container::svector<std::pair<Index, int>> const& levels) {
          if (vid >= rich.cells.size()) return;
          ValueCell const& c = rich.cells[vid];
          container::svector<std::pair<Index, int>> mode_depth;
          container::svector<bool> consumed(c.carried.size(), false);
          for (auto const& [axis, depth] : levels) {
            auto const space = axis.space().base_key();
            for (std::size_t j = 0; j < c.carried.size(); ++j) {
              if (consumed[j]) continue;
              if (std::wstring(c.carried[j].space().base_key()) != space)
                continue;
              consumed[j] = true;
              mode_depth.push_back({c.carried[j], depth});
              break;
            }
          }
          if (!mode_depth.empty())
            out.home_mode_depth[vid] = std::move(mode_depth);
        };
    std::function<void(ScopeBlock const&,
                       container::svector<std::pair<Index, int>> const&)>
        walk = [&](ScopeBlock const& b,
                   container::svector<std::pair<Index, int>> const& enc) {
          for (Step const& s : b.steps) {
            if (auto const* build = std::get_if<BuildStep>(&s.value)) {
              record(build->value_id, enc);
            } else if (auto const* child = std::get_if<ScopeBlock>(&s.value)) {
              container::svector<std::pair<Index, int>> inner = enc;
              // Store the FULL loop identity color (depth AND loop_slot), not
              // depth alone: one cache per LOOP (LoopKey::color), so the value-
              // id coloring and the home-scope filter distinguish same-group
              // sibling loops.
              inner.push_back(
                  {child->axis, static_cast<int>(child->level.key().color())});
              walk(*child, inner);
            }
          }
          // escape output homes one level OUT (enc minus this block's own
          // axis).
          container::svector<std::pair<Index, int>> out_scope = enc;
          if (!out_scope.empty()) out_scope.pop_back();
          for (auto const& [vid, kind] : b.outputs) {
            (void)kind;
            record(vid, out_scope);
          }
        };
    walk(out.root, {});
  }

  // Pillar 1 / B-full: persist the value/occurrence DAG edges (each value's
  // direct operand value_ids) the value-driven ordered executor consumes. `g`
  // is the same dep graph the topo-sort used above; its `depends_on` edges are
  // derived from every OccurrenceRec's consumer_point, so a split operand
  // resolves to the specific consumed value (not an ambiguous node hash).
  out.operand_vids = g.depends_on;

  SEQUANT_ASSERT(well_formed(out));
  return out;
}

namespace detail {

/// \brief One enclosing loop block on a value's realized eval scope: the
/// block's canonical representative \c axis and its \c DagScopeLevel (\c
/// depth / \c space / \c ordinal), read straight off \c ordered's block tree
/// (so the \c ordinal of a forced-split sibling is the REAL one the runtime
/// pushes, not a guessed one).
struct ScopeBlockAxisLevel {
  Index axis;
  DagScopeLevel level;
};

/// \brief \c compute_sliced_mode_assignment's block-tree walk: record, for
/// every value PRODUCED (a \c BuildStep) or ESCAPED (a block \c outputs
/// entry), the ordered list of enclosing loop blocks (\c axis + \c level) it
/// is EVALUATED under -- i.e. the loops whose batch it reads its operands
/// inside. This is exactly \c ordered_home_reads' \c build_scope, but keeping
/// each block's \c DagScopeLevel (not just its axis) so the map can name the
/// runtime \c BatchContextEntry::level.
///
/// \p enc is the root-to-\p block path INCLUDING \p block's own (axis, level)
/// -- a value built or escaping inside \p block reads its operands inside
/// \p block's own loop, so the block's own axis encloses that read (mirrors
/// \c ordered_home_reads' \c build_scope, which likewise includes the block
/// axis; only the HOME of an escape sits one level out, which is irrelevant
/// here -- we want the read/fetch scope, not the home).
inline void populate_build_scope_walk(
    ScopeBlock const& block, container::svector<ScopeBlockAxisLevel> const& enc,
    std::unordered_map<std::size_t, container::svector<ScopeBlockAxisLevel>>&
        build_scope) {
  for (Step const& step : block.steps) {
    if (auto const* build = std::get_if<BuildStep>(&step.value)) {
      build_scope[build->value_id] = enc;
    } else {
      ScopeBlock const& child = std::get<ScopeBlock>(step.value);
      container::svector<ScopeBlockAxisLevel> inner = enc;
      inner.push_back({child.axis, child.level});
      populate_build_scope_walk(child, inner, build_scope);
    }
  }
  for (auto const& [vid, kind] : block.outputs) {
    (void)kind;
    // escape is BUILT inside `block` (reads operands in this loop); its HOME is
    // one level out, but we want the read/fetch scope here. DEEPEST-escape
    // wins: a value carried on more than one nested axis escapes at EACH depth,
    // but the outer escapes only scatter an already-built value (read no
    // operands); the contraction runs at the innermost escape, so its operands
    // are sliced there. The child recursion above visits the deeper block
    // first, so a shallower outer escape must NOT clobber the deeper scope --
    // else the seam slices the contraction's operands on only the outer axis
    // while the inner axis is open (the is_range_set_congruent crash at the
    // multi-occ product).
    auto const it = build_scope.find(vid);
    if (it == build_scope.end() || enc.size() > it->second.size())
      build_scope[vid] = enc;
  }
}

/// \brief Debug safety net: \c compute_sliced_mode_assignment's canonical
/// level enumeration (\c enumerate_realized_levels, which folds every
/// realized \c ScopeBlock's \c DagScopeLevel into one \c LoopId per DISTINCT
/// level) leans on \c (level.depth, level.space, level.ordinal) being unique
/// GLOBALLY across the whole realized tree -- i.e. any two blocks sharing
/// that triple must also share the same representative \c axis, or the
/// enumeration could fold two structurally-different loops onto one \c
/// LoopId. \c well_formed only checks ordinal uniqueness among a block's OWN
/// same-axis DIRECT children (sibling-local, see \c
/// ordered_schedule_block_well_formed above) -- it says nothing about two
/// blocks at the same (depth, space, ordinal) that are NOT siblings (say,
/// nested under different parents). Walk every block in the tree and assert
/// the stronger, global invariant loudly: a violation means the scheduler
/// emitted two structurally-distinct loops the level-to-\c LoopId mapping
/// cannot tell apart. Debug-only insurance; expected to never fire on any
/// fixture.
inline void assert_global_level_axis_uniqueness(
    ScopeBlock const& block,
    std::map<std::tuple<std::size_t, std::wstring, int, int>, Index>& seen) {
  for (Step const& step : block.steps) {
    auto const* child = std::get_if<ScopeBlock>(&step.value);
    if (!child) continue;
    auto const key =
        std::make_tuple(child->level.depth, child->level.space,
                        child->level.loop_slot, child->level.latitude_ordinal);
    auto const [it, inserted] = seen.try_emplace(key, child->axis);
    SEQUANT_ASSERT(
        (inserted || it->second == child->axis) &&
        "assert_global_level_axis_uniqueness: two blocks at the same "
        "DAG-scope level (depth, space, loop_slot, latitude) realize DIFFERENT "
        "representative axes -- the canonical level enumeration assumes "
        "this triple names one axis globally, not just among a block's own "
        "direct siblings");
    assert_global_level_axis_uniqueness(*child, seen);
  }
}

}  // namespace detail

///
/// \brief Task 3 of the sliced-value canonical layout / loop-coloring design
/// (\c doc/dev/specs/2026-08-23-sliced-value-canonical-layout-loop-coloring-
/// design.md): a stable identifier for one DAG-scope loop realized by \c
/// build_ordered_schedule -- an index into \c SlicedModeAssignment::levels,
/// the schedule's own canonical (deterministic pre-order) enumeration of
/// every non-root \c ScopeBlock's \c DagScopeLevel. Two blocks that are
/// STRUCTURALLY the same loop (identical \c (depth, space, ordinal)) always
/// share one \c LoopId; two blocks that differ in ANY of those three --
/// including two forced-split SIBLING passes, which differ only in \c
/// ordinal -- get DISTINCT ids by construction (see the design's "Loop
/// identity is a slot color" section: producer/consumer passes must be
/// distinguishable colors, not folded).
///
/// \note Defined in \c dag_scope.hpp so the low-level \c LoopColoredSliceSeam
/// (the executor's runtime consumption view, consumed by \c CacheManager) can
/// name it without depending on this schedule header; re-exported here for the
/// schedule-side code that has always referred to \c eval::LoopId.
using sequant::LoopId;

///
/// \brief The per-(value, sliced-mode) -> DAG-scope-loop ASSIGNMENT: the
/// coloring input Task 5 feeds to \c canonicalize_slots's \c
/// NamedIndexColorMap. Unlike the per-cell \c ModeToLevel map this
/// deliberately superseded (removed by Task 8; see the design doc's
/// migration section), this is plain DATA keyed by the value's OWN physical
/// \c Index label for each mode it is sliced on -- not a per-cell position
/// map -- so a relabeled CSE participant (Task 4's regime-2 case) is keyed by
/// its own label, and a symmetric value's two occurrence-bound physical slots
/// are both recoverable (one entry per distinct (value, Index) pair actually
/// sliced).
struct SlicedModeAssignment {
  /// Every DISTINCT DAG-scope loop realized anywhere in the schedule, in
  /// canonical (deterministic pre-order over the block tree) order; a
  /// value's assigned \c LoopId is its position in this list. The root block
  /// itself (sentinel axis, outside every loop) is NEVER an entry here.
  container::vector<DagScopeLevel> levels;

  /// value_id -> the list of (this value's OWN sliced-mode Index, the LoopId
  /// that slices it) pairs actually recorded for that value. A value with no
  /// entry here (or whose list omits some Index it carries) is NOT sliced by
  /// any loop on that mode -- participation is respected via the same
  /// built-within gate and ectx-participant test
  /// \c compute_sliced_mode_assignment implements below (a value built
  /// INSIDE a loop, or one that does not itself carry that loop's occurrence
  /// label, is left unstamped).
  std::unordered_map<std::size_t, container::svector<std::pair<Index, LoopId>>>
      by_value;

  /// CONSUMER-attributed per-occurrence sliced-mode facts (sliced-value
  /// canonical-layout / loop-coloring design, PILLAR 2): each entry is
  /// (value_id, this occurrence's own sliced-mode Index, the slicing LoopId,
  /// the CONSUMER value_id -- the use-site whose fetch of \c value_id binds the
  /// loop to that Index). Recorded ONLY by the regime-2 (occurrence-driven)
  /// pass, the one pass that can attribute a stamp to a specific occurrence and
  /// hence to a specific consumer. This is the raw material the executor
  /// projects onto \c LoopColoredSliceSeam::by_hash_consumer to disambiguate
  /// the w8-symmetric case (one value, one loop, two free modes bound by two
  /// different consumers); \c by_value alone cannot express "pos0 here, pos1
  /// there" because it folds away which occurrence bound which mode.
  /// (value_id, this occurrence's own sliced-mode PHYSICAL POSITION -- its
  /// index in occ.carried, computed in THAT occurrence's own index-frame,
  /// LoopId, CONSUMER value_id). The position (not an Index label) is what the
  /// executor projects onto \c LoopColoredSliceSeam::by_hash_consumer, so the
  /// runtime never re-matches a label across index-frames.
  container::svector<std::tuple<std::size_t, std::size_t, LoopId, std::size_t>>
      occ_facts;

  /// \return the \c LoopId slicing \p value_id's own \p mode, or \c
  /// std::nullopt if \p mode is not one of \p value_id's sliced modes.
  [[nodiscard]] std::optional<LoopId> loop_of(std::size_t value_id,
                                              Index const& mode) const {
    auto const it = by_value.find(value_id);
    if (it == by_value.end()) return std::nullopt;
    for (auto const& [ix, lid] : it->second)
      if (ix == mode) return lid;
    return std::nullopt;
  }

  /// \return the realized \c DagScopeLevel a \p loop_id names (the inverse of
  /// the canonical enumeration \c levels holds).
  [[nodiscard]] DagScopeLevel const& level_of(LoopId loop_id) const {
    return levels.at(loop_id);
  }
};

namespace detail {

///
/// \brief Pre-order walk collecting every non-root \c ScopeBlock's \c
/// DagScopeLevel into \p out, in schedule (structural) order -- the canonical
/// enumeration \c SlicedModeAssignment::levels holds. Deterministic given a
/// fixed \p block (the same \c OrderedSchedule always yields the same
/// sequence), and duplicate-free by construction: distinct \c ScopeBlock
/// objects in the tree are, per \c assert_global_level_axis_uniqueness's
/// invariant (consulted by \c compute_sliced_mode_assignment before this
/// runs), each other's only witness for a given \c (depth, space, ordinal) --
/// i.e. no two DIFFERENT blocks visited here ever carry an equal \c level.
///
inline void enumerate_realized_levels(ScopeBlock const& block,
                                      container::vector<DagScopeLevel>& out) {
  for (Step const& step : block.steps) {
    auto const* child = std::get_if<ScopeBlock>(&step.value);
    if (!child) continue;
    out.push_back(child->level);
    enumerate_realized_levels(*child, out);
  }
}

}  // namespace detail

///
/// \brief Task 3 (sliced-value canonical layout / loop-coloring design):
/// build the per-(value, sliced-mode) -> DAG-scope-loop \c
/// SlicedModeAssignment from an already-built \p ordered (must be \c
/// build_ordered_schedule(rich, ...)'s return value for this SAME \p rich).
///
/// \details Uses a two-pass fetch-site walk (\c populate_build_scope_walk
/// for the enclosing-loop scope of every produced/escaped value, and \c
/// ordered_schedule_dep_graph for the operand edges, leaves included): an
/// EXACT pass matching a block's representative axis against a value's own
/// carried \c Index (step (3) below), then a REGIME-2 RELABEL pass (step (4)
/// below) recovering a CSE value's own label from its occurrences when it
/// was canonicalized independently of the block's representative axis (the
/// SAME physical loop, relabeled) -- but records the raw facts keyed by the
/// value's OWN \c Index label (not a \c ValueCell::carried POSITION),
/// consistency-checked the same way (two fetch sites disagreeing on one
/// value's mode's level is a scheduler bug, never resolved by averaging or
/// last-write-wins), then remapped through the canonical \c LoopId
/// enumeration (\c detail::enumerate_realized_levels) instead of storing the
/// \c DagScopeLevel directly -- this is what makes forced-split siblings
/// (same depth/space, different ordinal) come out as DISTINCT colors: they
/// are distinct entries in the canonical \c levels list by construction.
///
[[nodiscard]] inline SlicedModeAssignment compute_sliced_mode_assignment(
    OrderedSchedule const& ordered, RichSchedule const& rich) {
  // Debug safety net: the global (depth, space, ordinal) -> axis uniqueness
  // this construction (and the canonical level enumeration below) relies on.
  {
    std::map<std::tuple<std::size_t, std::wstring, int, int>, Index> seen;
    detail::assert_global_level_axis_uniqueness(ordered.root, seen);
  }

  SlicedModeAssignment result;
  detail::enumerate_realized_levels(ordered.root, result.levels);

  std::map<std::tuple<std::size_t, std::wstring, int, int>, LoopId> level_id;
  for (std::size_t i = 0; i < result.levels.size(); ++i) {
    DagScopeLevel const& L = result.levels[i];
    level_id.emplace(
        std::make_tuple(L.depth, L.space, L.loop_slot, L.latitude_ordinal), i);
  }
  auto const id_of = [&](DagScopeLevel const& L) -> LoopId {
    auto const it = level_id.find(
        std::make_tuple(L.depth, L.space, L.loop_slot, L.latitude_ordinal));
    SEQUANT_ASSERT(it != level_id.end());
    return it->second;
  };

  // (1) fetch/eval scope of every produced or escaped value: enclosing blocks
  // (axis + level), read off the realized block tree.
  std::unordered_map<std::size_t,
                     container::svector<detail::ScopeBlockAxisLevel>>
      build_scope;
  detail::populate_build_scope_walk(ordered.root, {}, build_scope);

  // (2) operand edges (leaves included).
  detail::OrderedScheduleDepGraph const g =
      detail::ordered_schedule_dep_graph(rich);

  // Consistency-checked write of one (value, OWN Index) -> level fact: two
  // fetch sites disagreeing on the level for the SAME (value, mode) pair is a
  // scheduler bug (the value should have been split), enforced by this
  // \c stamp_raw assert. Recorded as raw \c DagScopeLevel facts first (this
  // is the datum that must agree across fetch sites); the fold to \c LoopId
  // happens once, at the end.
  std::unordered_map<std::size_t,
                     container::svector<std::pair<Index, DagScopeLevel>>>
      raw_facts;
  auto const stamp_raw = [&raw_facts](std::size_t w_vid, Index const& mode,
                                      DagScopeLevel const& level) {
    auto& facts = raw_facts[w_vid];
    for (auto& [ix, lvl] : facts) {
      if (!(ix == mode)) continue;
      SEQUANT_ASSERT(
          lvl == level &&
          "compute_sliced_mode_assignment: divergent (value, mode) -> level "
          "across fetch sites -- the value should have been split");
      return;
    }
    facts.push_back({mode, level});
  };

  // (3) PER-OCCURRENCE POSITIONAL pass (2026-08-25 loop-open design; replaces
  // the former EXACT + REGIME-2 base_key passes). For each occurrence occ of a
  // value W consumed by C, the loops the runtime crosses when fetching W are
  // C's enclosing DAG blocks build_scope[owner(consumer_point)], outermost
  // first. occ.ectx -- now built from loop-OPENS (peak_profile), so one
  // physical loop appears exactly once -- names those same loops in W's OWN
  // frame, outermost first. Pairing them by NEST POSITION gives, per realized
  // enclosing loop, the occ-frame mode it slices; the physical slice position
  // is that mode's index in occ.carried. Everything is in occ's own frame: no
  // base_key, no cross-frame label match, no first-match guess. Divergent
  // (relabeled) occurrences and the symmetric case (one value sliced on
  // different positions by different consumers) are handled uniformly by the
  // consumer-keyed occ_facts. The old raw_facts / by_value path is left empty:
  // the ordered executor always fetches under a tracked consumer, so the
  // consumer-keyed facts cover every sliced fetch, and a value sliced on one
  // mode by two sibling loops via two consumers has no single consumer-blind
  // answer anyway.
  (void)g;          // dependency graph no longer consulted by this pass
  (void)stamp_raw;  // by_value fallback intentionally left empty (see above)
  std::unordered_map<std::size_t, std::size_t> point_owner;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences)
      point_owner[occ.point] = vc.value_id;

  // The VALUE-LEVEL component slot per carried position: the loop_slot the
  // union-find (Task 2) assigned to (value, position), recovered UNMASKED. Task
  // 2 stamps occ.loop_slot = -1 wherever a mode is not batched AT THAT
  // OCCURRENCE (i.e. the value was produced whole in that mode in that term).
  // But a value produced whole is still READ SLICED wherever a consumer sits
  // inside that mode's realized loop -- the atlas must key the slice fact off
  // which loops the consumer is IN and which modes the value CARRIES, not off
  // the per-occurrence production mask. So recover the value-level slot as any
  // non-(-1) loop_slot across the value's occurrences at each position.
  std::unordered_map<std::size_t, container::svector<int>> value_slot;
  for (ValueCell const& vc : rich.cells) {
    container::svector<int>& vs = value_slot[vc.value_id];
    vs.assign(vc.carried.size(), -1);
    for (OccurrenceRec const& occ : vc.occurrences)
      for (std::size_t p = 0; p < occ.loop_slot.size() && p < vs.size(); ++p)
        if (occ.loop_slot[p] >= 0) vs[p] = occ.loop_slot[p];
  }

  for (ValueCell const& w : rich.cells) {
    std::size_t const w_vid = w.value_id;
    for (OccurrenceRec const& occ : w.occurrences) {
      if (occ.consumer_point == occ.point) continue;  // forest root
      auto const oit = point_owner.find(occ.consumer_point);
      if (oit == point_owner.end()) continue;  // defensive
      auto const sit = build_scope.find(oit->second);
      if (sit == build_scope.end()) continue;  // consumer has no enclosing loop
      container::svector<detail::ScopeBlockAxisLevel> const& scope =
          sit->second;
      if (scope.empty()) continue;
      // Task 3 atlas: for each of the consumer's enclosing loops `scope[k]`
      // (space S, loop_slot L), the mode of W it slices is W's carried position
      // whose VALUE-LEVEL component slot is L and whose space is S. A LOOP-
      // IDENTITY match in W's own frame (occurrence-invariant), keyed off the
      // value-level slot -- so a value produced whole in a mode is still sliced
      // there when THIS consumer sits inside that mode's realized loop (the
      // measured i_3-contracted case). Over-nested loops the value carries no
      // mode of get no fact, correct (the value is invariant to them). Replaces
      // the old raw ectx<->scope POSITIONAL zip (which mis-sliced divergent
      // occurrences: the multi-occ collapse / ToTxToT deadlock).
      auto const& vs = value_slot[w_vid];
      for (std::size_t k = 0; k < scope.size(); ++k) {
        DagScopeLevel const& lvl = scope[k].level;
        for (std::size_t pos = 0; pos < occ.carried.size() && pos < vs.size();
             ++pos) {
          if (vs[pos] != lvl.loop_slot) continue;
          if (std::wstring(occ.carried[pos].space().base_key()) != lvl.space)
            continue;
          result.occ_facts.push_back(
              std::make_tuple(w_vid, pos, id_of(lvl), oit->second));
          break;  // one W-position per enclosing loop
        }
      }
    }
  }

  // Fold raw (value, Index) -> DagScopeLevel facts through the canonical
  // LoopId enumeration.
  for (auto const& [vid, facts] : raw_facts) {
    auto& out_list = result.by_value[vid];
    out_list.clear();
    for (auto const& [ix, level] : facts)
      out_list.push_back({ix, id_of(level)});
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
