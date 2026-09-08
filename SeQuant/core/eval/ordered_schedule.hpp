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
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstdlib>
#include <iostream>

#include <algorithm>
#include <climits>
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
  /// Per value_id, the value_ids of its DIRECT operands -- the value/
  /// occurrence DAG edges the value-driven ordered executor consumes to fetch
  /// each operand by its OWN cell id (see \c CellTable / \c CellRegistry).
  /// Recorded here from `ordered_schedule_dep_graph(rich).depends_on`, whose
  /// edges come from every `OccurrenceRec`'s `consumer_point` (so split
  /// operands resolve to the specific consumed value, not an ambiguous node
  /// hash). A leaf value (no operands) has no entry.
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
/// \brief Collect every \c BuildStep site (value_id + its root-to-block path)
/// into \p builds and every output escape site into \p sites. Both carry the
/// path: \c well_formed's built-and-escaped rule compares WHERE a value is
/// built against where it escapes, not merely whether it does both.
///
inline void collect_productions(ScopeBlock const& block,
                                container::svector<ScopeBlock const*>& path,
                                container::vector<OutputSite>& builds,
                                container::vector<OutputSite>& sites) {
  path.push_back(&block);
  for (Step const& step : block.steps) {
    if (auto const* b = std::get_if<BuildStep>(&step.value))
      builds.push_back(OutputSite{b->value_id, path});
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
///     no value_id is built (\c BuildStep) more than once; a built value_id may
///     ALSO escape only through its OWN chain -- every block listing it in \c
///     outputs either holds its \c BuildStep or is an ancestor of the block
///     that does (a member materialized across a forced loop split is built in
///     its home block for its in-nest readers and escapes from that same block
///     outward; see \c build_ordered_schedule's mixed-pass rule) -- and any
///     other combination of build and escape sites is duplicate production;
///     and a value_id may escape (\c outputs) at MORE than
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
  container::vector<detail::OutputSite> builds;
  container::vector<detail::OutputSite> sites;
  detail::collect_productions(sched.root, path, builds, sites);

  // (a) no value_id built more than once.
  {
    container::vector<std::size_t> b;
    for (auto const& s : builds) b.push_back(s.value_id);
    std::sort(b.begin(), b.end());
    if (auto const it = std::adjacent_find(b.begin(), b.end()); it != b.end()) {
      if (std::getenv("SEQUANT_DUMP_WF"))
        std::cerr << "[wf-fail] double-build vid=" << *it << std::endl;
      return false;
    }
  }
  // (b) a built value_id may ALSO escape, but only through its OWN chain:
  // every block listing it in `outputs` must either BE the block holding its
  // BuildStep (built and escaped in one block -- a mixed-pass member of a
  // forced split: its in-nest readers take the per-batch cell, the other pass
  // takes the assembled form) or an ANCESTOR of it (the chain levels above the
  // home block). An escape in a sibling or descendant of the home block is a
  // second, unrelated producer.
  {
    std::unordered_map<std::size_t, detail::OutputSite const*> build_of;
    for (auto const& b : builds) build_of.emplace(b.value_id, &b);
    for (auto const& s : sites) {
      auto const it = build_of.find(s.value_id);
      if (it == build_of.end()) continue;
      auto const& home = it->second->path;
      bool ok = s.path.size() <= home.size();
      for (std::size_t k = 0; ok && k < s.path.size(); ++k)
        ok = s.path[k] == home[k];
      if (!ok) {
        if (std::getenv("SEQUANT_DUMP_WF"))
          std::cerr << "[wf-fail] built-and-escaped vid=" << s.value_id
                    << " escape depth=" << s.path.size()
                    << " home depth=" << home.size() << std::endl;
        return false;
      }
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
  // DOWNWARD closure (2026-09-02): a LoopLocal member of this nest whose
  // consumers ALL sit in the consumer pass must move with them -- its
  // per-batch cell is pass-local, so left in the producer pass it is gone by
  // the time the consumer pass runs (w20 pVDZ-F12 strict walk: 59557, a
  // 4-occ member read only by two consumer-pass members at i#d1s0 lat 1 ->
  // "read-from-home value vanished"). A member with consumers in BOTH passes
  // is left where it is (it cannot be built twice); such a case must escape
  // instead and shows up in the strict walk.
  for (bool grew = true; grew;) {
    grew = false;
    container::svector<std::size_t> members(r.consumer_pass.begin(),
                                            r.consumer_pass.end());
    for (std::size_t v : members) {
      auto const dit = g.depends_on.find(v);
      if (dit == g.depends_on.end()) continue;
      for (std::size_t op : dit->second) {
        if (r.consumer_pass.count(op) || r.carried.count(op)) continue;
        auto const cit = g.consumers_of.find(op);
        if (cit == g.consumers_of.end() || cit->second.empty()) continue;
        bool all_in = true;
        for (std::size_t c : cit->second)
          if (!r.consumer_pass.count(c)) {
            all_in = false;
            break;
          }
        // All consumers on the consumer side: move it. Consumers on BOTH
        // sides: left where it is -- a per-batch cell is pass-local, so this
        // is a genuine scheduling conflict (recompute it in the consumer pass
        // vs materialize it across the pass boundary) that the schedule must
        // decide with a cost; it surfaces in the strict walk, by design.
        if (all_in && r.consumer_pass.insert(op).second) grew = true;
      }
    }
  }
  return r;
}

///
/// \brief Pass levels for the one forced-split axis space \p axis_key: every
/// value gets an integer pass such that a value reading a carried value's
/// completed (full) form sits in a later pass than the carried value.
///
/// Forward sweep (operands before consumers):
///   level(v) = max over direct operands o of (o carried ? pass(o) + 1
///              : level(o)), 0 with no operands;
///   a carried value's pass is its level; a non-carried value's BASE is its
///   level.
/// Reverse sweep (consumers before operands): a non-carried value with at
/// least one consumer is LIFTED to max(base, min over its direct consumers'
/// passes), so a value whose readers all sit later is built with them; a
/// value with readers in several passes keeps its base (and is materialized
/// by the builder's rule 4 when a later same-nest reader needs it).
///
/// Every dependency edge points to an equal or earlier pass. With only passes
/// 0 and 1 present this is exactly the former two-set partition.
///
struct ForcedSplitLevels {
  std::unordered_set<std::size_t> carried;       //!< LoopCarried-on-axis ids
  std::unordered_map<std::size_t, int> pass_of;  //!< value id -> pass
  int max_pass = 0;
  int pass(std::size_t vid) const {
    auto const it = pass_of.find(vid);
    return it == pass_of.end() ? 0 : it->second;
  }
};

inline ForcedSplitLevels forced_split_levels(std::wstring const& axis_key,
                                             RichSchedule const& rich,
                                             LegalitySchedule const& legality,
                                             OrderedScheduleDepGraph const& g) {
  ForcedSplitLevels r;
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

  std::size_t const n = rich.cells.size();
  // Topological order, operands before consumers (Kahn over depends_on).
  container::svector<std::size_t> topo;
  topo.reserve(n);
  {
    container::svector<std::size_t> indeg(n, 0);
    for (std::size_t v = 0; v < n; ++v) {
      auto const it = g.depends_on.find(v);
      if (it != g.depends_on.end()) indeg[v] = it->second.size();
    }
    container::svector<std::size_t> ready;
    for (std::size_t v = 0; v < n; ++v)
      if (indeg[v] == 0) ready.push_back(v);
    while (!ready.empty()) {
      std::size_t const v = ready.back();
      ready.pop_back();
      topo.push_back(v);
      auto const it = g.consumers_of.find(v);
      if (it == g.consumers_of.end()) continue;
      for (std::size_t u : it->second)
        if (--indeg[u] == 0) ready.push_back(u);
    }
    SEQUANT_ASSERT(topo.size() == n,
                   "forced_split_levels: the value dependency graph has a "
                   "cycle");
  }

  std::unordered_map<std::size_t, int> base;
  base.reserve(n);
  for (std::size_t v : topo) {
    int lv = 0;
    auto const it = g.depends_on.find(v);
    if (it != g.depends_on.end())
      for (std::size_t o : it->second)
        lv = std::max(lv, r.carried.count(o) ? base[o] + 1 : base[o]);
    base[v] = lv;
  }

  r.pass_of = base;
  for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
    std::size_t const v = *it;
    if (r.carried.count(v)) continue;
    auto const cit = g.consumers_of.find(v);
    if (cit == g.consumers_of.end() || cit->second.empty()) continue;
    int mn = INT_MAX;
    for (std::size_t u : cit->second) mn = std::min(mn, r.pass_of[u]);
    r.pass_of[v] = std::max(base[v], mn);
  }
  for (auto const& [v, p] : r.pass_of) r.max_pass = std::max(r.max_pass, p);
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
/// \note NOT LIVE on any schedule: \c analyze_legality takes its \c
/// DemotionSource as an optional argument and is inert without one (\c
/// derive_demotions returns immediately), and no caller passes one -- so this
/// function only runs in its own unit tests. The schedule builder applies the
/// same classification DIRECTLY, in \c build_ordered_schedule's per-value
/// placement loop (its "MIXED-PASS MEMBER" rule, section 2), where the pass
/// partition and the realized loop depths are both in hand; see there for the
/// live behavior.
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
  std::map<std::size_t, std::size_t> cluster_min;  // cluster -> outermost depth
  container::svector<std::size_t> order;
  {
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

  // 2a. PASS LEVELS (per-nest forced split design, section 3.1). The one
  // forced space is the space some cell is LoopCarried on; a second forced
  // space is rejected as before. Levels are global; each nest decides for
  // itself (step 4 below) whether it holds members of more than one pass.
  std::optional<detail::ForcedSplitLevels> levels;
  {
    container::svector<std::wstring> forced_spaces;
    for (CellLegality const& cl : legality.cells)
      for (Index const& ix : forced_split_types(cl)) {
        std::wstring const key{ix.space().base_key()};
        if (std::find(forced_spaces.begin(), forced_spaces.end(), key) ==
            forced_spaces.end())
          forced_spaces.push_back(key);
      }
    SEQUANT_ASSERT(forced_spaces.size() <= 1,
                   "build_ordered_schedule: more than one forced-split space "
                   "is not supported");
    if (!forced_spaces.empty())
      levels =
          detail::forced_split_levels(forced_spaces.front(), rich, legality, g);
  }
  auto const pass_of = [&](std::size_t vid) -> int {
    return levels ? levels->pass(vid) : 0;
  };
  // Legality record by value id, for reader lookups in rule 4.
  std::unordered_map<std::size_t, CellLegality const*> cl_by_vid;
  cl_by_vid.reserve(legality.cells.size());
  for (CellLegality const& c2 : legality.cells) {
    auto const it2 = value_id_of.find(c2.hash);
    if (it2 != value_id_of.end()) cl_by_vid.emplace(it2->second, &c2);
  }

  // 2. Per-value placement: home BuildStep (root-level bucket uses index n
  // as a sentinel "root" depth) vs. escape output at the deepest escaped
  // axis's depth.
  container::vector<detail::OrderedScheduleDepthBucket> buckets(n);
  container::svector<std::size_t> root_build_ids;
  // Values that keep their BuildStep AND escape it (the mixed-pass members of
  // a forced split, below): dump-only diagnostic (SEQUANT_DUMP_SCHEDULE) --
  // the live signal downstream is the per-value \c materialized_across_split
  // bool below, not this list.
  container::svector<std::size_t> materialized_across_split_ids;

  // The home depth of a value's plain BuildStep: the INNERMOST loop it is
  // LoopLocal on, resolved PER-INSTANCE by fusion loop_slot -- NOT by
  // shallowest same-space count. A value local to slots 2,3 (its own fusion
  // nest) homes in THAT nest, not the FIRST same-space nest a space-multiset
  // cover would pick; picking the wrong same-space nest homes the value where
  // its consumer's nest has not opened (or has already closed), so an
  // in-consumer-nest read misses and the value vanishes. Resolve each
  // LoopLocal mode's (space, fusion slot) to its realized depth (exactly as
  // the escape placement does), and home at the MAX such depth: within a
  // co-occurrence cluster larger depth nests inside smaller, so the innermost
  // of the value's own home slots is inside all of them. home_floor is the
  // LoopLocal subset, but it drops the pos->slot map, so walk per_axis
  // directly for the slot. Nullopt = no realized loop-local mode: root.
  auto const local_home_depth =
      [&](CellLegality const& cl) -> std::optional<std::size_t> {
    std::optional<std::size_t> target;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      if (cl.per_axis[pos].role != LoopRole::LoopLocal) continue;
      std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
      int const fs = fusion_slot(cl, pos);
      auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
      if (!d) continue;
      if (!target || *d > *target) target = *d;
    }
    return target;
  };

  // Dump-only diagnostic (SEQUANT_DUMP_SCHEDULE) for a rule-4 rejection:
  // print the value's axes with roles/slots/depths, its later-pass readers,
  // and the forced-split axis's carried set, each with home depth and nest.
  auto const dump_reject = [&](std::size_t v0, std::size_t nest,
                               container::svector<std::size_t> const& readers) {
    auto const role_str = [](LoopRole r) -> wchar_t const* {
      switch (r) {
        case LoopRole::LoopLocal:
          return L"L";
        case LoopRole::Reduction:
          return L"R";
        case LoopRole::LoopCarried:
          return L"C";
        default:
          return L"I";
      }
    };
    auto const describe = [&](std::size_t v) {
      auto const it2 = cl_by_vid.find(v);
      std::wcerr << L"v" << v << L"(pass=" << pass_of(v);
      if (it2 == cl_by_vid.end()) {
        std::wcerr << L" no legality)";
        return;
      }
      CellLegality const& c2 = *it2->second;
      auto const hd = local_home_depth(c2);
      std::wcerr << L" home=";
      if (hd)
        std::wcerr << *hd << L"/nest" << type_cluster[*hd];
      else
        std::wcerr << L"root";
      std::wcerr << L" axes={";
      for (std::size_t p2 = 0; p2 < c2.per_axis.size(); ++p2) {
        std::wstring const bk2{c2.per_axis[p2].axis.space().base_key()};
        int const fs2 = fusion_slot(c2, p2);
        auto const d2 = depth_of_instance(bk2, fs2 >= 0 ? fs2 : 0);
        std::wcerr << c2.per_axis[p2].axis.full_label() << L":"
                   << role_str(c2.per_axis[p2].role) << L"@slot" << fs2
                   << L"->d";
        if (d2)
          std::wcerr << *d2;
        else
          std::wcerr << L"?";
        std::wcerr << L" ";
      }
      std::wcerr << L"}" << (levels->carried.count(v) ? L" CARRIED" : L"")
                 << L")";
    };
    std::wcerr << L"[sched-reject] nest outermost depth "
               << cluster_min.at(nest) << L" value: ";
    describe(v0);
    std::wcerr << L"\n";
    for (std::size_t u : readers) {
      std::wcerr << L"[sched-reject]   later-pass reader ";
      describe(u);
      std::wcerr << L"\n";
    }
    for (std::size_t c : levels->carried) {
      std::wcerr << L"[sched-reject]   carried ";
      describe(c);
      std::wcerr << L"\n";
    }
  };

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

    // MIXED-PASS MEMBER (per-nest forced split design, section 3.3). A
    // same-nest, later-pass reader is the only reader that can see a
    // per-batch slice from another traversal of the nest: a reader outside
    // the nest or at root is never LoopLocal on the value and is served by
    // the value's role-driven escapes. Such a member keeps its BuildStep and
    // is scattered to full at every level of its nest it is LoopLocal on,
    // from its home out to the nest's outermost depth, so the later pass
    // reads the full form at root. The chain is measured in the value's OWN
    // nest; a chain that cannot reach root is rejected loudly (below).
    bool materialized_across_split = false;
    std::optional<std::size_t> const home_depth = local_home_depth(cl);
    // Direct readers homed in the same nest with a later pass.
    auto const later_same_nest_readers =
        [&](std::size_t nest) -> container::svector<std::size_t> {
      container::svector<std::size_t> out;
      auto const cons_it = g.consumers_of.find(vid);
      if (cons_it == g.consumers_of.end()) return out;
      for (std::size_t u : cons_it->second) {
        if (pass_of(u) <= pass_of(vid)) continue;
        auto const uit = cl_by_vid.find(u);
        if (uit == cl_by_vid.end()) continue;
        auto const uh = local_home_depth(*uit->second);
        if (uh && type_cluster[*uh] == nest) out.push_back(u);
      }
      return out;
    };
    if (escapes.empty() && levels && home_depth) {
      std::size_t const nest = type_cluster[*home_depth];
      auto const readers = later_same_nest_readers(nest);
      if (!readers.empty()) {
        for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
          if (cl.per_axis[pos].role != LoopRole::LoopLocal) continue;
          std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
          int const fs = fusion_slot(cl, pos);
          auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
          if (!d || type_cluster[*d] != nest || *d > *home_depth) continue;
          auto it = std::find_if(escapes.begin(), escapes.end(),
                                 [&](auto const& e) { return e.first == *d; });
          if (it == escapes.end())
            escapes.push_back({*d, OutputKind::AccumulateScatter});
          else
            it->second = OutputKind::AccumulateScatter;
          materialized_across_split = true;
        }
        // Every level of the nest from the home out to the outermost depth
        // must be escaped; a level the value is invariant to would end the
        // chain inside the nest, invisible to the later pass.
        std::size_t const outermost = cluster_min.at(nest);
        for (std::size_t d = outermost; d <= *home_depth; ++d) {
          if (type_cluster[d] != nest) continue;
          bool const escaped =
              std::any_of(escapes.begin(), escapes.end(),
                          [&](auto const& e) { return e.first == d; });
          if (!escaped) {
            if (std::getenv("SEQUANT_DUMP_SCHEDULE"))
              dump_reject(vid, nest, readers);
            throw Exception(
                "build_ordered_schedule: value " + std::to_string(vid) +
                " (pass " + std::to_string(pass_of(vid)) +
                ") is read by value " + std::to_string(readers.front()) +
                " (pass " + std::to_string(pass_of(readers.front())) +
                ") in a later pass of its nest (outermost depth " +
                std::to_string(outermost) +
                ") but is invariant to that nest's loop at depth " +
                std::to_string(d) +
                ": its escape chain cannot reach root (unsupported: hoisting "
                "an invariant assembled form)");
          }
        }
      }
    } else if (!escapes.empty() && levels && home_depth) {
      // A role-escaping value whose chain ends INSIDE its nest (its
      // outermost escaped level is not the nest's outermost depth) with a
      // same-nest later-pass reader: the assembled form is invisible to that
      // pass. Rejected loudly (parent design section 8.2's parked case).
      std::size_t const nest = type_cluster[*home_depth];
      std::size_t d_min = escapes.front().first;
      for (auto const& [d, kind] : escapes) {
        (void)kind;
        d_min = std::min(d_min, d);
      }
      if (type_cluster[d_min] == nest && d_min != cluster_min.at(nest)) {
        auto const readers = later_same_nest_readers(nest);
        if (!readers.empty()) {
          if (std::getenv("SEQUANT_DUMP_SCHEDULE"))
            dump_reject(vid, nest, readers);
          throw Exception(
              "build_ordered_schedule: value " + std::to_string(vid) +
              " (pass " + std::to_string(pass_of(vid)) +
              ") escapes only to depth " + std::to_string(d_min) +
              " inside its nest (outermost depth " +
              std::to_string(cluster_min.at(nest)) + ") but is read by value " +
              std::to_string(readers.front()) + " (pass " +
              std::to_string(pass_of(readers.front())) +
              ") in a later pass of that nest (unsupported: an inner "
              "escape read across passes)");
        }
      }
    }

    if (!escapes.empty()) {
      for (auto const& [d, kind] : escapes)
        buckets[d].outputs.push_back({vid, kind});
      // A value that escapes by its OWN per-axis roles has no BuildStep: its
      // production is the accumulation itself. One materialized across the
      // split by the rule above is different -- it is an ordinary LoopLocal
      // value that ALSO has producer-side consumers, which read it inside its
      // own nest, per batch. It keeps its BuildStep at its home block, so that
      // block both BUILDS it (for its producer-side readers, and as the
      // per-batch input of its own escape) and lists it as an output (for the
      // consumer pass). `well_formed` admits exactly this shape: every block
      // that lists it either holds its BuildStep or is an ancestor of the one
      // that does.
      if (!materialized_across_split) continue;
      materialized_across_split_ids.push_back(vid);
    }

    // Plain BuildStep: home at the INNERMOST loop the value is LoopLocal on
    // (\c local_home_depth, computed above as \c home_depth).
    if (home_depth)
      buckets[*home_depth].build_ids.push_back(vid);
    else
      root_build_ids.push_back(vid);
  }

  if (std::getenv("SEQUANT_DUMP_SCHEDULE") &&
      !materialized_across_split_ids.empty()) {
    std::wcerr << L"[sched-materialize] "
               << materialized_across_split_ids.size()
               << L" member(s) built AND escaped across the forced split:";
    for (std::size_t v : materialized_across_split_ids)
      std::wcerr << L" v" << v;
    std::wcerr << L"\n";
  }

  auto const union_into = [](container::svector<std::size_t>& acc,
                             container::svector<std::size_t> const& add) {
    for (std::size_t v : add)
      if (std::find(acc.begin(), acc.end(), v) == acc.end()) acc.push_back(v);
  };

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

  // Per-nest pass sets (design section 3.2): the passes of every build
  // homed at any of the nest's depths and of every escape output listed
  // there. A nest with one pass is one block (latitude = that pass); a nest
  // with several is one block per pass at its outermost depth.
  std::map<std::size_t, std::set<int>> nest_passes;
  for (std::size_t d = 0; d < n; ++d) {
    auto& ps = nest_passes[type_cluster[d]];
    for (std::size_t v : buckets[d].build_ids) ps.insert(pass_of(v));
    for (auto const& o : buckets[d].outputs) ps.insert(pass_of(o.first));
  }
  if (std::getenv("SEQUANT_DUMP_SCHEDULE"))
    for (auto const& [c, ps] : nest_passes) {
      std::wcerr << L"[sched-nest] outermost depth " << cluster_min.at(c)
                 << L" passes={";
      for (int k : ps) {
        std::size_t cnt = 0;
        for (std::size_t d = 0; d < n; ++d) {
          if (type_cluster[d] != c) continue;
          for (std::size_t v : buckets[d].build_ids)
            if (pass_of(v) == k) ++cnt;
        }
        std::wcerr << k << L":" << cnt << L" ";
      }
      std::wcerr << L"}\n";
    }

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

    std::size_t const nest = type_cluster[d];
    bool const outermost = d == cluster_min.at(nest);
    std::set<int> const& passes = nest_passes[nest];
    if (outermost && passes.size() > 1) {
      // Several passes at this nest's outermost depth (per-nest forced split
      // design, section 3.2): one block per pass, ascending, each holding the
      // builds and outputs of that pass at this depth plus the inner
      // sub-chain forked for that pass (fork_subchain applied once per pass;
      // the fork's predicate-true "consumer" side is the pass's own steps).
      // At the innermost axis pending is empty and fork_subchain returns an
      // empty inner list, reducing to a childless per-pass block emission.
      // The later pass's block reads the earlier pass's escaped (now-full)
      // outputs, so its `requires_` names them and the outer topo-sort orders
      // passes in ascending order.

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

      for (int k : passes) {
        auto const in_pass = [&](std::size_t v) { return pass_of(v) == k; };
        detail::ForkedSubchain forked =
            detail::fork_subchain(pending_steps, in_pass);
        // The "consumer" side of the two-way fork is the predicate-true side.
        container::svector<std::size_t> builds;
        for (std::size_t v : bucket.build_ids)
          if (in_pass(v)) builds.push_back(v);
        container::svector<std::pair<std::size_t, OutputKind>> outs;
        for (auto const& o : bucket.outputs)
          if (in_pass(o.first)) outs.push_back(o);
        auto metas = meta_for(forked.consumer);
        emit_pass(k, builds, outs, std::move(forked.consumer),
                  std::move(metas));
      }
    } else {
      ScopeBlock block = make_block(
          types[d], outermost && !passes.empty() ? *passes.begin() : 0, d + 1,
          type_slot[d], bucket.build_ids, bucket.outputs,
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

  // Persist the value/occurrence DAG edges (each value's direct operand
  // value_ids) the value-driven ordered executor consumes to fetch each
  // operand by its own cell id. `g` is the same dep graph the topo-sort used
  // above; its `depends_on` edges are derived from every OccurrenceRec's
  // consumer_point, so a split operand resolves to the specific consumed
  // value (not an ambiguous node hash).
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
/// inside. This is the build-scope walk itself, keeping each block's \c
/// DagScopeLevel (not just its axis) so the map can name the runtime \c
/// BatchContextEntry::level.
///
/// \p enc is the root-to-\p block path INCLUDING \p block's own (axis, level)
/// -- a value built or escaping inside \p block reads its operands inside
/// \p block's own loop, so the block's own axis encloses that read (this walk
/// always includes the block axis; only the HOME of an escape sits one level
/// out, which is irrelevant here -- we want the read/fetch scope, not the
/// home).
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
/// \note Defined in \c dag_scope.hpp so the low-level DAG-scope types can be
/// named without depending on this schedule header; re-exported here for the
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

  /// CONSUMER-attributed per-occurrence sliced-mode facts (sliced-value
  /// canonical-layout / loop-coloring design, PILLAR 2): each entry is
  /// (value_id, this occurrence's own sliced-mode Index, the slicing LoopId,
  /// the CONSUMER value_id -- the use-site whose fetch of \c value_id binds the
  /// loop to that Index). Recorded ONLY by the regime-2 (occurrence-driven)
  /// pass, the one pass that can attribute a stamp to a specific occurrence and
  /// hence to a specific consumer. This is the raw material the cell table
  /// builder (cell_table_builder.hpp) consumes to disambiguate the
  /// w8-symmetric case (one value, one loop, two free modes bound by two
  /// different consumers): a consumer-blind (value, mode) map cannot express
  /// "pos0 here, pos1 there" because it folds away which occurrence bound
  /// which mode.
  /// (value_id, this occurrence's own sliced-mode PHYSICAL POSITION -- its
  /// index in occ.carried, computed in THAT occurrence's own index-frame,
  /// LoopId, CONSUMER value_id). The position (not an Index label) is what the
  /// table builder consumes, so the runtime never re-matches a label across
  /// index-frames.
  container::svector<std::tuple<std::size_t, std::size_t, LoopId, std::size_t>>
      occ_facts;

  /// EXPLICIT per-occurrence INVARIANT facts: (value_id, LoopId, CONSUMER
  /// value_id) triples recording that this consumer's fetch of the value is
  /// correctly UNSLICED on this loop -- the consumer is building that loop's
  /// batch, but the value's occurrence IN THIS CONSUMER'S FRAME does not carry
  /// the loop's mode (a CSE-shared value can carry the mode in one frame and
  /// not another). Recording the negative decision lets the runtime tell
  /// "correctly invariant" apart from "no decision recorded" (a real gap), so
  /// the completeness guard fires only on the latter.
  container::svector<std::tuple<std::size_t, LoopId, std::size_t>>
      occ_invariant;

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
  // consumer-keyed occ_facts. There is no consumer-blind (value, mode) ->
  // loop fallback: the ordered executor always fetches under a tracked
  // consumer, so the consumer-keyed facts cover every sliced fetch, and a
  // value sliced on one mode by two sibling loops via two consumers has no
  // single consumer-blind answer anyway.
  (void)g;  // dependency graph no longer consulted by this pass
  std::unordered_map<std::size_t, std::size_t> point_owner;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences)
      point_owner[occ.point] = vc.value_id;

  // point -> that occurrence record (to reach the PARENT occurrence on an
  // edge: the consumer's occurrence in the SAME tree as the operand's).
  std::unordered_map<std::size_t, OccurrenceRec const*> point_occ;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences) point_occ[occ.point] = &occ;

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
      // PRODUCTION EDGE ONLY (2026-09-02): a consumer value is PRODUCED once,
      // from its canonical (front) occurrence's tree; its other occurrences
      // are READS of the finished value, and the W->C edges inside those
      // trees never execute. Slice facts recorded from such a non-production
      // edge are keyed by the same (W, C, loop) triple and MERGE with the
      // production edge's facts -- in a permuted frame that binds C's modes
      // to other loop instances, which sliced one operand position by two
      // loops (w20 strict dry-run walk: 20632 pos 1 by slots 1 AND 2).
      if (rich.cells[oit->second].occurrences.empty() ||
          rich.cells[oit->second].occurrences.front().point !=
              occ.consumer_point)
        continue;
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
        bool self_sliced = false;
        for (std::size_t pos = 0; pos < occ.carried.size() && pos < vs.size();
             ++pos) {
          // PER-OCCURRENCE slot (2026-09-02): value_slot folds every
          // occurrence's loop_slot into ONE per-value vector, so a CSE-shared
          // value whose mode is sliced by DIFFERENT loop instances in
          // different occurrences (w20: 51337's K position under K-loop slot
          // 0 for one consumer, slot 1 / depth 43 for a K-reduction
          // consumer) keeps only one slot and silently misses the other
          // occurrence (served WHOLE under a K batch -> TA sparse-gemm
          // out-of-bounds read). The occurrence's own loop_slot is the
          // truth; the per-value slot is only the fallback where it is
          // unstamped (-1).
          int const occ_slot =
              (pos < occ.loop_slot.size() && occ.loop_slot[pos] >= 0)
                  ? occ.loop_slot[pos]
                  : vs[pos];
          if (occ_slot != lvl.loop_slot) continue;
          if (std::wstring(occ.carried[pos].space().base_key()) != lvl.space)
            continue;
          result.occ_facts.push_back(
              std::make_tuple(w_vid, pos, id_of(lvl), oit->second));
          self_sliced = true;
          break;  // one W-position per enclosing loop
        }
        {
          static long const _kc = [] {
            char const* dh = std::getenv("SEQUANT_DUMP_FALLBACK_C");
            return dh ? std::strtol(dh, nullptr, 10) : -1L;
          }();
          if (_kc >= 0 && (rich.cells[oit->second].hash % 100000u) ==
                              static_cast<unsigned>(_kc))
            std::cerr << "[kloop] W=" << (rich.cells[w_vid].hash % 100000)
                      << " C=" << (rich.cells[oit->second].hash % 100000)
                      << " L(space=" << toUtf8(lvl.space)
                      << " slot=" << lvl.loop_slot
                      << ") self_sliced=" << self_sliced << "\n";
        }
        if (self_sliced) continue;
        // Use-induced slicing (Layer 2): W is WHOLE-produced on this loop's
        // space (no own sliced slot), but THIS CONSUMER C is sliced on a mode M
        // there. If M is a SHARED external W also carries, the C = ...*W
        // contraction binds W's M to C's sliced M -- so W must be sliced on M
        // at this loop for THIS fetch. Record it CONSUMER-KEYED in occ_facts,
        // never consumer-blind: W may be
        // CSE-shared between C (sliced here) and a DIFFERENT consumer that
        // reads it whole (invariant), and a blind fact would wrongly slice it
        // for both. Bounded to a mode C actually slices (conformability), NOT
        // every carried-mode coincidence, so a genuinely invariant shared
        // operand is not sliced.
        std::size_t const c_vid = oit->second;
        auto const& cvs = value_slot[c_vid];
        auto const& c_carried = rich.cells[c_vid].carried;
        static long const _fb_target = [] {
          char const* dh = std::getenv("SEQUANT_DUMP_FALLBACK");
          return dh ? std::strtol(dh, nullptr, 10) : -1L;
        }();
        static long const _fb_ctarget = [] {
          char const* dh = std::getenv("SEQUANT_DUMP_FALLBACK_C");
          return dh ? std::strtol(dh, nullptr, 10) : -1L;
        }();
        bool const _fb =
            (_fb_target >= 0 && (rich.cells[w_vid].hash % 100000u) ==
                                    static_cast<unsigned>(_fb_target)) ||
            (_fb_ctarget >= 0 && (rich.cells[c_vid].hash % 100000u) ==
                                     static_cast<unsigned>(_fb_ctarget));
        if (_fb) {
          std::cerr << "[fallback] W=" << (rich.cells[w_vid].hash % 100000)
                    << " C=" << (rich.cells[c_vid].hash % 100000)
                    << " L(space=" << toUtf8(lvl.space)
                    << " slot=" << lvl.loop_slot << " depth=" << lvl.depth
                    << ") W.occ.carried=[";
          for (auto const& x : occ.carried)
            std::cerr << toUtf8(x.full_label()) << " ";
          std::cerr << "] C.carried:vslot=[";
          for (std::size_t cp = 0; cp < c_carried.size(); ++cp)
            std::cerr << toUtf8(c_carried[cp].full_label()) << ":"
                      << (cp < cvs.size() ? cvs[cp] : -99) << " ";
          std::cerr << "]\n";
        }
        bool recorded = false;
        for (std::size_t cp = 0; cp < c_carried.size() && cp < cvs.size();
             ++cp) {
          // The consumer's slot for position cp IN THIS TREE: the parent
          // occurrence on this edge (point == occ.consumer_point). Mixing C's
          // OTHER occurrences (other trees, permuted frames) leaked a slot
          // binding from a different frame and sliced one operand position
          // by two loops (w20 strict dry-run walk: 20632 pos 1 by slot 1 AND
          // slot 2). The per-value value_slot is the fallback only when the
          // parent occurrence carries no stamp at all.
          bool c_sliced = false;
          if (auto const pit = point_occ.find(occ.consumer_point);
              pit != point_occ.end() && cp < pit->second->loop_slot.size() &&
              pit->second->loop_slot[cp] >= 0)
            c_sliced = pit->second->loop_slot[cp] == lvl.loop_slot;
          else
            c_sliced = cvs[cp] == lvl.loop_slot;
          if (!c_sliced) continue;
          if (std::wstring(c_carried[cp].space().base_key()) != lvl.space)
            continue;
          Index const& M = c_carried[cp];  // the mode C is sliced on here
          // W carries M at some own-occurrence position (shared external,
          // canonical-label match for a direct operand)?
          std::size_t pos_a = occ.carried.size();
          for (std::size_t pa = 0; pa < occ.carried.size(); ++pa)
            if (occ.carried[pa] == M) {
              pos_a = pa;
              break;
            }
          if (_fb)
            std::cerr << "[fallback]   C sliced on M=" << toUtf8(M.full_label())
                      << " -> W carries M? "
                      << (pos_a == occ.carried.size()
                              ? "NO (label mismatch across frames?)"
                              : "yes pos=" + std::to_string(pos_a))
                      << "\n";
          if (pos_a == occ.carried.size()) {
            // W does not carry C's sliced mode in THIS occurrence's frame: its
            // fetch by C is correctly unsliced on this loop. Record that
            // explicitly so the guard does not mistake it for a gap.
            result.occ_invariant.push_back(
                std::make_tuple(w_vid, id_of(lvl), oit->second));
            recorded = true;
            continue;
          }
          result.occ_facts.push_back(
              std::make_tuple(w_vid, pos_a, id_of(lvl), oit->second));
          recorded = true;
          if (std::getenv("SEQUANT_DUMP_USEINDUCED"))
            std::cerr << "[useinduced] value_h="
                      << (rich.cells[w_vid].hash % 100000)
                      << " M=" << toUtf8(M.full_label()) << " pos=" << pos_a
                      << " loop=" << id_of(lvl)
                      << " consumer_h=" << (rich.cells[c_vid].hash % 100000)
                      << "\n";
          break;  // one shared mode per enclosing loop
        }
        // REDUCTION consumer (2026-09-02): C REDUCES a mode under this loop
        // (its result does not carry it, so it is in no carried/value_slot
        // position) -- C is building a per-batch PARTIAL of that reduction,
        // so an operand occurrence carrying the reduced mode is sliced on
        // it; one that does not carry it is legitimately whole. Recorded
        // per occurrence from C's own reduced_slot (fusion union-find).
        if (!recorded) {
          for (auto const& cocc : rich.cells[c_vid].occurrences) {
            for (auto const& [rm, rs] : cocc.reduced_slot) {
              if (rs != lvl.loop_slot) continue;
              if (std::wstring(rm.space().base_key()) != lvl.space) continue;
              std::size_t pos_a = occ.carried.size();
              for (std::size_t pa = 0; pa < occ.carried.size(); ++pa)
                if (occ.carried[pa] == rm) {
                  pos_a = pa;
                  break;
                }
              if (_fb)
                std::cerr << "[fallback]   C REDUCES M="
                          << toUtf8(rm.full_label()) << " -> W carries M? "
                          << (pos_a == occ.carried.size()
                                  ? "NO"
                                  : "yes pos=" + std::to_string(pos_a))
                          << "\n";
              if (pos_a == occ.carried.size())
                result.occ_invariant.push_back(
                    std::make_tuple(w_vid, id_of(lvl), oit->second));
              else
                result.occ_facts.push_back(
                    std::make_tuple(w_vid, pos_a, id_of(lvl), oit->second));
              recorded = true;
              break;
            }
            if (recorded) break;
          }
        }
        // Nothing sliced this operand here and nothing declared it invariant:
        // record the invariant EXPLICITLY so the completeness guard (whose
        // consumer oracle is a per-value union over C's occurrences) does not
        // mistake a legitimately whole read for a gap. Range conformance is
        // enforced by the dry-run/wet backends themselves.
        if (!recorded)
          result.occ_invariant.push_back(
              std::make_tuple(w_vid, id_of(lvl), oit->second));
      }
    }
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
