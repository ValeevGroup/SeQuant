#ifndef SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
#define SEQUANT_EVAL_ORDERED_SCHEDULE_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/dag_scope.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/ordered_dump.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <map>
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
/// \brief The \c
/// OrderedSchedule IR -- an ordered tree of loop blocks and build steps.
/// Purely data + a structural \c well_formed check; the sequencer and
/// executor live elsewhere (\c ordered_executor.hpp).
///
/// \details \c OrderedSchedule threads builds and child loop blocks through a
/// single ordered sequence (\c ScopeBlock::steps): a value built after a
/// child block in that sequence may read the child block's accumulated
/// output, so relative order among steps is load-bearing, not incidental --
/// unlike an unordered per-scope bag of homed values, which cannot express
/// that dependence.
///

///
/// \brief What happens to a value at the close of its home \c ScopeBlock.
///
enum class OutputKind {
  Transient,          //!< the value's home is this block; nothing carries out.
                      //!< \note a Transient value never appears in
                      //!< any block's \c outputs list -- it is realized purely
                      //!< as a plain \c BuildStep, and an explicit Transient \c
                      //!< outputs entry would double-produce it (violating \c
                      //!< well_formed's single-producer check). Do not scan \c
                      //!< outputs for Transient; it is the absence of an escape
                      //!< output, not a recorded one.
  AccumulateSum,      //!< reduction: summed into an outer accumulator
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
/// root block, which sits outside every loop), containing an ordered
/// sequence of build steps and nested child blocks, plus the set of values
/// that leave this block (and how) when it closes.
///
struct ScopeBlock {
  Index axis{};  //!< the loop axis; default (sentinel) on the root block.
  int latitude_ordinal = 0;  //!< layout: the pass index --
                             //!< disambiguates recurring sibling blocks
                             //!< realizing the same axis (a forced-split nest
                             //!< emits one sibling block per pass it holds; see
                             //!< \c forced_split_levels).
  DagScopeLevel level{};     //!< this block's DAG-scope nest position (mirrors
                             //!< \c axis and \c latitude_ordinal; see \c
                             //!< DagScopeLevel's doc comment). Default-valued
                             //!< (depth 0, empty space, altitude/latitude
                             //!< ordinals 0) on the root block.
  BatchModeType kind =
      BatchModeType::Contracted;  //!< Contracted (accumulate on block exit)
                                  //!< or External (scatter on block exit);
                                  //!< meaningless on the root.
  container::vector<Step> steps;  //!< ordered: build-or-child-block,
                                  //!< interleaved (see \c Step's doc
                                  //!< comment for why \c container::vector,
                                  //!< not \c svector).
  container::svector<std::pair<std::size_t, OutputKind>>
      outputs{};  //!< value_id -> how it leaves this block on close.

  // Note: `steps` deliberately has no default member initializer: `{}` here
  // would instantiate std::vector<Step>'s default constructor in-class,
  // where Step is still incomplete (rejected by clang with libstdc++).
  // Declared here, defined below \c Step (see the out-of-line "= default"
  // definitions there). \c steps is a \c std::vector of the still-incomplete
  // \c Step, which the library tolerates only up to first use; defining
  // ScopeBlock's special members here -- including implicitly, by defaulting
  // them in-class or by letting the compiler generate them -- is that first
  // use, and instantiates std::vector<Step>'s own special members on an
  // incomplete element type. libstdc++ then does pointer arithmetic on \c
  // Step*, which clang rejects ("arithmetic on a pointer to an incomplete
  // type"); GCC happens to accept it. Moving the definitions past \c Step is
  // the standard way to close the ScopeBlock <-> Step cycle. Do not collapse
  // these back into in-class "= default".
  ScopeBlock();
  ScopeBlock(ScopeBlock const&);
  ScopeBlock(ScopeBlock&&);
  ScopeBlock& operator=(ScopeBlock const&);
  ScopeBlock& operator=(ScopeBlock&&);
  ~ScopeBlock();
};

///
/// \brief One ordered step of a \c ScopeBlock: build a value at this scope,
/// or enter a nested child loop block.
///
/// \details The design brief's target shape is a bare alias, \c using Step =
/// std::variant<BuildStep, ScopeBlock>. That is not directly expressible:
/// \c ScopeBlock::steps must hold a sequence of \c Step (so a build
/// interleaves with child blocks in one ordered list, per the class doc
/// above), which means \c Step has to be at least forward-declared before \c
/// ScopeBlock; but \c std::variant requires every alternative type complete
/// at the point the variant specialization is instantiated (unlike \c
/// container::vector / \c std::vector, which the C++17 library tolerates
/// holding an incomplete element type up to first use -- the same relaxation
/// \c ScopeBlock::steps itself relies on for its
/// self-reference). A bare \c using Step = std::variant<BuildStep,
/// ScopeBlock> declared before \c ScopeBlock therefore cannot compile
/// (\c ScopeBlock incomplete there), and a type alias cannot be
/// forward-declared separately from its definition (no "using Step;"
/// forward declaration exists in C++) to defer it to after \c ScopeBlock
/// either.
///
/// Hence \c Step is a real (forward-declarable) class wrapping the
/// variant, not a bare alias: \c ScopeBlock::steps holds \c
/// container::vector<Step> while \c Step is still only forward-declared
/// (legal, per the \c std::vector incomplete-type allowance above), and \c
/// Step's own definition (with the variant member) follows \c ScopeBlock,
/// where \c ScopeBlock is by then complete. This preserves the single
/// ordered sequence of interleaved build/child-block steps the design
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

// ScopeBlock's special members, defined here rather than in the class because
// \c Step -- std::vector<Step>'s element type -- is only complete from this
// point on. See the declarations in \c ScopeBlock.
inline ScopeBlock::ScopeBlock() = default;
inline ScopeBlock::ScopeBlock(ScopeBlock const&) = default;
inline ScopeBlock::ScopeBlock(ScopeBlock&&) = default;
inline ScopeBlock& ScopeBlock::operator=(ScopeBlock const&) = default;
inline ScopeBlock& ScopeBlock::operator=(ScopeBlock&&) = default;
inline ScopeBlock::~ScopeBlock() = default;

///
/// \brief The whole ordered schedule: the root block plus the total value
/// count (every \c BuildStep::value_id and every \c ScopeBlock::outputs
/// value_id is expected to be < \c num_values; see \c well_formed).
///
struct OrderedSchedule {
  ScopeBlock root{};
  std::size_t num_values = 0;
  /// Per value_id, the value_ids of its direct operands -- the value/
  /// occurrence DAG edges the value-driven ordered executor consumes to fetch
  /// each operand by its own cell id (see \c CellTable / \c CellRegistry).
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
  // of this block only; deeper levels are checked by the recursion above).
  for (std::size_t i = 0; i < block.steps.size(); ++i) {
    auto const* ci = std::get_if<ScopeBlock>(&block.steps[i].value);
    if (!ci) continue;
    for (std::size_t j = i + 1; j < block.steps.size(); ++j) {
      auto const* cj = std::get_if<ScopeBlock>(&block.steps[j].value);
      if (!cj) continue;
      // Two sibling blocks are the same realized loop only when their full
      // loop identity collides: (depth, loop_slot) and the latitude (pass
      // index). Keying on the axis space (a fusion color, not identity) would
      // wrongly reject two distinct same-space sibling loops that the un-fuse
      // legitimately emits at different (depth, loop_slot): two occupied
      // (space "i") nests at (1,0) and (2,1), same latitude 0, are different
      // loops, not a duplicate. (See LoopKey -- whose identity is the
      // (depth, loop_slot) pair -- and the same space-vs-identity correction
      // in the home-scope coloring.)
      if (ci->level.depth == cj->level.depth &&
          ci->level.loop_slot == cj->level.loop_slot &&
          ci->latitude_ordinal == cj->latitude_ordinal)
        return false;
    }
  }

  for (auto const& [value_id, kind] : block.outputs) {
    (void)kind;
    if (value_id >= num_values) return false;
  }

  return true;
}

///
/// \brief Append every value_id \p block produces -- as a \c BuildStep
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
/// \brief A single escape (output) site: a value_id and the root-to-block path
/// at which it escapes. Feeds \c well_formed's multi-level escape-chain check
/// (a value carried on an outer axis and reduced on an inner one escapes at
/// both -- see \c build_ordered_schedule's escape emission).
///
struct OutputSite {
  std::size_t value_id;
  container::svector<ScopeBlock const*>
      path;  //!< root .. this block, inclusive
};

///
/// \brief Collect every \c BuildStep site (value_id + its root-to-block path)
/// into \p builds and every output escape site into \p sites. Both carry the
/// path: \c well_formed's built-and-escaped rule compares where a value is
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
///   - single-producer (SSA-like), with the multi-level escape chain allowed:
///     no value_id is built (\c BuildStep) more than once; a built value_id may
///     also escape only through its own chain -- every block listing it in \c
///     outputs either holds its \c BuildStep or is an ancestor of the block
///     that does (a member materialized across a forced loop split is built
///     at its production site for its in-nest readers and escapes from that
///     same site outward; see \c build_ordered_schedule's mixed-pass rule) --
///     and any
///     other combination of build and escape sites is duplicate production;
///     and a value_id may escape (\c outputs) at more than
///     one block only when those blocks lie on a single root-to-node nesting
///     path (distinct depths, each shallower one an ancestor of the deepest) --
///     the inner-sum / outer-scatter escape chain of \c build_ordered_schedule.
///     Escapes in unrelated (sibling) blocks, or two escapes at one depth, are
///     rejected as duplicate production.
///
/// \note This checks single-producer (no duplicate production) but not
/// completeness (no value_id gaps -- that every id in `[0, num_values)` is
/// produced somewhere). Completeness holds by construction of \c
/// build_ordered_schedule and is covered by its own test, so a caller cannot
/// read \c well_formed as a guarantee that every value_id is present.
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
    if (auto const it = std::adjacent_find(b.begin(), b.end()); it != b.end())
      return false;
  }
  // (b) a built value_id may also escape, but only through its own chain:
  // every block listing it in `outputs` must either be the block holding its
  // BuildStep (built and escaped in one block -- a mixed-pass member of a
  // forced split: its in-nest readers take the per-batch cell, the other pass
  // takes the assembled form) or an ancestor of it (the chain levels above
  // the production site). An escape in a sibling or descendant of the
  // production site is a second, unrelated producer.
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
      if (!ok) return false;
    }
  }
  // (c) a value_id's escape sites (>1 => a multi-level chain) must lie on one
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
      for (auto const* s : ss) {
        if (!depths.insert(s->path.size()).second)  // same depth
          return false;
        if (s->path.size() > deepest->path.size()) return false;
        for (std::size_t k = 0; k < s->path.size(); ++k)
          if (s->path[k] != deepest->path[k])  // not an ancestor
            return false;
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
/// ScopeBlock, already built) makes visible to its own siblings at this
/// same block level (\c produced), which value_id's its content directly
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
/// paired index-for-index with \p meta) by the local dependency edges among
/// this block's own steps: step A must precede step B whenever B's \c
/// requires_ names a value_id that's in A's \c produced. Kahn's algorithm;
/// among simultaneously-ready steps, always picks the smallest \c tie_key
/// first, for a deterministic result when the true dependency order leaves
/// steps genuinely unordered relative to each other.
///
/// \details A single scalar per step can be made to sort a child block before
/// every value that reads its output (see \c build_ordered_schedule's own doc
/// comment, part 3), but it cannot also guarantee a step sorts after every
/// value its own content reads as an input -- those are two independent
/// constraints a single total order satisfies only when they happen to agree.
/// A topological sort over the actual per-step dependency edges satisfies both
/// directions by construction, and \c tie_key serves only to break ties among
/// steps with no dependency relation to each other at all.
///
/// \c SEQUANT_ASSERT's that every item is placed exactly once (a cycle in
/// this local edge set would be a bug -- these edges are a sub-relation of
/// the whole-forest DAG's edges, restricted to one block's own siblings, so
/// they inherit its acyclicity), then re-derives the local edges a second
/// time against the final order and \c SEQUANT_ASSERT's every one is
/// actually satisfied (a loud tripwire against any future violation of this
/// invariant, per the design review that requested it, rather than a silent
/// mis-order).
///
inline container::vector<Step> ordered_schedule_topo_sort_steps(
    container::vector<Step> items,
    container::vector<OrderedScheduleStepMeta> const& meta) {
  std::size_t const m = items.size();
  SEQUANT_ASSERT(meta.size() == m);

  // value_id -> which local item produces it. well_formed's whole-schedule
  // single-producer invariant guarantees at most one item at any level can
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
  // No cycle: see the function doc comment. Thrown rather than asserted so
  // this stays loud in a build with asserts disabled -- with a short \c
  // order, the code below would otherwise build \c out_steps from a
  // truncated \c order and silently drop the unplaced steps.
  if (order.size() != m) {
    std::size_t unsatisfied_edges = 0;
    for (std::size_t i = 0; i < m; ++i) unsatisfied_edges += indegree[i];
    throw Exception(
        "ordered_schedule_topo_sort_steps: cyclic step dependencies among " +
        std::to_string(m) + " sibling steps (" +
        std::to_string(unsatisfied_edges) +
        " prerequisite edges never satisfied)");
  }

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
/// \brief The global direct-dependency edges of a \c RichSchedule, recovered
/// from \c rich alone (no forest access): \c depends_on[p] lists every value_id
/// the value \c p directly reads, and \c consumers_of[c] lists every value_id
/// that directly reads \c c (the reverse). Same recovery \c
/// build_ordered_schedule uses inline (occurrence \c consumer_point ->
/// producing value_id via \c point_owner); factored here so \c
/// forced_split_levels and \c build_ordered_schedule agree edge-for-edge.
///
struct OrderedScheduleDepGraph {
  std::unordered_map<std::size_t, std::size_t> value_id_of;  //!< hash -> id
  std::unordered_map<std::size_t, container::svector<std::size_t>> depends_on;
  std::unordered_map<std::size_t, container::svector<std::size_t>> consumers_of;
};

/// Recovers the global direct-dependency edges of \p rich (see
/// \c OrderedScheduleDepGraph).
inline OrderedScheduleDepGraph ordered_schedule_dep_graph(
    RichSchedule const& rich) {
  OrderedScheduleDepGraph g;
  g.value_id_of.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells)
    g.value_id_of.emplace(value_key_of(vc), vc.value_id);

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
/// \brief Pass levels, global over every batched space (design section 9.2):
/// every value gets an integer pass such that a value that
/// needs another value's completed form sits in a later pass than that
/// other value.
///
/// Two kinds of dependency edge bump the reader's pass, both keyed on the
/// operand (the "source") rather than the axis space:
///   - the source is \c LoopCarried on some axis (any space): its full array
///     exists only after its own loop closes, so every direct reader is
///     bumped;
///   - the source is a \c Reduction on some instance (an \c AccumulateSum
///     escape) and the reader is produced inside that same instance (its
///     production site is at or below the reduction's depth, in the same
///     nest): such a reader would otherwise see the current batch's partial
///     sum rather than the completed reduction (the finding pinned by the
///     9.1 partial-sum check). A reader produced outside the reduced
///     instance reads the completed sum from the escape's residency scope
///     as always and needs no bump. \p inside(reader_vid, source_vid)
///     decides this per (reader, source) pair -- the caller's lambda
///     captures the loop chain (fusion_slot, depth_of_instance,
///     production_depth, type_cluster) needed to resolve the source's
///     reduced instance(s) and the reader's production site exactly as the
///     escape placement does.
///
/// Forward sweep (operands before consumers):
///   level(v) = max over direct operands o of (bump(o, v) ? pass(o) + 1
///              : level(o)), 0 with no operands;
///   a bumping-edge source's pass is its level (it is pinned, see below); a
///   value that is the source of no bumping edge has its level as its base.
/// Reverse sweep (consumers before operands): a non-pinned value with at
/// least one consumer is lifted to max(base, min over its direct consumers'
/// passes), so a value whose readers all sit later is built with them; a
/// value with readers in several passes keeps its base (and is materialized
/// by the builder's rule 4 when a later same-nest reader needs it).
///
/// Every dependency edge points to an equal or earlier pass. With only
/// LoopCarried bumps present (no Reduction-source bump fires) the passes
/// reduce to a two-set partition: carried sources and everything else.
///
struct ForcedSplitLevels {
  std::unordered_set<std::size_t> carried;  //!< LoopCarried (any space) ids
  std::unordered_set<std::size_t>
      pinned;  //!< sources of a bumping edge (carried values, and Reduction
               //!< sources with at least one in-loop reader): the reverse lift
               //!< skips every id in this set.
  std::unordered_map<std::size_t, int>
      pass_of;  //!< value id -> pass, for every value \c
                //!< ordered_schedule_dep_graph reached (has a legality cell
                //!< and takes part in the dependency graph); absent for a
                //!< value with no legality cell, which \c pass() below
                //!< reports as pass 0 rather than throwing.
  int max_pass = 0;
  int pass(std::size_t vid) const {
    auto const it = pass_of.find(vid);
    return it == pass_of.end() ? 0 : it->second;
  }
};

/// Assigns every value of \p rich its pass level (see \c ForcedSplitLevels),
/// reading the carried/reduction roles from \p legality, the dependency edges
/// from \p g, and deciding with \p inside whether a reader is produced inside
/// a source's own reduced loop instance.
inline ForcedSplitLevels forced_split_levels(
    RichSchedule const& rich, LegalitySchedule const& legality,
    OrderedScheduleDepGraph const& g,
    std::function<bool(std::size_t /*reader_vid*/,
                       std::size_t /*source_vid*/)> const& inside) {
  ForcedSplitLevels r;
  std::unordered_set<std::size_t> reduction_sources;
  for (CellLegality const& cl : legality.cells) {
    auto const it = g.value_id_of.find(cl.hash);
    if (it == g.value_id_of.end()) continue;
    bool const carried_here = std::any_of(
        cl.per_axis.begin(), cl.per_axis.end(),
        [&](AxisClass const& ac) { return ac.role == LoopRole::LoopCarried; });
    if (carried_here) r.carried.insert(it->second);
    bool const reduces_here = std::any_of(
        cl.per_axis.begin(), cl.per_axis.end(),
        [&](AxisClass const& ac) { return ac.role == LoopRole::Reduction; });
    if (reduces_here) reduction_sources.insert(it->second);
  }
  r.pinned = r.carried;

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
      for (std::size_t o : it->second) {
        bool const bump = r.carried.count(o) != 0 ||
                          (reduction_sources.count(o) != 0 && inside(v, o));
        if (bump) r.pinned.insert(o);
        if (bump && detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
          detail::dump_level_bump(v, o, r.carried.count(o) != 0,
                                  base.at(o) + 1);
        lv = std::max(lv, bump ? base.at(o) + 1 : base.at(o));
      }
    base[v] = lv;
  }

  r.pass_of = base;
  for (auto it = topo.rbegin(); it != topo.rend(); ++it) {
    std::size_t const v = *it;
    if (r.pinned.count(v)) continue;
    auto const cit = g.consumers_of.find(v);
    if (cit == g.consumers_of.end() || cit->second.empty()) continue;
    int mn = std::numeric_limits<int>::max();
    for (std::size_t u : cit->second) mn = std::min(mn, r.pass_of.at(u));
    r.pass_of[v] = std::max(base.at(v), mn);
  }
  for (auto const& [v, p] : r.pass_of) r.max_pass = std::max(r.max_pass, p);
  return r;
}

///
/// \brief The predicate-false and predicate-true copies of a forked inner
/// sub-chain (see \c fork_subchain).
///
/// \details The per-nest pass-split builder takes only \c consumer once per
/// pass (the predicate-true side) and never reads \c producer; the field is
/// kept for other callers (e.g. \c fork_subchain's own unit tests) that fork
/// on a two-way predicate rather than a per-pass one.
///
struct ForkedSubchain {
  container::vector<Step> producer;  //!< steps whose values are predicate-false
  container::vector<Step> consumer;  //!< steps whose values are predicate-true
};

///
/// \brief Fork an already-built inner sub-chain (an ordered list of \c Step)
/// into a predicate-false copy and a predicate-true copy, used once per pass
/// at a nest holding a forced-split axis (\c build_ordered_schedule): for pass
/// \p k, \p in_consumer(value_id) is `pass_of(value_id) == k`, so the
/// predicate-true side is exactly that pass's own steps out of the nest's
/// full pending sub-chain (its predicate-false side, everything else, is
/// picked up by a later call for a different \p k).
///
/// \details A \c BuildStep goes wholly to one side by \p in_consumer of its
/// value. A nested \c ScopeBlock (an inner loop) is recursively forked; each
/// side that has surviving steps or surviving escape \c outputs is rebuilt as a
/// per-side copy of the loop (same \c axis / \c ordinal / \c kind) carrying
/// only that side's steps and the escape \c outputs whose value lands on that
/// side; a side with neither is dropped (an empty loop is never emitted). Note:
/// a loop can be all-escape -- no \c BuildStep, only scatter \c outputs
/// contracted at the output step -- so "no surviving steps" does not imply "no
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
      // Emit the per-side loop when it has surviving steps *or* surviving
      // outputs. An all-escape loop -- one with no BuildStep, whose scatter
      // values are contracted at the output step itself -- is legitimate and
      // must not be dropped: doing so strands a whole nested loop (e.g. an
      // inner occ member loop or an aux loop) together with its escape outputs,
      // leaving the forced-split axis as the only realized loop (the
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
/// \brief The deterministic sequencer -- lowers the \c LegalitySchedule
/// (per-value \c home_floor / \c per_axis roles) plus the \c RichSchedule
/// (per-value \c first_use / \c last_use over the forest's single post-order
/// static-point timeline) into an \c OrderedSchedule.
///
/// The realized chain is per loop instance, not per axis type: one block per
/// distinct fusion \c loop_slot appearing on a space across all cells, so a
/// value carrying two same-space modes on different slots gets two nested
/// loops, not one. A nest holding members of more than one pass additionally
/// emits one sibling block per pass (latitude = pass), run in schedule order
/// (see step 2b and \c forced_split_levels). See the as-built design
/// section 6.1, \c
/// doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md.
///
/// \details Four-part algorithm, pure scheduling (no cost choice):
///
/// \par 1. The canonical chain
/// For each batch axis space (\c IndexSpace::base_key()) appearing in any
/// cell's \c CellLegality::per_axis (not just \c home_floor -- a \c
/// Reduction/\c LoopCarried-only axis must still get a block to host its
/// escape output, even though no cell is homed inside it in that role), one
/// realized depth per distinct fusion \c loop_slot on it (\c
/// slots_of_space / \c type_slot, resolved through \c fusion_slot). Spaces
/// are ordered by \p mode_order (most-significant/outermost first), ties/
/// unlisted spaces alphabetical, slot-ascending within a space -- but that is
/// only a stable tie-break: the nesting actually satisfies every \c
/// RichSchedule::loop_order pair, and a cycle among those is a builder error.
/// Loop members that never co-occur (no value is home-sliced on both) are then
/// cut into disjoint nests by a union-find (\c type_cluster / \c
/// cluster_min) and concatenated at root, rather than realized as one
/// over-deep chain.
///
/// \par 2. Per-value placement: home \c BuildStep vs. escape output
/// A value's \c CellLegality::per_axis has, by construction (see \c
/// CellLegality::per_axis's own doc comment), only \c LoopLocal, \c
/// Reduction, or \c LoopCarried entries (never the implicit \c
/// LoopInvariant). Two cases:
///   - every \c per_axis entry is \c LoopLocal (this includes the empty
///     case: no batch-axis dependence at all, e.g. water-20's \c
///     I(i,i;a,a)) -- the value is a plain \c BuildStep. Its home block is
///     the depth whose accumulated (root-to-depth) type set equals \c
///     home_floor's type set (root if \c home_floor is empty), by set
///     equality (an unmatched/non-prefix \c home_floor falls back to
///     root). This value
///     gets no \c outputs entry anywhere (see \c well_formed's single-
///     producer invariant): "Transient" (design point 4) is realized as
///     "produced by a \c BuildStep and nothing else", not as an explicit
///     \c OutputKind::Transient \c outputs record, since \c
///     well_formed::detail::collect_production_ids counts every \c outputs
///     entry (regardless of \c OutputKind) as an independent production
///     site -- a \c Transient \c outputs entry alongside the \c BuildStep
///     would be flagged as double-production.
///   - at least one \c per_axis entry is \c Reduction or \c LoopCarried
///     ("escapes" that axis, per design point 4: \c Reduction ->
///     accumulate-summed out, \c LoopCarried -> accumulate-scattered out)
///     -- the value has no \c BuildStep anywhere; instead it is recorded as
///     an \c outputs entry (kind \c AccumulateSum / \c AccumulateScatter)
///     of each escaped instance's block. With more than one escape this is a
///     multi-level escape chain: raw production at the deepest escape site,
///     pure forwarding at every shallower one. The innermost escaped loop is
///     where the accumulation the value's own node performs actually happens,
///     so that block is its true production site; an outer escape on a
///     shallower axis needs this value already complete before the outer loop
///     can close -- exactly consistent with an outer accumulator reading an
///     inner one. A chain may legitimately skip a level the value is
///     invariant on; that crossing is carried by residency plus \c
///     produce_if_absent, not by an escape (as-built section 6.2).
///
/// \par 3. Topological order within a block -- a real topological sort
/// Each block's own \c steps interleave its \c BuildStep's (one per value
/// homed there) with, if the chain continues, one nested child \c
/// ScopeBlock \c Step for the next-deeper axis. These are ordered by
/// \c detail::ordered_schedule_topo_sort_steps against a per-step dependency
/// graph, not a scalar key alone (see below for why a scalar key cannot
/// suffice), reconstructed from \p rich alone -- no forest access needed:
///   - global direct-dependency edges: for every \c OccurrenceRec of every
///     value, its \c consumer_point names the static point of its
///     structural parent node; resolving that point back to the value_id
///     whose own occurrence starts there (\c point_owner, built once up
///     front) recovers "this parent value directly reads that child value"
///     -- the exact same edges the forest itself encodes, without needing
///     the forest.
///   - Per level (one \c ScopeBlock's own \c steps list, including root),
///     each candidate step gets a \c detail::OrderedScheduleStepMeta:
///     - a \c BuildStep{v}'s \c produced = `{v}`; its \c requires_ = every
///       value_id \c v directly reads (raw, unfiltered -- irrelevant/
///       external entries are dropped by the topo-sort itself, since they
///       simply never match a local \c produced set).
///     - a nested child block's \c produced = that block's own top-level
///       \c outputs value_id's (what it makes visible to its own parent's
///       siblings -- its internal \c BuildStep's and any further-nested
///       child's content are never directly readable from outside it: by
///       construction, a value crossing a block boundary as an operand
///       must first have been resolved out of that axis, which is exactly
///       the escape/\c outputs case). Its \c requires_ is the full,
///       recursively bubbled external need of its whole subtree (built
///       bottom-up alongside the block itself: \c requires_all(level) =
///       (this level's own direct needs union its child's already-bubbled
///       \c requires_all) minus \c produced_all(level), where
///       \c produced_all is everything ever produced anywhere in the
///       subtree, recursively) -- so a need that is only satisfiable
///       several levels further out (e.g. a root-homed common factor
///       consumed by a value nested two axes deep) still surfaces at
///       whichever level can actually satisfy it.
///   - Ties (two ready steps with no dependency relation to each other) are
///     broken by \c tie_key ascending: a \c BuildStep's is its value's own
///     \c ValueCell::first_use; a child block's is the min \c first_use
///     over its own \c produced_all (deterministic, and -- though correctness
///     rests on the real edges, which enforce both directions -- it places a
///     block as early as its own true dependency slack allows).
///
/// A single scalar key alone cannot express both directions of this at once.
/// Taking the min \c first_use as the sort key itself is sound for "the block
/// sorts before every true consumer", but carries no corresponding guarantee
/// for "the block sorts after every true input it reads": a value produced by
/// a same-level sibling \c BuildStep (e.g. a root-homed operand consumed by
/// content nested inside a child block) can land, by raw point value, after
/// the block's min-derived key, mis-ordering the schedule with no structural
/// check to catch it. The topological sort above satisfies both directions by
/// construction and is
/// checked twice (no-cycle placement count, then a second pass confirming
/// every edge survived the final order) -- see \c
/// ordered_schedule_topo_sort_steps's own doc comment.
///
/// \p policy is accepted for interface symmetry with the rest of the pipeline
/// (every stage from \c analyze_legality onward threads it)
/// and as a hook for a split threshold; the logic here
/// only consults \p rich and \p legality; the batchable-axis filtering
/// \p policy would otherwise provide is already baked into \c
/// CellLegality::per_axis by \c analyze_legality.
///
namespace detail {

///
/// \brief True iff \p mode's index type (\c IndexSpace::base_key()) ever
/// survives, un-summed, into a forest-root value's own carried slots.
///
/// \details \c RichSchedule does not carry \c BatchModeType directly: the
/// cross-occurrence meet in \c stamp_lifetime_masks folds every kind
/// (External or Contracted) into one \c Index set (see \c
/// stamp_residency_impl's doc comment -- "any BatchModeType"), so a mode's
/// kind is not a field anywhere on \c ValueCell. It is still recoverable from
/// \p rich alone: a batch mode realized in the forest is either Contracted
/// (summed away below every root -- an accumulate loop, so it never appears
/// on a root's own result) or External (a free/spectator index of some final
/// output -- a scatter loop, one output slice per block, so it does appear on
/// a root's own result). A forest-root occurrence is identified by \c
/// OccurrenceRec::point == \c consumer_point: \c compute_dag_boulevard's
/// post-order walk only ever overwrites a child's \c consumer_point (to its
/// parent's point, always later in the walk), so a node that is never
/// anyone's child -- a top-level tree of the forest -- keeps its
/// default-seeded \c consumer_point == its own \c point.
///
/// Matches by type (\c base_key()) -- not by exact \c Index identity (\c
/// Index::operator== compares space and ordinal). More than one physical \c
/// Index label of the same type can appear across \p rich's cells (e.g. one
/// occurrence's K_1 vs another's K_2), so an exact-Index match could miss a
/// root occurrence that uses a different physical label of the very type
/// being classified, silently misclassifying an External mode as Contracted.
///
inline bool mode_is_external(RichSchedule const& rich, Index const& mode) {
  auto const& mode_bk = mode.space().base_key();
  for (auto const& cell : rich.cells) {
    bool const is_root =
        std::any_of(cell.occurrences.begin(), cell.occurrences.end(),
                    [](OccurrenceRec const& occ) {
                      return occ.point == occ.consumer_point;
                    });
    if (!is_root) continue;
    // A batch axis (occ/aux) is external iff its space appears on a root's own
    // result slots. Slots are matched as-is: batch axes are plain, and a root
    // carrying a composite PAO leg a<i,j> always carries its occ pair i,j
    // plainly too, so no space is reachable only through a proto.
    bool const has_type = std::any_of(
        cell.carried.begin(), cell.carried.end(),
        [&](Index const& ix) { return ix.space().base_key() == mode_bk; });
    if (has_type) return true;
  }
  return false;
}

}  // namespace detail

/// \brief Builds the \c OrderedSchedule of \p rich: the loop-block tree, the
///        build steps inside each block, and each value's escape chain.
///
/// \param rich the schedule's values, their occurrences and loop instances
/// \param legality the per-value, per-axis roles \c analyze_legality derived
/// \param policy accepted for interface symmetry with the rest of the
///        pipeline; the axis filtering it would provide is already baked into
///        \c CellLegality::per_axis
/// \param mode_order optional outer-to-inner order of index-space base keys,
///        used to order the loop chain
/// \return a schedule satisfying \c well_formed
/// \throw Exception when a batched mode has no loop identity to place its
///        escape on
[[nodiscard]] inline OrderedSchedule build_ordered_schedule(
    RichSchedule const& rich, LegalitySchedule const& legality,
    [[maybe_unused]] BatchPolicy const& policy,
    std::initializer_list<std::wstring> mode_order = {}) {
  OrderedSchedule out;
  out.num_values = rich.cells.size();

  // hash -> value_id and the global direct-dependency edges, recovered from
  // rich alone (see the function doc comment's part 3, and \c
  // ordered_schedule_dep_graph): for every occurrence of every value, its
  // consumer_point names its structural parent's own production point,
  // resolved back to the parent value_id.
  auto const g = detail::ordered_schedule_dep_graph(rich);
  auto const& value_id_of = g.value_id_of;
  static container::svector<std::size_t> const kNoDeps{};
  auto const requires_of =
      [&](std::size_t vid) -> container::svector<std::size_t> const& {
    auto const it = g.depends_on.find(vid);
    return it == g.depends_on.end() ? kNoDeps : it->second;
  };

  // The fusion-assigned loop_slot (compute_dag_boulevard) of a
  // legality cell's per_axis[pos]: which member of its same-space loop group
  // slices it -- an occurrence-invariant identity, established by producer->
  // consumer connectivity, not a within-cell position. Read off a
  // representative occurrence in the value's own frame (per_axis modes match
  // the occurrence's `carried` by label). Returns -1 if unavailable (a leaf, a
  // not-batched mode, or a divergent occurrence): the LoopCarried caller
  // falls back to slot 0, but the Reduction caller (the escape-placement
  // loop below) treats -1 as a hard error and throws -- a Reduction mode
  // with no slot has no loop identity at all (compute_dag_boulevard's
  // union-find never numbered a component for it), so guessing would silently
  // place the escape in the wrong nest.
  std::unordered_map<std::size_t, ValueCell const*> hash_to_rich;
  hash_to_rich.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells)
    hash_to_rich.emplace(value_key_of(vc), &vc);
  auto const fusion_slot = [&](CellLegality const& cl, std::size_t pos) -> int {
    auto const hit = hash_to_rich.find(cl.hash);
    if (hit == hash_to_rich.end() || hit->second->occurrences.empty())
      return -1;
    OccurrenceRec const& occ = hit->second->occurrences.front();
    Index const& m = cl.per_axis[pos].axis;
    auto const cit = std::find(occ.carried.begin(), occ.carried.end(), m);
    if (cit == occ.carried.end()) {
      // Not a carried mode: a Reduction mode is contracted at this value and
      // has no carried position. Its slot is stamped by compute_dag_boulevard
      // either via carried->reduced propagation (uniting the home-sliced
      // operand's carried-mode node with a synthetic reduction node) or,
      // absent any home-sliced operand, by seeding that synthetic node
      // directly; read that slot so the reduction escape lands in the same
      // same-space nest as the operand it reduces, instead of a different
      // nest (the operand vanishes before the reduction reaches it) or,
      // absent a slot altogether, the hard-error throw above.
      for (auto const& [rm, rs] : occ.reduced_slot)
        if (rm == m) return rs;
      return -1;
    }
    std::size_t const p = static_cast<std::size_t>(cit - occ.carried.begin());
    return p < occ.loop_slot.size() ? occ.loop_slot[p] : -1;
  };

  // 1. The canonical chain: one representative Index per distinct axis type
  // present in any cell's per_axis (LoopLocal, Reduction, or LoopCarried --
  // not just home_floor; see the function doc comment's part 1).
  // Per-instance loop chain (position-based): for each space, m_s = the max
  // over cells of the number of same-space per_axis modes in one cell; emit
  // m_s consecutive depths, one per within-space slot (loop_slot 0..m_s-1). A
  // value carrying two same-space batched modes (a doubles amplitude's two occ
  // externals) thus gets two distinct loops rather than one. `types[d]` is a
  // representative Index of the depth's space; `type_slot[d]` is its
  // within-space slot (the DagScopeLevel loop_slot). Same-space slots occupy
  // distinct depths here (the assembly is one loop per depth); depth carries
  // both group and member nesting.
  // Members present per space = the distinct fusion loop_slots that
  // appear on that space across all cells. Each distinct slot becomes one
  // realized loop. The fusion slot is the occurrence-invariant identity: a
  // count of same-space per_axis modes with local 0..m-1 slots would be a
  // per-position numbering the atlas cannot match to a value's own frame.
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
  container::svector<int> type_slot;  // the fusion loop_slot of this loop
  for (auto const& bk : spaces)
    for (int s : slots_of_space.at(bk)) {  // ascending (std::set)
      types.push_back(rep.at(bk));
      type_slot.push_back(s);
    }

  // Nest the chain as the DP realized it: RichSchedule::loop_order holds the
  // (outer, inner) instance pairs read off the occurrences' enclosing
  // contexts. The space-major, slot-ascending order above is only the
  // tie-break (a stable topological order: among the ready instances the
  // earliest in that order goes first). A cycle means the loop identity fused
  // two physical loops that nest in opposite orders -- a builder error, not a
  // silent choice.
  if (!rich.loop_order.empty()) {
    std::size_t const nn = types.size();
    std::map<std::pair<std::wstring, int>, std::size_t> item_of;
    for (std::size_t d = 0; d < nn; ++d)
      item_of[{std::wstring{types[d].space().base_key()}, type_slot[d]}] = d;
    container::svector<container::svector<std::size_t>> succ(nn);
    container::svector<std::size_t> indeg(nn, 0);
    for (auto const& [pair, witness] : rich.loop_order) {
      auto const& [outer, inner] = pair;
      (void)witness;
      auto const a = item_of.find(outer);
      auto const b = item_of.find(inner);
      if (a == item_of.end() || b == item_of.end() || a->second == b->second)
        continue;
      succ[a->second].push_back(b->second);
      ++indeg[b->second];
    }
    std::set<std::size_t> ready;
    for (std::size_t d = 0; d < nn; ++d)
      if (indeg[d] == 0) ready.insert(d);
    container::svector<std::size_t> perm;
    while (!ready.empty()) {
      std::size_t const d = *ready.begin();
      ready.erase(ready.begin());
      perm.push_back(d);
      for (std::size_t s : succ[d])
        if (--indeg[s] == 0) ready.insert(s);
    }
    if (perm.size() != nn) {
      // The instances left unplaced are on (or downstream of) the cycle;
      // list the constraints among them with their witnesses.
      std::string msg =
          "build_ordered_schedule: the loop nesting constraints read off the "
          "batched realization are contradictory (a cycle among loop "
          "instances): the loop identity fused two physical loops that nest "
          "in opposite orders; unplaced constraints:";
      auto const narrow = [](std::wstring const& w) {
        return std::string(w.begin(), w.end());
      };
      for (auto const& [pair, witness] : rich.loop_order) {
        auto const a = item_of.find(pair.first);
        auto const b = item_of.find(pair.second);
        if (a == item_of.end() || b == item_of.end()) continue;
        if (indeg[a->second] == 0 && indeg[b->second] == 0) continue;
        msg += " [" + narrow(pair.first.first) + "#" +
               std::to_string(pair.first.second) + " > " +
               narrow(pair.second.first) + "#" +
               std::to_string(pair.second.second) + " by v" +
               std::to_string(witness) + "]";
      }
      throw Exception(msg);
    }
    container::svector<Index> types2;
    container::svector<int> type_slot2;
    for (std::size_t d : perm) {
      types2.push_back(types[d]);
      type_slot2.push_back(type_slot[d]);
    }
    types = std::move(types2);
    type_slot = std::move(type_slot2);
  }

  if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
    detail::dump_loop_chain(rich.loop_order, types, type_slot);

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
  // is home-sliced on both (carries both in one term, so its fusion loop_slots
  // -- home-based -- name both). Members that never co-occur live
  // in disjoint nests: e.g. two residual sub-DAGs that both batch occ i,j but
  // are connected only through a full symmetric intermediate (home meet empties
  // it, so it seeds no loop) -- they are genuinely separate loop groups.
  // Realizing them as one over-deep nested chain (slot0 superset ... superset
  // slotK) is the structure that deadlocks; each cluster must be a separate
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
  // (larger d) first, so the loop wraps a cluster's members into one nest and
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

  // Legality record by value id, for reader lookups in rule 4 and (below)
  // for the pass levels' \c inside predicate.
  std::unordered_map<std::size_t, CellLegality const*> cl_by_vid;
  cl_by_vid.reserve(legality.cells.size());
  for (CellLegality const& c2 : legality.cells) {
    auto const it2 = value_id_of.find(c2.hash);
    if (it2 != value_id_of.end()) cl_by_vid.emplace(it2->second, &c2);
  }

  // 2. Per-value placement: BuildStep at its production site (root-level
  // bucket uses index n as a sentinel "root" depth; see \c build_depth below
  // for a value materialized across a forced split, whose production site
  // can sit deeper than its LoopLocal home) vs. escape output at the deepest
  // escaped axis's depth.
  container::vector<detail::OrderedScheduleDepthBucket> buckets(n);
  container::svector<std::size_t> root_build_ids;
  // Values that keep their BuildStep and escape it (the mixed-pass members of
  // a forced split, below): dump-only diagnostic (SEQUANT_DUMP_SCHEDULE) --
  // the live signal downstream is the per-value \c materialized_across_split
  // bool below, not this list.
  container::svector<std::size_t> materialized_across_split_ids;

  // The LoopLocal home depth of a value: the innermost loop it is LoopLocal
  // on, resolved per-instance by fusion loop_slot -- not by shallowest
  // same-space count. A value local to slots 2,3 (its own fusion nest) homes
  // in that nest, not the first same-space nest a space-multiset cover would
  // pick; picking the wrong same-space nest homes the value where its
  // consumer's nest has not opened (or has already closed), so an
  // in-consumer-nest read misses and the value vanishes. Resolve each
  // LoopLocal mode's (space, fusion slot) to its realized depth (exactly as
  // the escape placement does), and home at the max such depth: within a
  // co-occurrence cluster larger depth nests inside smaller, so the innermost
  // of the value's own home slots is inside all of them. home_floor is the
  // LoopLocal subset, but it drops the pos->slot map, so walk per_axis
  // directly for the slot. Nullopt = no realized loop-local mode: root. This
  // is a value's plain BuildStep production site unless it is materialized
  // across a forced split and a role escape nests deeper than this depth, in
  // which case \c build_depth (below, in the placement loop) seeds from this
  // and deepens it to the true production site.
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

  // The DAG-scope nest a value is produced inside, for rule-4 reader
  // classification: the deepest depth any of its own per_axis modes (any
  // role, not just LoopLocal) resolves to. A value with only carried/
  // reduction roles -- a forest root delivered in full, or a carried value
  // of a later pass -- is still produced per batch inside its own nest (its
  // production is the accumulation folded into its escape bucket), so
  // testing only LoopLocal modes (local_home_depth) would report such a
  // value as homed at root: rule 4 would then neither fire the mixed-pass
  // materialization for a value it reads, nor guard the tripwire against
  // it. production_depth instead considers every per_axis mode regardless
  // of role. A mode whose fusion slot does not resolve (\c fusion_slot
  // returns -1) is skipped rather than guessed at slot 0 -- a guessed slot
  // can land in the wrong nest (fusion_slot's own doc comment), and this
  // result feeds the outside-its-nest tripwire below, where a wrong nest
  // decides whether to throw. Nullopt = no mode resolves at all, whether
  // because the value is genuinely unbatched (root) or because every one of
  // its modes has an unresolvable fusion slot -- the two are
  // indistinguishable here; the table validator's visibility rule is the
  // net for whichever of those a later-pass reader turns out to be.
  auto const production_depth =
      [&](CellLegality const& cl) -> std::optional<std::size_t> {
    std::optional<std::size_t> target;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
      int const fs = fusion_slot(cl, pos);
      if (fs < 0) continue;  // unresolved: never guess slot 0
      auto const d = depth_of_instance(bk, fs);
      if (!d) continue;
      if (!target || *d > *target) target = *d;
    }
    return target;
  };

  // 2a. Pass levels (per-nest forced split design, sections 3.1 and 9.2).
  // Global over every batched space: a bumping edge is either
  // a LoopCarried source (any space, every direct reader bumped) or a
  // Reduction source with a direct reader produced inside that same
  // reduced instance. \c inside resolves the second kind exactly as the
  // escape placement does: for each of the source's Reduction axes,
  // fusion_slot (which already reads \c reduced_slot for a Reduction mode)
  // gives the instance's slot, depth_of_instance gives its depth, and the
  // reader is "inside" when its own production_depth is at or below that
  // depth in the same nest (type_cluster equal). Defined here because it
  // needs production_depth and cl_by_vid. Computed unconditionally, no
  // per-space loop, no "more than one forced space" throw: an empty
  // carried and reduction-source set yields all-zero passes, so a schedule
  // with no forced-split axis at all is unaffected.
  auto const inside = [&](std::size_t reader_vid,
                          std::size_t source_vid) -> bool {
    auto const sit = cl_by_vid.find(source_vid);
    if (sit == cl_by_vid.end()) return false;
    CellLegality const& scl = *sit->second;
    auto const rit = cl_by_vid.find(reader_vid);
    if (rit == cl_by_vid.end()) return false;
    auto const rd = production_depth(*rit->second);
    if (!rd) return false;
    for (std::size_t pos = 0; pos < scl.per_axis.size(); ++pos) {
      if (scl.per_axis[pos].role != LoopRole::Reduction) continue;
      std::wstring const bk{scl.per_axis[pos].axis.space().base_key()};
      int const fs = fusion_slot(scl, pos);
      if (fs < 0) continue;
      auto const d = depth_of_instance(bk, fs);
      if (!d) continue;
      if (*rd >= *d && type_cluster[*rd] == type_cluster[*d]) return true;
    }
    return false;
  };
  detail::ForcedSplitLevels const levels =
      detail::forced_split_levels(rich, legality, g, inside);
  auto const pass_of = [&](std::size_t vid) -> int { return levels.pass(vid); };

  // Dump-only diagnostic (SEQUANT_DUMP_SCHEDULE) for the outside-nest
  // tripwire rejection: print the value's axes with roles/slots/depths, its
  // later-pass readers, and the carried set, each with home depth and nest.
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
      std::wcerr << L"}" << (levels.carried.count(v) ? L" CARRIED" : L"")
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
    for (std::size_t c : levels.carried) {
      std::wcerr << L"[sched-reject]   carried ";
      describe(c);
      std::wcerr << L"\n";
    }
  };

  for (CellLegality const& cl : legality.cells) {
    auto const vid_it = value_id_of.find(cl.hash);
    SEQUANT_ASSERT(vid_it != value_id_of.end());
    std::size_t const vid = vid_it->second;

    // A forest leaf is an input fetched on demand by its consumers, never a
    // computed value: it must not be emitted as a BuildStep (doing so makes the
    // executor evaluate it as a standalone root -- one wasted leaf fetch per
    // iteration -- which the forest descent never does). Its consumers reach it
    // through the leaf evaluator exactly as in forest descent.
    if (rich.cells[vid].is_leaf) continue;

    // Emit an escape at every non-local axis depth (the multi-level escape
    // chain of a non-innermost split): a value that reduces an inner axis and
    // is carried on an outer one escapes at both -- AccumulateSum at the inner
    // (into the accumulator one level out) then AccumulateScatter at the outer
    // (that accumulator to full). The bottom-up assembly materializes them
    // inner -> outer. A value non-local on a single axis keeps exactly one
    // escape, unchanged. Same-depth axis-classes (e.g. two carried occ indices,
    // or a same-type reduce+carry pair) collapse to one escape at that depth,
    // with LoopCarried (AccumulateScatter) dominating Reduction: a carried axis
    // must materialize to full even if a same-type index reduces.
    container::svector<std::pair<std::size_t, OutputKind>> escapes;
    for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
      AxisClass const& ac = cl.per_axis[pos];
      if (ac.role == LoopRole::LoopLocal) continue;
      std::wstring const bk{ac.axis.space().base_key()};
      int const fs = fusion_slot(cl, pos);
      // A Reduction mode with no fusion slot has no loop identity at all
      // (peak_profile's union-find never numbered a component for it) --
      // guessing slot 0 places the AccumulateSum escape in whatever nest
      // happens to own slot 0, not the loop the operands were actually
      // sliced on, and every batch then silently contracts the full
      // operands (the sum over n batches overcounts by n). A LoopCarried
      // mode keeps the existing slot-0 fallback: it always carries a real
      // position, so home_scope + compute_dag_boulevard resolve it in the
      // ordinary case and the fallback is purely defensive there.
      if (ac.role == LoopRole::Reduction && fs < 0)
        throw Exception(
            "build_ordered_schedule: value " + std::to_string(vid) +
            " reduces a batched mode with no loop identity (no reduced_slot "
            "stamped); the escape cannot be placed");
      auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
      SEQUANT_ASSERT(d.has_value());  // types was built from this same union
      OutputKind const kind =
          (ac.role == LoopRole::Reduction)
              ? OutputKind::AccumulateSum
              : OutputKind::AccumulateScatter;  // LoopCarried
      // Each per-instance mode escapes to its own depth (distinct loop_slot),
      // so same-space modes keep separate escapes; the dedup below
      // only merges a genuine reduce+carry pair that lands at one depth.
      auto it = std::find_if(escapes.begin(), escapes.end(),
                             [&](auto const& e) { return e.first == *d; });
      if (it == escapes.end())
        escapes.push_back({*d, kind});
      else if (kind == OutputKind::AccumulateScatter)
        it->second = OutputKind::AccumulateScatter;
    }

    // A value with >1 non-local mode of one space (a doubles amplitude / PPL
    // product carrying two occ externals) has its distinct per-instance
    // escapes collapsed to fewer escapes (one per depth == one per space).
    // Dump those cells (SEQUANT_DUMP_SCHEDULE) to expose the collapse.
    if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE")) {
      container::svector<std::pair<Index, wchar_t const*>> nonlocal;
      for (AxisClass const& ac : cl.per_axis)
        if (ac.role != LoopRole::LoopLocal)
          nonlocal.push_back(
              {ac.axis, ac.role == LoopRole::Reduction ? L"R" : L"C"});
      bool collapse = false;
      for (std::size_t a = 0; a < nonlocal.size() && !collapse; ++a)
        for (std::size_t b = a + 1; b < nonlocal.size(); ++b)
          if (nonlocal[a].first.space().base_key() ==
              nonlocal[b].first.space().base_key()) {
            collapse = true;
            break;
          }
      if (collapse) detail::dump_sched_collapse(cl.hash, nonlocal, escapes);
    }

    // Mixed-pass member (per-nest forced split design, sections 3.3 and 7).
    // A reader of a later pass produced inside this value's nest is the only
    // reader that can see a per-batch slice from another traversal; it must
    // read a full form that resides at root. The value keeps its Build step
    // at its production site and gains an AccumulateScatter escape at every
    // instance of its nest it is loop-local on and no role already escapes
    // it as a Scatter (a reduced instance sums, a carried one scatters; a
    // loop-local mode sharing a depth with a Reduction dominates it to a
    // Scatter too -- its own batches are disjoint, so summing them would be
    // silently wrong). Every instance the value is sliced by is then
    // escaped, so the outermost assembled form is bound to no enclosing
    // instance and the table's residency rule places it at root; a level the
    // value is invariant to is simply skipped.
    bool materialized_across_split = false;
    std::optional<std::size_t> const home_depth = local_home_depth(cl);
    // Whether this value has a LoopLocal instance of its own nest that is
    // not covered by a role-driven escape (before rule 4 adds anything) --
    // the invariant the outside-nest tripwire below actually needs. A
    // per-batch-only instance like that is never delivered to root, so a
    // later-pass reader outside the nest contradicts legality regardless of
    // whether some other instance of this same value happens to be
    // role-escaped elsewhere, at a different depth (a coarser gate on
    // "any role escape at all" would miss exactly this two-different-depth
    // case).
    bool unescaped_local_instance = false;
    if (home_depth) {
      std::size_t const home_nest = type_cluster[*home_depth];
      for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
        if (cl.per_axis[pos].role != LoopRole::LoopLocal) continue;
        std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
        int const fs = fusion_slot(cl, pos);
        auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
        if (!d || type_cluster[*d] != home_nest) continue;
        bool const escaped =
            std::any_of(escapes.begin(), escapes.end(),
                        [&](auto const& e) { return e.first == *d; });
        if (!escaped) {
          unescaped_local_instance = true;
          break;
        }
      }
    }
    // Direct readers produced in the same nest with a later pass. Nest
    // membership is decided by production_depth, not by local_home_depth: a
    // reader with only carried/reduction roles is still produced per batch
    // inside its own nest even though it reports no LoopLocal home.
    auto const later_same_nest_readers =
        [&](std::size_t nest) -> container::svector<std::size_t> {
      container::svector<std::size_t> out;
      auto const cons_it = g.consumers_of.find(vid);
      if (cons_it == g.consumers_of.end()) return out;
      for (std::size_t u : cons_it->second) {
        if (pass_of(u) <= pass_of(vid)) continue;
        auto const uit = cl_by_vid.find(u);
        if (uit == cl_by_vid.end()) continue;
        auto const uh = production_depth(*uit->second);
        if (uh && type_cluster[*uh] == nest) out.push_back(u);
      }
      return out;
    };
    if (home_depth) {
      std::size_t const nest = type_cluster[*home_depth];
      auto const readers = later_same_nest_readers(nest);
      if (!readers.empty()) {
        // Rule 4 (section 7.3): escape every instance of this nest the
        // value is loop-local on and not already escaped by a role as a
        // Scatter; a role escape already scattering (LoopCarried) that
        // instance is left as-is, a role escape summing it (Reduction) is
        // upgraded to a Scatter (a loop-local mode's batches are disjoint,
        // so summing them across a depth it also shares with a reduced mode
        // would silently combine values that must stay separate -- the same
        // dominance the role loop above applies to a reduce+carry pair), and
        // an instance the value is invariant to (its depth does not
        // resolve, or resolves outside this nest) is simply skipped --
        // completeness then falls out for free, since the value is sliced
        // by exactly the instances this loop escapes.
        for (std::size_t pos = 0; pos < cl.per_axis.size(); ++pos) {
          std::wstring const bk{cl.per_axis[pos].axis.space().base_key()};
          int const fs = fusion_slot(cl, pos);
          auto const d = depth_of_instance(bk, fs >= 0 ? fs : 0);
          if (!d || type_cluster[*d] != nest) continue;
          auto it = std::find_if(escapes.begin(), escapes.end(),
                                 [&](auto const& e) { return e.first == *d; });
          if (cl.per_axis[pos].role == LoopRole::LoopLocal) {
            if (it == escapes.end()) {
              escapes.push_back({*d, OutputKind::AccumulateScatter});
              materialized_across_split = true;
            } else if (it->second == OutputKind::AccumulateSum) {
              it->second = OutputKind::AccumulateScatter;
            }
          } else if (it == escapes.end()) {
            // Should be impossible: every non-LoopLocal per_axis mode was
            // already pushed into escapes, for every depth it resolves to,
            // by the role loop above.
            if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
              dump_reject(vid, nest, readers);
            throw Exception(
                "build_ordered_schedule: value " + std::to_string(vid) +
                " is read by value " + std::to_string(readers.front()) +
                " in a later pass of its nest, but its instance at depth " +
                std::to_string(*d) +
                " is neither loop-local nor escaped by a role");
          }
        }
      } else if (unescaped_local_instance) {
        // Tripwire (controller ruling I3; reader test corrected by ruling
        // I4): a value with an unescaped LoopLocal instance of its own nest
        // (checked before rule 4 above ran, via unescaped_local_instance)
        // has a per-batch-only form of that instance that is never
        // delivered to root -- and there is no same-nest later-pass reader
        // to trigger rule 4 above and fix it (this is the `readers.empty()`
        // branch). A direct later-pass reader whose production site
        // resolves to a nest other than this one is the reader's location,
        // not this value's escapes: it cannot see the per-batch home form,
        // and legality and the schedule disagree, regardless of whether
        // some other instance of this value happens to be role-escaped
        // elsewhere. A value whose every LoopLocal instance of its own nest
        // is already role-escaped is exempt: each such escape already
        // assembles a full form with root residency (rule 4's own
        // scatter-dominance above ensures no instance is left half-summed),
        // which any later-pass reader, in any nest, can see.
        // "Produced outside this nest" is decided by production_depth, not
        // by local_home_depth: a reader with only carried/reduction roles --
        // a forest root delivered in full, or a carried value of a later
        // pass -- is still produced per batch inside its own nest, so it
        // must not spuriously trip this guard.
        //
        // production_depth never guesses: a reader whose production depth
        // does not resolve at all (every one of its modes has an
        // unresolvable fusion slot, or it has no per_axis modes) is neither
        // confidently inside this nest nor confidently outside it, so it
        // neither trips this guard nor counts as an in-nest reader -- the
        // guard fires only for a confidently resolved different nest. An
        // unresolved reader is not silently accepted, either: the table
        // validator's visibility rule is the net that catches a reader the
        // builder could not locate. Thrown loudly, for the case this guard
        // does catch, rather than producing a silent mis-schedule.
        auto const cons_it = g.consumers_of.find(vid);
        if (cons_it != g.consumers_of.end())
          for (std::size_t u : cons_it->second) {
            if (pass_of(u) <= pass_of(vid)) continue;
            auto const uit = cl_by_vid.find(u);
            std::optional<std::size_t> const uh =
                uit == cl_by_vid.end() ? std::nullopt
                                       : production_depth(*uit->second);
            if (!uh)
              continue;  // unresolved: neither trips nor counts as
                         // inside the nest
            if (type_cluster[*uh] == nest) continue;  // same nest: n/a
            if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
              dump_reject(vid, nest, container::svector<std::size_t>{u});
            throw Exception("build_ordered_schedule: value " +
                            std::to_string(vid) + " (pass " +
                            std::to_string(pass_of(vid)) +
                            ") is read by value " + std::to_string(u) +
                            " (pass " + std::to_string(pass_of(u)) +
                            ") outside its nest (outermost depth " +
                            std::to_string(cluster_min.at(nest)) +
                            "): legality and the schedule disagree");
          }
      }
    }

    // Where the value's plain BuildStep sits: \c home_depth (the innermost
    // loop it is LoopLocal on), unless it is materialized across the split
    // (below) and its escape chain reaches deeper than that -- a role
    // escape nested inside the LoopLocal home (a Reduction axis, say) --
    // in which case the deepest site on that chain is the true production
    // site and home's rule-4 escape is pure forwarding, like any other link
    // in the chain.
    // Completeness invariant (design section 9.2): a value
    // reduced over a loop instance (an AccumulateSum escape at depth d) is
    // complete only after that loop closes; a reader produced inside that
    // instance (production depth at or below d in the same nest) in the
    // same pass would be served the current batch's partial sum. The pass
    // levels above (2a) now bump exactly such a reader to a later pass via
    // the Reduction-source bumping edge, so this shape should never survive
    // to here; this check stays as a loud tripwire on the levels
    // themselves -- its firing means the levels failed to bump a read they
    // should have, a builder defect, not an expected outcome.
    for (auto const& [d, kind] : escapes) {
      if (kind != OutputKind::AccumulateSum) continue;
      auto const cons_it = g.consumers_of.find(vid);
      if (cons_it == g.consumers_of.end()) continue;
      for (std::size_t u : cons_it->second) {
        auto const uit = cl_by_vid.find(u);
        if (uit == cl_by_vid.end()) continue;
        auto const pd = production_depth(*uit->second);
        if (!pd || *pd < d || type_cluster[*pd] != type_cluster[d]) continue;
        if (pass_of(u) != pass_of(vid)) continue;
        throw Exception(
            "build_ordered_schedule: value " + std::to_string(vid) +
            " is reduced over the loop at depth " + std::to_string(d) +
            " but value " + std::to_string(u) +
            " reads it inside that loop in the same pass (pass " +
            std::to_string(pass_of(u)) +
            "): the read would see a partial sum; a reader of a reduction "
            "inside its own loop needs a later pass of that loop "
            "(unsupported)");
      }
    }
    std::optional<std::size_t> build_depth = home_depth;
    if (!escapes.empty()) {
      for (auto const& [d, kind] : escapes)
        buckets[d].outputs.push_back({vid, kind});
      // A value that escapes by its own per-axis roles has no BuildStep: its
      // production is the accumulation itself, at the deepest escape site
      // (the multi-level chain's bottom-up assembly: raw production at the
      // deepest site, pure forwarding at every shallower one). One
      // materialized across the split by rule 4 above is different: it is
      // produced, at the deepest site of its full chain (role escapes and
      // the rule-4 scatter together), which same-pass consumers read inside
      // its own nest, per batch. It keeps its BuildStep there, so that
      // block both builds it (for its same-pass in-nest readers, and as the
      // per-batch input of the rest of its chain) and lists it as an output
      // (for the later-pass reader). `well_formed` admits exactly this
      // shape: every block that lists a value in `outputs` either holds its
      // BuildStep or is an ancestor of the one that does -- always true
      // here since the BuildStep sits at the chain's deepest site and every
      // other site is, by construction, an ancestor of it.
      if (!materialized_across_split) continue;
      for (auto const& [d, kind] : escapes) {
        (void)kind;
        if (!build_depth || d > *build_depth) build_depth = d;
      }
      materialized_across_split_ids.push_back(vid);
    }

    if (build_depth)
      buckets[*build_depth].build_ids.push_back(vid);
    else
      root_build_ids.push_back(vid);
  }

  if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE") &&
      !materialized_across_split_ids.empty())
    detail::dump_sched_materialize(materialized_across_split_ids);

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
    // The block's kind is the kind of its loop instance (the open that
    // created the (space, slot) component -- RichSchedule::loop_kind); a
    // space can hold both a contracted-in-batches instance and an external
    // one, so the per-space test (does the space appear on a root result?)
    // is only the fallback for an instance no open recorded.
    if (auto const kit = rich.loop_kind.find(
            {std::wstring{axis.space().base_key()}, loop_slot});
        kit != rich.loop_kind.end())
      block.kind = kit->second;
    else
      block.kind = detail::mode_is_external(rich, axis)
                       ? BatchModeType::External
                       : BatchModeType::Contracted;
    block.steps =
        detail::ordered_schedule_topo_sort_steps(std::move(items), meta);
    block.outputs.assign(outputs.begin(), outputs.end());
    return block;
  };

  // 3. Assemble the chain bottom-up (innermost first). Thread up a list of
  // child block Steps to embed one level out -- normally one (the single block
  // for the deeper axis), but one per pass at a nest's
  // outermost depth when that nest holds a forced-split axis (its pass
  // blocks) -- each paired with the meta the outer topo-sort needs (its own
  // escape outputs as `produced`, its whole subtree's external need as
  // `requires_`, and a deterministic `tie_key`). child_produced_all /
  // child_requires_all carry the full recursive produced/external-need sets of
  // everything built at this depth (identical whichever way -- the outer
  // level sees the same production/need set either way) to grow the next
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
  if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
    detail::dump_sched_nests(nest_passes, cluster_min,
                             [&](std::size_t c, int k) {
                               std::size_t cnt = 0;
                               for (std::size_t d = 0; d < n; ++d) {
                                 if (type_cluster[d] != c) continue;
                                 for (std::size_t v : buckets[d].build_ids)
                                   if (pass_of(v) == k) ++cnt;
                               }
                               return cnt;
                             });

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
            // Recursive produced (not just this pass's own outputs): a pass
            // block is now always a root-level sibling of every other nest
            // (the split always lands at the nest's outermost depth), so an
            // under-reported `produced` orders a sibling nest that requires a
            // value this pass builds (no escape of its own) before this pass
            // runs -- the same read-before-build hazard the single-block
            // path's `produced_all` comment already explains.
            m.produced.assign(pass_produced.begin(), pass_produced.end());
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
      // Recursive produced (not just this block's own outputs): a nest
      // advertises everything it produces so the topo sort can order a sibling
      // nest that consumes a value produced by an inner block of this one. With
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
  // Finalize the last cluster's nest.
  for (std::size_t k = 0; k < pending_steps.size(); ++k) {
    finished_steps.push_back(std::move(pending_steps[k]));
    finished_metas.push_back(std::move(pending_metas[k]));
  }

  // Root assembly: root-level BuildStep's plus the separate cluster nests as
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

  if (detail::dump_enabled("SEQUANT_DUMP_SCHEDULE"))
    detail::dump_schedule_tree(out.root, 0);
  // Refusal: the builder's own structural self-check. A throw, not a
  // SEQUANT_ASSERT -- the latter is compiled out in non-Debug builds, which
  // would hand an ill-formed schedule to the executor in exactly the builds
  // that run production work.
  if (!well_formed(out))
    throw Exception(
        "build_ordered_schedule: the schedule it built is not well-formed "
        "(see OrderedSchedule::well_formed)");
  return out;
}

namespace detail {

/// \brief One enclosing loop block on a value's realized eval scope: the
/// block's canonical representative \c axis and its \c DagScopeLevel (\c
/// depth / \c space / \c ordinal), read straight off \c ordered's block tree
/// (so the \c ordinal of a forced-split sibling is the real one the runtime
/// pushes, not a guessed one).
struct ScopeBlockAxisLevel {
  Index axis;
  DagScopeLevel level;
};

/// \brief \c compute_sliced_mode_assignment's block-tree walk: record, for
/// every value produced (a \c BuildStep) or escaped (a block \c outputs
/// entry), the ordered list of enclosing loop blocks (\c axis + \c level) it
/// is evaluated under -- i.e. the loops whose batch it reads its operands
/// inside. This is the build-scope walk itself, keeping each block's \c
/// DagScopeLevel (not just its axis) so the map can name the runtime \c
/// BatchContextEntry::level.
///
/// \p enc is the root-to-\p block path including \p block's own (axis, level)
/// -- a value built or escaping inside \p block reads its operands inside
/// \p block's own loop, so the block's own axis encloses that read (this walk
/// always includes the block axis; only the home of an escape sits one level
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
    // escape is built inside `block` (reads operands in this loop); its home is
    // one level out, but we want the read/fetch scope here. Deepest-escape
    // wins: a value carried on more than one nested axis escapes at each depth,
    // but the outer escapes only scatter an already-built value (read no
    // operands); the contraction runs at the innermost escape, so its operands
    // are sliced there. The child recursion above visits the deeper block
    // first, so a shallower outer escape must not clobber the deeper scope --
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
/// realized \c ScopeBlock's \c DagScopeLevel into one \c LoopId per distinct
/// level) leans on \c (level.depth, level.space, level.ordinal) being unique
/// globally across the whole realized tree -- i.e. any two blocks sharing
/// that triple must also share the same representative \c axis, or the
/// enumeration could fold two structurally-different loops onto one \c
/// LoopId. \c well_formed only checks ordinal uniqueness among a block's own
/// same-axis direct children (sibling-local, see \c
/// ordered_schedule_block_well_formed above) -- it says nothing about two
/// blocks at the same (depth, space, ordinal) that are not siblings (say,
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
/// \brief A stable identifier for one DAG-scope loop realized by \c
/// build_ordered_schedule -- an index into \c SlicedModeAssignment::levels,
/// the schedule's own canonical (deterministic pre-order) enumeration of
/// every non-root \c ScopeBlock's \c DagScopeLevel. Two blocks that are
/// structurally the same loop (identical \c (depth, space, ordinal)) always
/// share one \c LoopId; two blocks that differ in any of those three --
/// including two of a nest's own pass blocks (one per pass, latitude =
/// pass), which differ only in \c latitude_ordinal -- get distinct ids by
/// construction: a nest's pass blocks must be distinguishable colors, not
/// folded (as-built design section 5.1, \c
/// doc/dev/specs/2026-09-12-batched-array-dag-eval-as-built.md).
///
/// \note Defined in \c dag_scope.hpp so the low-level DAG-scope types can be
/// named without depending on this schedule header; re-exported here so the
/// schedule-side code can name it as \c eval::LoopId.
using sequant::LoopId;

///
/// \brief The per-(value, sliced-mode) -> DAG-scope-loop assignment: the
/// coloring input fed to \c canonicalize_slots's \c NamedIndexColorMap.
/// This is plain data keyed by the value's own physical \c Index
/// label for each mode it is sliced on -- not a per-cell position map -- so a
/// relabeled CSE participant is keyed by
/// its own label, and a symmetric value's two occurrence-bound physical slots
/// are both recoverable (one entry per distinct (value, Index) pair actually
/// sliced).
struct SlicedModeAssignment {
  /// Every distinct DAG-scope loop realized anywhere in the schedule, in
  /// canonical (deterministic pre-order over the block tree) order; a
  /// value's assigned \c LoopId is its position in this list. The root block
  /// itself (sentinel axis, outside every loop) is never an entry here.
  container::vector<DagScopeLevel> levels;

  /// Consumer-attributed per-occurrence sliced-mode facts (sliced-value
  /// canonical-layout / loop-coloring design, pillar 2): each entry is
  /// (value_id, this occurrence's own sliced-mode Index, the slicing LoopId,
  /// the consumer value_id -- the use-site whose fetch of \c value_id binds the
  /// loop to that Index). Recorded only by the regime-2 (occurrence-driven)
  /// pass, the one pass that can attribute a stamp to a specific occurrence and
  /// hence to a specific consumer. This is the raw material the cell table
  /// builder (cell_table_builder.hpp) consumes to disambiguate the
  /// symmetric case (one value, one loop, two free modes bound by two
  /// different consumers): a consumer-blind (value, mode) map cannot express
  /// "pos0 here, pos1 there" because it folds away which occurrence bound
  /// which mode.
  /// (value_id, this occurrence's own sliced-mode physical position -- its
  /// index in occ.carried, computed in that occurrence's own index-frame,
  /// LoopId, consumer value_id). The position (not an Index label) is what the
  /// table builder consumes, so the runtime never re-matches a label across
  /// index-frames.
  /// The fifth element is the operand's leg at the consumer (its index in
  /// the consumer occurrence's \c operand_points; \c npos when unknown): a
  /// value read on both legs of one consumer under different labels has
  /// facts per leg, and the table builder pairs each leg's read with its own.
  container::svector<
      std::tuple<std::size_t, std::size_t, LoopId, std::size_t, std::size_t>>
      occ_facts;

  /// Explicit per-occurrence invariant facts: (value_id, LoopId, consumer
  /// value_id) triples recording that this consumer's fetch of the value is
  /// correctly unsliced on this loop -- the consumer is building that loop's
  /// batch, but the value's occurrence in this consumer'S frame does not carry
  /// the loop's mode (a CSE-shared value can carry the mode in one frame and
  /// not another). Recording the negative decision lets the runtime tell
  /// "correctly invariant" apart from "no decision recorded" (a real gap), so
  /// the completeness guard fires only on the latter.
  /// Fourth element: the operand's leg at the consumer, as in \c occ_facts.
  container::svector<std::tuple<std::size_t, LoopId, std::size_t, std::size_t>>
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
/// i.e. no two different blocks visited here ever carry an equal \c level.
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
/// \brief
/// build the per-(value, sliced-mode) -> DAG-scope-loop \c
/// SlicedModeAssignment from an already-built \p ordered (must be \c
/// build_ordered_schedule(rich, ...)'s return value for this same \p rich).
///
/// \details Uses a two-pass fetch-site walk (\c populate_build_scope_walk
/// for the enclosing-loop scope of every produced/escaped value, and \c
/// ordered_schedule_dep_graph for the operand edges, leaves included): an
/// exact pass matching a block's representative axis against a value's own
/// carried \c Index (step (3) below), then a regime-2 relabel pass (step (4)
/// below) recovering a CSE value's own label from its occurrences when it
/// was canonicalized independently of the block's representative axis (the
/// same physical loop, relabeled) -- but records the raw facts keyed by the
/// value's own \c Index label (not a \c ValueCell::carried position),
/// consistency-checked the same way (two fetch sites disagreeing on one
/// value's mode's level is a scheduler bug, never resolved by averaging or
/// last-write-wins), then remapped through the canonical \c LoopId
/// enumeration (\c detail::enumerate_realized_levels) instead of storing the
/// \c DagScopeLevel directly -- this is what makes forced-split siblings
/// (same depth/space, different ordinal) come out as distinct colors: they
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

  // (3) per-occurrence positional pass (the loop-open design). For each
  // occurrence occ of a value W consumed by C, the loops the runtime crosses
  // when fetching W are C's enclosing DAG blocks
  // build_scope[owner(consumer_point)], outermost first. occ.ectx -- built
  // from loop-opens (peak_profile), so one physical loop appears exactly once
  // -- names those same loops in W's own
  // frame, outermost first. Pairing them by nest position gives, per realized
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
  (void)g;  // the dependency graph is not consulted by this pass
  std::unordered_map<std::size_t, std::size_t> point_owner;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences)
      point_owner[occ.point] = vc.value_id;

  // point -> that occurrence record (to reach the parent occurrence on an
  // edge: the consumer's occurrence in the same tree as the operand's).
  std::unordered_map<std::size_t, OccurrenceRec const*> point_occ;
  for (ValueCell const& vc : rich.cells)
    for (OccurrenceRec const& occ : vc.occurrences) point_occ[occ.point] = &occ;

  // The value-level component slot per carried position: the loop_slot the
  // union-find assigned to (value, position), recovered unmasked.
  // compute_dag_boulevard stamps occ.loop_slot = -1 wherever a mode is not
  // batched at that occurrence (i.e. the value was produced whole in that mode
  // in that term). But a value produced whole is still read sliced wherever a
  // consumer sits inside that mode's realized loop -- the atlas must key the
  // slice fact off which loops the consumer is in and which modes the value
  // carries, not off the per-occurrence production mask. So recover the
  // value-level slot as any non-(-1) loop_slot across the value's occurrences
  // at each position.
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
      // production edge only: a consumer value is produced once,
      // from its canonical (front) occurrence's tree; its other occurrences
      // are reads of the finished value, and the W->C edges inside those
      // trees never execute. Slice facts recorded from such a non-production
      // edge are keyed by the same (W, C, loop) triple and merge with the
      // production edge's facts -- in a permuted frame that binds C's modes
      // to other loop instances, which sliced one operand position by two
      // loops (observed on the water-20 strict dry-run walk: one value's
      // operand position 1 was sliced by two different loop slots at once).
      //
      // This guard is also what keeps the order-independent Product identity
      // (canonical_children, eval_node_compare.hpp) safe. Two spellings of one
      // contraction -- (X,Y) in one term, (Y,X) in another -- now fold into a
      // single value, so a consumer cell can own occurrences whose trees carry
      // X and Y on opposite legs. Occurrence facts are recorded only from the
      // consumer's canonical (front) production occurrence, and
      // cell_table_builder derives an operand's `my_leg` from those
      // (operand, consumer) facts, so every leg comes from one tree and the
      // legs cannot be crossed. Relaxing the guard to admit non-production
      // edges would mix the two spellings' legs and mis-pair them.
      if (rich.cells[oit->second].occurrences.empty() ||
          rich.cells[oit->second].occurrences.front().point !=
              occ.consumer_point)
        continue;
      auto const sit = build_scope.find(oit->second);
      if (sit == build_scope.end()) continue;  // consumer has no enclosing loop
      container::svector<detail::ScopeBlockAxisLevel> const& scope =
          sit->second;
      if (scope.empty()) continue;
      // The atlas: for each of the consumer's enclosing loops `scope[k]`
      // (space S, loop_slot L), the mode of W it slices is W's carried position
      // whose value-level component slot is L and whose space is S. A loop-
      // identity match in W's own frame (occurrence-invariant), keyed off the
      // value-level slot -- so a value produced whole in a mode is still sliced
      // there when this consumer sits inside that mode's realized loop (the
      // measured i_3-contracted case). Over-nested loops the value carries no
      // mode of get no fact, correct (the value is invariant to them).
      auto const& vs = value_slot[w_vid];
      // This occurrence's leg at its consumer (left 0 / right 1).
      std::size_t leg = static_cast<std::size_t>(-1);
      {
        auto const& cpts =
            rich.cells[oit->second].occurrences.front().operand_points;
        for (std::size_t li = 0; li < cpts.size(); ++li)
          if (cpts[li] == occ.point) leg = li;
      }
      for (std::size_t k = 0; k < scope.size(); ++k) {
        DagScopeLevel const& lvl = scope[k].level;
        bool self_sliced = false;
        for (std::size_t pos = 0; pos < occ.carried.size() && pos < vs.size();
             ++pos) {
          // Per-occurrence slot: value_slot folds every occurrence's
          // loop_slot into one per-value vector, so a CSE-shared value whose
          // mode is sliced by different loop instances in different
          // occurrences (one consumer slicing a mode under a K-loop slot,
          // another under a second slot for a K-reduction) keeps only one slot
          // and silently misses the other occurrence (served whole under a K
          // batch -> TA sparse-gemm out-of-bounds read). The occurrence's own
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
              std::make_tuple(w_vid, pos, id_of(lvl), oit->second, leg));
          self_sliced = true;
          break;  // one W-position per enclosing loop
        }
        if (self_sliced) continue;
        // Use-induced slicing (Layer 2): W is whole-produced on this loop's
        // space (no own sliced slot), but this consumer C is sliced on a mode M
        // there. If M is a shared external W also carries, the C = ...*W
        // contraction binds W's M to C's sliced M -- so W must be sliced on M
        // at this loop for this fetch. Record it consumer-keyed in occ_facts,
        // never consumer-blind: W may be
        // CSE-shared between C (sliced here) and a different consumer that
        // reads it whole (invariant), and a blind fact would wrongly slice it
        // for both. Bounded to a mode C actually slices (conformability), not
        // every carried-mode coincidence, so a genuinely invariant shared
        // operand is not sliced.
        std::size_t const c_vid = oit->second;
        auto const& cvs = value_slot[c_vid];
        auto const& c_carried = rich.cells[c_vid].carried;
        bool recorded = false;
        for (std::size_t cp = 0; cp < c_carried.size() && cp < cvs.size();
             ++cp) {
          // The consumer's slot for position cp in this tree: the parent
          // occurrence on this edge (point == occ.consumer_point). Mixing C's
          // other occurrences (other trees, permuted frames) leaked a slot
          // binding from a different frame and sliced one operand position
          // by two loops (observed on the water-20 strict dry-run walk: one
          // value's operand position 1 was sliced by two different loop slots
          // at once). The per-value value_slot is the fallback only when the
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
          if (pos_a == occ.carried.size()) {
            // W does not carry C's sliced mode in this occurrence's frame: its
            // fetch by C is correctly unsliced on this loop. Record that
            // explicitly so the guard does not mistake it for a gap.
            result.occ_invariant.push_back(
                std::make_tuple(w_vid, id_of(lvl), oit->second, leg));
            recorded = true;
            continue;
          }
          result.occ_facts.push_back(
              std::make_tuple(w_vid, pos_a, id_of(lvl), oit->second, leg));
          recorded = true;
          if (detail::dump_enabled("SEQUANT_DUMP_USEINDUCED"))
            detail::dump_use_induced(rich.cells[w_vid].hash % 100000, M, pos_a,
                                     id_of(lvl),
                                     rich.cells[c_vid].hash % 100000);
          break;  // one shared mode per enclosing loop
        }
        // Reduction consumer: C reduces a mode under this loop
        // (its result does not carry it, so it is in no carried/value_slot
        // position) -- C is building a per-batch partial of that reduction,
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
              if (pos_a == occ.carried.size())
                result.occ_invariant.push_back(
                    std::make_tuple(w_vid, id_of(lvl), oit->second, leg));
              else
                result.occ_facts.push_back(std::make_tuple(
                    w_vid, pos_a, id_of(lvl), oit->second, leg));
              recorded = true;
              break;
            }
            if (recorded) break;
          }
        }
        // Nothing sliced this operand here and nothing declared it invariant:
        // record the invariant explicitly so the completeness guard (whose
        // consumer oracle is a per-value union over C's occurrences) does not
        // mistake a legitimately whole read for a gap. Range conformance is
        // enforced by the dry-run/wet backends themselves.
        if (!recorded)
          result.occ_invariant.push_back(
              std::make_tuple(w_vid, id_of(lvl), oit->second, leg));
      }
    }
  }

  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
