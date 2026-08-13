#ifndef SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
#define SEQUANT_EVAL_ORDERED_SCHEDULE_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <limits>
#include <memory>
#include <optional>
#include <string>
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
  Index axis{};     //!< the loop axis; default (sentinel) on the root block.
  int ordinal = 0;  //!< disambiguates recurring sibling blocks realizing the
                    //!< same axis (e.g. two disjoint passes over the same
                    //!< axis TYPE at one nesting level).
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
      if (ci->axis.space().base_key() == cj->axis.space().base_key() &&
          ci->ordinal == cj->ordinal)
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

}  // namespace detail

///
/// \brief Structural sanity check on \p sched (no sequencer logic):
///   - every \c BuildStep::value_id is < \c sched.num_values;
///   - ordinals are unique among same-axis (\c IndexSpace::base_key())
///     sibling blocks within a parent;
///   - every \c ScopeBlock::outputs value_id is < \c sched.num_values;
///   - SINGLE-PRODUCER (SSA-like): across the WHOLE schedule, no value_id
///     appears more than once total among every \c BuildStep and every
///     block's \c outputs entries combined -- a value built directly must
///     not also be a loop output, a loop-accumulated value must not also be
///     built directly, and no value may be produced twice.
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

  container::vector<std::size_t> production_ids;
  detail::collect_production_ids(sched.root, production_ids);
  std::sort(production_ids.begin(), production_ids.end());
  return std::adjacent_find(production_ids.begin(), production_ids.end()) ==
         production_ids.end();
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
template <typename node_t>
[[nodiscard]] inline std::function<std::size_t(std::size_t)> ordered_home_reads(
    OrderedSchedule const& ordered, RichSchedule const& rich,
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
  auto build_scope = std::make_shared<
      std::unordered_map<std::size_t, container::svector<Index>>>();
  auto home_scope = std::make_shared<
      std::unordered_map<std::size_t, container::svector<Index>>>();
  auto const walk = [&build_scope, &home_scope](
                        auto&& self, ScopeBlock const& b,
                        container::svector<Index> const& enc) -> void {
    for (Step const& s : b.steps) {
      if (auto const* build = std::get_if<BuildStep>(&s.value)) {
        (*build_scope)[build->value_id] = enc;
        (*home_scope)[build->value_id] = enc;
      } else if (auto const* child = std::get_if<ScopeBlock>(&s.value)) {
        container::svector<Index> inner = enc;
        inner.push_back(child->axis);  // child's steps are inside child's loop
        self(self, *child, inner);
      }
    }
    // An escape output is BUILT inside b (build_scope == enc, reading its
    // operands once per batch) but HOMES one level OUT (home_scope == enc minus
    // b's own axis). Root (sentinel axis, enc=={}) has no outputs.
    for (auto const& [vid, kind] : b.outputs) {
      (void)kind;
      (*build_scope)[vid] = enc;
      container::svector<Index> out_scope = enc;
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
             std::size_t, container::svector<Index>>> const& m,
         std::size_t vid) -> container::svector<Index> {
    auto it = m->find(vid);
    return it == m->end() ? container::svector<Index>{} : it->second;
  };
  auto const type_in = [](container::svector<Index> const& s, Index const& m) {
    for (Index const& x : s)
      if (x.space().base_key() == m.space().base_key()) return true;
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
    container::svector<Index> const scopeW = scope_of(build_scope, wc.value_id);
    auto const add_edge = [&](node_t const& child) {
      auto const cvit = hash_to_vid.find(child->hash_value());
      if (cvit == hash_to_vid.end()) return;
      std::size_t const cvid = cvit->second;
      container::svector<Index> const homeC = scope_of(home_scope, cvid);
      std::size_t prod = 1;
      for (Index const& L : scopeW)
        if (!type_in(homeC, L)) prod *= n_blocks(L);
      (*reads)[cvid] += prod;
    };
    add_edge(wit->second.left());
    add_edge(wit->second.right());
  }

  return [reads](std::size_t vid) -> std::size_t {
    auto it = reads->find(vid);
    return (it == reads->end() ? std::size_t{0} : it->second) + 1;
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

  // 1. The canonical chain: one representative Index per distinct axis TYPE
  // present in ANY cell's per_axis (LoopLocal, Reduction, OR LoopCarried --
  // NOT just home_floor; see the function doc comment's part 1).
  container::svector<Index> types;
  for (CellLegality const& cl : legality.cells)
    for (AxisClass const& ac : cl.per_axis) {
      auto const& bk = ac.axis.space().base_key();
      bool const seen = std::any_of(
          types.begin(), types.end(),
          [&](Index const& ix) { return ix.space().base_key() == bk; });
      if (!seen) types.push_back(ac.axis);
    }

  auto const rank_of = [&](Index const& ix) -> std::size_t {
    std::size_t i = 0;
    for (auto const& key : mode_order) {
      if (ix.space().base_key() == key) return i;
      ++i;
    }
    return static_cast<std::size_t>(-1);
  };
  std::sort(types.begin(), types.end(), [&](Index const& a, Index const& b) {
    auto const ra = rank_of(a), rb = rank_of(b);
    if (ra != rb) return ra < rb;
    return a.space().base_key() < b.space().base_key();
  });

  std::size_t const n = types.size();
  auto const depth_of_type =
      [&](Index const& ix) -> std::optional<std::size_t> {
    auto const& bk = ix.space().base_key();
    for (std::size_t d = 0; d < n; ++d)
      if (types[d].space().base_key() == bk) return d;
    return std::nullopt;
  };

  // 2. Per-value placement: home BuildStep (root-level bucket uses index n
  // as a sentinel "root" depth) vs. escape output at the deepest escaped
  // axis's depth.
  container::vector<detail::OrderedScheduleDepthBucket> buckets(n);
  container::svector<std::size_t> root_build_ids;

  for (CellLegality const& cl : legality.cells) {
    auto const vid_it = value_id_of.find(cl.hash);
    SEQUANT_ASSERT(vid_it != value_id_of.end());
    std::size_t const vid = vid_it->second;

    std::optional<std::size_t> escape_depth;
    OutputKind escape_kind = OutputKind::Transient;  // overwritten if found
    for (AxisClass const& ac : cl.per_axis) {
      if (ac.role == LoopRole::LoopLocal) continue;
      auto const d = depth_of_type(ac.axis);
      SEQUANT_ASSERT(d.has_value());  // types was built from this same union
      if (!escape_depth || *d > *escape_depth) {
        escape_depth = d;
        escape_kind = (ac.role == LoopRole::Reduction)
                          ? OutputKind::AccumulateSum
                          : OutputKind::AccumulateScatter;  // LoopCarried
      }
    }

    if (escape_depth) {
      buckets[*escape_depth].outputs.push_back({vid, escape_kind});
      continue;
    }

    // Plain BuildStep: home depth by SET equality of the accumulated
    // (root-to-depth) TYPE set against home_floor's TYPE set (root if
    // home_floor is empty), mirroring build_scope_schedule's step 4.
    container::svector<std::wstring> want;
    for (Index const& m : cl.home_floor) {
      auto const& bk = m.space().base_key();
      if (std::find(want.begin(), want.end(), bk) == want.end())
        want.push_back(bk);
    }

    std::optional<std::size_t> target;
    if (!want.empty()) {
      container::svector<std::wstring> enclosing;
      for (std::size_t d = 0; d < n; ++d) {
        enclosing.push_back(types[d].space().base_key());
        if (enclosing.size() == want.size() &&
            std::is_permutation(enclosing.begin(), enclosing.end(),
                                want.begin()))
          target = d;
      }
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
  std::optional<detail::ForcedSplitPasses> split_passes;
  if (n > 0) {
    std::wstring const inner_key{types[n - 1].space().base_key()};
    bool const inner_forced =
        std::any_of(legality.cells.begin(), legality.cells.end(),
                    [&](CellLegality const& cl) {
                      for (Index const& ix : forced_split_types(cl))
                        if (ix.space().base_key() == inner_key) return true;
                      return false;
                    });
    if (inner_forced) {
      auto passes = detail::forced_split_passes(inner_key, legality, g);
      bool has_prod = false, has_cons = false;
      auto const tally = [&](std::size_t vid) {
        (passes.consumer_pass.count(vid) ? has_cons : has_prod) = true;
      };
      for (std::size_t v : buckets[n - 1].build_ids) tally(v);
      for (auto const& o : buckets[n - 1].outputs) tally(o.first);
      if (has_prod && has_cons) split_passes = std::move(passes);
    }
  }

  // Fail-safe for the non-innermost deferral above. The split fires ONLY at the
  // innermost realized axis (d == n-1); an OUTER realized axis (d < n-1) that
  // is genuinely forced to split would fall through to a single unbroken block
  // -- NOT a safe conservative choice but a SILENTLY WRONG cross-iteration
  // schedule (the loop is never actually split). No current input reaches here
  // (no multi-axis forced split exists), so this assert is inert today; it
  // turns that future silent mis-schedule into a loud failure the moment SP3 or
  // a later change introduces one. (n <= 1 skips the loop: the sole axis is by
  // definition the innermost.)
  for (std::size_t d = 0; d + 1 < n; ++d) {
    std::wstring const outer_key{types[d].space().base_key()};
    bool const outer_forced =
        std::any_of(legality.cells.begin(), legality.cells.end(),
                    [&](CellLegality const& cl) {
                      for (Index const& ix : forced_split_types(cl))
                        if (ix.space().base_key() == outer_key) return true;
                      return false;
                    });
    SEQUANT_ASSERT(!outer_forced,
                   "build_ordered_schedule: non-innermost forced split not yet "
                   "supported (SP2 deferral)");
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
      [&](Index const& axis, int ordinal,
          container::svector<std::size_t> const& build_ids,
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
    block.ordinal = ordinal;
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

  for (std::size_t i = 0; i < n; ++i) {
    std::size_t const d = n - 1 - i;  // innermost (n-1) to outermost (0)
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

    if (split_passes && d == n - 1) {
      // Innermost split: partition this bucket's homed BuildSteps and escape
      // outputs into the producer (ordinal 0) and consumer (ordinal 1) passes
      // by consumer_pass membership. d == n-1 => no pending child blocks to
      // fork. The consumer block reads the producer block's escaped (now-full)
      // outputs, so its `requires_` names them and the outer topo-sort orders
      // the producer pass first.
      auto const in_consumer = [&](std::size_t vid) {
        return split_passes->consumer_pass.count(vid) != 0;
      };
      auto const emit_pass =
          [&](int ordinal, container::svector<std::size_t> const& builds,
              container::svector<std::pair<std::size_t, OutputKind>> const&
                  outs) {
            container::svector<std::size_t> pass_produced = builds;
            for (auto const& o : outs) pass_produced.push_back(o.first);
            next_steps.push_back(
                Step{make_block(types[d], ordinal, builds, outs, {}, {})});
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

      emit_pass(0, prod_builds, prod_outs);
      emit_pass(1, cons_builds, cons_outs);
    } else {
      ScopeBlock block =
          make_block(types[d], 0, bucket.build_ids, bucket.outputs,
                     std::move(pending_steps), std::move(pending_metas));
      next_steps.push_back(Step{std::move(block)});
      detail::OrderedScheduleStepMeta m;
      for (auto const& out_entry : bucket.outputs)
        m.produced.push_back(out_entry.first);
      m.requires_ = requires_all;
      m.tie_key = min_first_use(produced_all);
      next_metas.push_back(std::move(m));
    }

    pending_steps = std::move(next_steps);
    pending_metas = std::move(next_metas);
    child_produced_all = std::move(produced_all);
    child_requires_all = std::move(requires_all);
  }

  // Root assembly: root-level BuildStep's plus the outermost block(s) as
  // Step(s) (one in the usual chain; two if the outermost axis itself was the
  // innermost -- i.e. n == 1 -- and got split).
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
  for (std::size_t k = 0; k < pending_steps.size(); ++k) {
    root_items.push_back(std::move(pending_steps[k]));
    root_meta.push_back(std::move(pending_metas[k]));
  }
  out.root.steps = detail::ordered_schedule_topo_sort_steps(
      std::move(root_items), root_meta);

  SEQUANT_ASSERT(well_formed(out));
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
