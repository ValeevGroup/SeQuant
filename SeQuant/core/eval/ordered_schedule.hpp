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
#include <initializer_list>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
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
  Transient,          //!< the value's home IS this block; nothing carries out
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

}  // namespace detail

///
/// \brief Task 3 of the ordered-scope batched-eval design (SP2): the
/// deterministic sequencer -- lowers SP1's \c LegalitySchedule (per-value
/// \c home_floor / \c per_axis roles) plus the \c RichSchedule (per-value
/// \c first_use / \c last_use over the forest's single post-order static-
/// point timeline) into an \c OrderedSchedule. NON-SPLIT case only: every
/// batch axis TYPE realizes exactly ONE loop block, chained (not branched)
/// exactly as \c build_scope_schedule's single canonical chain -- splitting
/// one axis TYPE into several disjoint passes is Task 4.
///
/// \details Three-part algorithm, pure scheduling (no cost choice):
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
/// \par 3. Topological order within a block
/// Each block's own \c steps interleave its \c BuildStep's (one per value
/// homed there) with, if the chain continues, ONE nested child \c
/// ScopeBlock \c Step for the next-deeper axis. Ordered by a
/// \c (representative_point, category) key, category 0 (child block) before
/// category 1 (\c BuildStep), stable-sorted ascending:
///   - a \c BuildStep's key is \c (value's \c ValueCell::first_use, 1) --
///     the point the value is first produced.
///   - a child block's key is \c (representative_point, 0), where \c
///     representative_point is the MIN \c ValueCell::first_use over every
///     value produced ANYWHERE in that block's subtree (its own \c
///     BuildStep's and \c outputs entries, folded recursively through any
///     further-nested child) -- NOT the max. \c compute_dag_boulevard's
///     \c point counter is monotone WITHIN one tree (post-order: every
///     operand's point precedes its consumer's) but carries NO cross-tree
///     meaning -- the forest's later trees (terms) simply continue
///     counting, so a single loop block realized ONCE across the WHOLE
///     forest (the non-split case) aggregates content from EVERY term that
///     touches its axis, and an EARLY term's own consumer can have a much
///     SMALLER point than a LATER term's content also homed in this same
///     block. Using the MAX would place the block after only the latest
///     term's content and thus, wrongly, before some earlier term's true
///     consumer (this was tried and failed the water-20 test below: an
///     early-term consumer sorted BEFORE the block whose output it reads).
///     The MIN is provably safe in the "before every true consumer"
///     direction instead: for ANY value V that directly reads a value O
///     the block produces, \c V.first_use == O's own \c consumer_point
///     (V is O's structural parent) which is >= O's own \c point, which is
///     >= the block's \c representative_point by construction (O is one of
///     the values the MIN ranges over) -- so \c representative_point <=
///     V.first_use always, and the tie case (equality, exactly water-20's
///     shape: I(i,i;a,a)'s \c first_use IS I(μ̃,μ̃)'s own point) is broken
///     correctly by the category tie-break (0 < 1, block first). The
///     symmetric risk (scheduling the block too EARLY, before some OTHER
///     value it genuinely depends on as an external input) is not derived
///     as rigorously -- see the report for this JUDGMENT CALL's scope.
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

  // hash -> value_id (== the value's slot in rich.cells, per ValueCell's own
  // doc comment).
  std::unordered_map<std::size_t, std::size_t> value_id_of;
  value_id_of.reserve(rich.cells.size());
  for (ValueCell const& vc : rich.cells)
    value_id_of.emplace(vc.hash, vc.value_id);

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

  // 3. Assemble the chain bottom-up (innermost first), threading each
  // block's representative point up to its parent depth.
  std::optional<ScopeBlock> child_block;
  std::size_t child_rep_point = 0;

  for (std::size_t i = 0; i < n; ++i) {
    std::size_t const d = n - 1 - i;  // innermost (n-1) to outermost (0)
    detail::OrderedScheduleDepthBucket const& bucket = buckets[d];

    container::vector<std::pair<std::pair<std::size_t, int>, Step>> items;
    for (std::size_t v : bucket.build_ids)
      items.push_back({{rich.cells[v].first_use, 1}, Step{BuildStep{v}}});
    if (child_block)
      items.push_back({{child_rep_point, 0}, Step{std::move(*child_block)}});
    std::stable_sort(
        items.begin(), items.end(),
        [](auto const& a, auto const& b) { return a.first < b.first; });

    ScopeBlock block;
    block.axis = types[d];
    block.ordinal = 0;
    block.kind = detail::mode_is_external(rich, types[d])
                     ? BatchModeType::External
                     : BatchModeType::Contracted;
    for (auto& [key, step] : items) block.steps.push_back(std::move(step));
    block.outputs.assign(bucket.outputs.begin(), bucket.outputs.end());

    std::size_t rep = std::numeric_limits<std::size_t>::max();
    for (std::size_t v : bucket.build_ids)
      rep = std::min(rep, rich.cells[v].first_use);
    for (auto const& [vid, kind] : bucket.outputs)
      rep = std::min(rep, rich.cells[vid].first_use);
    if (child_block) rep = std::min(rep, child_rep_point);
    if (rep == std::numeric_limits<std::size_t>::max())
      rep = 0;  // empty block (no direct content, no child): harmless
                // fallback, see the function doc comment's part 3.

    child_block = std::move(block);
    child_rep_point = rep;
  }

  // Root assembly: root-level BuildStep's plus, if the chain is non-empty,
  // the outermost block as one Step.
  container::vector<std::pair<std::pair<std::size_t, int>, Step>> root_items;
  for (std::size_t v : root_build_ids)
    root_items.push_back({{rich.cells[v].first_use, 1}, Step{BuildStep{v}}});
  if (child_block)
    root_items.push_back({{child_rep_point, 0}, Step{std::move(*child_block)}});
  std::stable_sort(
      root_items.begin(), root_items.end(),
      [](auto const& a, auto const& b) { return a.first < b.first; });
  for (auto& [key, step] : root_items)
    out.root.steps.push_back(std::move(step));

  SEQUANT_ASSERT(well_formed(out));
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_ORDERED_SCHEDULE_HPP
