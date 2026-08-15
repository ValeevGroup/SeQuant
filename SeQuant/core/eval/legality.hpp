#ifndef SEQUANT_EVAL_LEGALITY_HPP
#define SEQUANT_EVAL_LEGALITY_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>

namespace sequant::eval {

///
/// \brief Legality analysis output types (SP1 of the ordered-scope batched-
/// eval design): per-value classification of which batch-loop axes a value's
/// computation depends on, and (later tasks) whether it is legal to home the
/// value at each such axis.
///

///
/// \brief The role a batch-loop axis plays at one value (Task 2 fills \c
/// per_axis; this enum is declared now so \c AxisClass / \c CellLegality have
/// a complete type in Task 1).
///
enum class LoopRole {
  LoopLocal,      //!< the axis is fully consumed within the value's own build
  Reduction,      //!< the axis is summed (contracted) at or below the value
  LoopCarried,    //!< the axis survives into the value's own result indices
  LoopInvariant,  //!< the axis encloses the value but the value does not
                  //!< depend on it (a pure replication context)
};

///
/// \brief The classification of one (value, build-site axis) pair.
///
struct AxisClass {
  Index axis;
  LoopRole role;
};

///
/// \brief The legality record of one value (one \c ValueCell::hash): its
/// build-site (Task 1), per-axis role (Task 2), and home-placement legality
/// (Task 3/4).
///
struct CellLegality {
  std::size_t hash = 0;  //!< == ValueCell::hash

  //!< The batch-loop axes this value depends on AT ITS OWN NODE -- carried in
  //!< its own result OR contracted at its own node (NOT recursively over its
  //!< operand subtrees; each operand is itself a value with its own,
  //!< separately-analyzed \c CellLegality). Filled by \c build_site_of.
  container::svector<Index> build_site;

  //!< One entry per \c build_site axis, classifying its role. Filled in
  //!< Task 2. NOTE for SP2: \c LoopRole::LoopInvariant is NEVER stored here --
  //!< every \c build_site axis is by construction carried or contracted, so
  //!< \c classify_axis returns only LoopLocal/Reduction/LoopCarried at these
  //!< axes. A LoopInvariant axis is the IMPLICIT case (batchable + enclosing
  //!< but absent from \c build_site); recover it from that absence, not from a
  //!< \c per_axis entry (see \c home_floor).
  container::svector<AxisClass> per_axis;

  //!< The shallowest legal home mode-set for this value: the \c per_axis
  //!< axes whose role is \c LoopLocal (the axes the value stays sliced on --
  //!< the value is homed INSIDE these). Every other \c build_site axis --
  //!< \c Reduction, \c LoopCarried -- is lifted OUT (the value is homed
  //!< above it), as is any axis absent from \c build_site altogether (the
  //!< IMPLICIT \c LoopInvariant case: the value never depended on it, so it
  //!< was trivially hoisted over it). Filled in Task 3 by \c
  //!< analyze_legality.
  container::svector<Index> home_floor;

  //!< Loops that must re-enter (the value cannot be homed above them without
  //!< a re-materializing split). Filled in Task 3/4. Entries are per-\c
  //!< Index-INSTANCE, not per-space-type -- an outer product like \c
  //!< A{;i_3}*A{;i_4} lists BOTH i_3 and i_4 (both LoopCarried on space i),
  //!< yet they name a SINGLE occ loop to split. Consumers that want the
  //!< distinct loops (axes) to split, not the raw per-instance list, should
  //!< use \c forced_split_types below.
  container::svector<Index> forced_split_axes;
};

///
/// \brief A forest-wide legality analysis: one \c CellLegality per distinct
/// value. Assembled in Task 2 (\c analyze_legality); this header only
/// declares the container shape.
///
struct LegalitySchedule {
  container::svector<CellLegality> cells;
};

///
/// \brief A sequencer-supplied source of monotone demotions for \c
/// analyze_legality's fixpoint (SP2, Task 4).
///
/// \details Given the CURRENT \c LegalitySchedule (one round's classification),
/// returns the set of \c (ValueCell::hash, axis-SPACE-\c base_key()) pairs that
/// must be demoted \c LoopLocal -> \c LoopCarried on that axis in the NEXT
/// round
/// -- e.g. a value that is \c LoopLocal on a forced-split axis \c L but is used
/// in BOTH resulting \c L-passes, which cannot keep a single per-\c L copy
/// alive across the split (see \c ordered_schedule.hpp's \c
/// forced_split_demotions, the concrete SP2 source). The demotion is a property
/// of SP2's ORDERED split pass structure, not of the DAG SP1 sees, which is why
/// SP1 supplies NONE (the default-constructed, empty source): with no source \c
/// analyze_legality is BYTE-IDENTICAL to its pre-SP2 behavior (the fixpoint
/// runs exactly one inert round). The source is consulted afresh each round, so
/// a value already demoted in a prior round is simply not re-reported (it is no
/// longer \c LoopLocal), which -- together with the caller de-duplicating
/// against the running demotion set -- keeps the fixpoint monotone and bounded
/// by \c Sum|per_axis|+1.
///
using DemotionSource =
    std::function<container::svector<std::pair<std::size_t, std::wstring>>(
        LegalitySchedule const&)>;

///
/// \brief Group \p cell's \c forced_split_axes by axis SPACE TYPE (\c
/// base_key()), collapsing multiple same-type \c Index INSTANCES into the
/// single loop (axis) they jointly force to split.
///
/// \details \c forced_split_axes is recorded per \c Index instance (see its
/// field doc): an outer product like \c A{;i_3}*A{;i_4} lists BOTH \c i_3 and
/// \c i_4, both \c LoopCarried on the occ space, but they name only ONE occ
/// loop that must re-enter. This returns one representative \c Index per
/// distinct \c base_key(), in \c forced_split_axes's discovery order, so a
/// consumer that needs "which loops must split" (not "which indices are
/// carried") gets a de-duplicated-by-type answer.
///
[[nodiscard]] inline container::svector<Index> forced_split_types(
    CellLegality const& cell) {
  container::svector<Index> out;
  for (Index const& ix : cell.forced_split_axes) {
    auto const same_type = [&](Index const& o) {
      return o.space().base_key() == ix.space().base_key();
    };
    if (std::find_if(out.begin(), out.end(), same_type) == out.end())
      out.push_back(ix);
  }
  return out;
}

///
/// \brief The AT-NODE build-site of the value rooted at \p node: the
/// batch-loop axes (per \p policy.is_batchable_index()) that appear AT this
/// node -- carried in \p node's own result indices, OR contracted AT \p node.
///
/// \details Deliberately NOT a subtree union: in the ordered-scope model
/// every node is itself a VALUE, and a value is one contraction of already
/// -cached OPERAND values (its children), so its build-site is only the axes
/// it carries or contracts at its own node -- the axes below it belong to
/// the operand VALUES' own (separately-analyzed) build-sites, not to this
/// one. Concretely, the union of
///   - \p node's own \c canon_indices() (its result/free indices), and
///   - \p node's own \c contracted_indices(node) (see eval.hpp; empty for a
///     leaf or a non-product node),
/// filtered to the batchable subset. The result is de-duplicated by \c Index
/// identity (space + ordinal + proto-indices); order follows first
/// discovery (carried indices, then contracted indices).
///
[[nodiscard]] inline container::svector<Index> build_site_of(
    meta::eval_node auto const& node, BatchPolicy const& policy) {
  container::svector<Index> result;
  auto const add_if_new = [&](Index const& ix) {
    if (std::find(result.begin(), result.end(), ix) == result.end())
      result.push_back(ix);
  };

  // Source the build-site axes from the cost model's ACTUAL per-node decision,
  // NOT from policy.is_batchable_index() (what COULD be batched). Using the
  // predicate made the ordered schedule loop over every batchable axis (e.g.
  // aux whenever batch:aux_target_size>0) at every node carrying/contracting
  // it, regardless of whether the peak-constrained optimizer decided to slice
  // it -- always-on batching independent of peak_threshold, and over-scoping
  // relative to the DP even under a finite budget. The two authoritative
  // sources are:
  //   - a RESULT (carried) index is a build-site axis iff it slices this node's
  //     own result slots per the cross-occurrence meet (\c sliced_modes, which
  //     already folds in ancestor batch loops the value is variant to); and
  //   - a CONTRACTED-at-node index iff the DP batched it here (a \c Contracted
  //     \c batched_here stamp == the node's chosen aprime).
  // With an infinite budget (or no batchable axis) the DP emits no stamps, so
  // both are empty and the schedule is flat -- matching forest descent.
  auto const& sliced = node->sliced_modes();
  auto const& stamps = node->batched_here();
  for (Index const& ix : node->canon_indices())
    if (std::find(sliced.begin(), sliced.end(), ix) != sliced.end())
      add_if_new(ix);
  for (Index const& ix : contracted_indices(node))
    if (std::any_of(stamps.begin(), stamps.end(), [&](auto const& p) {
          return p.second == BatchModeType::Contracted && p.first == ix;
        }))
      add_if_new(ix);
  (void)policy;
  return result;
}

///
/// \brief Classify one batch-loop axis \p axis against one value, given the
/// value's own result indices \p carried, the axes it contracts AT its own
/// node \p contracted_below (see \c contracted_indices in eval.hpp), and
/// every use-site \p occurrences of the value in the forest.
///
/// \details The four-way decision tree (SP1 design, Task 2):
///   - Q1: does \p carried hold an index of \p axis's \c IndexSpace (compared
///     by \c base_key(), i.e. TYPE not identity)?
///     - NO  -> Q2a: does \p contracted_below hold an index of that type
///       (the value reduces the axis at its own node)? -> \c Reduction;
///       otherwise the axis merely encloses the value without touching it
///       -> \c LoopInvariant.
///     - YES -> Q2b: for every occurrence that has an ENCLOSING loop of that
///       axis type in its \c OccurrenceRec::ectx, compare (via the
///       ordinal-and-proto-aware \c Index::operator==) that loop's own
///       \c Index against EVERY one of the occurrence's own carried indices
///       of the same type (an occurrence may carry more than one same-space
///       index, e.g. an outer product carrying both a lockstep and a free
///       slot of the same space). Bound to the SAME \c Index at EVERY such
///       slot, at every such occurrence (lockstep with the enclosing loop)
///       -> \c LoopLocal; bound to a DIFFERENT \c Index at ANY slot, at any
///       occurrence (a free / cross-iteration read) -> \c LoopCarried. If no
///       occurrence has an enclosing loop of that type at all, the axis is a
///       plain free result index with nothing to lock it to a loop iteration
///       -> \c LoopCarried.
///
[[nodiscard]] inline LoopRole classify_axis(
    container::svector<Index> const& carried,
    container::svector<Index> const& contracted_below, Index const& axis,
    container::svector<OccurrenceRec> const& occurrences) {
  auto const same_type = [&](Index const& ix) {
    return ix.space().base_key() == axis.space().base_key();
  };

  bool const carries_axis =
      std::any_of(carried.begin(), carried.end(), same_type);

  if (!carries_axis) {
    bool const reduces_axis = std::any_of(contracted_below.begin(),
                                          contracted_below.end(), same_type);
    return reduces_axis ? LoopRole::Reduction : LoopRole::LoopInvariant;
  }

  bool found_enclosing = false;
  for (OccurrenceRec const& occ : occurrences) {
    auto const ectx_it =
        std::find_if(occ.ectx.begin(), occ.ectx.end(),
                     [&](auto const& e) { return same_type(e.first); });
    if (ectx_it == occ.ectx.end())
      continue;  // no enclosing loop of this
                 // type at this occurrence
    found_enclosing = true;

    // Compare EVERY same-type carried slot at this occurrence against the
    // enclosing loop's own Index, not just the first: a value can carry TWO
    // indices of L's space (e.g. i_3 lockstep, i_4 free) under a realized L
    // loop, and the free slot must not be shadowed by an earlier lockstep
    // slot. Only an occurrence whose EVERY same-type carried slot equals the
    // enclosing loop's own Index is truly lockstep; a single mismatched slot
    // makes the whole occurrence a free / cross-iteration read.
    bool matched_any_same_type = false;
    for (Index const& c : occ.carried) {
      if (!same_type(c)) continue;
      matched_any_same_type = true;
      if (!(c == ectx_it->first))
        return LoopRole::LoopCarried;  // free / cross-iteration binding
    }
    if (!matched_any_same_type)
      return LoopRole::LoopCarried;  // no same-type carried slot at all
  }
  return found_enclosing ? LoopRole::LoopLocal : LoopRole::LoopCarried;
}

///
/// \brief Fold a \c RichSchedule into the first real \c LegalitySchedule: one
/// \c CellLegality per \c ValueCell, with its at-node \c build_site (\c
/// build_site_of), \c per_axis classification (\c classify_axis), \c home_floor
/// (the \c LoopLocal subset of \c per_axis), and \c forced_split_axes (the \c
/// LoopCarried subset of \c per_axis, Task 4) filled in.
///
/// \details \c ValueCell does not itself retain the forest node its \c hash
/// resolves to (only \c carried, a copy of the node's \c canon_indices()), so
/// \p forest is re-walked once up front into a hash -> representative-node
/// map (same pattern as \c placement_remat.hpp's router build and \c
/// scope_executor.hpp's \c build_value_node_map) to recover each cell's own
/// \c contracted_indices(node). Under perfect CSE every occurrence of a value
/// shares one node shape, so the first-visited representative suffices.
///
/// \par Forced splits (Task 4, SOUND SP1 core)
/// A value classified \c LoopCarried on axis \c L survives the axis into its
/// own result indices, so its producing \c L-loop must CLOSE before any
/// cross-iteration consumer can read it (the consumer, at a different \c L
/// iteration, reads a value bound to a foreign \c L index). \c L is therefore
/// recorded in that value's \c forced_split_axes: the axis cannot be a single
/// unbroken loop around this value. This is directly derivable from the
/// per-axis roles alone and is the deliverable this task ships as sound.
///
/// \par The monotone fixpoint (Task 4)
/// The analysis is wrapped in a bounded fixpoint. The design's intended body
/// is a DEMOTION: once \c L is forced to split, a value used in BOTH resulting
/// \c L-passes cannot keep a single per-\c L copy alive across the split, so it
/// is demoted \c LoopLocal -> \c LoopCarried on \c L, which lifts its home
/// floor; classification/floors are re-derived and the round repeats until no
/// floor rises. Because a home can only ever LIFT (a role only moves toward
/// \c LoopCarried, never back), the fixpoint is monotone and terminates; the
/// hard cap \c Sum_over_cells|per_axis|+1 (each productive round demotes at
/// least one per-axis entry) guards against a non-terminating bug via \c
/// SEQUANT_ASSERT.
///
/// \par SP2 wiring of the demotion (Task 4)
/// "Used in both \c L-passes" is a property of the PASS STRUCTURE -- the
/// ordered, split schedule that SP2 builds -- not of the DAG SP1 sees. It is
/// therefore supplied FROM OUTSIDE, via the optional \p demotion_source
/// callback: the SP2 sequencer passes \c forced_split_demotions (\c
/// ordered_schedule.hpp), which knows the producer/consumer pass partition and
/// reports every \c LoopLocal value read from BOTH passes of a forced split.
/// Each fixpoint round consults it and grows the monotone demotion set; a
/// productive round demotes at least one \c per_axis entry, so the tightened
/// cap
/// \c Sum|per_axis|+1 still bounds termination. When \p demotion_source is
/// EMPTY (the SP1 default), \c derive_demotions is inert (the loop runs exactly
/// one round) and this function is BYTE-IDENTICAL to its pre-SP2 behavior --
/// deriving a demotion from a DAG guess is thereby avoided, not approximated.
///
template <meta::eval_node_range R>
[[nodiscard]] inline LegalitySchedule analyze_legality(
    RichSchedule const& rich, R const& forest, BatchPolicy const& policy,
    DemotionSource demotion_source = {}) {
  using Node = std::ranges::range_value_t<R>;

  std::unordered_map<std::size_t, Node> node_of;
  auto visit = [&](auto&& self, Node const& n) -> void {
    node_of.emplace(n->hash_value(), n);
    if (!n.leaf()) {
      self(self, n.left());
      self(self, n.right());
    }
  };
  for (auto const& tree : forest) visit(visit, tree);

  (void)policy;  // build-site now sourced from the DP decision, not the policy

  // Monotone demotion set (the fixpoint state): hash -> the axis SPACE
  // base_keys forced from LoopLocal to LoopCarried by a prior round. Grows
  // only; empty in SP1 (see derive_demotions).
  std::unordered_map<std::size_t, container::svector<std::wstring>> demotions;
  auto const is_demoted = [&](std::size_t hash, Index const& axis) {
    auto const it = demotions.find(hash);
    if (it == demotions.end()) return false;
    std::wstring const key{axis.space().base_key()};
    return std::find(it->second.begin(), it->second.end(), key) !=
           it->second.end();
  };

  // One classification round over every cell, honoring the current demotion
  // set. Pure in `demotions` (and the fixed inputs), so re-running it is what
  // makes the fixpoint well-defined.
  auto const build_cells = [&]() -> LegalitySchedule {
    LegalitySchedule out;
    out.cells.reserve(rich.cells.size());
    for (ValueCell const& vc : rich.cells) {
      CellLegality cl;
      cl.hash = vc.hash;

      // Every RichSchedule cell was produced by walking THIS SAME forest (see
      // compute_dag_boulevard), so its hash must resolve here; a miss would
      // silently leave contracted_below empty and misclassify the cell rather
      // than surface the underlying bug.
      auto const it = node_of.find(vc.hash);
      SEQUANT_ASSERT(it != node_of.end());
      container::svector<Index> contracted_below;
      for (Index const& ix : contracted_indices(it->second))
        contracted_below.push_back(ix);

      // Build-site axes come from the cost model's ACTUAL per-node decision,
      // NOT policy.is_batchable_index() (what COULD be batched). Sourcing from
      // the predicate made the ordered schedule loop over every batchable axis
      // (e.g. aux whenever batch:aux_target_size>0) regardless of whether the
      // peak-constrained optimizer sliced it -- always-on batching independent
      // of peak_threshold, and over-scoping relative to the DP under a finite
      // budget. A RESULT (carried) axis is a build-site axis iff it slices this
      // node's own result slots per the cross-occurrence meet (sliced_modes,
      // which folds in the ancestor batch loops the value is variant to); a
      // CONTRACTED-at-node axis iff the DP batched it here (a Contracted
      // batched_here stamp == the node's chosen aprime). With an infinite
      // budget (or no batchable axis) the DP emits no stamps, so the site is
      // empty and the schedule is flat -- identical to forest descent.
      container::svector<Index> site;
      auto const add_if_new = [&](Index const& ix) {
        if (std::find(site.begin(), site.end(), ix) == site.end())
          site.push_back(ix);
      };
      auto const& dp_sliced = it->second->sliced_modes();
      auto const& dp_stamps = it->second->batched_here();
      for (Index const& ix : vc.carried)
        if (std::find(dp_sliced.begin(), dp_sliced.end(), ix) !=
            dp_sliced.end())
          add_if_new(ix);
      for (Index const& ix : contracted_below)
        if (std::any_of(dp_stamps.begin(), dp_stamps.end(), [&](auto const& p) {
              return p.second == BatchModeType::Contracted && p.first == ix;
            }))
          add_if_new(ix);

      for (Index const& axis : site) {
        AxisClass ac;
        ac.axis = axis;
        ac.role =
            classify_axis(vc.carried, contracted_below, axis, vc.occurrences);
        // Monotone demotion (SP2 hook): a prior fixpoint round can force an
        // axis LoopLocal -> LoopCarried; roles only ever move toward
        // LoopCarried, never back, which is what makes the fixpoint monotone.
        if (ac.role == LoopRole::LoopLocal && is_demoted(vc.hash, axis))
          ac.role = LoopRole::LoopCarried;
        cl.per_axis.push_back(std::move(ac));
      }

      for (AxisClass const& ac : cl.per_axis) {
        // home_floor: the LoopLocal subset of per_axis -- the axes the value
        // stays sliced on (homed inside). Every other build-site axis
        // (Reduction, LoopCarried) is lifted out, as is any axis outside
        // build_site altogether (the implicit LoopInvariant case).
        if (ac.role == LoopRole::LoopLocal) cl.home_floor.push_back(ac.axis);
        // forced_split_axes: the LoopCarried subset. The value survives the
        // axis into its own result, so its producing loop must close before
        // any cross-iteration consumer -- the loop cannot stay a single
        // unbroken pass around this value.
        if (ac.role == LoopRole::LoopCarried)
          cl.forced_split_axes.push_back(ac.axis);
      }

      cl.build_site = std::move(site);
      out.cells.push_back(std::move(cl));
    }
    return out;
  };

  // The monotone demotion step. Given the current schedule (its per_axis roles
  // and forced_split_axes), enlarge `demotions` with every (value, axis) the
  // sequencer-supplied `demotion_source` reports as used in BOTH L-passes of a
  // forced L-split (and which must therefore lift its floor). Returns true iff
  // it added anything not already recorded.
  //
  // With NO source (the SP1 default) this is DELIBERATELY INERT (returns false
  // without touching `demotions`): "used in both L-passes" is a property of
  // SP2's ordered, split pass structure, which is supplied from outside, not
  // guessed from the DAG here. The de-duplication against the running set below
  // is what keeps the fixpoint monotone (a pair already demoted is never
  // counted as growth) and guarantees convergence within the cap even if the
  // source re-reports a stale pair.
  auto const derive_demotions = [&](LegalitySchedule const& current) -> bool {
    if (!demotion_source) return false;  // SP1: inert, byte-identical
    bool grew = false;
    for (auto const& [hash, key] : demotion_source(current)) {
      auto& keys = demotions[hash];
      if (std::find(keys.begin(), keys.end(), key) == keys.end()) {
        keys.push_back(key);
        grew = true;
      }
    }
    return grew;
  };

  LegalitySchedule out = build_cells();

  // Monotone fixpoint. `derive_demotions` demotes a single (cell, axis) pair
  // per round -- SP1's is inert (returns false immediately, so the loop below
  // runs exactly one round), but once SP2 wires it live, a cell can be lifted
  // on several axes over several rounds, so the tight non-termination bound
  // is Sum_over_cells |per_axis| + 1 (every productive round demotes at
  // least one per_axis entry, and roles only ever move toward LoopCarried,
  // never back), NOT cells.size()+1 (which only holds for CELL-granular
  // progress, i.e. >= one whole cell demoted per round). Asserted below as a
  // hard non-termination tripwire.
  [[maybe_unused]] std::size_t const cap = [&] {
    std::size_t sum = 0;
    for (CellLegality const& cl : out.cells) sum += cl.per_axis.size();
    return sum + 1;
  }();
  for (std::size_t iter = 0;; ++iter) {
    SEQUANT_ASSERT(iter <= cap,
                   "analyze_legality: forced-split fixpoint did not converge");
    if (!derive_demotions(out)) break;
    out = build_cells();
  }
  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_LEGALITY_HPP
