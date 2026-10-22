#ifndef SEQUANT_EVAL_LEGALITY_HPP
#define SEQUANT_EVAL_LEGALITY_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sequant::eval {

///
/// \brief Legality analysis output types of the ordered scheduler: per-value
/// classification of which batch-loop axes a value's computation depends on,
/// and whether it is legal to home the value at each such axis.
///

///
/// \brief The role a batch-loop axis plays at one value (\c analyze_legality
/// fills \c CellLegality::per_axis with these).
///
enum class LoopRole {
  LoopLocal,      //!< the axis is fully consumed within the value's own build
  Reduction,      //!< the axis is summed (contracted) at or below the value
  LoopCarried,    //!< the axis survives into the value's own result indices
  LoopInvariant,  //!< the axis encloses the value but the value does not
                  //!< depend on it (a pure replication context)
};

///
/// \brief The role of one occurrence of a value with respect to one of the
/// value's batched positions: \c LoopLocal when the loop instance that slices
/// the position encloses the occurrence (a lockstep, per-batch read), \c
/// LoopCarried when it does not (the reader needs the assembled value).
///
struct OccurrenceAxisRole {
  std::size_t point;  //!< the occurrence (\c OccurrenceRec::point)
  LoopRole role;      //!< \c LoopLocal or \c LoopCarried
};

///
/// \brief The classification of one (value, build-site axis) pair.
///
struct AxisClass {
  Index axis;
  /// The value's summary role: \c Reduction or \c LoopInvariant from the
  /// value's own node, else \c LoopLocal if any occurrence is lockstep with
  /// the axis's loop instance and \c LoopCarried if none is.
  LoopRole role;
  /// The role at each occurrence (empty for \c Reduction / \c
  /// LoopInvariant). A \c LoopLocal axis may still list \c LoopCarried
  /// occurrences: the value is built inside the loop for its lockstep
  /// readers and also escapes assembled for the others.
  container::svector<OccurrenceAxisRole> per_occurrence = {};
};

///
/// \brief The legality record of one value (one \c ValueCell::hash): its
/// build-site, per-axis role, and home-placement legality.
///
struct CellLegality {
  std::size_t hash = 0;  //!< == value_key_of(ValueCell): the value id (node
                         //!< id + home-sliced positions), the join key with
                         //!< the rich schedule

  /// The batch-loop axes this value depends on at its own node -- carried in
  /// its own result or contracted at its own node (not recursively over its
  /// operand subtrees; each operand is itself a value with its own,
  /// separately-analyzed \c CellLegality). Filled by \c build_site_of.
  container::svector<Index> build_site;

  /// One entry per \c build_site axis, classifying its role. Filled by
  /// \c analyze_legality.
  /// \note \c LoopRole::LoopInvariant is never stored here --
  /// every \c build_site axis is by construction carried or contracted, so
  /// \c classify_axis returns only LoopLocal/Reduction/LoopCarried at these
  /// axes. A LoopInvariant axis is the implicit case (batchable + enclosing
  /// but absent from \c build_site); recover it from that absence, not from a
  /// \c per_axis entry (see \c home_floor).
  container::svector<AxisClass> per_axis;

  /// The shallowest legal home mode-set for this value: the \c per_axis
  /// axes whose role is \c LoopLocal (the axes the value stays sliced on --
  /// the value is homed inside these). Every other \c build_site axis --
  /// \c Reduction, \c LoopCarried -- is lifted out (the value is homed
  /// above it), as is any axis absent from \c build_site altogether (the
  /// implicit \c LoopInvariant case: the value never depended on it, so it
  /// was trivially hoisted over it). Filled by \c
  /// analyze_legality.
  container::svector<Index> home_floor;

  /// Loops that must re-enter (the value cannot be homed above them without
  /// a re-materializing split). Entries are per-\c
  /// Index-instance, not per-basis-key -- an outer product like \c
  /// A{;i_3}*A{;i_4} lists both i_3 and i_4 (both LoopCarried on space i),
  /// yet they name a single occ loop to split. Consumers that want the
  /// distinct loops (axes) to split, not the raw per-instance list, should
  /// use \c forced_split_types below.
  container::svector<Index> forced_split_axes;
};

///
/// \brief A forest-wide legality analysis: one \c CellLegality per distinct
/// value. Assembled by \c analyze_legality; this header only
/// declares the container shape.
///
struct LegalitySchedule {
  container::svector<CellLegality> cells;
};

///
/// \brief Group \p cell's \c forced_split_axes by axis type (\c
/// IndexBasis::base_key()), collapsing multiple same-type \c Index instances
/// into the single loop (axis) they jointly force to split.
///
/// \details \c forced_split_axes is recorded per \c Index instance (see its
/// field doc): an outer product like \c A{;i_3}*A{;i_4} lists both \c i_3 and
/// \c i_4, both \c LoopCarried on the occ space, but they name only one occ
/// loop that must re-enter. This returns one representative \c Index per
/// distinct \c IndexBasis::base_key(), in \c forced_split_axes's discovery
/// order, so a consumer that needs "which loops must split" (not "which indices
/// are carried") gets a de-duplicated-by-type answer.
///
[[nodiscard]] inline container::svector<Index> forced_split_types(
    CellLegality const& cell) {
  container::svector<Index> out;
  for (Index const& ix : cell.forced_split_axes) {
    auto const same_type = [&](Index const& o) {
      return o.basis().base_key() == ix.basis().base_key();
    };
    if (std::find_if(out.begin(), out.end(), same_type) == out.end())
      out.push_back(ix);
  }
  return out;
}

///
/// \brief The at-node build-site of the value rooted at \p node: the
/// batch-loop axes (per \p policy.is_batchable_index()) that appear at this
/// node -- carried in \p node's own result indices, or contracted at \p node.
///
/// \details Deliberately not a subtree union: in the ordered-scope model
/// every node is itself a value, and a value is one contraction of already
/// -cached operand values (its children), so its build-site is only the axes
/// it carries or contracts at its own node -- the axes below it belong to
/// the operand values' own (separately-analyzed) build-sites, not to this
/// one. Concretely, the union of
///   - \p node's own \c canon_indices() (its result/free indices), and
///   - \p node's own \c contracted_indices(node) (see eval_expr.hpp;
///     empty for a leaf or a non-product node),
/// filtered to the batchable subset. The result is de-duplicated by \c Index
/// identity (basis + ordinal + proto-indices); order follows first
/// discovery (carried indices, then contracted indices).
///
[[nodiscard]] inline container::svector<Index> build_site_of(
    meta::eval_node auto const& node, BatchPolicy const& policy) {
  container::svector<Index> result;
  auto const add_if_new = [&](Index const& ix) {
    if (std::find(result.begin(), result.end(), ix) == result.end())
      result.push_back(ix);
  };

  // Source the build-site axes from the cost model's actual per-node decision,
  // not from policy.is_batchable_index() (what could be batched). Using the
  // predicate made the ordered schedule loop over every batchable axis (e.g.
  // aux whenever batch:aux_target_size>0) at every node carrying/contracting
  // it, regardless of whether the peak-constrained optimizer decided to slice
  // it -- always-on batching independent of peak_threshold, and over-scoping
  // relative to the DP even under a finite budget. The two authoritative
  // sources are:
  //   - a result (carried) index is a build-site axis iff it slices this node's
  //     own result slots per the cross-occurrence meet (\c sliced_modes, which
  //     already folds in ancestor batch loops the value is variant to); and
  //   - a contracted-at-node index iff the DP batched it here (a \c Contracted
  //     \c node_slice_mask stamp == the node's chosen aprime).
  // With an infinite budget (or no batchable axis) the DP emits no stamps, so
  // both are empty and the schedule is flat -- matching forest descent.
  auto const& sliced = sequant::home_scope(node);
  auto const& stamps = node->node_slice_mask();
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
/// \brief A batch-loop instance: (IndexBasis::base_key(), fusion loop_slot),
/// the key of \c RichSchedule::loop_kind. Loop identity is established once,
/// by \c compute_dag_boulevard's union-find, and stamped on every occurrence
/// (\c OccurrenceRec::loop_slot / \c reduced_slot); nothing here rederives it.
///
using LoopInstance = std::pair<std::wstring, int>;

///
/// \brief The loop instance an occurrence opens on mode \p ix with kind \p
/// kind, read off the opener's own stamps: a Contracted open is the \c
/// reduced_slot entry for \p ix, an External open the \c loop_slot of \p ix's
/// carried position. Nullopt when the opener carries no stamp for it.
///
[[nodiscard]] inline std::optional<LoopInstance> opened_instance(
    OccurrenceRec const& opener, Index const& ix, BatchModeType kind) {
  std::wstring const bk{ix.basis().base_key()};
  if (kind == BatchModeType::Contracted) {
    for (auto const& [rm, rs] : opener.reduced_slot)
      if (rm == ix) return LoopInstance{bk, rs};
    return std::nullopt;
  }
  for (std::size_t p = 0; p < opener.carried.size(); ++p)
    if (opener.carried[p] == ix)
      return p < opener.loop_slot.size() && opener.loop_slot[p] >= 0
                 ? std::optional<LoopInstance>{{bk, opener.loop_slot[p]}}
                 : std::nullopt;
  return std::nullopt;
}

///
/// \brief The loop instances enclosing \p occ, one per \c OccurrenceRec::ectx
/// entry, each resolved through the occurrence that opened it (\c
/// ectx_opener_point) and that opener's stamps (\c opened_instance). \p
/// point_occ is \c point_occurrences of the forest's rich schedule.
///
/// Throws \c Exception when an entry cannot be resolved: the record names no
/// opener for it, the opener is not in \p point_occ, opens no such mode, or
/// carries no stamp for it. Each of those is a defect in the boulevard's
/// record, and a guess here (a label-only comparison, say) would turn it into
/// a silently wrong schedule.
///
[[nodiscard]] inline container::svector<LoopInstance> enclosing_instances(
    OccurrenceRec const& occ,
    std::unordered_map<std::size_t, PointOccurrence> const& point_occ) {
  container::svector<LoopInstance> out;
  if (occ.ectx.size() != occ.ectx_opener_point.size())
    throw Exception(
        "enclosing_instances: the occurrence at point " +
        std::to_string(occ.point) + " records " +
        std::to_string(occ.ectx.size()) + " enclosing loop(s) but " +
        std::to_string(occ.ectx_opener_point.size()) + " opener(s)");
  for (std::size_t k = 0; k < occ.ectx.size(); ++k) {
    auto const pit = point_occ.find(occ.ectx_opener_point[k]);
    if (pit == point_occ.end())
      throw Exception(
          "enclosing_instances: the opener of an enclosing loop (point " +
          std::to_string(occ.ectx_opener_point[k]) + ") is not an occurrence");
    OccurrenceRec const& opener = *pit->second.occ;
    Index const& L = occ.ectx[k].first;
    std::optional<LoopInstance> inst;
    for (auto const& [oix, okind] : opener.opens)
      if (oix == L) {
        inst = opened_instance(opener, oix, okind);
        break;
      }
    if (!inst)
      throw Exception("enclosing_instances: the enclosing loop over " +
                      toUtf8(L.label()) + " at point " +
                      std::to_string(occ.point) +
                      " has no loop identity at its opener (point " +
                      std::to_string(occ.ectx_opener_point[k]) + ")");
    out.push_back(*inst);
  }
  return out;
}

///
/// \brief The \c classify_axis resolver of enclosing loop instances over every
/// occurrence of \p rich: \c enclosing_instances over \c point_occurrences of
/// the forest, memoized per occurrence (one value's axes all ask about the
/// same occurrences). \p rich must outlive the returned function.
///
[[nodiscard]] inline std::function<
    container::svector<LoopInstance>(OccurrenceRec const&)>
enclosing_instances_of(RichSchedule const& rich) {
  struct State {
    std::unordered_map<std::size_t, PointOccurrence> point_occ;
    std::unordered_map<std::size_t, container::svector<LoopInstance>> memo;
  };
  auto st = std::make_shared<State>();
  st->point_occ = point_occurrences(rich);
  return [st](OccurrenceRec const& occ) {
    auto const it = st->memo.find(occ.point);
    if (it != st->memo.end()) return it->second;
    auto inst = enclosing_instances(occ, st->point_occ);
    st->memo.emplace(occ.point, inst);
    return inst;
  };
}

///
/// \brief Classify one batch-loop axis \p axis against one value, given the
/// value's own result indices \p carried, the axes it contracts at its own
/// node \p contracted_below (see \c contracted_indices in eval_expr.hpp), and
/// every use-site \p occurrences of the value in the forest.
///
/// \details Two questions, the second answered per occurrence:
///   - Q1: is \p axis one of \p carried (by \c Index identity)?
///     - no  -> does \p contracted_below hold it (the value reduces the axis
///       at its own node)? -> \c Reduction; otherwise the axis merely
///       encloses the value without touching it -> \c LoopInvariant.
///     - yes -> Q2: for every occurrence, is the loop instance that slices
///       \p axis's carried position (the occurrence's own \c loop_slot at that
///       position) among the instances enclosing the occurrence (\p
///       enclosing_instances)? Enclosed -> the occurrence reads the value
///       inside its own loop, per batch (\c LoopLocal at that occurrence);
///       not enclosed -> the occurrence reads it where that loop is not
///       running: the value opens the loop itself, or the reader sits after
///       it, or in a sibling nest -- in every case the reader needs the
///       assembled value (\c LoopCarried at that occurrence).
///     The value's role is \c LoopLocal if any occurrence is enclosed and \c
///     LoopCarried otherwise (also with no occurrences at all). The
///     per-occurrence verdicts are appended to \p per_occurrence when given;
///     a \c LoopLocal value with \c LoopCarried occurrences is built inside
///     its loop for the former and additionally escapes assembled for the
///     latter (\c build_ordered_schedule).
///
/// Positions, not labels, align the occurrences: every occurrence's \c
/// carried list is positionally aligned with \p carried (each occurrence's
/// frame relabels the same slots), so \p axis's position in \p carried is
/// the one read at every occurrence. Loop identity is the recorded one
/// (\c OccurrenceRec::loop_slot and the openers' stamps). An occurrence that
/// is not aligned, or a batched position with no stamp, is a defect of the
/// record and throws.
///
[[nodiscard]] inline LoopRole classify_axis(
    container::svector<Index> const& carried,
    container::svector<Index> const& contracted_below, Index const& axis,
    container::svector<OccurrenceRec> const& occurrences,
    std::function<container::svector<LoopInstance>(OccurrenceRec const&)> const&
        enclosing_instances,
    container::svector<OccurrenceAxisRole>* per_occurrence = nullptr) {
  // The role is decided for this axis, by identity, not for its basis key: a
  // value may carry one index of a basis key (loop-local on that instance) and
  // contract another index of the same key in batches (a reduction over a
  // different instance) -- e.g. a residual-pair intermediate carrying the
  // external pair ij while reducing over a contracted pair kl, once the
  // occupied space is batchable in both roles.
  auto const carried_pos = std::find(carried.begin(), carried.end(), axis);
  if (carried_pos == carried.end()) {
    bool const reduces_axis =
        std::find(contracted_below.begin(), contracted_below.end(), axis) !=
        contracted_below.end();
    return reduces_axis ? LoopRole::Reduction : LoopRole::LoopInvariant;
  }
  std::size_t const axis_slot =
      static_cast<std::size_t>(carried_pos - carried.begin());
  std::wstring const bk{axis.basis().base_key()};

  bool any_lockstep = false;
  for (OccurrenceRec const& occ : occurrences) {
    if (axis_slot >= occ.carried.size() ||
        occ.carried[axis_slot].basis().base_key() != bk)
      throw Exception("classify_axis: the occurrence at point " +
                      std::to_string(occ.point) +
                      " is not positionally aligned with the value's carried "
                      "list at position " +
                      std::to_string(axis_slot));
    int const own_slot =
        axis_slot < occ.loop_slot.size() ? occ.loop_slot[axis_slot] : -1;
    if (own_slot < 0)
      throw Exception("classify_axis: batched position " +
                      std::to_string(axis_slot) + " of the value at point " +
                      std::to_string(occ.point) +
                      " has no loop identity (loop_slot unstamped)");
    LoopInstance const own{bk, own_slot};
    auto const encl = enclosing_instances(occ);
    bool const lockstep =
        std::find(encl.begin(), encl.end(), own) != encl.end();
    any_lockstep = any_lockstep || lockstep;
    if (per_occurrence)
      per_occurrence->push_back(
          {occ.point, lockstep ? LoopRole::LoopLocal : LoopRole::LoopCarried});
  }
  return any_lockstep ? LoopRole::LoopLocal : LoopRole::LoopCarried;
}

///
/// \brief Fold a \c RichSchedule into the first real \c LegalitySchedule: one
/// \c CellLegality per \c ValueCell, with its at-node \c build_site (\c
/// build_site_of), \c per_axis classification (\c classify_axis), \c home_floor
/// (the \c LoopLocal subset of \c per_axis), and \c forced_split_axes (the \c
/// LoopCarried subset of \c per_axis) filled in.
///
/// \details \c ValueCell does not itself retain the forest node its \c hash
/// resolves to (only \c carried, a copy of the node's \c canon_indices()), so
/// \p forest is re-walked once up front into a hash -> representative-node
/// map (the same pattern \c build_value_node_map, value_node_map.hpp, uses)
/// to recover each cell's own
/// \c contracted_indices(node). Under perfect CSE every occurrence of a value
/// shares one node shape, so the first-visited representative suffices.
///
/// \par Forced splits
/// A value classified \c LoopCarried on axis \c L survives the axis into its
/// own result indices, so its producing \c L-loop must close before any
/// cross-iteration consumer can read it (the consumer, at a different \c L
/// iteration, reads a value bound to a foreign \c L index). \c L is therefore
/// recorded in that value's \c forced_split_axes: the axis cannot be a single
/// unbroken loop around this value. This is directly derivable from the
/// per-axis roles alone.
///
/// \par No fixpoint: a value read across a forced split is not demoted here
/// A value \c LoopLocal on a forced-split axis \c L that is read by a
/// later-pass reader is handled entirely by the sequencer (\c
/// build_ordered_schedule's per-nest pass placement, materialization rule),
/// which materializes such a value across the pass boundary at schedule-build
/// time -- "used across the split" is a property of the ordered pass structure
/// the sequencer builds, not of the DAG this function sees, so classification
/// here needs to react to it. \c analyze_legality is therefore a single
/// deterministic round: one classification per cell, with no re-derivation
/// loop.
///
template <meta::eval_node_range R>
[[nodiscard]] inline LegalitySchedule analyze_legality(
    RichSchedule const& rich, R const& forest, BatchPolicy const& policy) {
  using Node = std::ranges::range_value_t<R>;

  auto const enclosing = enclosing_instances_of(rich);

  // One representative forest node per value (keyed by value id, not node
  // id): a value's roles are read off one of its own occurrences' nodes -- a
  // node of the same hash home-sliced elsewhere carries that frame's slice
  // annotations, not this value's.
  //
  // Held by pointer, not by value: \c Node's copy constructor deep-copies the
  // whole subtree, so one entry per node costs O(subtree) and a by-value map is
  // O(n^2) in forest size -- gigabytes on a thousand-summand spine, which is
  // OOM long before any stack limit matters. \p forest outlives
  // this call (it is the caller's), so the pointees stay valid; the entries are
  // read-only here.
  //
  // Taking addresses requires the range to yield references to the caller's
  // nodes; a range of prvalues (a transform view, say) would hand us dangling
  // pointers. Same guard, same reason, as cache_manager()'s pointer-keyed DAG
  // walk (cache_manager.hpp) and the value->node bridges
  // (value_node_map.hpp).
  // `R const`, not `R`: the walk below iterates `forest`, which is bound as
  // `R const&`, so the const-qualified range's reference type is the one that
  // has to be a reference. A range that yields references when mutable and
  // prvalues when const would slip past the unqualified form.
  static_assert(
      std::is_reference_v<std::ranges::range_reference_t<R const>>,
      "analyze_legality(): the forest range must yield references to nodes "
      "that outlive the call (the value->node map keys on their addresses)");
  std::unordered_map<std::size_t, Node const*> node_of;
  {
    // Iterative pre-order: the residual's Sum spine is as deep as the number
    // of terms, so a recursive descent would overflow the call stack. Right
    // child pushed before left, so the pop order is the recursion's.
    std::vector<Node const*> stack;
    for (auto const& tree : forest) {
      stack.push_back(&tree);
      while (!stack.empty()) {
        Node const& n = *stack.back();
        stack.pop_back();
        node_of.emplace(value_key_of(n), &n);
        if (!n.leaf()) {
          stack.push_back(&n.right());
          stack.push_back(&n.left());
        }
      }
    }
  }

  (void)policy;  // build-site now sourced from the DP decision, not the policy

  // The single classification round over every cell.
  auto const build_cells = [&]() -> LegalitySchedule {
    LegalitySchedule out;
    out.cells.reserve(rich.cells.size());
    for (ValueCell const& vc : rich.cells) {
      CellLegality cl;
      cl.hash = value_key_of(vc);

      // Every RichSchedule cell was produced by walking this same forest (see
      // compute_dag_boulevard), so its hash must resolve here; a miss would
      // silently leave contracted_below empty and misclassify the cell rather
      // than surface the underlying bug.
      auto const it = node_of.find(value_key_of(vc));
      SEQUANT_ASSERT(it != node_of.end());
      container::svector<Index> contracted_below;
      for (Index const& ix : contracted_indices(*it->second))
        contracted_below.push_back(ix);

      // Build-site axes come from the cost model's actual per-node decision,
      // not policy.is_batchable_index() (what could be batched). Sourcing from
      // the predicate made the ordered schedule loop over every batchable axis
      // (e.g. aux whenever batch:aux_target_size>0) regardless of whether the
      // peak-constrained optimizer sliced it -- always-on batching independent
      // of peak_threshold, and over-scoping relative to the DP under a finite
      // budget. A result (carried) axis is a build-site axis iff it slices this
      // node's own result slots per the cross-occurrence meet (sliced_modes,
      // which folds in the ancestor batch loops the value is variant to); a
      // contracted-at-node axis iff the DP batched it here (a Contracted
      // node_slice_mask stamp == the node's chosen aprime). With an infinite
      // budget (or no batchable axis) the DP emits no stamps, so the site is
      // empty and the schedule is flat -- identical to forest descent.
      container::svector<Index> site;
      auto const add_if_new = [&](Index const& ix) {
        if (std::find(site.begin(), site.end(), ix) == site.end())
          site.push_back(ix);
      };
      auto const& dp_sliced = sequant::home_scope(*it->second);
      auto const& dp_stamps = (*it->second)->node_slice_mask();
      // Positional: the representative node's home is labeled in its tree's
      // frame, vc.carried in the first occurrence's; positions are canonical
      // across occurrences ("Loop identity and value identity"), labels are
      // not.
      {
        auto const& rep_carried = (*it->second)->canon_indices();
        for (std::size_t p = 0; p < vc.carried.size() && p < rep_carried.size();
             ++p)
          if (std::find(dp_sliced.begin(), dp_sliced.end(), rep_carried[p]) !=
              dp_sliced.end())
            add_if_new(vc.carried[p]);
      }
      for (Index const& ix : contracted_below)
        if (std::any_of(dp_stamps.begin(), dp_stamps.end(), [&](auto const& p) {
              return p.second == BatchModeType::Contracted && p.first == ix;
            }))
          add_if_new(ix);

      for (Index const& axis : site) {
        AxisClass ac;
        ac.axis = axis;
        ac.role = classify_axis(vc.carried, contracted_below, axis,
                                vc.occurrences, enclosing, &ac.per_occurrence);
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

  return build_cells();
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_LEGALITY_HPP
