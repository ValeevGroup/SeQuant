#ifndef SEQUANT_EVAL_SCOPE_SCHEDULE_HPP
#define SEQUANT_EVAL_SCOPE_SCHEDULE_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/index.hpp>

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <string>

namespace sequant::eval {

///
/// \brief Task 1 of the whole-scope batched DAG execution design (see
/// `doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md`):
/// the scope-tree schedule data structure and its builder. Purely a
/// placement PROJECTION of the existing DAG-scope placement (\c
/// compute_dag_boulevard's \c RichSchedule) -- no execution here; the Task-2+
/// executor walks this tree.
///

///
/// \brief One node of the (narrow, single-chain) scope tree: one batch loop,
/// entered at its parent scope. The root -> node path is the loop nest a
/// value homed at that node is built under.
///
struct ScopeNode {
  Index mode{};  //!< the batch loop mode realized at this node; default
                 //!< (sentinel) on the root, which sits OUTSIDE every loop.
  BatchModeType kind =
      BatchModeType::Contracted;  //!< Contracted (accumulate on block exit)
                                  //!< or External (scatter on block exit);
                                  //!< meaningless on the root (never
                                  //!< consulted there).
  container::svector<std::size_t>
      homed_values{};  //!< \c ValueCell::value_id's homed AT this node --
                       //!< built once per block here, not deeper.
  container::vector<ScopeNode>
      children{};  //!< deeper loops, nested under this one. \c
                   //!< container::vector (not \c svector): \c ScopeNode is
                   //!< self-referential through this member, and \c
                   //!< boost::container::small_vector requires \c T complete
                   //!< at the point of instantiation (unlike \c std::vector,
                   //!< which tolerates an incomplete element type here).
};

///
/// \brief The scope tree (see \c ScopeNode) for a whole forest's \c
/// RichSchedule, plus the total value count.
///
struct ScopeSchedule {
  ScopeNode root{};            //!< the top scope (outside every loop).
  std::size_t num_values = 0;  //!< total value count (== rich.cells.size()).
};

namespace detail {

///
/// \brief True iff \p mode ever survives, un-summed, into a forest-ROOT
/// value's own (proto-expanded) carried slots.
///
/// \details \c RichSchedule does not carry \c BatchModeType directly: the
/// cross-occurrence meet in \c stamp_lifetime_masks folds EVERY kind
/// (External or Contracted) into one \c Index set (see \c
/// stamp_residency_impl's doc comment -- "any BatchModeType"), so a mode's
/// kind is not a field anywhere on \c ValueCell. It is still recoverable from
/// \p rich alone: a batch mode realized in the forest is either Contracted
/// (summed away below every root -- an ACCUMULATE loop, so it never appears
/// on a root's own result) or External (a free/spectator index of some final
/// output -- a SCATTER loop, one output slice per block, so it DOES appear on
/// a root's own result). A forest-root occurrence is identified by \c
/// OccurrenceRec::point == \c consumer_point: \c compute_dag_boulevard's
/// post-order walk only ever overwrites a child's \c consumer_point (to its
/// parent's point, always later in the walk), so a node that is never
/// anyone's child -- a top-level tree of the forest -- keeps its
/// default-seeded \c consumer_point == its own \c point.
///
inline bool mode_is_external(RichSchedule const& rich, Index const& mode) {
  for (auto const& cell : rich.cells) {
    bool const is_root =
        std::any_of(cell.occurrences.begin(), cell.occurrences.end(),
                    [](OccurrenceRec const& occ) {
                      return occ.point == occ.consumer_point;
                    });
    if (!is_root) continue;
    container::svector<Index> slots;
    for (auto const& s : cell.carried)
      sequant::detail::proto_expand_into(slots, s);
    if (std::find(slots.begin(), slots.end(), mode) != slots.end()) return true;
  }
  return false;
}

}  // namespace detail

///
/// \brief Build the NARROW scope tree from an existing DAG-scope placement
/// (design section "Narrow (current) vs general (future) trees"): one \c
/// ScopeNode per batch-mode INDEX TYPE present across \p rich's cells (keyed
/// by \c IndexSpace::base_key(), not raw \c Index identity -- two physical
/// labels of the same space fold into one node), chained in a single
/// canonical order shared by the whole forest -- effectively one linear loop
/// nest every value's home is a prefix of. Each \c ValueCell is assigned to
/// the node at the depth equal to the number of distinct index TYPES in its
/// \c home_modes (root when empty): the deepest level whose accumulated
/// (root-to-node) type set equals the cell's \c home_modes type set.
///
/// No execution -- purely a placement projection ready for the Task-2+
/// executor to walk.
///
/// \p mode_order ranks the canonical chain order: a brace-init list of
/// base-key strings (e.g. `{L"Κ"}`), most-significant (outermost) first.
/// (Taken as \c std::initializer_list<Key> -- rather than a generic \c auto
/// range -- specifically so a literal `{L"Κ"}` argument deduces: a bare
/// generic template parameter is a non-deduced context for a braced-init-list
/// argument; \c std::initializer_list<Key> is the standard exception.) A
/// present type absent from \p mode_order sorts after every listed type, in
/// alphabetical order. Per the design note ("Order does not affect
/// correctness here"), \p mode_order only pins a DETERMINISTIC order; the
/// per-value node assignment is by SET equality against the accumulated type
/// set at each depth, not positional.
///
template <typename Key>
ScopeSchedule build_scope_schedule(RichSchedule const& rich,
                                   std::initializer_list<Key> mode_order) {
  // 1. One representative Index per distinct index TYPE present across every
  // cell's home_modes.
  container::svector<Index> present;
  for (auto const& cell : rich.cells)
    for (auto const& m : cell.home_modes) {
      auto const& bk = m.space().base_key();
      bool const seen = std::any_of(
          present.begin(), present.end(),
          [&](Index const& ix) { return ix.space().base_key() == bk; });
      if (!seen) present.push_back(m);
    }

  // 2. Canonical order: rank by position in mode_order (matched by
  // base_key()); ties / unlisted types sort alphabetically after every
  // listed type.
  auto const rank_of = [&](Index const& ix) -> std::size_t {
    std::size_t i = 0;
    for (auto const& key : mode_order) {
      if (ix.space().base_key() == key) return i;
      ++i;
    }
    return static_cast<std::size_t>(-1);
  };
  std::sort(present.begin(), present.end(),
            [&](Index const& a, Index const& b) {
              auto const ra = rank_of(a), rb = rank_of(b);
              if (ra != rb) return ra < rb;
              return a.space().base_key() < b.space().base_key();
            });

  // 3. Build the chain, one ScopeNode per present type, deepest last.
  ScopeSchedule out;
  out.num_values = rich.cells.size();
  container::svector<ScopeNode*> chain;
  ScopeNode* cur = &out.root;
  for (auto const& m : present) {
    ScopeNode child;
    child.mode = m;
    child.kind = detail::mode_is_external(rich, m) ? BatchModeType::External
                                                   : BatchModeType::Contracted;
    cur->children.push_back(std::move(child));
    cur = &cur->children.back();
    chain.push_back(cur);
  }

  // 4. Assign each value to the node whose accumulated (root-to-node) type
  // set equals its home_modes type set (root when empty).
  for (auto const& cell : rich.cells) {
    container::svector<std::wstring> want;
    for (auto const& m : cell.home_modes) {
      auto const& bk = m.space().base_key();
      if (std::find(want.begin(), want.end(), bk) == want.end())
        want.push_back(bk);
    }

    ScopeNode* target = &out.root;
    if (!want.empty()) {
      container::svector<std::wstring> enclosing;
      for (auto* node : chain) {
        enclosing.push_back(node->mode.space().base_key());
        if (enclosing.size() == want.size() &&
            std::is_permutation(enclosing.begin(), enclosing.end(),
                                want.begin()))
          target = node;
      }
    }
    target->homed_values.push_back(cell.value_id);
  }

  return out;
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_SCOPE_SCHEDULE_HPP
