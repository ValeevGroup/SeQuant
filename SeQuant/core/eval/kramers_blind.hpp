//
// Kramers-blind eval-node identity (design spec:
// doc/dev/specs/2026-09-18-block-shared-t-independent-intermediates.md in
// MPQC). A caller declares which leaf slots are served independently of the
// Kramers flavour of the index occupying them; nodes whose only flavoured
// indices sit in such slots are then one value for every flavour and share
// one eval slot.
//

#ifndef SEQUANT_CORE_EVAL_KRAMERS_BLIND_HPP
#define SEQUANT_CORE_EVAL_KRAMERS_BLIND_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/index_space_registry.hpp>
#include <SeQuant/core/space.hpp>

#include <algorithm>
#include <functional>
#include <ranges>
#include <stdexcept>

namespace sequant::eval {

/// \brief Caller-declared Kramers blindness of leaf slots.
///
/// \details A slot is blind if the value served for the leaf does not depend
/// on the Kramers flavour of the index occupying it (a plain slot) or of the
/// proto indices of the composite occupying it (a composite slot): e.g. the
/// PNS composite slot of a Kramers-union CSV projector C{a~; a<ij>}, served
/// from one array for every flavour of the pair labels i, j. Empty functions
/// mean no erasure: identities are bit-identical to a hook-less
/// binarization.
struct KramersBlindness {
  /// true if slot \p slot (position in Tensor::const_slots()) of the LEAF
  /// \p t is blind
  std::function<bool(Tensor const&, std::size_t slot)> blind_slot;
  /// the flavour-erased image of a flavoured index space: a canonical
  /// representative such as the ↑ partner. It must return the space itself
  /// for a space that carries no flavour, and it must NOT map onto a space
  /// whose leaves are served differently (e.g. the spin-free union space, an
  /// index of which selects both Kramers halves: a genuine union-dummy leaf
  /// and a flavour-erased pair leaf would then share one identity but denote
  /// different arrays)
  std::function<IndexSpace(IndexSpace const&)> erase_space;
  [[nodiscard]] bool active() const noexcept {
    return static_cast<bool>(blind_slot) && static_cast<bool>(erase_space);
  }
};

/// \brief The indices erasable over a network of factors: every flavoured
///        plain index all of whose occurrences in LEAF factors are blind.
///
/// \details Occurrences in a leaf decide: a plain slot nominates (blind) or
/// pins (non-blind) its index; a composite slot nominates (blind) or pins
/// (non-blind) its flavoured proto indices, except a proto that also occupies
/// a plain slot of the SAME leaf, which is a reference to that slot and
/// follows it. A composite index itself is never erased (its own flavour is
/// value-distinctive). Non-leaf factors (Sum-rooted
/// intermediates entering a product) are neutral: whether their value depends
/// on a flavour is already encoded in their own identity hash, which the
/// product combines. An index that occurs in no leaf is never erased.
///
/// \param tensors the factors (Tensor expressions; others are skipped)
/// \param leaf_flags parallel to \p tensors: whether each is a leaf; empty
///        means every factor is a leaf
template <std::ranges::input_range Rng>
  requires std::convertible_to<std::ranges::range_value_t<Rng>, ExprPtr>
container::set<Index> erasable_indices(
    Rng const& tensors, KramersBlindness const& kb,
    container::svector<bool> const& leaf_flags = {}) {
  container::set<Index> candidates, pinned;
  if (!kb.active()) return candidates;
  auto const flavoured = [&kb](Index const& ix) {
    return ix.space() != kb.erase_space(ix.space());
  };
  // design guards (a violation is a caller bug, never a runtime condition):
  // a blind slot holds a pure-occupied plain index or a composite whose
  // protos are pure occupied, and one leaf never reports an index blind in
  // one slot and non-blind in another
  auto const isr = get_default_context().index_space_registry();
  auto const pure_occ = [&isr](Index const& ix) {
    return isr && isr->is_pure_occupied(ix.space());
  };
  std::size_t k = 0;
  for (ExprPtr const& e : tensors) {
    std::size_t const pos = k++;
    if (!e->is<Tensor>()) continue;
    if (!leaf_flags.empty() && !leaf_flags[pos]) continue;  // neutral
    auto const& t = e->as<Tensor>();
    container::set<Index> blind_here, pinned_here, plain_slots;
    for (auto const& ix : t.const_slots())
      if (!ix.has_proto_indices()) plain_slots.emplace(ix);
    std::size_t slot = 0;
    for (auto const& ix : t.const_slots()) {
      bool const blind = kb.blind_slot(t, slot++);
      auto note = [&](Index const& p) {
        if (!flavoured(p)) return;
        if (blind) {
          if (!pure_occ(p))
            throw std::invalid_argument(
                "KramersBlindness: a blind slot must hold (or, for a "
                "composite, be indexed by) pure-occupied indices");
          blind_here.emplace(p);
          candidates.emplace(p);
        } else {
          pinned_here.emplace(p);
          pinned.emplace(p);
        }
      };
      if (ix.has_proto_indices()) {
        for (auto const& p : ix.proto_indices())
          if (!plain_slots.contains(p)) note(p);  // else: follows its slot
      } else {
        note(ix);
      }
    }
    for (auto const& ix : blind_here)
      if (pinned_here.contains(ix))
        throw std::invalid_argument(
            "KramersBlindness: an index is blind in one slot and not in "
            "another slot of the same leaf");
  }
  for (auto const& p : pinned) candidates.erase(p);
  return candidates;
}

/// \brief The erasure map of a network: every erasable index (see
///        erasable_indices) mapped to its flavour-erased placeholder.
///
/// \details Placeholders are fresh temporary indices (Index::make_tmp_index)
/// of the erased space, minted in FIRST-OCCURRENCE order among the slots of
/// \p tensors (tensor order, slot order, a composite's protos in proto
/// order). Identity is label-blind, so only their relative order matters:
/// being the newest temporaries they sort after every index already in the
/// network and among themselves in occurrence order, hence two networks that
/// differ only in the flavours of erasable indices get the same erased
/// spelling up to a renaming that keeps the occurrence order (what the layout
/// fingerprint keys on), and a placeholder never collides with an index
/// already present. A composite's proto list is rewritten through the same
/// map (Index::transform).
using ErasureMap = container::map<Index, Index>;

template <std::ranges::input_range Rng>
  requires std::convertible_to<std::ranges::range_value_t<Rng>, ExprPtr>
ErasureMap erasure_map(Rng const& tensors, KramersBlindness const& kb,
                       container::svector<bool> const& leaf_flags = {}) {
  ErasureMap result;
  auto const erasable = erasable_indices(tensors, kb, leaf_flags);
  if (erasable.empty()) return result;
  auto place = [&](Index const& ix) {
    if (!erasable.contains(ix) || result.contains(ix)) return;
    result.emplace(ix, Index::make_tmp_index(kb.erase_space(ix.space())));
  };
  for (ExprPtr const& e : tensors) {
    if (!e->is<Tensor>()) continue;
    for (auto const& ix : e->as<Tensor>().const_slots()) {
      if (ix.has_proto_indices())
        for (auto const& p : ix.proto_indices()) place(p);
      else
        place(ix);
    }
  }
  return result;
}

/// \return a clone of \p t with the erasure map applied to every slot (a
///         composite's protos included); slot order and labels' occurrence
///         order are kept
inline Tensor erase_indices(Tensor const& t, ErasureMap const& m) {
  if (m.empty()) return t;
  Tensor result = t;
  result.transform_indices(m);
  result.reset_tags();
  return result;
}

/// \return \p ix with the erasure map applied (itself, or its protos)
inline Index erase_index(Index const& ix, ErasureMap const& m) {
  if (m.empty()) return ix;
  Index result = ix;
  result.transform(m);
  result.reset_tag();
  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_KRAMERS_BLIND_HPP
