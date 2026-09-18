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
/// on the Kramers flavour of the index occupying it (e.g. the outer pair slots
/// of a Kramers-union CSV projector, served from one array for every pair
/// flavour). Empty functions mean no erasure: identities are bit-identical to
/// a hook-less binarization.
struct KramersBlindness {
  /// true if slot \p slot (position in Tensor::const_slots()) of \p t is blind
  std::function<bool(Tensor const&, std::size_t slot)> blind_slot;
  /// the flavour-erased (spin-free) image of a flavoured index space; must
  /// return the space itself for a space that carries no flavour
  std::function<IndexSpace(IndexSpace const&)> erase_space;
  [[nodiscard]] bool active() const noexcept {
    return static_cast<bool>(blind_slot) && static_cast<bool>(erase_space);
  }
};

/// \brief The indices erasable over a set of tensors: every flavoured plain
///        index all of whose plain-slot occurrences among \p tensors are
///        blind slots.
///
/// \details Only plain-slot occurrences decide: an index that occupies any
/// non-blind slot is pinned everywhere; one with no plain occurrence at all is
/// never erased. Proto occurrences neither nominate nor pin (a composite's
/// own flavour is value-distinctive and is never erased; its proto list is
/// rewritten to follow whatever the plain slots decided).
template <std::ranges::input_range Rng>
  requires std::convertible_to<std::ranges::range_value_t<Rng>, ExprPtr>
container::set<Index> erasable_indices(Rng const& tensors,
                                       KramersBlindness const& kb) {
  container::set<Index> candidates, pinned;
  if (!kb.active()) return candidates;
  auto const flavoured = [&kb](Index const& ix) {
    return ix.space() != kb.erase_space(ix.space());
  };
  auto note = [&](Index const& ix, bool blind) {
    if (ix.has_proto_indices()) return;  // protos follow the plain slots
    if (flavoured(ix)) (blind ? candidates : pinned).emplace(ix);
  };
  // design guards (a violation is a caller bug, never a runtime condition):
  // a blind slot holds a plain pure-occupied index, and one tensor never
  // reports the same index blind in one slot and non-blind in another
  auto const isr = get_default_context().index_space_registry();
  for (ExprPtr const& e : tensors) {
    if (!e->is<Tensor>()) continue;
    auto const& t = e->as<Tensor>();
    container::set<Index> blind_here, plain_here;
    std::size_t slot = 0;
    for (auto const& ix : t.const_slots()) {
      bool const blind = kb.blind_slot(t, slot++);
      if (blind) {
        if (ix.has_proto_indices())
          throw std::invalid_argument(
              "KramersBlindness: a blind slot must hold a plain (proto-free) "
              "index");
        if (!isr || !isr->is_pure_occupied(ix.space()))
          throw std::invalid_argument(
              "KramersBlindness: a blind slot must be pure occupied");
        blind_here.emplace(ix);
      } else if (!ix.has_proto_indices()) {
        plain_here.emplace(ix);
      }
      note(ix, blind);
    }
    for (auto const& ix : blind_here)
      if (plain_here.contains(ix))
        throw std::invalid_argument(
            "KramersBlindness: an index is blind in one slot and not in "
            "another slot of the same tensor");
  }
  for (auto const& p : pinned) candidates.erase(p);
  return candidates;
}

/// \brief The erasure map of a network: every erasable plain index (see
///        erasable_indices) mapped to its flavour-erased placeholder.
///
/// \details Placeholders are numbered by FIRST OCCURRENCE among the plain
/// slots of \p tensors (in tensor order, slot order), starting past the
/// largest ordinal any index of the network carries, so that (a) two
/// networks that differ only in the flavours of erasable indices get the
/// same erased spelling -- the same spelling up to a renaming that keeps the
/// occurrence order, which is what the layout fingerprint keys on -- and (b)
/// a placeholder can never collide with a spin-free index already present
/// (e.g. a union-contracted dummy). A composite's proto list is rewritten
/// through the same map (Index::transform).
using ErasureMap = container::map<Index, Index>;

template <std::ranges::input_range Rng>
  requires std::convertible_to<std::ranges::range_value_t<Rng>, ExprPtr>
ErasureMap erasure_map(Rng const& tensors, KramersBlindness const& kb) {
  ErasureMap result;
  auto const erasable = erasable_indices(tensors, kb);
  if (erasable.empty()) return result;
  // largest ordinal in the network (plain slots and protos)
  std::size_t max_ord = 0;
  auto note_ord = [&max_ord](Index const& ix) {
    if (ix.ordinal()) max_ord = std::max<std::size_t>(max_ord, *ix.ordinal());
  };
  for (ExprPtr const& e : tensors) {
    if (!e->is<Tensor>()) continue;
    for (auto const& ix : e->as<Tensor>().const_slots()) {
      note_ord(ix);
      for (auto const& p : ix.proto_indices()) note_ord(p);
    }
  }
  std::size_t next = max_ord + 1;
  for (ExprPtr const& e : tensors) {
    if (!e->is<Tensor>()) continue;
    for (auto const& ix : e->as<Tensor>().const_slots()) {
      if (ix.has_proto_indices() || !erasable.contains(ix)) continue;
      if (result.contains(ix)) continue;
      result.emplace(ix, Index{kb.erase_space(ix.space()), next++});
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
