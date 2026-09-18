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
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <functional>
#include <ranges>

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
  for (ExprPtr const& e : tensors) {
    if (!e->is<Tensor>()) continue;
    auto const& t = e->as<Tensor>();
    std::size_t slot = 0;
    for (auto const& ix : t.const_slots()) note(ix, kb.blind_slot(t, slot++));
  }
  for (auto const& p : pinned) candidates.erase(p);
  return candidates;
}

/// \return the flavour-erased image of a plain erasable index: same ordinal,
///         erased space
inline Index erased_image(Index const& ix, KramersBlindness const& kb) {
  SEQUANT_ASSERT(!ix.has_proto_indices());
  return Index{kb.erase_space(ix.space()), ix.ordinal()};
}

/// \return a clone of \p t in which every index of \p erasable -- as a slot
///         or inside a composite's proto list -- is replaced by its erased
///         image; labels' ordinals, slot order and proto order are kept
inline Tensor erase_indices(Tensor const& t,
                            container::set<Index> const& erasable,
                            KramersBlindness const& kb) {
  if (erasable.empty()) return t;
  container::map<Index, Index> repl;
  for (auto const& ix : erasable) repl.emplace(ix, erased_image(ix, kb));
  Tensor result = t;
  // Index::transform rewrites a composite's protos through the map
  result.transform_indices(repl);
  result.reset_tags();
  return result;
}

/// \return \p ix with erase_indices' replacement applied (plain index or
///         composite whose protos are rewritten)
inline Index erase_index(Index const& ix, container::set<Index> const& erasable,
                         KramersBlindness const& kb) {
  if (erasable.empty()) return ix;
  container::map<Index, Index> repl;
  for (auto const& e : erasable) repl.emplace(e, erased_image(e, kb));
  Index result = ix;
  result.transform(repl);
  result.reset_tag();
  return result;
}

}  // namespace sequant::eval

#endif  // SEQUANT_CORE_EVAL_KRAMERS_BLIND_HPP
