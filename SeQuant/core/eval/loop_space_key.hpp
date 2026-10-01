#ifndef SEQUANT_EVAL_LOOP_SPACE_KEY_HPP
#define SEQUANT_EVAL_LOOP_SPACE_KEY_HPP

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <functional>
#include <string>
#include <utility>

namespace sequant::eval {

/// \brief The key under which the batching layer identifies a loop's index
///        space.
///
/// \details Every structure that names a batch loop by the space of its
/// axis -- the boulevard's loop groups and fusion slots, the legality
/// classification, the ordered schedule's chain and clusters, the cell
/// table's scope paths, the dag-scope strings -- goes through this key. By
/// default it is \c IndexSpace::base_key(). Under a \em loop-space erasure
/// (\c set_loop_space_erasure) it is the key of the erased image, so loops
/// over spaces that erase to one space are ONE loop: the opt-in fusion of the
/// Kramers-flavoured occupied loops of a Kramers-union evaluation, whose ↑ and
/// ↓ batches run over the same doublet range and whose t-independent
/// intermediates are one value under Kramers-blind identity. Without the
/// fusion those intermediates are rebuilt once per batch of every
/// other-flavour loop they do not depend on (the flavour blocks' loops chain
/// into one nest).
///
/// The erasure is process-global and meant as a temporary opt-in: the
/// caller (MPQC) installs it before an evaluation and clears it after. It
/// must map a space onto a space of the same extent and batching.
inline std::function<IndexSpace(IndexSpace const&)>& loop_space_erasure() {
  static std::function<IndexSpace(IndexSpace const&)> f;
  return f;
}

/// installs (or, with an empty function, clears) the loop-space erasure
inline void set_loop_space_erasure(
    std::function<IndexSpace(IndexSpace const&)> f) {
  loop_space_erasure() = std::move(f);
}

[[nodiscard]] inline std::wstring loop_space_key(IndexSpace const& s) {
  auto const& f = loop_space_erasure();
  return std::wstring(f ? f(s).base_key() : s.base_key());
}

[[nodiscard]] inline std::wstring loop_space_key(Index const& ix) {
  return loop_space_key(ix.space());
}

/// true iff \p a and \p b name the same loop space (see loop_space_key)
[[nodiscard]] inline bool same_loop_space(IndexSpace const& a,
                                          IndexSpace const& b) {
  if (a == b) return true;
  auto const& f = loop_space_erasure();
  return f ? f(a).base_key() == f(b).base_key() : false;
}

[[nodiscard]] inline bool same_loop_space(Index const& a, Index const& b) {
  return same_loop_space(a.space(), b.space());
}

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_LOOP_SPACE_KEY_HPP
