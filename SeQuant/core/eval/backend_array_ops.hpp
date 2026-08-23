#ifndef SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP
#define SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <functional>
#include <utility>

namespace sequant {

/// \brief Backend-provided realizations of the two operations external-axis
/// batching needs but that are backend-specific: constructing a zero
/// destination array and chunking an axis into batches.
///
/// \details The neutral eval layer names only INDICES (which carry their
/// spaces); the backend (the "user", e.g. mpqc) supplies these closures, so no
/// backend artifact -- a TiledArray tiling has no meaning for, say, an on-disk
/// backend -- ever leaks into the eval layer.
///
/// This replaces the old "carrier" model, in which the batched executor
/// borrowed an axis's tiling from whichever array in the DAG happened to carry
/// it (\c Result::pre_sized_zeros_over_mode / \c Result::mode_batches) and then
/// had to reconcile that array's Result TYPE and mode ordinal against the
/// scatter destination's. Tiling is a property of the space, not of any one
/// array, so it is sourced once, backend-side, from the index alone.
struct BackendArrayOps {
  /// Construct a sufficiently-initialized ZERO result shaped by \p descriptor
  /// -- a FULL (unsliced) index list, e.g. a node's \c canon_indices(). The
  /// backend maps each index's space to its own artifact and applies its own
  /// outer/inner split for proto-bearing (nested) indices, so flat-vs-nested
  /// is decided by the descriptor, not by any type reconciliation here.
  /// "Sufficiently initialized" is backend-defined (TA: a World + TiledRange,
  /// zero-filled; a nested result gets empty inner tiles, filled by the
  /// subsequent scatter writes).
  std::function<ResultPtr(container::vector<Index> const& descriptor)>
      make_zeros;

  /// Enumerate the half-open [lo,hi) element ranges chunking \p axis at
  /// ~\p target_batch_size. The backend owns the chunking rule (TA lands on
  /// tile boundaries). Per-space: two indices of one space chunk identically.
  std::function<container::svector<std::pair<std::size_t, std::size_t>>(
      Index const& axis, std::size_t target_batch_size)>
      axis_batches;

  /// True iff both closures are installed (a batched run requires them).
  explicit operator bool() const noexcept {
    return static_cast<bool>(make_zeros) && static_cast<bool>(axis_batches);
  }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP
