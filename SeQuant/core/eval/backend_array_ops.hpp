#ifndef SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP
#define SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>

namespace sequant {

/// \brief Backend-provided realizations of the two operations external-axis
/// batching needs but that are backend-specific: constructing a zero
/// destination array and chunking an axis into batches.
///
/// \details The neutral eval layer names only indices (which carry their
/// spaces); the backend (the "user", e.g. mpqc) supplies these closures, so no
/// backend artifact -- a TiledArray tiling has no meaning for, say, an on-disk
/// backend -- ever leaks into the eval layer.
///
/// Tiling is a property of the space, not of any one array, so it is sourced
/// once, backend-side, from the index alone rather than borrowed from whichever
/// array in the DAG carries that axis.
struct BackendArrayOps {
  /// Construct a sufficiently-initialized zero result shaped by \p descriptor
  /// -- a full (unsliced) index list, e.g. a node's \c canon_indices(). The
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

  /// The realization of an \c EvalOp::KramersFlip node: \p phase times the
  /// time-reversal flip of \p operand over its OUTER modes \p modes, each a
  /// Kramers-union axis blocked [⇑ half | ⇓ half]: per mode out[⇑] = +conj
  /// in[⇓], out[⇓] = −conj in[⇑] (the sign is that of the target half; the
  /// signs multiply over modes, the conjugation is applied once; no union
  /// mode = the conjugation alone). What a union axis IS (its blocking,
  /// tiling and half-swap) is the backend's ("user's") knowledge, so the
  /// operation lives here and not in the eval layer's array backends.
  /// Required by any run whose trees hold a KramersFlip node
  /// (\c BinarizationOptions::kramers_fold_intermediates).
  std::function<ResultPtr(Result const& operand,
                          container::svector<std::size_t> const& modes,
                          std::int8_t phase)>
      kramers_flip;

  /// True iff the two batching closures are installed (a batched run requires
  /// them); \c kramers_flip is independent of them.
  explicit operator bool() const noexcept {
    return static_cast<bool>(make_zeros) && static_cast<bool>(axis_batches);
  }
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_BACKEND_ARRAY_OPS_HPP
