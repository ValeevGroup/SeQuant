#ifndef SEQUANT_OPTIMIZE_OPTIMIZE_HPP
#define SEQUANT_OPTIMIZE_OPTIMIZE_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/options.hpp>

#include <cstddef>
#include <optional>

namespace sequant {

/// A forest-level batching signal: the single spectator external index (e.g. a
/// CC residual's own external OCCUPIED index, carried only as a composite
/// protoindex on the output PNOs and therefore contracted at NO node) over
/// which the whole optimized term should be evaluated in blocks of
/// \c block_size, so every intermediate carrying that axis shrinks by
/// ~block_size/extent while total work is unchanged (partitioning a spectator
/// is work-neutral). Because a never-contracted axis has no eval node to hang
/// on, this is reported per optimized term rather than as a per-node
/// \c EvalExpr::batch_axes entry; the runtime (P3) and MPQC (P4) read it here.
struct ExternalBatchAxis {
  /// The spectator external axis to block over.
  Index axis;
  /// The per-block extent (== \c batch_target_size(axis), snapped down to the
  /// axis extent).
  std::size_t block_size = 0;
};

/// Result of \ref optimize_result: the optimized expression plus an optional
/// \ref ExternalBatchAxis (empty when no spectator-index batching was
/// selected).
struct OptimizeResult {
  ExprPtr expr;
  std::optional<ExternalBatchAxis> external_batch_axis;
};

/// Optimize the expression for lower evaluation cost.
///
/// \param expr  Expression to be optimized.
/// \param opts  Optimization parameters; see \c OptimizeOptions. By default:
///              the cost metric is flop count, index extents are taken from
///              \c IndexSpace::approximate_size(), and the summands of a sum
///              are reordered to cluster terms that share intermediates.
/// \return Optimized expression.
ExprPtr optimize(ExprPtr const& expr, OptimizeOptions opts = {});

/// Forest-batching entry point: optimize \p expr, and additionally report an
/// optional \ref ExternalBatchAxis --- the spectator external axis over which
/// the whole term should be evaluated in blocks of
/// \c opts.batch_policy.batch_target_size(axis) (snapped down to the axis
/// extent).
///
/// The axis is emitted only when ALL of the following hold:
///   (i)   \p opts.batch_policy.batch_spectator_indices == true (false =
///         disabled; then \c external_batch_axis is empty and the optimized
///         expression is byte-identical to \ref optimize);
///   (ii)  \p opts.objective_function == \c DenseTimeSpaceBatched (the
///         perf-first batched objective);
///   (iii) \p expr is a single tensor Product with a genuine spectator external
///         axis (open on the root, contracted at no node); and
///   (iv)  the perf-first factorization's UNSEEDED peak exceeds
///         \p opts.batch_policy.peak_threshold --- i.e. the term actually needs
///         batching to fit the budget. (perf-first ignores \c peak_threshold
///         for FACTORIZATION selection; it is consulted here only to decide
///         WHEN to emit the axis.)
/// Otherwise \c external_batch_axis is empty.
///
/// The chosen factorization (\c OptimizeResult::expr) is byte-identical to what
/// \ref optimize produces under the same \p opts; only the extra signal is
/// computed.
///
/// \param expr Expression to be optimized.
/// \param opts Optimization parameters; see \ref OptimizeOptions. The spectator
///             per-block extent comes from
///             \c opts.batch_policy.batch_target_size(axis).
/// \return The optimized expression paired with the optional forest signal.
OptimizeResult optimize_result(ExprPtr const& expr, OptimizeOptions opts = {});

/// \copydoc optimize(ExprPtr const&, OptimizeOptions)
ResultExpr& optimize(ResultExpr& expr, OptimizeOptions opts = {});

/// \copydoc optimize(ExprPtr const&, OptimizeOptions)
[[nodiscard]] ResultExpr& optimize(ResultExpr&& expr,
                                   OptimizeOptions opts = {});

// Overloads for backwards compatibility

/// \deprecated Use the \c OptimizeOptions overload instead.
///
/// Equivalent to calling the primary overload with default \c OptimizeOptions
/// and \c OptimizeOptions::reorder_sum set to \p reorder_sum.
///
/// \param expr         Expression to be optimized.
/// \param reorder_sum  If true, reorder the summands of a sum to cluster terms
///                     that share intermediates.
/// \return Optimized expression.
[[deprecated(
    "use the OptimizeOptions"
    " overload of optimize() instead")]] ExprPtr
optimize(ExprPtr const& expr, bool reorder_sum);

/// \copydoc optimize(ExprPtr const&, bool)
[[deprecated(
    "use the OptimizeOptions"
    " overload of optimize() instead")]] ResultExpr&
optimize(ResultExpr& expr, bool reorder_sum);

/// \copydoc optimize(ExprPtr const&, bool)
[[nodiscard, deprecated("use the OptimizeOptions"
                        " overload of optimize() instead")]] ResultExpr&
optimize(ResultExpr&& expr, bool reorder_sum);

}  // namespace sequant

#endif  // SEQUANT_OPTIMIZE_OPTIMIZE_HPP
