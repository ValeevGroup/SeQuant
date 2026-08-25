//
// Per-contraction-node batch annotation emitted by the order-aware batched
// cost model and threaded to the eval tree.
//

#ifndef SEQUANT_EVAL_NODE_BATCH_ANNOTATION_HPP
#define SEQUANT_EVAL_NODE_BATCH_ANNOTATION_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <utility>

namespace sequant {

/// \brief Per-contraction-node batch annotation emitted by the order-aware
///        batched cost model (\c PeakBatchedModel::reconstruct_batched_modes)
///        and threaded to the eval tree via
///        \c BinarizationOptions::node_batch_axes.
///
/// One instance per contraction node, in the optimizer's left-first post-order.
/// Beyond \c axes, this is the annotation bridge to a runtime batched
/// evaluator that hoists loop-invariant intermediates: \c order_aware (the
/// gate) feeds per-level placement together with the runtime residency
/// \c EvalExpr::sliced_modes (the all-batched-modes cross-occurrence meet,
/// which now carries the enclosing-contracted residency directly); \c
/// effective_count is the node's effective use count -- how many times its
/// value is (re)referenced across the enclosing batch loops it does NOT
/// carry, the product of per-mode batch counts over its escaped-outer set.
///
/// The defaults (\c order_aware == false, \c effective_count == 1) are the
/// order-blind / OFF-path values, so an unannotated node rides along inert and
/// leaves \c order_aware_recompute == false byte-identical.
struct NodeBatchAnnotation {
  /// Batchable indices sliced AT this node, each tagged with its
  /// \c BatchModeType (see \c EvalExpr::node_slice_mask).
  container::svector<std::pair<Index, BatchModeType>> axes{};
  /// Subset of \c axes for which THIS node is the loop-OPEN site (the outermost
  /// node at which the physical batch loop over the index is introduced), as
  /// opposed to a deeper node that merely carries the sliced mode. An external
  /// batch loop opens once, at the term root (an external mode is on the final
  /// result, so the root is its outermost carrier); a contracted batch loop
  /// opens at its (unique) contraction node. Distinguished from \c axes so a
  /// consumer that needs the loop NEST (\c peak_profile's ectx) counts each
  /// physical loop once instead of once-per-carrying-node. Default empty (OFF
  /// path and every non-open node). See \c EvalExpr::batch_loops_opened_here.
  container::svector<std::pair<Index, BatchModeType>> opened_here{};
  /// Effective use count; see the class doc.
  std::size_t effective_count = 1;
  /// True iff the order-aware cost model emitted this node (set for EVERY node
  /// on the order-aware path, including a whole-nest invariant whose residency
  /// union is empty). The per-level placement order-aware GATE: a positive
  /// signal an empty union cannot provide, distinguishing an OFF-path all-full
  /// node (do not hoist) from an order-aware whole-nest invariant (hoist to the
  /// run/term-scope root). Default false keeps the OFF path byte-identical. See
  /// \c EvalExpr::batch_order_aware.
  bool order_aware = false;
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_BATCH_ANNOTATION_HPP
