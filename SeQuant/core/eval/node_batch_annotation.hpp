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
/// evaluator that hoists loop-invariant intermediates: \c contracted_modes
/// (this node's enclosing-contracted residency) and \c order_aware (the gate)
/// feed per-level placement together with \c EvalExpr::sliced_modes; \c
/// effective_count is the node's effective use count -- how many times its
/// value is (re)referenced across the enclosing batch loops it does NOT
/// carry, the product of per-mode batch counts over its escaped-outer set.
///
/// The defaults (\c order_aware == false, \c effective_count == 1,
/// \c contracted_modes empty) are the order-blind / OFF-path values, so an
/// unannotated node rides along inert and leaves \c order_aware_recompute ==
/// false byte-identical.
struct NodeBatchAnnotation {
  /// Batchable indices sliced AT this node, each tagged with its
  /// \c BatchModeType (see \c EvalExpr::batched_here).
  container::svector<std::pair<Index, BatchModeType>> axes{};
  /// The enclosing CONTRACTED (aux) batch modes this node carries open on its
  /// result -- the contracted-residency signal per-level placement unions with
  /// the external lifetime mask (\c EvalExpr::sliced_modes) to decide hoist
  /// placement. Emitted only on the order-aware path; empty by default (OFF
  /// path, byte-identical) and empty for a node invariant to every enclosing
  /// contracted loop. See \c EvalExpr::contracted_modes.
  container::svector<Index> contracted_modes{};
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
