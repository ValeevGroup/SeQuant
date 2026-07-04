#ifndef SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_EVAL_EXPR_HPP
#define SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_EVAL_EXPR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/conversion.hpp>

#include <memory>
#include <type_traits>
#include <utility>

namespace sequant::eval::dryrun {

///
/// \brief Extends EvalExpr with an annot() method so DryRun eval nodes can be
///        evaluated.
///
/// Unlike \c EvalExprTAPP (opaque \c int64_t hashes of index labels -- see
/// \c backends/tapp/eval_expr.hpp), DryRun's annotation IS the plain literal
/// (canon-order) index list itself: \c Result::prod/sum/permute need each
/// index's actual space/extent (via \c CostModel), not just its identity, to
/// compute a modeled size.
///
class EvalExprDryRun final : public EvalExpr {
 public:
  using annot_t = dryrun::annot_t;  // container::svector<Index>

  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  explicit EvalExprDryRun(Args&&... args)
      : EvalExpr{std::forward<Args>(args)...} {
    annot_ = canon_indices() | ranges::to<annot_t>;
  }

  ///
  /// \return Annotation (container::svector<Index>) for DryRun tensors.
  ///
  [[nodiscard]] annot_t const& annot() const noexcept { return annot_; }

 private:
  annot_t annot_;
};

/// Type alias for DryRun evaluation nodes
using EvalNodeDryRun = EvalNode<EvalExprDryRun>;

static_assert(meta::eval_node<EvalNodeDryRun>);
static_assert(meta::can_evaluate<EvalNodeDryRun>);

///
/// \brief Leaf yielder: turns each IR leaf (a tensor/constant/variable node)
///        into a zero-data DryRun Result. This is the `F` in
///        \c evaluate<Trace>(node, layout, F, cache).
///
/// A tensor leaf's literal (canon-order) index list decides flat vs nested:
/// \c make_dryrun_result builds a flat \c ResultDryRun if none of the leaf's
/// indices are proto-indexed, or a nested \c ResultDryRunNested (a CSV/PNO
/// amplitude or coefficient) if any are -- and threads that SAME literal list
/// through as the nested result's canon-order position map, so a later
/// \c slice_mode()/\c mode_batches() call (which the batched runtime only
/// ever issues against a LEAF's result) resolves its positional `mode`
/// argument correctly regardless of the leaf's flat/nested-ness.
///
struct DryRunLeafEvaluator {
  std::shared_ptr<CostModel const> cm;

  [[nodiscard]] ResultPtr operator()(EvalNodeDryRun const& leaf) const {
    SEQUANT_ASSERT(leaf.leaf());
    if (!leaf->is_tensor()) {
      // Constant / Variable leaf: a bare scalar. No real numeric value is
      // ever tracked by this zero-data backend (only sizes/costs), so 1.0 is
      // a placeholder never meant to be read as a physical result.
      return eval_result<ResultScalar<double>>(1.0);
    }
    container::svector<Index> idx = leaf->canon_indices() | ranges::to<annot_t>;
    return make_dryrun_result(std::move(idx), cm);
  }
};

}  // namespace sequant::eval::dryrun

#endif  // SEQUANT_CORE_EVAL_BACKENDS_DRYRUN_EVAL_EXPR_HPP
