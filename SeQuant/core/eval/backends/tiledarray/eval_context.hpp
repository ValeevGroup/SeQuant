#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>

#include <tiledarray.h>

#include <any>
#include <functional>
#include <optional>

namespace sequant {

/// \brief Backend context for the TiledArray eval backend.
///
/// Holds TiledArray-specific evaluation state that travels with the
/// CacheManager (via its shaped-product hook) to the binary-product evaluation
/// site in sequant::evaluate().  All TA-specific types are kept here and inside
/// the hook closure so the generic eval and CacheManager remain TA-free.
///
/// Usage (consumer side):
/// \code
///   TAEvalContext ta_ctx;
///   ta_ctx.result_shape_provider = [&](auto const& node, auto const& tr)
///       -> std::optional<TA::SparseShape<float>> {
///     // inspect node, build a TA::SparseShape<float>, or return nullopt
///     return std::nullopt;
///   };
///   cache.set_shaped_product_hook(
///       TAEvalContext::make_hook<double, TA::SparsePolicy>(ta_ctx));
/// \endcode
struct TAEvalContext {
  /// Called by the eval engine at each internal binary-Product node.
  ///
  /// \param node  The Product node being evaluated.  Full IR access is
  ///              available: indices, op type, etc.
  /// \param result_outer_trange  The TiledRange of the Product's outer result,
  ///              over which the returned SparseShape is built.
  /// \return a \c TA::SparseShape<float> to impose on the result, or
  ///         \c std::nullopt to leave the node unshaped (eval then takes the
  ///         normal einsum path).
  std::function<std::optional<TA::SparseShape<float>>(
      FullBinaryNode<EvalExprTA> const& node,
      TA::TiledRange const& result_outer_trange)>
      result_shape_provider;

  /// Build a shaped-product hook for CacheManager::set_shaped_product_hook()
  /// from a TAEvalContext.
  ///
  /// The hook, when consulted at a binary-Product node:
  ///   1. any_casts the node to FullBinaryNode<EvalExprTA>;
  ///   2. computes the result's outer TiledRange from the operands and
  ///      annotations;
  ///   3. calls result_shape_provider(node, trange);
  ///   4. if it returns nullopt, returns a null ResultPtr (eval falls through
  ///      to the unshaped prod()); else emits the shaped product (via
  ///      apply_shaped_product) and returns it.
  ///
  /// The provider \c std::function is captured BY VALUE, so the hook does not
  /// dangle on \p ctx after this call returns.
  ///
  /// \tparam NumericT  The numeric type of the operand arrays (e.g. double).
  /// \tparam PolicyT   The TA policy of the operand arrays (SparsePolicy for a
  ///                   shape to be meaningful).
  /// \param ctx  The TAEvalContext whose result_shape_provider is captured.
  template <typename NumericT, typename PolicyT>
  static CacheManager<FullBinaryNode<EvalExprTA>>::shaped_product_hook_type
  make_hook(TAEvalContext const& ctx) {
    // Capture the provider BY VALUE so the hook owns its copy and does not
    // dangle on ctx.
    auto provider = ctx.result_shape_provider;
    return [provider = std::move(provider)](
               std::any const& node_any, Result const& left,
               Result const& right,
               std::array<std::any, 3> const& annot) -> ResultPtr {
      if (!provider) return nullptr;

      auto const& node =
          std::any_cast<
              std::reference_wrapper<FullBinaryNode<EvalExprTA> const>>(
              node_any)
              .get();

      // The result's outer TiledRange, over which the provider builds a shape.
      auto const trange = result_outer_trange_from_results<NumericT, PolicyT>(
          left, right, annot);

      auto shape = provider(node, trange);
      if (!shape) return nullptr;  // decline => unshaped prod()

      // de_nest: ToT * ToT -> flat (both operands nested, result is not). Read
      // from the IR node, matching the eval site's computation.
      bool const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();

      // shape must outlive the assignment inside apply_shaped_product (TA holds
      // it by pointer); it does (local here, passed by const&, used fully
      // within the call which fences before returning).
      return apply_shaped_product<NumericT, PolicyT>(left, right, annot, *shape,
                                                     de_nest);
    };
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
