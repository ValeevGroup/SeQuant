#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>

#include <tiledarray.h>

#include <functional>
#include <optional>

namespace sequant {

/// \brief Backend context for the TiledArray eval backend.
///
/// Holds TiledArray-specific evaluation state that travels with the
/// CacheManager (via its \c backend_context field) to the binary-product
/// evaluation site in sequant::evaluate().  All TA-specific types are kept
/// here so the generic eval and CacheManager remain TA-free.
///
/// Usage (consumer side):
/// \code
///   TAEvalContext ta_ctx;
///   ta_ctx.result_shape_provider = [&](auto const& node, auto const& tr) {
///     // inspect node, build a TA::SparseShape<float>, or return nullopt
///     return std::nullopt;
///   };
///   cache.set_product_node_visitor(
///       TAEvalContext::make_visitor(ta_ctx));
/// \endcode
struct TAEvalContext {
  /// Called by the eval engine at each internal binary-Product node.
  ///
  /// \param node  The Product node being evaluated.  Full IR access is
  ///              available: indices, op type, etc.
  /// \param result_outer_trange  The TiledRange of the Product's outer
  ///              result (used in Task 2 to build the SparseShape).  Passed
  ///              as \c TA::TiledRange{} in Task 1.
  /// \return a \c TA::SparseShape<float> to impose on the result (Task 2),
  ///         or \c std::nullopt to leave the node unshaped (Task 1 always
  ///         returns nullopt).
  std::function<std::optional<TA::SparseShape<float>>(
      FullBinaryNode<EvalExprTA> const& node,
      TA::TiledRange const& result_outer_trange)>
      result_shape_provider;

  /// Build a type-erased visitor suitable for
  /// CacheManager::set_product_node_visitor() from a TAEvalContext.
  ///
  /// The visitor \c any_cast's the node from the std::any argument (which
  /// holds a std::reference_wrapper<FullBinaryNode<EvalExprTA> const>) and
  /// forwards it to \p ctx.result_shape_provider if set.
  ///
  /// \param ctx  The TAEvalContext whose result_shape_provider is invoked.
  ///             The caller must ensure \p ctx outlives the visitor.
  static std::function<void(std::any const&)> make_visitor(
      TAEvalContext const& ctx) {
    return [&ctx](std::any const& node_any) {
      if (!ctx.result_shape_provider) return;
      auto const& node =
          std::any_cast<
              std::reference_wrapper<FullBinaryNode<EvalExprTA> const>>(
              node_any)
              .get();
      ctx.result_shape_provider(node, TA::TiledRange{});
    };
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
