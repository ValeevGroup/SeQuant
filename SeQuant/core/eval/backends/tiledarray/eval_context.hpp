#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/backends/tiledarray/eval_expr.hpp>
#include <SeQuant/core/eval/backends/tiledarray/result.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/logger.hpp>

#include <tiledarray.h>

#include <algorithm>
#include <any>
#include <cstddef>
#include <functional>
#include <optional>

namespace sequant {

namespace detail {

/// Average nested-scalar count per outer element of a ToT operand result
/// (total scalar count / dense outer-element count); 1.0 for a flat or scalar
/// operand. Used only on the gated predicted-footprint trace path, so the
/// collective size reduction it performs (size_in_bytes sums over ranks) is
/// paid only when tracing. For a ToT * T contraction the surviving inner
/// (e.g. the PNO index) comes from the ToT operand, so this factor lifts the
/// result's outer-element prediction to a full nested-size prediction.
template <typename NumericT, typename PolicyT, typename InnerTileT>
[[nodiscard]] double tot_avg_inner_factor(Result const& r) {
  using ToTArray = TA::DistArray<TA::Tensor<InnerTileT>, PolicyT>;
  using ToTResult = ResultTensorOfTensorTA<ToTArray>;
  if (!r.is<ToTResult>()) return 1.0;
  auto const& a = r.template get<ToTArray>();
  std::size_t outer = 1;
  for (std::size_t d = 0; d < a.trange().rank(); ++d)
    outer *= a.trange().dim(d).extent();
  if (outer == 0) return 1.0;
  double const total_scalars = static_cast<double>(r.size_in_bytes()) /
                               static_cast<double>(sizeof(NumericT));
  return total_scalars / static_cast<double>(outer);
}

/// Nonzero outer-element count over a result TiledRange, honoring an optional
/// SparseShape (sum of kept-tile element volumes); the full element volume when
/// the shape is absent (dense / unconstrained). O(total tiles).
[[nodiscard]] inline std::size_t pred_outer_elems(
    TA::TiledRange const& tr,
    std::optional<TA::SparseShape<float>> const& shape) {
  if (!shape) return tr.elements_range().volume();
  std::size_t n = 0;
  for (auto const& idx : tr.tiles_range())
    if (!shape->is_zero(idx)) n += tr.make_tile_range(idx).volume();
  return n;
}

}  // namespace detail

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
  /// \tparam InnerTileT  The inner tile type of a nested (ToT) operand;
  /// defaults
  ///                   to the plain \c TA::Tensor<NumericT>. The CSV/PNO path
  ///                   produces arena-pinned ToT operands
  ///                   (\c TA::ArenaTensor<NumericT>), a distinct Result kind;
  ///                   instantiate with that inner tile so the hook recognizes
  ///                   (and shapes) CSV intermediates instead of declining
  ///                   them.
  /// \param ctx  The TAEvalContext whose result_shape_provider is captured.
  template <typename NumericT, typename PolicyT,
            typename InnerTileT = TA::Tensor<NumericT>>
  static CacheManager<FullBinaryNode<EvalExprTA>>::shaped_product_hook_type
  make_hook(TAEvalContext const& ctx) {
    // Capture the provider BY VALUE so the hook owns its copy and does not
    // dangle on ctx.
    auto provider = ctx.result_shape_provider;
    return
        [provider = std::move(provider)](
            std::any const& node_any, Result const& left, Result const& right,
            std::array<std::any, 3> const& annot) -> ResultPtr {
          if (!provider) return nullptr;

          auto const& node =
              std::any_cast<
                  std::reference_wrapper<FullBinaryNode<EvalExprTA> const>>(
                  node_any)
                  .get();

          // The result's outer TiledRange computation (below) reads each
          // operand's array via a kind-dispatched get<>(); it handles only the
          // flat (ResultTensorTA) and nested (ResultTensorOfTensorTA) kinds.  A
          // product may legitimately have a SCALAR operand (ResultScalar),
          // which carries no TiledRange and would mis-cast.  Such a product has
          // no outer tensor result to shape, so decline early before computing
          // the trange.
          using FlatArray = TA::DistArray<TA::Tensor<NumericT>, PolicyT>;
          using ToTArray = TA::DistArray<TA::Tensor<InnerTileT>, PolicyT>;
          using FlatResult = ResultTensorTA<FlatArray>;
          using ToTResult = ResultTensorOfTensorTA<ToTArray>;
          auto is_tensor_like = [](Result const& r) {
            return r.is<FlatResult>() || r.is<ToTResult>();
          };
          if (!is_tensor_like(left) || !is_tensor_like(right)) return nullptr;

          // The result's outer TiledRange, over which the provider builds a
          // shape.
          auto const trange =
              result_outer_trange_from_results<NumericT, PolicyT, InnerTileT>(
                  left, right, annot);

          auto shape = provider(node, trange);
          if (!shape) return nullptr;  // decline => unshaped prod()

          // de_nest: ToT * ToT -> flat (both operands nested, result is not).
          // Read from the IR node, matching the eval site's computation.
          bool const de_nest =
              node.left()->tot() && node.right()->tot() && !node->tot();

          // shape must outlive the assignment inside apply_shaped_product (TA
          // holds it by pointer); it does (local here, passed by const&, used
          // fully within the call which fences before returning).
          return apply_shaped_product<NumericT, PolicyT, InnerTileT>(
              left, right, annot, *shape, de_nest);
        };
  }

  /// Build an observe-only predicted-footprint hook for
  /// CacheManager::set_predict_hook() from a TAEvalContext.
  ///
  /// The hook, consulted at a binary-Product node *before* materialization
  /// (only when tracing), computes the result's outer TiledRange from the
  /// operands, mirrors the shaping decision WITHOUT materializing (consulting
  /// the result_shape_provider only when \p shapeable -- i.e. when the same
  /// cache actually carries a shaped-product hook; the batched scratch does
  /// not, so its products are predicted dense), estimates the result footprint
  /// (outer-element count from the trange/shape times the ToT operand's average
  /// inner extent), and emits a "Predict" trace line.  It never alters the
  /// result.  Template parameters match make_hook().
  template <typename NumericT, typename PolicyT,
            typename InnerTileT = TA::Tensor<NumericT>>
  static CacheManager<FullBinaryNode<EvalExprTA>>::predict_hook_type
  make_predict_hook(TAEvalContext const& ctx) {
    auto provider = ctx.result_shape_provider;
    return [provider = std::move(provider)](
               std::any const& node_any, Result const& left,
               Result const& right, std::array<std::any, 3> const& annot,
               bool shapeable) -> void {
      if (Logger::instance().eval.level < 2) return;

      auto const& node =
          std::any_cast<
              std::reference_wrapper<FullBinaryNode<EvalExprTA> const>>(
              node_any)
              .get();

      using FlatArray = TA::DistArray<TA::Tensor<NumericT>, PolicyT>;
      using ToTArray = TA::DistArray<TA::Tensor<InnerTileT>, PolicyT>;
      using FlatResult = ResultTensorTA<FlatArray>;
      using ToTResult = ResultTensorOfTensorTA<ToTArray>;
      auto is_tensor_like = [](Result const& r) {
        return r.is<FlatResult>() || r.is<ToTResult>();
      };
      if (!is_tensor_like(left) || !is_tensor_like(right)) return;

      auto const trange =
          result_outer_trange_from_results<NumericT, PolicyT, InnerTileT>(
              left, right, annot);

      // Mirror the shaping decision without materializing: a shape is consulted
      // (and reported) only when this cache will actually try to shape it.
      std::optional<TA::SparseShape<float>> shape =
          (shapeable && provider) ? provider(node, trange) : std::nullopt;

      double const inner =
          node->tot()
              ? std::max(
                    detail::tot_avg_inner_factor<NumericT, PolicyT, InnerTileT>(
                        left),
                    detail::tot_avg_inner_factor<NumericT, PolicyT, InnerTileT>(
                        right))
              : 1.0;
      std::size_t const outer = detail::pred_outer_elems(trange, shape);
      auto const& tiles = trange.tiles_range();
      std::size_t total_tiles = tiles.volume(), nz_tiles = total_tiles;
      if (shape) {
        nz_tiles = 0;
        for (auto const& idx : tiles)
          if (!shape->is_zero(idx)) ++nz_tiles;
      }
      std::size_t const pred_bytes =
          static_cast<std::size_t>(static_cast<double>(outer) * inner *
                                   static_cast<double>(sizeof(NumericT)));

      char const* const shape_str = !shapeable ? "unshaped(batched)"
                                    : shape    ? "SHAPED"
                                               : "plain";
      write_log(Logger::instance(), "Predict", " | ", node->label(),
                " | shape=", shape_str, " | nnz_tiles=", nz_tiles, "/",
                total_tiles, " | pred_result=", pred_bytes, "B", '\n');
    };
  }
};

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_CONTEXT_HPP
