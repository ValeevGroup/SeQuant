#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_EXPR_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_EXPR_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/eval_expr.hpp>

#include <range/v3/range/conversion.hpp>

#include <string>

namespace sequant {

///
/// \brief This class extends the EvalExpr class by adding an annot() method so
///        that it can be used to evaluate using TiledArray.
///
class EvalExprTA final : public EvalExpr {
 public:
  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  EvalExprTA(Args&&... args) : EvalExpr{std::forward<Args>(args)...} {
    annot_ = indices_annot();
  }

  [[nodiscard]] inline auto const& annot() const noexcept { return annot_; }

  ///
  /// \param head The head of a result: the tensor that an expression defines.
  /// \pre The indices of \p head are those of the result, as for the head
  ///      of the ResultExpr that the evaluated node was binarized from;
  ///      neither this function nor sequant::evaluate checks it.
  /// \return The layout of a result with its modes in the order that the
  ///         indices of \p head are written in: bra, then ket, then aux;
  ///         those with protoindices follow the others, after a semicolon,
  ///         as the modes of the inner tensors.
  ///         Give it to sequant::evaluate to lay a result out as its
  ///         head.
  ///
  [[nodiscard]] static std::string layout(Tensor const& head) {
    return indices_annot(head.const_braketaux_indices() |
                         ranges::to<index_vector>);
  }

 private:
  std::string annot_;
};

/// Type alias for TiledArray evaluation nodes
using EvalNodeTA = EvalNode<EvalExprTA>;

static_assert(meta::eval_node<EvalNodeTA>);
static_assert(meta::can_evaluate<EvalNodeTA>);

}  // namespace sequant

#endif  // SEQUANT_HAS_TILEDARRAY

#endif  // SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_EXPR_HPP
