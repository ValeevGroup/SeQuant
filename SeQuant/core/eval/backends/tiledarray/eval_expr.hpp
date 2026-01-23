#ifndef SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_EXPR_HPP
#define SEQUANT_EVAL_BACKENDS_TILEDARRAY_EVAL_EXPR_HPP

#ifdef SEQUANT_HAS_TILEDARRAY

#include <SeQuant/core/eval/eval_expr.hpp>

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
