#ifndef SEQUANT_EVAL_BACKENDS_TAPP_EVAL_EXPR_HPP
#define SEQUANT_EVAL_BACKENDS_TAPP_EVAL_EXPR_HPP

#ifdef SEQUANT_HAS_TAPP

#include <SeQuant/core/eval/eval_expr.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>

#include <range/v3/view.hpp>

#include <cstdint>

namespace sequant {

///
/// \brief This class extends the EvalExpr class by adding an annot() method so
///        that it can be used to evaluate using TAPP.
///
class EvalExprTAPP final : public EvalExpr {
 public:
  using annot_t = container::svector<int64_t>;

  ///
  /// \param bk iterable of Index objects.
  /// \return vector of int64_t-type hash values
  ///         of the labels of indices in \c bk
  ///
  template <typename Iterable>
  static auto index_hash(Iterable&& bk) {
    return ranges::views::transform(
        std::forward<Iterable>(bk), [](auto const& idx) {
          return static_cast<int64_t>(sequant::hash::value(Index{idx}.label()));
        });
  }

  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  EvalExprTAPP(Args&&... args) : EvalExpr{std::forward<Args>(args)...} {
    annot_ = index_hash(canon_indices()) | ranges::to<annot_t>;
  }

  ///
  /// \return Annotation (container::svector<int64_t>) for TAPP tensors.
  ///
  [[nodiscard]] inline annot_t const& annot() const noexcept { return annot_; }

 private:
  annot_t annot_;
};

/// Type alias for TAPP evaluation nodes
using EvalNodeTAPP = EvalNode<EvalExprTAPP>;

static_assert(meta::eval_node<EvalNodeTAPP>);
static_assert(meta::can_evaluate<EvalNodeTAPP>);

}  // namespace sequant

#endif  // SEQUANT_HAS_TAPP

#endif  // SEQUANT_EVAL_BACKENDS_TAPP_EVAL_EXPR_HPP
