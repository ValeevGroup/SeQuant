#ifndef SEQUANT_EVAL_BACKENDS_BTAS_EVAL_EXPR_HPP
#define SEQUANT_EVAL_BACKENDS_BTAS_EVAL_EXPR_HPP

#ifdef SEQUANT_HAS_BTAS

#include <SeQuant/core/eval/eval_expr.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/index.hpp>

#include <range/v3/view.hpp>

namespace sequant {

///
/// \brief This class extends the EvalExpr class by adding an annot() method so
///        that it can be used to evaluate using BTAS.
///
class EvalExprBTAS final : public EvalExpr {
 public:
  using annot_t = container::svector<long>;

  ///
  /// \param bk iterable of Index objects.
  /// \return vector of long-type hash values
  ///         of the labels of indices in \c bk
  ///
  template <typename Iterable>
  static auto index_hash(Iterable&& bk) {
    return ranges::views::transform(
        std::forward<Iterable>(bk), [](auto const& idx) {
          //
          // WARNING!
          // The BTAS uses long for scalar indexing by default.
          // Hence, here we explicitly cast the size_t values to long
          // Which is a potentially narrowing conversion leading to
          // integral overflow. Hence, the values in the returned
          // container are mixed negative and positive integers (long type)
          //
          return static_cast<long>(sequant::hash::value(Index{idx}.label()));
        });
  }

  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  EvalExprBTAS(Args&&... args) : EvalExpr{std::forward<Args>(args)...} {
    annot_ = index_hash(canon_indices()) | ranges::to<annot_t>;
  }

  ///
  /// \return Annotation (container::svector<long>) for BTAS::Tensor.
  ///
  [[nodiscard]] inline annot_t const& annot() const noexcept { return annot_; }

 private:
  annot_t annot_;
};

/// Type alias for BTAS evaluation nodes
using EvalNodeBTAS = EvalNode<EvalExprBTAS>;

static_assert(meta::eval_node<EvalNodeBTAS>);
static_assert(meta::can_evaluate<EvalNodeBTAS>);

}  // namespace sequant

#endif  // SEQUANT_HAS_BTAS

#endif  // SEQUANT_EVAL_BACKENDS_BTAS_EVAL_EXPR_HPP
