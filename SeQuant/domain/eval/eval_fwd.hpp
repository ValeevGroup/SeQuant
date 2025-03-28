//
// Created by Bimal Gaudel on 3/27/25.
//

#ifndef SEQUANT_EVAL_FWD_HPP
#define SEQUANT_EVAL_FWD_HPP

#include <SeQuant/core/eval_expr.hpp>

#include <any>

namespace sequant {

struct CacheManager;
struct EvalResult;

using ERPtr = std::shared_ptr<EvalResult>;

namespace meta {

template <typename T>
concept has_annot = requires(T t) {
  t.annot();
  requires !std::is_void_v<decltype(t.annot())>;
};

template <typename T>
concept can_evaluate = eval_node<T> && requires(T n) {
  { *n } -> has_annot;
};

template <typename Rng>
concept can_evaluate_range =
    std::ranges::range<Rng> && can_evaluate<std::ranges::range_value_t<Rng>>;

template <typename Node, typename F>
concept leaf_node_evaluator =
    can_evaluate<Node> && requires(F f, Node const& n) {
      { f(n) } -> std::same_as<ERPtr>;
    };

template <typename C>
concept can_manage_cache =
    std::is_lvalue_reference_v<C> &&
    std::is_same_v<std::remove_reference_t<C>, CacheManager>;

}  // namespace meta

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
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

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
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
  static auto index_hash(Iterable const& bk) {
    return ranges::views::transform(bk, [](auto const& idx) {
      //
      // WARNING!
      // The BTAS expects index types to be long by default.
      // There is no straight-forward way to turn the default.
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
};  // EvalExprBTAS

static_assert(meta::eval_node<EvalNode<EvalExpr>>);
static_assert(meta::eval_node<EvalNode<EvalExprTA>>);
static_assert(meta::eval_node<EvalNode<EvalExprBTAS>>);

static_assert(!meta::can_evaluate<EvalNode<EvalExpr>>);
static_assert(meta::can_evaluate<EvalNode<EvalExprTA>>);
static_assert(meta::can_evaluate<EvalNode<EvalExprBTAS>>);

///
/// \brief This class represents a triplet of annotations used in a tensor
///        contraction or summation.
///
/// \tparam T Annotation type eg. TA::DistArray takes std::string.
///
template <typename T>
struct Annot {
  explicit Annot(std::array<std::any, 3> const& a)
      : lannot(std::any_cast<T>(a[0])),
        rannot(std::any_cast<T>(a[1])),
        this_annot(std::any_cast<T>(a[2])) {}

  /// Annotation of the left operand.
  T const lannot;

  /// Annotation of the right operand.
  T const rannot;

  /// Annotation of the result.
  T const this_annot;
};

// namespace meta {

// namespace {
// template <typename, typename = void>
// constexpr bool is_annot{true};

// // template <typename T>
// // constexpr bool is_annot<Annot<T>>{true};

// }  // namespace
// template <typename T>
// concept annot = is_annot<T>;

// }  // namespace meta

}  // namespace sequant
#endif  // SEQUANT_EVAL_FWD_HPP
