#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <string>

namespace sequant {

class Tensor;

///
/// \brief defines types of binary
///        operations that can be performed on two EvalExpr objects.
enum class EvalOp {

  ///
  /// \brief Represents the sum of two EvalExpr objects. Such as a tensor plus
  ///        a tensor and a constant plus a constant.
  Sum,

  ///
  /// \brief Represents the product of two EvalExpr objects. Such as tensor
  ///        times a tensor, a constant times a constant, or a tensor times a
  ///        constant
  ///        (in either order).
  Product
};

///
/// \brief The ResultType enum.
///
/// \details Represents the type of the result of @c EvalOp on two EvalExpr
///          objects (or, single EvalExpr object if `this->is_primary()`).
///
enum class ResultType { Tensor, Scalar };

///
/// \brief The EvalExpr is a building block of binary trees used to evaluate
/// expressions.
///
/// \details The EvalExpr class itself is not a proper node of the binary tree.
/// It is rather a data that a node should hold.
///
class EvalExpr {
 public:
  using index_vector = Index::index_vector;

  ///
  /// \brief Construct an EvalExpr object from a tensor.
  ///
  explicit EvalExpr(Tensor const& tnsr);

  ///
  /// \brief Construct an EvalExpr object from a Constant.
  ///
  explicit EvalExpr(Constant const& c);

  ///
  /// \brief Construct an EvalExpr object from a Variable.
  ///
  explicit EvalExpr(Variable const& v);

  ///
  /// @param op Evaluation operation resulting to this object.
  /// @param res Evaluation result type that will be produced.
  /// @param expr A SeQuant expression corresponding to @c res.
  /// @param ixs Canonical indices used for annotating the result's modes if @c
  ///            res is tensor type. Possibly empty for non-tensor @c res type.
  /// @param phase Phase that was part of the tensor network canonicalization.
  ///              Considered for reusing sub-expressions.
  /// @param hash A hash value that is equal for two EvalExpr objects that
  ///             produce the same evaluated result modulo the @c phase.
  ///
  EvalExpr(EvalOp op, ResultType res, ExprPtr const& expr, index_vector ixs,
           std::int8_t phase, size_t hash);

  ///
  /// \return Operation type of this expression, or null if this is a primary
  /// expression.
  ///
  [[nodiscard]] const std::optional<EvalOp>& op_type() const noexcept;

  ///
  /// \return The ResultType of the evaluation performed on this node.
  ///
  [[nodiscard]] ResultType result_type() const noexcept;

  ///
  /// \brief Compute the hash value of this EvalExpr object.
  ///
  /// \details The hash value is computed during construction of the object
  ///          by also looking at the hash values of the EvalExpr objects if
  ///          passed.
  ///
  /// \return The hash value of this EvalExpr object.
  ///
  [[nodiscard]] size_t hash_value() const noexcept;

  ///
  /// \return The ExprPtr object that this EvalExpr object holds.
  ///
  [[nodiscard]] ExprPtr expr() const noexcept;

  ///
  /// \return True if this EvalExpr object contains a SeQuant tensor with
  ///         proto-indices, false otherwise.
  ///
  [[nodiscard]] bool tot() const noexcept;

  ///
  /// \return Returns the result of calling to_latex() on the ExprPtr object
  ///         contained by this object.
  ///
  [[nodiscard]] std::wstring to_latex() const noexcept;

  ///
  /// \return True if the ExprPtr held by this object is Tensor and equivalently
  ///         the result of evaluation is tensor.
  ///
  [[nodiscard]] bool is_tensor() const noexcept;

  ///
  /// \return True if the ExprPtr held by this object is scalar(Constant, or
  ///         Variable) and equivalently the result of evaluation is scalar.
  ///
  [[nodiscard]] bool is_scalar() const noexcept;

  ///
  /// \return True if ExprPtr held by this object is Constant.
  ///
  [[nodiscard]] bool is_constant() const noexcept;

  ///
  /// \return True if ExprPtr held by this object is Variable.
  ///
  [[nodiscard]] bool is_variable() const noexcept;

  ///
  /// \return True if this is a primary expression (i.e. a leaf on expression
  /// tree)
  ///
  [[nodiscard]] bool is_primary() const noexcept;

  ///
  /// \return True if this expression is a product.
  ///
  [[nodiscard]] bool is_product() const noexcept;

  ///
  /// \return True if this expression is a sum.
  ///
  [[nodiscard]] bool is_sum() const noexcept;

  ///
  /// \brief Calls to<Tensor>() on ExprPtr held by this object.
  ///
  /// \return Tensor const&
  ///
  [[nodiscard]] Tensor const& as_tensor() const;

  ///
  /// \brief Calls to<Constant>() on ExprPtr held by this object.
  ///
  /// \return Constant const&.
  ///
  [[nodiscard]] Constant const& as_constant() const;

  ///
  /// \brief Calls to<Variable>() on ExprPtr held by this object.
  ///
  /// \return Variable const&
  ///
  [[nodiscard]] Variable const& as_variable() const;

  ///
  /// \brief Get the label for this object useful for logging.
  ///
  [[nodiscard]] std::string label() const noexcept;

  ///
  /// \return A string usable as TiledArray annotation if is_tensor() true,
  ///         empty string otherwise.
  ///
  [[nodiscard]] std::string indices_annot() const noexcept;

  ///
  /// \return Canonically ordered indices -- non-empty if this object represents
  /// a tensor result.
  ///
  [[nodiscard]] index_vector const& canon_indices() const noexcept;

  ///
  /// \return The canonicalization phase (+1 or -1).
  ///
  [[nodiscard]] std::int8_t canon_phase() const noexcept;

 private:
  std::optional<EvalOp> op_type_ = std::nullopt;

  ResultType result_type_;

  size_t hash_value_;

  index_vector canon_indices_;

  std::int8_t canon_phase_{1};

  ExprPtr expr_;
};

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
};

namespace meta {

namespace detail {
template <typename, typename = void>
constexpr bool is_eval_expr{};

template <typename T>
constexpr bool
    is_eval_expr<T, std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>{true};

template <typename, typename = void>
constexpr bool is_eval_node{};

template <typename T>
constexpr bool
    is_eval_node<FullBinaryNode<T>, std::enable_if_t<is_eval_expr<T>>>{true};

}  // namespace detail

///
/// \brief A type satisfies eval_expr if it is convertible to EvalExpr.
/// \see   EvalExpr
///
template <typename T>
concept eval_expr = detail::is_eval_expr<T>;

///
/// \brief A type satisfies eval_node if it is a FullBinaryNode of a type
///        that satisfies the eval_expr concept.
///
template <typename T>
concept eval_node = detail::is_eval_node<std::remove_cvref_t<T>>;

///
/// \brief Satisfied by a range type of eval_node objects.
///
template <typename Rng>
concept eval_node_range =
    std::ranges::range<Rng> && eval_node<std::ranges::range_value_t<Rng>>;

}  // namespace meta

namespace impl {
FullBinaryNode<EvalExpr> binarize(ExprPtr const&);
}

///
/// \brief A type alias for the types that satisfy the eval_node concept.
///
template <meta::eval_expr T>
using EvalNode = FullBinaryNode<T>;

///
/// Creates a binary tree of evaluation.
///
template <typename ExprT = EvalExpr,
          typename = std::enable_if_t<std::is_constructible_v<ExprT, EvalExpr>>>
FullBinaryNode<ExprT> binarize(ExprPtr const& expr) {
  if constexpr (std::is_same_v<ExprT, EvalExpr>) return impl::binarize(expr);
  return transform_node(impl::binarize(expr),
                        [](auto&& val) { return ExprT{val}; });
}

///
/// Converts an `EvalExpr` to `ExprPtr`.
///
ExprPtr to_expr(meta::eval_node auto const& node) {
  auto const op = node->op_type();
  auto const& evxpr = *node;

  if (node.leaf()) return evxpr.expr();

  if (op == EvalOp::Product) {
    auto prod = Product{};

    ExprPtr lexpr = to_expr(node.left());
    ExprPtr rexpr = to_expr(node.right());

    prod.append(1, lexpr, Product::Flatten::No);
    prod.append(1, rexpr, Product::Flatten::No);

    assert(!prod.empty());

    if (prod.size() == 1 && !prod.factor(0)->is<Tensor>()) {
      return ex<Product>(Product{prod.scalar(), prod.factor(0)->begin(),
                                 prod.factor(0)->end(), Product::Flatten::No});
    } else {
      return ex<Product>(std::move(prod));
    }

  } else {
    assert(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
