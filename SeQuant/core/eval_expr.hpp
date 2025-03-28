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

namespace detail {
std::size_t next_eval_expr_id();
}

///
/// \brief The EvalOp enum
///
/// \details The EvalOp enum is used to distinguish between the different binary
///          operations that can be performed on two EvalExpr objects.
enum class EvalOp {

  ///
  /// \brief Represents the identity evaluation. It is not a binary operation
  ///        per se. An atomic tensor or an atomic constant is always evaluated
  ///        this way.
  Id,

  ///
  /// \brief Represents the sum of two EvalExpr objects. Such as a tensor plus
  ///        a tensor and a constant plus a constant.
  Sum,

  ///
  /// \brief Represents the product of two EvalExpr objects. Such as tensor
  ///        times a tensor, a constant times a constant, or a tensor times a
  ///        constant
  ///        (in either order).
  Prod
};

///
/// \brief The ResultType enum.
///
/// \details Represents the type of the result of @c EvalOp on two EvalExpr
///          objects (or, single EvalExpr object if the EvalOp is Id).
///
enum class ResultType { Tensor, Scalar };

///
/// \brief The EvalExpr class represents the object that go into the nodes of
///        the binary tree that is used to evaluate the sequant expressions.
///
/// \details The EvalExpr class itself is not a proper node. It is rather a data
///          that a node should hold.
///
class EvalExpr {
 public:
  using index_vector = Index::index_vector;

  ///
  /// \brief Construct an EvalExpr object from a tensor. The EvalOp is Id.
  ///
  explicit EvalExpr(Tensor const& tnsr);

  ///
  /// \brief Construct an EvalExpr object from a Constant. The EvalOp is Id.
  ///
  explicit EvalExpr(Constant const& c);

  ///
  /// \brief Construct an EvalExpr object from a Variable. The EvalOp is Id.
  ///
  explicit EvalExpr(Variable const& v);

  ///
  /// @param op Evaluation operation resulting to this object.
  /// @param res Evaluation result type that will be produced.
  /// @param expr A sequant expression corresponding to @c res.
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
  /// \return The EvalOp resulting into this EvalExpr object.
  ///
  [[nodiscard]] EvalOp op_type() const noexcept;

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

  /// \ return The ID of this object
  [[nodiscard]] std::size_t id() const noexcept;

  ///
  /// \return The ExprPtr object that this EvalExpr object holds.
  ///
  [[nodiscard]] ExprPtr expr() const noexcept;

  void set_expr(ExprPtr expr);

  ///
  /// \return True if this EvalExpr object contains a sequant tensor with
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
  /// \brief Calls to<Tensor>() on ExprPtr held by this object.
  ///
  /// \return Tensor const&
  ///
  [[nodiscard]] Tensor const& as_tensor() const noexcept;

  ///
  /// \brief Calls to<Constant>() on ExprPtr held by this object.
  ///
  /// \return Constant const&.
  ///
  [[nodiscard]] Constant const& as_constant() const noexcept;

  ///
  /// \brief Calls to<Variable>() on ExprPtr held by this object.
  ///
  /// \return Variable const&
  ///
  [[nodiscard]] Variable const& as_variable() const noexcept;

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
  EvalOp op_type_;

  ResultType result_type_;

  size_t hash_value_;

  std::size_t id_ = detail::next_eval_expr_id();

  index_vector canon_indices_;

  std::int8_t canon_phase_{1};

  ExprPtr expr_;
};

namespace meta {
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

template <typename T>
constexpr bool is_eval_node<const FullBinaryNode<T>,
                            std::enable_if_t<is_eval_expr<T>>>{true};

template <typename, typename = void>
constexpr bool is_eval_node_range{};

template <typename Rng>
constexpr bool is_eval_node_range<
    Rng, std::enable_if_t<is_eval_node<meta::range_value_t<Rng>>>> = true;

}  // namespace meta

namespace impl {
FullBinaryNode<EvalExpr> binarize(ExprPtr const&);
}

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
template <typename NodeT,
          typename = std::enable_if_t<meta::is_eval_node<NodeT>>>
ExprPtr to_expr(NodeT const& node) {
  auto const op = node->op_type();
  auto const& evxpr = *node;

  if (node.leaf()) return evxpr.expr();

  if (op == EvalOp::Prod) {
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
