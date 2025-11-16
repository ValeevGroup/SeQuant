#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <cstddef>
#include <memory>
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

struct EvalOpSetter;

///
/// \brief The EvalExpr is a building block of binary trees used to evaluate
/// expressions.
///
/// \details The EvalExpr class itself is not a proper node of the binary tree.
/// It is rather a data that a node should hold.
///
class EvalExpr {
 public:
  friend struct EvalOpSetter;
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
  /// @param connectivity The graph representing the connectivity. May be null
  ///                     to indicate that no graph is present/necessary.
  ///
  EvalExpr(EvalOp op, ResultType res, ExprPtr const& expr, index_vector ixs,
           std::int8_t phase, size_t hash,
           std::shared_ptr<bliss::Graph> connectivity);

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

  ///
  /// \return Whether this expression has a connectivity graph
  /// \see connectivity_graph
  ///
  [[nodiscard]] bool has_connectivity_graph() const noexcept;

  ///
  /// \return The graph representing the connectivity of two factors in a
  /// product \note If has_connectivity_graph returns false, this function must
  /// not be called \see has_connectivity_graph
  ///
  [[nodiscard]] const bliss::Graph& connectivity_graph() const noexcept;

  ///
  /// \return A copy of the graph representing the connectivity of two factors
  /// in a product
  ///
  [[nodiscard]] std::shared_ptr<bliss::Graph> copy_connectivity_graph()
      const noexcept;

 protected:
  std::optional<EvalOp> op_type_ = std::nullopt;

  ResultType result_type_;

  ExprPtr expr_;

  index_vector canon_indices_;

  std::int8_t canon_phase_{1};

  size_t hash_value_;

  std::shared_ptr<bliss::Graph> connectivity_;
};

struct EvalOpSetter {
  void set(EvalExpr& expr, EvalOp op) { expr.op_type_ = op; }
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
/// Creates a binary tree for evaluation.
///
template <typename ExprT = EvalExpr>
  requires std::is_constructible_v<ExprT, EvalExpr>
FullBinaryNode<ExprT> binarize(ExprPtr const& expr) {
  if constexpr (std::is_same_v<ExprT, EvalExpr>) return impl::binarize(expr);
  return transform_node(impl::binarize(expr),
                        [](auto&& val) { return ExprT{val}; });
}

///
/// Creates a binary tree for evaluation.
///
template <typename ExprT = EvalExpr>
  requires std::is_constructible_v<ExprT, EvalExpr>
FullBinaryNode<ExprT> binarize(ResultExpr const& res) {
  FullBinaryNode<ExprT> tree = binarize<ExprT>(res.expression());

  const bool is_scalar =
      res.bra().empty() && res.ket().empty() && res.aux().empty();

  if (tree.size() < 2) {
    // We want to have a result node with the result from the ResultExpr.
    // In order for that to work, we need a dedicated result node in the first
    // place. Hence, we adapt the represented expression for terminals to be
    // that terminal multiplied by 1.
    ExprT result = [&]() {
      if (is_scalar) {
        return *binarize<ExprT>(ex<Variable>(res.result_as_variable()));
      }
      return *binarize<ExprT>(ex<Tensor>(res.result_as_tensor()));
    }();
    EvalOpSetter{}.set(result, EvalOp::Product);

    tree = FullBinaryNode<ExprT>(std::move(result), std::move(tree),
                                 binarize<ExprT>(ex<Constant>(1)));
  }
  SEQUANT_ASSERT(tree.size() > 1);

  if (is_scalar) {
    if (res.has_label()) {
      tree->expr().template as<Variable>().set_label(res.label());
    }
  } else {
    Tensor& tensor = tree->expr().template as<Tensor>();

    tensor = res.result_as_tensor();

    // if (res.has_label()) {
    //   tensor.set_label(res.label());
    // }

    // SEQUANT_ASSERT(tensor.num_slots() ==
    //        res.bra().size() + res.ket().size() + res.aux().size());
    // tensor.set_bra(res.bra());
    // tensor.set_ket(res.ket());
    // tensor.set_aux(res.aux());
  }

  return tree;
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

    SEQUANT_ASSERT(!prod.empty());

    if (prod.size() == 1 && !prod.factor(0)->is<Tensor>()) {
      return ex<Product>(Product{prod.scalar(), prod.factor(0)->begin(),
                                 prod.factor(0)->end(), Product::Flatten::No});
    } else {
      return ex<Product>(std::move(prod));
    }

  } else {
    SEQUANT_ASSERT(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
