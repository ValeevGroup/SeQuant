#ifndef SEQUANT_EVAL_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EVAL_EXPR_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/slot_symmetry.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <range/v3/algorithm/all_of.hpp>

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
  Product,

  ///
  /// \brief Represents the adjoint (conjugate transpose) of one EvalExpr
  ///        object. The result equals the bra/ket-swapped, complex-conjugated
  ///        operand. The IR representation is "structurally binary, lexically
  ///        unary": an Adjoint node holds the bare-label operand as its left
  ///        child and a Constant(1) sentinel as its right child (so the
  ///        FullBinaryNode invariant — every non-leaf has both children —
  ///        is preserved). Evaluate dispatches on this op_type and only uses
  ///        the left operand; the right is ignored.
  Adjoint
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
  /// \brief Construct an EvalExpr object from a Power.
  ///
  explicit EvalExpr(Power const& p);

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
  /// \return The type id of the Expr held by this object.
  ///
  [[nodiscard]] Expr::type_id_type type_id() const noexcept;

  ///
  /// \return True if the ExprPtr held by this object is Tensor and equivalently
  ///         the result of evaluation is tensor.
  ///
  [[nodiscard]] bool is_tensor() const noexcept;

  ///
  /// \return True if the ExprPtr held by this object is scalar (Constant,
  ///         Variable, or Power) and equivalently the result of evaluation
  ///         is scalar.
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
  /// \return True if ExprPtr held by this object is Power.
  ///
  [[nodiscard]] bool is_power() const noexcept;

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
  /// \return True if this expression is an adjoint (unary) node.
  ///
  [[nodiscard]] bool is_adjoint() const noexcept;

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
  /// \brief Calls to<Power>() on ExprPtr held by this object.
  ///
  /// \return Power const&
  ///
  [[nodiscard]] Power const& as_power() const;

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

  ///
  /// \return The out-of-band permutational-symmetry descriptor for this node's
  ///         result tensor.  Default-constructed (empty) until a later
  ///         deduction pass writes it.  This field is intentionally NOT
  ///         referenced by any hashing, canonicalization, bliss-graph, or
  ///         export code path.
  ///
  [[nodiscard]] SlotSymmetry const& slot_symmetry() const noexcept;

 protected:
  std::optional<EvalOp> op_type_ = std::nullopt;

  ResultType result_type_;

  ExprPtr expr_;

  index_vector canon_indices_;

  std::int8_t canon_phase_{1};

  size_t hash_value_;

  std::shared_ptr<bliss::Graph> connectivity_;

  /// Out-of-band permutational-symmetry descriptor.
  /// NOT referenced by hashing, canonicalization, bliss-graph, or export.
  SlotSymmetry slot_symmetry_{};
};

struct EvalOpSetter {
  void set(EvalExpr& expr, EvalOp op) { expr.op_type_ = op; }
  void reset(EvalExpr& expr) { expr.op_type_ = std::nullopt; }
  /// Set the out-of-band slot-symmetry descriptor (used by the binarize
  /// deduction pass; the descriptor is not part of the node's identity).
  void set_slot_symmetry(EvalExpr& expr, SlotSymmetry sym) {
    expr.slot_symmetry_ = std::move(sym);
  }
};

struct BinarizationOptions {
  /// Whether to merge indices of intermediate tensors into a single list
  /// (stored as aux indices) instead of retaining the bra, ket and aux
  /// separation
  bool merge_indices = false;
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

///
/// \brief Satisfied by a type with a method named `annot` that returns
///        a non-void type.
///
template <typename T>
concept has_annot = requires(T t) {
  t.annot();
  requires !std::is_void_v<decltype(t.annot())>;
};

///
/// \brief Satisfied by an eval_node whose dereferenced type satisfies the
///        has_annot method.
/// \example
///          * `static_assert(!meta::can_evaluate<EvalNode<EvalExpr>>)`
///          * `static_assert(meta::can_evaluate<EvalNodeTA>)` (where EvalNodeTA
///            is defined in backends/tiledarray/eval_expr.hpp)
///
template <typename T>
concept can_evaluate = eval_node<T> && requires(T n) {
  { *n } -> has_annot;
};

///
/// \brief Satisfied by a range type of objects satisfying can_evaluate.
///
template <typename Rng>
concept can_evaluate_range =
    std::ranges::range<Rng> && can_evaluate<std::ranges::range_value_t<Rng>>;

///
/// \brief \tparam F is a leaf node evaluator of type \tparam Node if
///        an object (a function object) of type \tparam F returns ResultPtr
///        when called with the single argument of const ref type to
///        \tparam Node and the \tparam Node satisfies can_evaluate.
///
template <typename Node, typename F>
concept leaf_node_evaluator =
    can_evaluate<Node> && requires(F f, Node const& n) {
      { f(n) } -> std::same_as<ResultPtr>;
    };

}  // namespace meta

namespace impl {

FullBinaryNode<EvalExpr> binarize(ExprPtr const&, IndexSet const& uncontract,
                                  const BinarizationOptions& opts);
}  // namespace impl

///
/// \brief A type alias for the types that satisfy the eval_node concept.
///
template <meta::eval_expr T>
using EvalNode = FullBinaryNode<T>;

static_assert(meta::eval_node<EvalNode<EvalExpr>>);
static_assert(!meta::can_evaluate<EvalNode<EvalExpr>>);

/// Creates a binary tree for evaluation.
/// @param expr an expression to binarize
/// @param external additional external (uncontracted) indices; only needed for
///        indices that appear more than once in @p expr (i.e., hyperindices)
///        but should not be contracted. Indices appearing once are
///        automatically treated as external.
///
/// @deprecated The root EvalExpr's tensor has a positional bra/ket split: each
///        surviving external ends up in whichever bra/ket slot it occupied in
///        its source factor (see eval_expr.cpp ~ "target_indices" lambda).
///        For terms equivalent under bra<->ket-swap, this can yield a head
///        with an unconventional bra/ket layout (e.g., 3:1 for a CCSD T2
///        residual) that breaks downstream code reading the result tensor by
///        slot position. Prefer the ResultExpr overload, which lets the
///        caller declare the head's bra/ket layout explicitly via the LHS
///        and is space/convention-agnostic in the IR.
///
/// @sa binarize(ResultExpr const&)
template <typename ExprT = EvalExpr>
  requires std::is_constructible_v<ExprT, EvalExpr>
[[deprecated(
    "binarize(ExprPtr) builds the eval-tree head's bra/ket from external "
    "slot positions in the source factors and can yield an unconventional "
    "layout (e.g. T2 residual head 3:1 instead of 2:2). Use "
    "binarize(ResultExpr) and supply the desired head as the "
    "LHS.")]] FullBinaryNode<ExprT>
binarize(ExprPtr const& expr, IndexSet const& external = {},
         const BinarizationOptions& opts = {}) {
  SEQUANT_ASSERT(
      ranges::all_of(external, [](const auto& idx) { return idx.nonnull(); }));
  auto tree = impl::binarize(expr, external, opts);
  if constexpr (std::is_same_v<ExprT, EvalExpr>)
    return tree;
  else
    return transform_node(std::move(tree),
                          [](auto&& val) { return ExprT{val}; });
}

/// Creates a binary tree for evaluation from a ResultExpr.
/// External indices are deduced from @p res (its non-null slots).
template <typename ExprT = EvalExpr>
  requires std::is_constructible_v<ExprT, EvalExpr>
FullBinaryNode<ExprT> binarize(ResultExpr const& res,
                               const BinarizationOptions& opts = {}) {
  // collect result indices so they are not contracted in the expression tree
  IndexSet uncontract;
  for (auto&& ix : res.indices()) uncontract.emplace(ix);

  // The recursive calls below dispatch to the (now deprecated)
  // binarize(ExprPtr) overload — the head's bra/ket would be positional, but
  // we overwrite it below with res.result_as_tensor() so the layout is
  // caller-determined and the deprecation does not apply.
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  FullBinaryNode<ExprT> tree =
      binarize<ExprT>(res.expression(), uncontract, opts);

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
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  SEQUANT_ASSERT(tree.size() > 1);

  if (is_scalar) {
    if (res.has_label()) {
      tree->expr().template as<Variable>().set_label(res.label());
    }
  } else {
    Tensor& tensor = tree->expr().template as<Tensor>();

    tensor = res.result_as_tensor();
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

  // Adjoint is unary and stores the marker-bearing tensor directly in its
  // own ExprPtr; the bare-leaf left child and Constant(1) right child are
  // structural plumbing for the IR, not part of the symbolic form.
  if (op == EvalOp::Adjoint) return evxpr.expr();

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

#endif  // SEQUANT_EVAL_EVAL_EXPR_HPP
