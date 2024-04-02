#ifndef SEQUANT_EVAL_EXPR_HPP
#define SEQUANT_EVAL_EXPR_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <string>

namespace sequant {

class Tensor;

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
/// \brief Represents the outer indices and the inner indices of a nested
/// tensor.
///
/// \note The nested tensor is a concept that generalizes the sequant::Tensor
/// with and without proto indices. sequant::Tensors with proto indices have
/// outer and inner indices, whereas, those without proto indices only have
/// outer indices.
///
struct NestedTensorIndices {
  container::svector<Index> outer, inner;

  explicit NestedTensorIndices(Tensor const&);
};

///
/// \brief The EvalExpr class represents the object that go into the nodes of
///        the binary tree that is used to evaluate the sequant expressions.
///
/// \details The EvalExpr class itself is not a proper node. It is rather a data
///          that a node should hold.
///
class EvalExpr {
 public:
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
  /// \brief Construct an EvalExpr object from two EvalExpr objects and an
  ///        EvalOp. The EvalOp is either Sum or Prod.
  ///
  EvalExpr(EvalExpr const& left, EvalExpr const& right, EvalOp op);

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

  ///
  /// \brief Get the unique id of this EvalExpr object. Useful for tracing. Not
  ///        used by evaluation.
  ///
  /// \return The unique id of this EvalExpr object.
  ///
  [[nodiscard]] size_t id() const noexcept;

  ///
  /// \return The ExprPtr object that this EvalExpr object holds.
  ///
  [[nodiscard]] ExprPtr expr() const noexcept;

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
  [[nodiscard]] std::string braket_annot() const noexcept;

 private:
  EvalOp op_type_;

  ResultType result_type_;

  size_t hash_value_;

  size_t id_;

  ExprPtr expr_;

  bool tot_;

  static size_t global_id_;
};

}  // namespace sequant

#endif  // SEQUANT_EVAL_EXPR_HPP
