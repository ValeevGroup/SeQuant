#ifndef SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP
#define SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/rational.hpp>

#include <string>
#include <type_traits>

namespace sequant {

/// @brief Represents the element-wise rational power of a c-number tensor
/// leaf, e.g. `D^{-1}` for a perturbation-theory denominator.
///
/// The power is taken element by element, hence the indices, and so the
/// bra/ket/aux structure, are those of the base.
class TensorPower : public Expr, public AbstractTensor {
 public:
  using exponent_type = rational;

  TensorPower() = delete;
  virtual ~TensorPower() = default;
  /// deep-clones the base
  TensorPower(const TensorPower& other);
  TensorPower(TensorPower&&) = default;
  /// deep-clones the base
  TensorPower& operator=(const TensorPower& other);
  TensorPower& operator=(TensorPower&&) = default;

  /// @param[in] base a c-number AbstractTensor Expr that is not a TensorPower;
  ///            cloned
  /// @param[in] exponent any rational; a non-integer exponent requires the
  ///            symmetry of @p base to not be Symmetry::Antisymm
  /// @throw Exception if a precondition is violated (via SEQUANT_ENFORCE)
  TensorPower(ExprPtr base, exponent_type exponent);

  /// @return the base expression
  const ExprPtr& base() const;

  /// @return the base as an AbstractTensor
  const AbstractTensor& base_tensor() const;

  /// @return the rational exponent
  const exponent_type& exponent() const;

  static constexpr type_rank_type type_rank = expr_type_rank::tensor_power;
  static constexpr std::string static_type_name(
      std::type_identity<TensorPower> = {}) {
    return "sequant::TensorPower";
  }

  type_id_type type_id() const override;

  bool is_scalar() const override;

  ExprPtr clone() const override;

  /// @brief adjoint of the power is the power of the adjoint
  void adjoint() override;

  std::wstring to_latex() const override;

  /// @return the symmetry of the base, except that for an even integer
  /// exponent an Symmetry::Antisymm becomes Symmetry::Symm
  Symmetry _symmetry() const override;
  /// @return the braket symmetry of the base
  BraKetSymmetry _braket_symmetry() const override;
  /// @return the column symmetry of the base
  ColumnSymmetry _column_symmetry() const override;

 private:
  ExprPtr base_;
  exponent_type exponent_;

  hash_type memoizing_hash() const override;

  bool static_equal(const Expr& that) const override;

  bool static_less_than(const Expr& that) const override;
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP
