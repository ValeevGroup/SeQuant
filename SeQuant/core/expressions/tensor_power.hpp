#ifndef SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP
#define SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/rational.hpp>

#include <memory>
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
  /// exponent an antisymmetric base becomes Symmetry::Symm
  Symmetry _symmetry() const override;
  /// @return the braket symmetry of the base
  BraKetSymmetry _braket_symmetry() const override;
  /// @return the column symmetry of the base
  ColumnSymmetry _column_symmetry() const override;

  // the rest of the AbstractTensor interface forwards to the base; mutating
  // members flush the memoized hash
  AbstractTensor* _clone() const override;
  std::shared_ptr<AbstractTensor> _clone_shared() const override;
  const_any_view_randsz _bra() const override;
  const_any_view_randsz _ket() const override;
  const_any_view_randsz _aux() const override;
  const_any_view_rand _braket() const override;
  const_any_view_rand _braketaux() const override;
  const_any_view_rand _slots() const override;
  std::size_t _bra_rank() const override;
  std::size_t _bra_net_rank() const override;
  std::size_t _ket_rank() const override;
  std::size_t _ket_net_rank() const override;
  std::size_t _aux_rank() const override;
  std::size_t _num_slots() const override;
  std::size_t _num_indices() const override;
  Hermiticity _hermiticity() const override;
  Field _base_field() const override;
  std::size_t _color() const override;
  bool _is_cnumber() const override;
  std::wstring_view _label() const override;
  std::wstring _to_latex() const override;
  std::size_t _hash_value() const override;
  bool operator<(const AbstractTensor& that) const override;
  bool _transform_indices(
      const container::map<Index, Index>& index_map) override;
  void _reset_tags() override;
  void _swap_bra_ket() override;
  any_view_randsz _bra_mutable() override;
  any_view_randsz _ket_mutable() override;
  any_view_randsz _aux_mutable() override;

 private:
  ExprPtr base_;
  exponent_type exponent_;

  AbstractTensor& base_tensor_mut();

  hash_type memoizing_hash() const override;

  bool static_equal(const Expr& that) const override;

  bool static_less_than(const Expr& that) const override;
};

/// @brief raises a scalar or a c-number tensor leaf to a rational power
///
/// A Constant or Variable base gives a Power, a c-number tensor base a
/// TensorPower. Exponents of nested powers are composed when that is exact:
/// for a TensorPower `T^a` and exponent `b`, `(T^a)^b = T^{ab}` for real `T`
/// only if `b` is an integer or `a` has an odd numerator.
/// @param[in] base a Constant, Variable, Power, c-number AbstractTensor or
///            TensorPower
/// @param[in] exponent any rational
/// @throw Exception if @p base is not supported or the composition is inexact
ExprPtr pow(ExprPtr base, rational exponent);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_TENSOR_POWER_HPP
