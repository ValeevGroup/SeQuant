#include <SeQuant/core/expressions/tensor_power.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <utility>

namespace sequant {

TensorPower::TensorPower(ExprPtr base, exponent_type exponent)
    : base_{}, exponent_{std::move(exponent)} {
  SEQUANT_ENFORCE(base && base->is<AbstractTensor>() &&
                      base->as<AbstractTensor>()._is_cnumber(),
                  "TensorPower: base must be a c-number tensor; name a "
                  "compound tensor-valued expression as a ResultExpr and "
                  "raise the resulting leaf");
  SEQUANT_ENFORCE(!base->is<TensorPower>(),
                  "TensorPower: use pow() to compose exponents");
  SEQUANT_ENFORCE(
      denominator(exponent_) == 1 ||
          base->as<AbstractTensor>()._symmetry() != Symmetry::Antisymm,
      "TensorPower: a non-integer power of an antisymmetric "
      "tensor is not defined");
  base_ = base->clone();
}

TensorPower::TensorPower(const TensorPower& other)
    : Expr(other),
      AbstractTensor(other),
      base_(other.base_->clone()),
      exponent_(other.exponent_) {}

TensorPower& TensorPower::operator=(const TensorPower& other) {
  if (this != &other) {
    Expr::operator=(other);
    AbstractTensor::operator=(other);
    base_ = other.base_->clone();
    exponent_ = other.exponent_;
  }
  return *this;
}

const ExprPtr& TensorPower::base() const { return base_; }

const AbstractTensor& TensorPower::base_tensor() const {
  return base_->as<AbstractTensor>();
}

const TensorPower::exponent_type& TensorPower::exponent() const {
  return exponent_;
}

Expr::type_id_type TensorPower::type_id() const {
  return get_type_id<TensorPower>();
}

bool TensorPower::is_scalar() const { return false; }

ExprPtr TensorPower::clone() const {
  auto cloned = ex<TensorPower>(base_, exponent_);
  cloned->as<TensorPower>().copy_canonical_mark(*this);
  return cloned;
}

void TensorPower::adjoint() {
  base_->adjoint();
  reset_hash_value();
}

ExprPtr TensorPower::canonicalize(CanonicalizeOptions opts) {
  return canonicalize_tensor_leaf(
      *this, opts, [this](const AbstractTensor& canonical) {
        *this = static_cast<const TensorPower&>(canonical);
      });
}

std::wstring TensorPower::to_latex() const {
  return io::latex::to_string(*this);
}

Symmetry TensorPower::_symmetry() const {
  const auto base_symmetry = base_tensor()._symmetry();
  if (denominator(exponent_) == 1 && numerator(exponent_) % 2 == 0 &&
      base_symmetry == Symmetry::Antisymm)
    return Symmetry::Symm;
  return base_symmetry;
}

BraKetSymmetry TensorPower::_braket_symmetry() const {
  return base_tensor()._braket_symmetry();
}

ColumnSymmetry TensorPower::_column_symmetry() const {
  return base_tensor()._column_symmetry();
}

AbstractTensor& TensorPower::base_tensor_mut() {
  return base_->as<AbstractTensor>();
}

AbstractTensor* TensorPower::_clone() const { return new TensorPower(*this); }

std::shared_ptr<AbstractTensor> TensorPower::_clone_shared() const {
  return std::make_shared<TensorPower>(*this);
}

AbstractTensor::const_any_view_randsz TensorPower::_bra() const {
  return base_tensor()._bra();
}

AbstractTensor::const_any_view_randsz TensorPower::_ket() const {
  return base_tensor()._ket();
}

AbstractTensor::const_any_view_randsz TensorPower::_aux() const {
  return base_tensor()._aux();
}

AbstractTensor::const_any_view_rand TensorPower::_braket() const {
  return base_tensor()._braket();
}

AbstractTensor::const_any_view_rand TensorPower::_braketaux() const {
  return base_tensor()._braketaux();
}

AbstractTensor::const_any_view_rand TensorPower::_slots() const {
  return base_tensor()._slots();
}

std::size_t TensorPower::_bra_rank() const { return base_tensor()._bra_rank(); }

std::size_t TensorPower::_bra_net_rank() const {
  return base_tensor()._bra_net_rank();
}

std::size_t TensorPower::_ket_rank() const { return base_tensor()._ket_rank(); }

std::size_t TensorPower::_ket_net_rank() const {
  return base_tensor()._ket_net_rank();
}

std::size_t TensorPower::_aux_rank() const { return base_tensor()._aux_rank(); }

std::size_t TensorPower::_num_slots() const {
  return base_tensor()._num_slots();
}

std::size_t TensorPower::_num_indices() const {
  return base_tensor()._num_indices();
}

Hermiticity TensorPower::_hermiticity() const {
  return base_tensor()._hermiticity();
}

Field TensorPower::_base_field() const { return base_tensor()._base_field(); }

std::size_t TensorPower::_color() const { return base_tensor()._color(); }

std::size_t TensorPower::_color_salt() const {
  auto salt = hash::value(exponent_);
  hash::combine(salt, static_cast<std::size_t>(get_type_id<TensorPower>()));
  return salt ? salt : 1;
}

bool TensorPower::_is_cnumber() const { return base_tensor()._is_cnumber(); }

std::wstring_view TensorPower::_label() const { return base_tensor()._label(); }

std::wstring TensorPower::_to_latex() const { return to_latex(); }

std::size_t TensorPower::_hash_value() const { return hash_value(); }

bool TensorPower::operator<(const AbstractTensor& that) const {
  if (const auto* other = dynamic_cast<const TensorPower*>(&that))
    return static_less_than(*other);
  return cross_type_less(*this, that);
}

bool TensorPower::_transform_indices(
    const container::map<Index, Index>& index_map) {
  const auto result = base_tensor_mut()._transform_indices(index_map);
  reset_hash_value();
  return result;
}

void TensorPower::_reset_tags() {
  base_tensor_mut()._reset_tags();
  reset_hash_value();
}

void TensorPower::_swap_bra_ket() {
  base_tensor_mut()._swap_bra_ket();
  reset_hash_value();
}

AbstractTensor::any_view_randsz TensorPower::_bra_mutable() {
  reset_hash_value();
  return base_tensor_mut()._bra_mutable();
}

AbstractTensor::any_view_randsz TensorPower::_ket_mutable() {
  reset_hash_value();
  return base_tensor_mut()._ket_mutable();
}

AbstractTensor::any_view_randsz TensorPower::_aux_mutable() {
  reset_hash_value();
  return base_tensor_mut()._aux_mutable();
}

Expr::hash_type TensorPower::memoizing_hash() const {
  auto compute_hash = [this]() {
    auto val = hash::value(*base_);
    hash::combine(val, hash::value(exponent_));
    return val;
  };

  if (!hash_value_) {
    hash_value_ = compute_hash();
  } else {
    SEQUANT_ASSERT(*hash_value_ == compute_hash());
  }
  return *hash_value_;
}

bool TensorPower::static_equal(const Expr& that) const {
  const auto& other = static_cast<const TensorPower&>(that);
  return exponent_ == other.exponent_ && *base_ == *other.base_;
}

bool TensorPower::static_less_than(const Expr& that) const {
  const auto& other = static_cast<const TensorPower&>(that);
  if (*base_ != *other.base_) return *base_ < *other.base_;
  return exponent_ < other.exponent_;
}

namespace {

ExprPtr make_scalar_power(const ExprPtr& base, rational exponent) {
  if (base->is<Constant>() && base->as<Constant>().is_zero() && exponent < 0)
    throw Exception("pow: 0 raised to a negative exponent is undefined");
  return ex<Power>(base, std::move(exponent));
}

/// rejects the inputs that the TensorPower constructor enforces against
ExprPtr make_tensor_power(const ExprPtr& base, rational exponent) {
  const auto& tensor = base->as<AbstractTensor>();
  if (!tensor._is_cnumber())
    throw Exception("pow: base must be a c-number tensor");
  if (denominator(exponent) != 1 && tensor._symmetry() == Symmetry::Antisymm)
    throw Exception(
        "pow: a non-integer power of an antisymmetric tensor is not defined");
  return ex<TensorPower>(base, std::move(exponent));
}

}  // namespace

ExprPtr pow(ExprPtr base, rational exponent) {
  if (exponent == 1) return base->clone();

  if (base->is<Power>()) {
    const auto& inner = base->as<Power>();
    if (denominator(exponent) != 1 && numerator(inner.exponent()) % 2 == 0)
      throw Exception(
          "pow: (x^a)^b with even a and non-integer b is |x|^(ab), not "
          "x^(ab)");
    const auto composed = inner.exponent() * exponent;
    if (composed == 1) {
      auto result = inner.base()->clone();
      if (inner.conjugated()) result->adjoint();
      return result;
    }
    auto result = make_scalar_power(inner.base(), composed);
    if (inner.conjugated()) result->as<Power>().conjugate();
    return result;
  }

  if (base->is_scalar()) {
    if (!base->is<Constant>() && !base->is<Variable>())
      throw Exception(
          "pow: a scalar base must be a Constant, Variable or Power");
    return make_scalar_power(base, std::move(exponent));
  }

  if (base->is<TensorPower>()) {
    const auto& inner = base->as<TensorPower>();
    if (denominator(exponent) != 1 && numerator(inner.exponent()) % 2 == 0)
      throw Exception(
          "pow: (T^a)^b with even a and non-integer b is |T|^(ab), not "
          "T^(ab); name T^a as a ResultExpr and raise the resulting leaf");
    const auto composed = inner.exponent() * exponent;
    if (composed == 1) return inner.base()->clone();
    return make_tensor_power(inner.base(), composed);
  }

  if (!base->is<AbstractTensor>())
    throw Exception(
        "pow: base must be a scalar or a c-number tensor; name a compound "
        "tensor-valued expression as a ResultExpr and raise the resulting "
        "leaf");
  return make_tensor_power(base, std::move(exponent));
}

}  // namespace sequant
