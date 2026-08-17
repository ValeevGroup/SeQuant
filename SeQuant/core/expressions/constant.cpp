#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>

namespace sequant {

std::wstring Constant::to_latex() const {
  return L"{" + io::latex::to_string(value()) + L"}";
}

Expr::type_id_type Constant::type_id() const { return get_type_id<Constant>(); }

bool Constant::is_scalar() const { return true; }

void Constant::adjoint() {
  value_ = conj(value_);
  reset_hash_value();
}

Constant &Constant::operator*=(const Constant &that) {
  value_ *= that.value();

  reset_hash_value();

  return *this;
}

Constant &Constant::operator*=(const Expr &that) {
  if (!that.is<Constant>()) {
    throw Exception("Constant::operator*=(that): not valid for that");
  }

  return *this *= that.as<Constant>();
}

Constant &Constant::operator+=(const Constant &that) {
  value_ += that.value();

  reset_hash_value();

  return *this;
}

Constant &Constant::operator+=(const Expr &that) {
  if (!that.is<Constant>()) {
    throw Exception("Constant::operator+=(that): not valid for that");
  }

  return *this += that.as<Constant>();
}

Constant &Constant::operator-=(const Constant &that) {
  value_ -= that.value();

  reset_hash_value();

  return *this;
}

Constant &Constant::operator-=(const Expr &that) {
  if (!that.is<Constant>()) {
    throw Exception("Constant::operator-=(that): not valid for that");
  }

  return *this -= that.as<Constant>();
}

bool Constant::is_zero(scalar_type v) { return v.is_zero(); }

bool Constant::is_zero() const { return is_zero(this->value()); }

std::unique_ptr<Expr> Constant::unique_copy() const {
  auto result = std::make_unique<Constant>(this->value());
  result->copy_canonical_mark(*this);
  return result;
}

Expr::hash_type Constant::memoizing_hash() const {
  if (!hash_value_) {
    hash_value_ = hash::value(value_);
  } else {
    SEQUANT_ASSERT(*hash_value_ == hash::value(value_));
  }
  return *hash_value_;
}

bool Constant::static_equal(const Expr &that) const {
  return value() == static_cast<const Constant &>(that).value();
}

bool Constant::static_less_than(const Expr &that) const {
  const auto &that_value = static_cast<const Constant &>(that).value_;
  if (value_.real() != that_value.real())
    return value_.real() < that_value.real();
  return value_.imag() < that_value.imag();
}

std::weak_ordering Constant::static_lexicographic_compare(
    const Expr &that, const LexicographicCompareOptions &) const {
  const auto &that_value = static_cast<const Constant &>(that).value_;
  if (auto cmp =
          std::compare_weak_order_fallback(value_.real(), that_value.real());
      cmp != 0)
    return cmp;
  return std::compare_weak_order_fallback(value_.imag(), that_value.imag());
}

Constant operator*(const Constant &lhs, const Constant &rhs) {
  Constant result(lhs);

  result *= rhs;

  return result;
}

Constant operator+(const Constant &lhs, const Constant &rhs) {
  Constant result(lhs);

  result += rhs;

  return result;
}

Constant operator-(const Constant &lhs, const Constant &rhs) {
  Constant result(lhs);

  result -= rhs;

  return result;
}

}  // namespace sequant
