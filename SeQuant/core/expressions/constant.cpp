#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

namespace sequant {

std::wstring Constant::to_latex() const {
  return L"{" + io::latex::to_string(value()) + L"}";
}

Expr::type_id_type Constant::type_id() const { return get_type_id<Constant>(); }

bool Constant::is_scalar() const { return true; }

ExprPtr Constant::clone() const { return ex<Constant>(this->value()); }

void Constant::adjoint() {
  value_ = conj(value_);
  reset_hash_value();
}

Constant &Constant::operator*=(const Expr &that) {
  if (that.is<Constant>()) {
    value_ *= that.as<Constant>().value();
  } else {
    throw Exception("Constant::operator*=(that): not valid for that");
  }

  reset_hash_value();

  return *this;
}

Constant &Constant::operator+=(const Expr &that) {
  if (that.is<Constant>()) {
    value_ += that.as<Constant>().value();
  } else {
    throw Exception("Constant::operator+=(that): not valid for that");
  }

  reset_hash_value();

  return *this;
}

Constant &Constant::operator-=(const Expr &that) {
  if (that.is<Constant>()) {
    value_ -= that.as<Constant>().value();
  } else {
    throw Exception("Constant::operator-=(that): not valid for that");
  }

  reset_hash_value();

  return *this;
}

bool Constant::is_zero(scalar_type v) { return v.is_zero(); }

bool Constant::is_zero() const { return is_zero(this->value()); }

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

}  // namespace sequant
