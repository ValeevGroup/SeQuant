#ifndef SEQUANT_EXPRESSIONS_CONSTANT_HPP
#define SEQUANT_EXPRESSIONS_CONSTANT_HPP

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include <string>

namespace sequant {

namespace {
template <typename X>
X numeric_cast(const sequant::rational &r) {
  if constexpr (std::is_integral_v<X>) {
    SEQUANT_ASSERT(denominator(r) == 1);
    return boost::numeric_cast<X>(numerator(r));
  } else {
    return boost::numeric_cast<X>(numerator(r)) /
           boost::numeric_cast<X>(denominator(r));
  }
};
}  // namespace

/// @brief a constant number

/// This is represented as a "compile-time" complex rational number
class Constant : public Expr {
 public:
  using scalar_type = Complex<sequant::rational>;

 public:
  Constant() = delete;
  virtual ~Constant() = default;
  Constant(const Constant &) = default;
  Constant(Constant &&) = default;
  Constant &operator=(const Constant &) = default;
  Constant &operator=(Constant &&) = default;
  template <typename U, typename = std::enable_if_t<!is_constant_v<U>>>
  explicit Constant(U &&value) : value_(std::forward<U>(value)) {}

  /// @tparam T the result type; default to the type of value_
  /// @return the value cast to ResultType
  /// @throw std::invalid_argument if conversion to T is not possible
  /// @throw boost::numeric::positive_overflow or
  /// boost::numeric::negative_overflow if cast fails
  template <typename T = scalar_type>
  auto value() const {
    if constexpr (std::is_arithmetic_v<T>) {
      SEQUANT_ASSERT(value_.imag() == 0);
      return numeric_cast<T>(value_.real());
    } else if constexpr (meta::is_complex_v<T>) {
      return T(numeric_cast<typename T::value_type>(value_.real()),
               numeric_cast<typename T::value_type>(value_.imag()));
    } else
      throw std::invalid_argument(
          "Constant::value<T>: cannot convert value to type T");
  }

  std::wstring to_latex() const override {
    return L"{" + sequant::to_latex(value()) + L"}";
  }

  std::wstring to_wolfram() const override {
    return sequant::to_wolfram(value());
  }

  type_id_type type_id() const override { return get_type_id<Constant>(); }

  ExprPtr clone() const override { return ex<Constant>(this->value()); }

  /// @brief adjoint of a Constant is its complex conjugate
  virtual void adjoint() override;

  virtual Expr &operator*=(const Expr &that) override {
    if (that.is<Constant>()) {
      value_ *= that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator*=(that): not valid for that");
    }
    return *this;
  }

  virtual Expr &operator+=(const Expr &that) override {
    if (that.is<Constant>()) {
      value_ += that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator+=(that): not valid for that");
    }
    return *this;
  }

  virtual Expr &operator-=(const Expr &that) override {
    if (that.is<Constant>()) {
      value_ -= that.as<Constant>().value();
    } else {
      throw std::logic_error("Constant::operator-=(that): not valid for that");
    }
    return *this;
  }

  /// @param[in] v a scalar
  /// @return true if this is a soft zero, i.e. its magnitude is less than
  /// `std::sqrt(std::numeric_limits<float>::epsilon())`
  static bool is_zero(scalar_type v) { return v.is_zero(); }

  /// @return `Constant::is_zero(this->value())`
  bool is_zero() const { return is_zero(this->value()); }

 private:
  scalar_type value_;

  hash_type memoizing_hash() const override {
    if (!hash_value_) {
      hash_value_ = hash::value(value_);
    } else {
      SEQUANT_ASSERT(*hash_value_ == hash::value(value_));
    }
    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return value() == static_cast<const Constant &>(that).value();
  }
};  // class Constant

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_CONSTANT_HPP
