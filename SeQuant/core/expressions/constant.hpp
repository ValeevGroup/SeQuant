#ifndef SEQUANT_EXPRESSIONS_CONSTANT_HPP
#define SEQUANT_EXPRESSIONS_CONSTANT_HPP

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include <memory>
#include <string>

namespace sequant {

class ExprPtr;

// implementation details of Constant; prefer sequant::detail over an unnamed
// namespace in a header (see CppCoreGuidelines SF.21)
namespace detail {
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
}  // namespace detail

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
  template <typename U>
    requires(!is_constant_v<U> && !is_an_expr_v<std::remove_reference_t<U>> &&
             !Expr::is_shared_ptr_of_expr_or_derived<
                 std::remove_reference_t<U>>::value &&
             std::constructible_from<scalar_type, U>)
  explicit Constant(U &&value) : value_(std::forward<U>(value)) {}

  /// @tparam T the result type; default to the type of value_
  /// @return the value cast to ResultType
  /// @throw Exception if conversion to T is not possible
  /// @throw boost::numeric::positive_overflow or
  /// boost::numeric::negative_overflow if cast fails
  template <typename T = scalar_type>
  auto value() const {
    if constexpr (std::is_arithmetic_v<T>) {
      SEQUANT_ASSERT(value_.imag() == 0);
      return detail::numeric_cast<T>(value_.real());
    } else if constexpr (meta::is_complex_v<T>) {
      return T(detail::numeric_cast<typename T::value_type>(value_.real()),
               detail::numeric_cast<typename T::value_type>(value_.imag()));
    } else
      throw Exception("Constant::value<T>: cannot convert value to type T");
  }

  std::wstring to_latex() const override;

  type_id_type type_id() const override;

  bool is_scalar() const override;

  /// @brief adjoint of a Constant is its complex conjugate
  virtual void adjoint() override;

  Constant &operator*=(const Expr &that);

  Constant &operator+=(const Expr &that);

  Constant &operator-=(const Expr &that);

  /// @param[in] v a scalar
  /// @return true if this is zero
  static bool is_zero(scalar_type v);

  /// @return `Constant::is_zero(this->value())`
  bool is_zero() const final;

 protected:
  std::unique_ptr<Expr> unique_copy() const override;

 private:
  scalar_type value_;

  hash_type memoizing_hash() const override;

  bool static_equal(const Expr &that) const override;

};  // class Constant

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_CONSTANT_HPP
