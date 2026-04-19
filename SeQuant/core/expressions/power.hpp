#ifndef SEQUANT_EXPRESSIONS_POWER_HPP
#define SEQUANT_EXPRESSIONS_POWER_HPP

#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/macros.hpp>

namespace sequant {

/// @brief Represents base^exponent where base is a scalar (Constant or
/// Variable) and exponent is a complex rational number.
class Power : public Expr {
 public:
  using exponent_type = Complex<rational>;

  Power() = delete;
  virtual ~Power() = default;
  Power(const Power&) = default;
  Power(Power&&) = default;
  Power& operator=(const Power&) = default;
  Power& operator=(Power&&) = default;

  /// @param[in] base the base expression; must be a Constant, Variable, or
  ///   Power (in the last case the result is flattened:
  ///   `Power(Power(b,e1), e2)` → `Power(b, e1*e2)`)
  /// @param[in] exponent complex rational exponent
  Power(ExprPtr base, exponent_type exponent)
      : base_{}, exponent_{std::move(exponent)} {
    SEQUANT_ASSERT(base);
    if (base->is<Power>()) {
      auto& inner = base->as<Power>();
      base_ = inner.base();
      exponent_ = inner.exponent() * exponent_;
    } else {
      SEQUANT_ASSERT(base->is<Constant>() || base->is<Variable>());
      base_ = std::move(base);
    }
  }

  /// @overload constructs a `Variable` base from @p label
  template <typename L>
    requires std::constructible_from<std::wstring, L> &&
             (!std::convertible_to<L, ExprPtr>)
  Power(L&& label, exponent_type exponent)
      : Power(ex<Variable>(std::forward<L>(label)), std::move(exponent)) {}

  /// @overload constructs a `Constant` base from scalar @p value
  template <typename V>
    requires(!std::constructible_from<std::wstring, V> &&
             !std::convertible_to<V, ExprPtr> &&
             std::constructible_from<Constant::scalar_type, V>)
  Power(V&& value, exponent_type exponent)
      : Power(ex<Constant>(std::forward<V>(value)), std::move(exponent)) {}

  /// @return the base expression
  const ExprPtr& base() const { return base_; }

  /// @return the complex rational exponent
  const exponent_type& exponent() const { return exponent_; }

  /// @return true if the base is zero and the exponent has positive real part
  bool is_zero() const override {
    return exponent_.real() > 0 && base_->is<Constant>() &&
           base_->as<Constant>().is_zero();
  }

  /// @brief Attempts to flatten a Power into a Constant, mutating @p expr
  /// in place. No-op unless @p expr holds a Power whose:
  ///   - base is a Constant
  ///   - exponent is a real integer
  /// On success, @p expr is rebound to the folded Constant; otherwise it is
  /// left unchanged.
  static void flatten(ExprPtr& expr) {
    if (!expr || !expr->is<Power>()) return;
    const auto& self = expr->as<Power>();
    if (!self.base_->is<Constant>()) return;
    if (self.exponent_.imag() != 0) return;
    const auto& expr_real = self.exponent_.real();
    if (denominator(expr_real) != 1) return;

    auto exp_numerator = numerator(expr_real);
    const auto& base_val = self.base_->as<Constant>().value();
    using scalar_type = Constant::scalar_type;

    // zero power
    if (exp_numerator == 0) {
      expr = ex<Constant>(scalar_type{1});
      return;
    }
    // negative power
    const bool negate = exp_numerator < 0;
    if (negate) exp_numerator = -exp_numerator;

    scalar_type value{1};
    for (auto i = exp_numerator; i > 0; --i) value *= base_val;
    if (negate) value = scalar_type{1} / value;

    expr = ex<Constant>(std::move(value));
  }

  std::wstring to_latex() const override {
    if (exponent_ == 1) {
      return base_->to_latex();
    }
    if (exponent_.imag() == 0) {
      return base_->to_latex() + L"^" + io::latex::to_string(exponent_.real());
    }
    return base_->to_latex() + L"^" + exponent_.to_latex();
  }

  type_id_type type_id() const override { return get_type_id<Power>(); }

  ExprPtr clone() const override {
    return ex<Power>(base_->clone(), exponent_);
  }

  /// @brief adjoint of Power
  void adjoint() override {
    base_ = ::sequant::adjoint(base_);
    exponent_ = conj(exponent_);
    reset_hash_value();
  }

  /// @brief Combines exponents when bases match:
  ///   - `b^e1 *= b^e2` → `b^(e1+e2)`
  ///   - `b^e *= b` → `b^(e+1)` (treats bare base as base^1)
  /// @throw Exception if bases differ or @p that is not combinable
  Expr& operator*=(const Expr& that) override {
    if (that.is<Power>()) {
      const auto& other = that.as<Power>();
      if (*base_ == *other.base_) {
        exponent_ += other.exponent_;
        reset_hash_value();
        return *this;
      }
    } else if (*base_ == that) {
      exponent_ += rational{1};
      reset_hash_value();
      return *this;
    }
    throw Exception("Power::operator*=(that): not valid for that");
  }

 private:
  ExprPtr base_;
  exponent_type exponent_;

  hash_type memoizing_hash() const override {
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

  bool static_equal(const Expr& that) const override {
    const auto& other = static_cast<const Power&>(that);
    return exponent_ == other.exponent_ && *base_ == *other.base_;
  }

  bool static_less_than(const Expr& that) const override {
    const auto& other = static_cast<const Power&>(that);
    if (*base_ != *other.base_) return *base_ < *other.base_;
    if (exponent_.real() != other.exponent_.real())
      return exponent_.real() < other.exponent_.real();
    return exponent_.imag() < other.exponent_.imag();
  }
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_POWER_HPP
