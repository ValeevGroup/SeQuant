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

  /// @overload
  /// @param[in] exponent real rational exponent (imaginary part is zero)
  Power(ExprPtr base, rational exponent)
      : Power(std::move(base), exponent_type{std::move(exponent), 0}) {}

  /// @return the base expression
  const ExprPtr& base() const { return base_; }

  /// @return the complex rational exponent
  const exponent_type& exponent() const { return exponent_; }

  /// @return true if the base is zero and the exponent has positive real part
  bool is_zero() const override {
    return exponent_.real() > 0 && base_->is<Constant>() &&
           base_->as<Constant>().is_zero();
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

  /// @brief adjoint of Power: (base^exp)† = (base†)^(exp*)
  void adjoint() override {
    base_->adjoint();
    exponent_ = conj(exponent_);
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
};

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_POWER_HPP
