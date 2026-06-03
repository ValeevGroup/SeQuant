#ifndef SEQUANT_EXPRESSIONS_POWER_HPP
#define SEQUANT_EXPRESSIONS_POWER_HPP

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
/// Variable) and exponent is a rational number.
class Power : public Expr {
 public:
  using exponent_type = rational;

  Power() = delete;
  virtual ~Power() = default;
  Power(const Power&) = default;
  Power(Power&&) = default;
  Power& operator=(const Power&) = default;
  Power& operator=(Power&&) = default;

  /// @param[in] base the base expression; must be a Constant or Variable.
  /// @param[in] exponent rational exponent
  Power(ExprPtr base, exponent_type exponent)
      : base_{}, exponent_{std::move(exponent)} {
    SEQUANT_ASSERT(base);
    SEQUANT_ASSERT(base->is<Constant>() || base->is<Variable>());
    // clone on construction so that external
    // mutations of the input cannot invalidate our memoized hash
    base_ = base->clone();
    // 0^n is defined only for n >= 0 (0^0 = 1 by convention)
    SEQUANT_ASSERT(!base_->is<Constant>() || !base_->as<Constant>().is_zero() ||
                   exponent_ >= 0);
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

  /// @return the rational exponent
  const exponent_type& exponent() const { return exponent_; }

  /// @return whether this Power has been complex-conjugated via adjoint()
  /// @note Conjugation is tracked as a flag because, in general,
  /// `conj(base^exponent) != conj(base)^exponent`
  bool conjugated() const { return conjugated_; }

  /// @brief toggles the conjugation flag
  void conjugate() {
    conjugated_ = !conjugated_;
    reset_hash_value();
  }

  /// @return true if the base is zero and the exponent is positive
  /// @note Construction rejects all undefined 0^n cases; 0^0 is legal and
  /// treated as 1.
  bool is_zero() const override {
    return exponent_ > 0 && base_->is<Constant>() &&
           base_->as<Constant>().is_zero();
  }

  /// @brief Attempts to flatten a Power, mutating @p expr in place. Folds
  /// when @p expr holds a Power and any of:
  ///   - the exponent is 1 (then `b^1 = b` and conjugate if needed);
  ///   - the exponent is 0 (then `b^0 = 1` for any base);
  ///   - the base is the `Constant` 1 (then `1^k = 1` for any rational @c k);
  ///   - the base is a Constant and the exponent is a real integer;
  ///   - the base is a Constant and the exponent has the form `m/2` and the
  ///     base is a non-negative real rational `p/q` with both @c p and @c q
  ///     perfect squares (e.g. `4^{1/2} -> 2`, `(1/4)^{1/2} -> 1/2`).
  /// On success, @p expr is rebound to the folded expression; otherwise it
  /// is left unchanged.
  /// @note Only square-root exponents are folded (that is the only
  /// case needed in practice right now). Extending to general n-th roots only
  /// requires replacing the integer-square-root step with an integer n-th-root.
  static void flatten(ExprPtr& expr) {
    if (!expr || !expr->is<Power>()) return;
    const auto& pw = expr->as<Power>();

    // b^1 = b and conjugate if needed
    if (pw.exponent_ == 1) {
      auto lifted = pw.base_->clone();
      if (pw.conjugated_) lifted->adjoint();
      expr = std::move(lifted);
      return;
    }
    // b^0 = 1 for any base (the ctor rejects 0^(negative)
    if (pw.exponent_ == 0) {
      expr = ex<Constant>(Constant::scalar_type{1});
      return;
    }
    if (!pw.base_->is<Constant>()) return;

    using scalar_type = Constant::scalar_type;
    const auto& base_val = pw.base_->as<Constant>().value();

    // 1^k = 1 for any rational k.
    if (base_val == scalar_type{1}) {
      expr = ex<Constant>(scalar_type{1});
      return;
    }

    // Both remaining fold cases share one shape — `rational base raised to
    // an integer exponent` — so we normalize to that shape and run a single
    // exp-by-squaring loop. Anything else is left untouched.
    //
    // Case A: integer exponent (any Constant base, real or complex).
    //   `base` is just `base_val`.
    // Case B: half-integer exponent on a non-negative real rational base
    //   `p/q` with both `p` and `q` perfect squares. Then
    //     (p/q)^(m/2) = (sqrt(p)/sqrt(q))^m,
    //   so we replace `base` with `sqrt(p)/sqrt(q)` (still a rational) and
    //   keep `exp_int = m`.

    // initialize the base
    scalar_type base{0};
    auto exp_nr = numerator(pw.exponent_);  // numerator of exponent

    if (denominator(pw.exponent_) == 1) {
      base = base_val;
    } else if (denominator(pw.exponent_) == 2 && base_val.imag() == 0 &&
               base_val.real() >= 0) {
      intmax_t p = numerator(base_val.real());
      intmax_t q = denominator(base_val.real());  // > 0 by Boost's convention,
                                                  // sign is with the numerator

      // check for perfect squares
      intmax_t p_rem{0}, q_rem{0};
      intmax_t p_root = boost::multiprecision::sqrt(p, p_rem);
      intmax_t q_root = boost::multiprecision::sqrt(q, q_rem);
      // fold if p and q are perfect squares, else return
      if (p_rem != 0 || q_rem != 0) return;
      base = scalar_type{rational(p_root) / rational(q_root)};
    } else {
      return;
    }

    // Standard exp-by-squaring; for negative exponents we power the
    // magnitude and invert at the end.
    const bool negate = exp_nr < 0;
    if (negate) exp_nr = -exp_nr;
    scalar_type value{1};
    scalar_type b = base;
    while (exp_nr > 0) {
      if (exp_nr % 2 != 0) value *= b;
      exp_nr /= 2;
      if (exp_nr > 0) b *= b;
    }
    if (negate) value = scalar_type{1} / value;

    if (pw.conjugated_) value = conj(value);
    expr = ex<Constant>(std::move(value));
  }

  type_id_type type_id() const override { return get_type_id<Power>(); }

  bool is_scalar() const override { return true; }

  ExprPtr clone() const override {
    auto cloned = ex<Power>(base_, exponent_);
    if (conjugated_) cloned->as<Power>().conjugate();
    return cloned;
  }

  /// @brief adjoint of Power: flips the conjugation flag.
  void adjoint() override { conjugate(); }

  /// @brief Combines exponents when effective bases match:
  ///   - `b^e1 *= b^e2` → `b^(e1+e2)` when this and @p that share the same
  ///      Power-level conjugation flag.
  ///   - `b^e *= B` → `b^(e+1)` for a bare scalar @p B equal to this
  ///      Power's effective base. For a Variable base the labels must match
  ///      and the effective conjugation parities align. For a Constant base
  ///      only the fully unconjugated case combines.
  /// @throw Exception if @p that is not combinable.
  Expr& operator*=(const Expr& that) override {
    // b^e1 *= b^e2  ->  b^(e1+e2)
    if (that.is<Power>()) {
      const auto& other = that.as<Power>();
      if (conjugated_ == other.conjugated_ && *base_ == *other.base_) {
        exponent_ += other.exponent_;
        reset_hash_value();
        return *this;
      }
    }
    // (b^e)* *= b*  ->  (b^(e+1))*
    else if (base_->is<Variable>() && that.is<Variable>()) {
      // check effective conjugation of Variable in this and that, if valid
      // operation iff they match
      const auto& base_var = base_->as<Variable>();
      const auto& that_var = that.as<Variable>();
      if (base_var.label() == that_var.label() &&
          (base_var.conjugated() ^ conjugated_) == that_var.conjugated()) {
        exponent_ += rational{1};
        reset_hash_value();
        return *this;
      }
    }
    // C^e *= C  ->  C^(e+1)
    else if (!conjugated_ && *base_ == that) {
      exponent_ += rational{1};
      reset_hash_value();
      return *this;
    }
    throw Exception("Power::operator*=(that): not valid for that");
  }

 private:
  ExprPtr base_;
  exponent_type exponent_;
  bool conjugated_ = false;

  /// @return hash of this Power
  /// @note when exponent is 1 and not conjugated the hash matches the base's
  hash_type memoizing_hash() const override {
    auto compute_hash = [this]() {
      if (exponent_ == 1 && !conjugated_) return hash::value(*base_);
      auto val = hash::value(*base_);
      hash::combine(val, hash::value(exponent_));
      hash::combine(val, conjugated_);
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
    return exponent_ == other.exponent_ && conjugated_ == other.conjugated_ &&
           *base_ == *other.base_;
  }

  bool static_less_than(const Expr& that) const override {
    const auto& other = static_cast<const Power&>(that);
    if (*base_ != *other.base_) return *base_ < *other.base_;
    if (exponent_ != other.exponent_) return exponent_ < other.exponent_;
    return conjugated_ < other.conjugated_;
  }
};
}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_POWER_HPP
