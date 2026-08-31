#ifndef SEQUANT_EXPRESSIONS_POWER_HPP
#define SEQUANT_EXPRESSIONS_POWER_HPP

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/rational.hpp>

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
  Power(ExprPtr base, exponent_type exponent);

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
  const ExprPtr& base() const;

  /// @return the rational exponent
  const exponent_type& exponent() const;

  /// @return whether this Power has been complex-conjugated via adjoint()
  /// @note Conjugation is tracked as a flag because, in general,
  /// `conj(base^exponent) != conj(base)^exponent`
  bool conjugated() const;

  /// @brief toggles the conjugation flag
  void conjugate();

  /// @return true if the base is zero and the exponent is positive
  /// @note Construction rejects all undefined 0^n cases; 0^0 is legal and
  /// treated as 1.
  bool is_zero() const override;

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
  static void flatten(ExprPtr& expr);

  type_id_type type_id() const override;

  bool is_scalar() const override;

  ExprPtr clone() const override;

  /// @brief adjoint of Power: flips the conjugation flag.
  void adjoint() override;

  /// @brief Combines exponents when effective bases match:
  ///   - `b^e1 *= b^e2` → `b^(e1+e2)` when this and @p that share the same
  ///      Power-level conjugation flag.
  ///   - `b^e *= B` → `b^(e+1)` for a bare scalar @p B equal to this
  ///      Power's effective base. For a Variable base the labels must match
  ///      and the effective conjugation parities align. For a Constant base
  ///      only the fully unconjugated case combines.
  /// @throw Exception if @p that is not combinable.
  Power& operator*=(const Expr& that);

 private:
  ExprPtr base_;
  exponent_type exponent_;
  bool conjugated_ = false;

  /// @return hash of this Power
  /// @note when exponent is 1 and not conjugated the hash matches the base's
  hash_type memoizing_hash() const override;

  bool static_equal(const Expr& that) const override;

  bool static_less_than(const Expr& that) const override;
};
}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_POWER_HPP
