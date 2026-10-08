#ifndef SEQUANT_EXPRESSIONS_VARIABLE_HPP
#define SEQUANT_EXPRESSIONS_VARIABLE_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/labeled.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string>
#include <string_view>
#include <type_traits>

namespace sequant {

class ExprPtr;

/// This is represented as a "run-time" complex rational number
class Variable : public Expr, public MutatableLabeled {
 public:
  Variable() = delete;
  virtual ~Variable() = default;
  Variable(const Variable &) = default;
  Variable(Variable &&) = default;
  Variable &operator=(const Variable &) = default;
  Variable &operator=(Variable &&) = default;
  template <typename U>
    requires(!is_variable_v<U> && !is_an_expr_v<std::remove_reference_t<U>> &&
             !Expr::is_shared_ptr_of_expr_or_derived<
                 std::remove_reference_t<U>>::value &&
             std::constructible_from<std::wstring, U>)
  explicit Variable(U &&label) : label_(std::forward<U>(label)) {
    adopt_marks();
  }

  /// @param label the name; a trailing `꙳` (sequant::conjugate_label) is
  ///        adopted as the conjugated state, so that a printed name
  ///        (decorated_label()) reads back as the variable it printed
  /// @throw Exception if @p label ends in an adjoint mark `⁺`, which a
  ///        variable has no state for, or repeats a mark
  Variable(std::wstring label);

  /// @copydoc Variable(std::wstring)
  Variable(const std::string &label);

  /// @return variable label
  /// @warning conjugation does not change it
  std::wstring_view label() const override;

  /// @param label the new label; a trailing `꙳` sets the conjugated state and
  ///        a plain label keeps it, as in the constructors
  /// @throw Exception as the constructors do, in which case the variable is
  ///        left as the call found it
  void set_label(std::wstring label) override;

  /// complex-conjugates this
  void conjugate();

  /// @return whether this object has been conjugated
  bool conjugated() const;

  /// @return label() followed by the `꙳` of a conjugated variable: the printed
  ///         name, the counterpart of Tensor::decorated_label()
  std::wstring decorated_label() const;

  std::wstring to_latex() const override;

  static constexpr type_rank_type type_rank = expr_type_rank::variable;
  static constexpr std::string static_type_name(
      std::type_identity<Variable> = {}) {
    return "sequant::Variable";
  }

  type_id_type type_id() const override;

  bool is_scalar() const override;

  ExprPtr clone() const override;

  /// @brief adjoint of a Variable is its complex conjugate
  [[nodiscard]] virtual std::int8_t adjoint() override;

  /// @brief K-conjugate of a Variable is its complex conjugate
  [[nodiscard]] std::int8_t kconjugate() override;

 private:
  /// adopts a trailing `꙳` of label_ into conjugated_; refuses `⁺`
  void adopt_marks();

  std::wstring label_;
  bool conjugated_ = false;

  hash_type memoizing_hash() const override;

  bool static_equal(const Expr &that) const override;

};  // class Variable

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_VARIABLE_HPP
