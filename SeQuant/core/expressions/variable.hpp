#ifndef SEQUANT_EXPRESSIONS_VARIABLE_HPP
#define SEQUANT_EXPRESSIONS_VARIABLE_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/labeled.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <memory>
#include <string>
#include <string_view>

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
  explicit Variable(U &&label) : label_(std::forward<U>(label)) {}

  Variable(std::wstring label);

  Variable(const std::string &label);

  /// @return variable label
  /// @warning conjugation does not change it
  std::wstring_view label() const override;

  void set_label(std::wstring label) override;

  /// complex-conjugates this
  void conjugate();

  /// @return whether this object has been conjugated
  bool conjugated() const;

  std::wstring to_latex() const override;

  type_id_type type_id() const override;

  bool is_scalar() const override;

  /// @brief adjoint of a Variable is its complex conjugate
  virtual void adjoint() override;

 protected:
  std::unique_ptr<Expr> unique_copy() const override;

 private:
  std::wstring label_;
  bool conjugated_ = false;

  hash_type memoizing_hash() const override;

  bool static_equal(const Expr &that) const override;

};  // class Variable

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_VARIABLE_HPP
