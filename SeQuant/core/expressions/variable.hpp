#ifndef SEQUANT_EXPRESSIONS_VARIABLE_HPP
#define SEQUANT_EXPRESSIONS_VARIABLE_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/labeled.hpp>
#include <SeQuant/core/hash.hpp>

#include <string>
#include <string_view>
#include <cassert>

namespace sequant {

/// This is represented as a "run-time" complex rational number
class Variable : public Expr, public MutatableLabeled {
 public:
  Variable() = delete;
  virtual ~Variable() = default;
  Variable(const Variable &) = default;
  Variable(Variable &&) = default;
  Variable &operator=(const Variable &) = default;
  Variable &operator=(Variable &&) = default;
  template <typename U, typename = std::enable_if_t<!is_variable_v<U>>>
  explicit Variable(U &&label) : label_(std::forward<U>(label)) {}

  Variable(std::wstring label) : label_(std::move(label)), conjugated_(false) {}

  /// @return variable label
  /// @warning conjugation does not change it
  std::wstring_view label() const override;

  void set_label(std::wstring label) override;

  /// complex-conjugates this
  void conjugate();

  /// @return whether this object has been conjugated
  bool conjugated() const;

  std::wstring to_latex() const override;

  type_id_type type_id() const override { return get_type_id<Variable>(); }

  ExprPtr clone() const override;

  /// @brief adjoint of a Variable is its complex conjugate
  virtual void adjoint() override;

 private:
  std::wstring label_;
  bool conjugated_ = false;

  hash_type memoizing_hash() const override {
    auto compute_hash = [this]() {
      auto val = hash::value(label_);
      hash::combine(val, conjugated_);
      return val;
    };

    if (!hash_value_) {
      hash_value_ = compute_hash();
    } else {
      assert(*hash_value_ == compute_hash());
    }

    return *hash_value_;
  }

  bool static_equal(const Expr &that) const override {
    return label_ == static_cast<const Variable &>(that).label_ &&
           conjugated_ == static_cast<const Variable &>(that).conjugated_;
  }
};  // class Variable

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_VARIABLE_HPP
