#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string>
#include <string_view>

namespace sequant {

Variable::Variable(std::string label)
    : label_(std::move(label)), conjugated_(false) {
  SEQUANT_ASSERT(is_valid_utf8(label_), "Variable label must be UTF-8 encoded");
}

Expr::type_id_type Variable::type_id() const { return get_type_id<Variable>(); }

bool Variable::is_scalar() const { return true; }

Expr::hash_type Variable::memoizing_hash() const {
  auto compute_hash = [this]() {
    auto val = hash::value(label_);
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

bool Variable::static_equal(const Expr &that) const {
  return label_ == static_cast<const Variable &>(that).label_ &&
         conjugated_ == static_cast<const Variable &>(that).conjugated_;
}

std::string_view Variable::label() const { return label_; }

void Variable::set_label(std::string label) {
  SEQUANT_ASSERT(is_valid_utf8(label), "Variable label must be UTF-8 encoded");
  label_ = std::move(label);
  reset_hash_value();
}

void Variable::conjugate() {
  conjugated_ = !conjugated_;
  reset_hash_value();
}

bool Variable::conjugated() const { return conjugated_; }

std::string Variable::to_latex() const {
  std::string result = "{" + io::latex::utf_to_string(label_) + "}";
  if (conjugated_) result = "{" + result + "^*" + "}";
  return result;
}

ExprPtr Variable::clone() const { return ex<Variable>(*this); }

void Variable::adjoint() { conjugate(); }

}  // namespace sequant
