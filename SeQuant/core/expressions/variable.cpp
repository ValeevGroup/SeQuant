#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string>
#include <string_view>
#include <tuple>

namespace sequant {

Variable::Variable(std::wstring label)
    : label_(std::move(label)), conjugated_(false) {
  adopt_marks();
}

Variable::Variable(const std::string &label)
    : label_(sequant::toUtf16(label)), conjugated_(false) {
  adopt_marks();
}

void Variable::adopt_marks() {
  const auto [adjointed, kconjugated] = split_state_marks(label_);
  if (adjointed)
    throw Exception("Variable: a variable name carries no adjoint mark");
  if (kconjugated) conjugated_ = true;
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

std::wstring_view Variable::label() const { return label_; }

void Variable::set_label(std::wstring label) {
  auto before = std::tuple{std::move(label_), conjugated_};
  label_ = std::move(label);
  try {
    adopt_marks();
  } catch (...) {
    std::tie(label_, conjugated_) = std::move(before);
    throw;
  }
  reset_hash_value();
}

void Variable::conjugate() {
  conjugated_ = !conjugated_;
  reset_hash_value();
}

bool Variable::conjugated() const { return conjugated_; }

std::wstring Variable::decorated_label() const {
  std::wstring result(label_);
  if (conjugated_) result.push_back(sequant::conjugate_label);
  return result;
}

std::wstring Variable::to_latex() const {
  std::wstring result = L"{" + io::latex::utf_to_string(label_) + L"}";
  if (conjugated_) result = L"{" + result + L"^{*}}";
  return result;
}

ExprPtr Variable::clone() const { return ex<Variable>(*this); }

std::int8_t Variable::adjoint() {
  conjugate();
  return 1;
}

std::int8_t Variable::kconjugate() { return adjoint(); }

}  // namespace sequant
