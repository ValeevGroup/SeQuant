//
// Created by Ajay Melekamburath on 4/27/26.
//

#include <SeQuant/core/export/utils.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cctype>
#include <cstdint>
#include <sstream>
#include <string>
#include <utility>

namespace sequant::detail {

std::string basis_instance_tag(const Index &idx) {
  if (!idx.basis().has_basis_instance()) return {};
  // widened so that the magnitude of the most negative instance fits
  const std::int64_t instance = *idx.basis().basis_instance();
  return instance < 0 ? "_" + std::to_string(-instance)
                      : std::to_string(instance);
}

std::string block_tag(std::string space_tag, const Index &idx) {
  const auto instance = basis_instance_tag(idx);
  if (!instance.empty() && !space_tag.empty() &&
      (std::isdigit(static_cast<unsigned char>(space_tag.back())) ||
       space_tag.back() == '_'))
    throw Exception("export: the tag '" + space_tag + "' of the space of " +
                    toUtf8(idx.full_label()) +
                    " ends in a digit or '_', so a basis instance appended to "
                    "it cannot be told from the tag; tag the space otherwise");
  return space_tag += instance;
}

std::string dim_name(std::string space_dim, const Index &idx) {
  const auto instance = basis_instance_tag(idx);
  if (instance.empty()) return space_dim;
  if (!space_dim.empty() && space_dim.back() == '_')
    throw Exception("export: the dimension name '" + space_dim +
                    "' of the space of " + toUtf8(idx.full_label()) +
                    " ends in '_', so the name of a basis instance appended "
                    "to it is ambiguous; name the dimension otherwise");
  return space_dim += "_" + instance;
}

std::string format_power_exponent(const Power::exponent_type &exponent,
                                  bool double_slash) {
  std::stringstream ss;
  if (denominator(exponent) == 1) {
    const auto n = numerator(exponent);
    if (n < 0) {
      ss << "(" << n << ")";
    } else {
      ss << n;
    }
  } else {
    ss << "(" << numerator(exponent) << (double_slash ? "//" : "/")
       << denominator(exponent) << ")";
  }
  return ss.str();
}

std::string format_power_base(const ExprPtr &base, std::string base_str) {
  if (base->is<Constant>()) {
    const auto &v = base->as<Constant>().value();
    if (v.imag() == 0 &&
        (denominator(v.real()) != 1 || numerator(v.real()) < 0)) {
      return "(" + std::move(base_str) + ")";
    }
  }
  return base_str;
}

std::string format_power(
    const Power &power, std::string base_str, std::string_view pow_op,
    bool double_slash,
    const std::function<std::string(std::string)> &wrap_conj) {
  const ExprPtr &base = power.base();
  if (base->is<Variable>() && base->as<Variable>().conjugated()) {
    base_str = wrap_conj(std::move(base_str));
  }
  auto s = format_power_base(base, std::move(base_str)) + std::string(pow_op) +
           format_power_exponent(power.exponent(), double_slash);
  if (power.conjugated()) s = wrap_conj(std::move(s));
  return s;
}

}  // namespace sequant::detail
