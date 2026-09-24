//
// Created by Ajay Melekamburath on 4/27/26.
//

#include <SeQuant/core/export/utils.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/rational.hpp>

#include <sstream>
#include <utility>

namespace sequant::detail {

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
