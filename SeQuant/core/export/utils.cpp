//
// Created by Ajay Melekamburath on 4/27/26.
//

#include <SeQuant/core/export/utils.hpp>

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
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

}  // namespace sequant::detail
