//
// Created by Eduard Valeyev on 2/8/19.
//

#ifndef SEQUANT_WOLFRAM_HPP
#define SEQUANT_WOLFRAM_HPP

#include <type_traits>
#include "meta.hpp"

namespace sequant {

template <typename T>
std::enable_if_t<meta::has_memfn_to_wolfram_v<std::decay_t<T>>, std::wstring>
to_wolfram(T &&t) {
  return t.to_wolfram();
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, std::wstring>
to_wolfram(T &&t) {
  using ::sequant::to_wstring;
  return to_wstring(t);
}

template <typename T>
std::wstring to_wolfram(const std::complex<T> &t) {
  using ::sequant::to_wstring;
  if (t.imag() == 0)
    return to_wolfram(t.real());
  else
    return std::wstring(L"Complex[") + to_wstring(t.real()) + L"," +
           to_wstring(t.imag()) + L"]";
}

template <typename T>
std::wstring to_wolfram(const boost::rational<T> &t) {
  using ::sequant::to_wstring;
  if (t.denominator() == 1) {
    // n.b. use to_string to skip extra braces so that output agrees with code
    // that used scalars return to_wolfram(t.numerator());
    return to_wstring(t.numerator());
  } else
    return std::wstring(L"Rational[") + to_wolfram(t.numerator()) + L"," +
           to_wolfram(t.denominator()) + L"]";
}

}  // namespace sequant

#endif  // SEQUANT_WOLFRAM_HPP
