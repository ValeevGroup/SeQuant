//
// Created by Eduard Valeyev on 2/8/19.
//

#ifndef SEQUANT_WOLFRAM_HPP
#define SEQUANT_WOLFRAM_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/wstring.hpp>

#include <string>
#include <type_traits>

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

}  // namespace sequant

#endif  // SEQUANT_WOLFRAM_HPP
