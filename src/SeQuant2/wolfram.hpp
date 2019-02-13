//
// Created by Eduard Valeyev on 2/8/19.
//

#ifndef SEQUANT2_WOLFRAM_HPP
#define SEQUANT2_WOLFRAM_HPP

#include <type_traits>

namespace sequant2 {

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, std::wstring>
to_wolfram(T &&t) {
  using ::sequant2::to_wstring;
  return to_wstring(t);
}

template <typename T>
std::wstring to_wolfram(const std::complex<T> &t) {
  using ::sequant2::to_wstring;
  if (t.imag() == 0)
    return to_wolfram(t.real());
  else
    return std::wstring(L"Complex[") + to_wstring(t.real()) + L"," +
           to_wstring(t.imag()) + L"]";
}

}  // namespace sequant2

#endif  // SEQUANT2_WOLFRAM_HPP
