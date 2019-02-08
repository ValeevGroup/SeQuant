//
// Created by Eduard Valeyev on 2/8/19.
//

#ifndef SEQUANT2_WOLFRAM_HPP
#define SEQUANT2_WOLFRAM_HPP

#include <type_traits>

namespace sequant2 {

template<typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, std::wstring> to_wolfram(T &&t) {
  return std::to_wstring(t);
}

template<typename T>
std::wstring to_wolfram(const std::complex<T> &t) {
  if (t.imag() == 0)
    return to_wolfram(t.real());
  else
    return std::wstring(L"Complex[") + std::to_wstring(t.real()) + L"," + std::to_wstring(t.imag()) + L"]";
}

}

#endif //SEQUANT2_WOLFRAM_HPP
