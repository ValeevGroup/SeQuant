//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT2_LATEX_HPP
#define SEQUANT2_LATEX_HPP

#include <type_traits>

namespace sequant2 {

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>>, std::wstring> to_latex(T&& t) {
  std::wstring result = L"{";
  result += std::to_wstring(t) + L"}";
  return result;
}

template <typename T>
std::wstring to_latex(const std::complex<T>& t) {
  std::wstring result = L"{";
  result += std::to_wstring(t.real());
  if (t.imag() > 0) {
    result += L" + i " + std::to_wstring(t.imag());
  } else if (t.imag() < 0)
    result += L" - i " + std::to_wstring(-t.imag());
  result += L"}";
  return result;
}

}

#endif //SEQUANT2_LATEX_HPP
