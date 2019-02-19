//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT2_LATEX_HPP
#define SEQUANT2_LATEX_HPP

#include <type_traits>
#include "wstring.hpp"

namespace sequant2 {

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     !std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_latex(T&& t) {
  std::wstring result = L"{";
  using ::sequant2::to_wstring;
  result += to_wstring(t) + L"}";
  return result;
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_latex(T&& t) {
  using Real = std::decay_t<T>;
  static const auto eps_sqrt = std::sqrt(std::numeric_limits<Real>::epsilon());

  std::wstring result = L"{";
  using ::sequant2::to_wstring;

  if (std::floor(t) == t)  // exact integer
    result += to_wstring(t) + L"}";
  else {  // TODO detect rationals
    const auto inv_t = Real(1) / t;
    const auto round_inv_t = round(inv_t);
    if (std::abs(round_inv_t - inv_t) < eps_sqrt) {  // exact inverse of an
                                                     // integer
      long denom = long(round_inv_t);
      using namespace std::literals;
      result += (std::signbit(t) ? L"-"s : L""s) + L"\\frac{1}{"s +
                std::to_wstring(abs(denom)) + L"}}"s;
    } else
      result += to_wstring(t) + L"}";
  }
  return result;
}

template <typename T>
std::wstring to_latex(const std::complex<T>& t) {
  std::wstring result = L"{";
  result += to_latex(t.real());
  if (t.imag() > 0) {
    result += L" + i " + to_latex(t.imag());
  } else if (t.imag() < 0)
    result += L" - i " + to_latex(-t.imag());
  result += L"}";
  return result;
}

}  // namespace sequant2

#endif //SEQUANT2_LATEX_HPP
