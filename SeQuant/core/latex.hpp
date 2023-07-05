//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_LATEX_HPP
#define SEQUANT_LATEX_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/wstring.hpp>
#include <type_traits>

#include <boost/rational.hpp>

namespace sequant {

template <typename T>
std::enable_if_t<meta::has_memfn_to_latex_v<std::decay_t<T>>, std::wstring>
to_latex(T&& t) {
  return t.to_latex();
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     !std::is_floating_point_v<std::decay_t<T>>,
                 std::wstring>
to_latex(T&& t) {
  std::wstring result = L"{";
  using ::sequant::to_wstring;
  result += to_wstring(t) + L"}";
  return result;
}
/// for multiprecison integer
template <typename T>
std::enable_if_t<
    std::is_same_v<std::decay_t<T>, sequant::mp_int_type>,
    std::wstring>
to_latex(T&& t) {
  std::wstring result = L"{";
  using ::sequant::to_wstring;
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
  using ::sequant::to_wstring;

  const long round_t = std::lround(t);
  if (std::abs(round_t - t) < eps_sqrt)  // exact integer
    result += to_wstring(round_t) + L"}";
  else {  // TODO detect rationals
    const auto inv_t = Real(1) / t;
    const long round_inv_t = std::lround(inv_t);
    if (std::abs(round_inv_t - inv_t) < eps_sqrt) {  // exact inverse of an
                                                     // integer
      long denom = round_inv_t;
      using namespace std::literals;
      result += (std::signbit(t) ? L"-"s : L""s) + L"\\frac{1}{"s +
                std::to_wstring(std::abs(denom)) + L"}}"s;
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

template <typename T>
std::wstring to_latex(const boost::rational<T>& t) {
  // n.b. skip enclosing braces to make Constant::to_latex to produce same
  // output as before std::wstring result = L"{";
  std::wstring result;
  if (t.denominator() == 1)
    result += to_latex(t.numerator());
  else {
    const auto num = t.numerator();
    // n.b. extra braces around \frac and use of to_wstring instead of to_latex
    // to avoid extra braces around args to \frac
    if (num > 0) {
      result += L"{\\frac{" + to_wstring(t.numerator()) + L"}{" +
                to_wstring(t.denominator()) + L"}}";
    } else if (num < 0) {
      result += L"{-\\frac{" + to_wstring(-num) + L"}{" +
                to_wstring(t.denominator()) + L"}}";
    } else
      result += L"0";
  }
  // n.b.
  // result += L"}";
  return result;
}

}  // namespace sequant

#endif  // SEQUANT_LATEX_HPP
