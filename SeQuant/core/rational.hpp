#ifndef SEQUANT_CORE_RATIONAL_H
#define SEQUANT_CORE_RATIONAL_H

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/wstring.hpp>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

namespace boost {

inline auto hash_value(const boost::multiprecision::cpp_rational& i) {
  auto val = sequant::hash::value(numerator(i));
  sequant::hash::combine(val, denominator(i));
  return val;
}

}  // namespace boost

namespace sequant {

using rational = boost::multiprecision::cpp_rational;

/// shorter=sweeter? sometimes
using ratio = rational;

// clang-format off
///@note: boost::multiprecision::cpp_int only has limited support for constexpr
/// see: https://www.boost.org/doc/libs/1_82_0/libs/multiprecision/doc/html/boost_multiprecision/tut/lits.html
// clang-format on

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
inline rational to_rational(
    T t, T eps = std::sqrt(std::numeric_limits<T>::epsilon()),
    std::size_t max_niter = 1000) {
  if (std::isnan(t) || std::isinf(t)) {
    throw std::invalid_argument(
        "sequant::to_rational: cannot make a rational out of " +
        std::to_string(t));
  }
  return rational(t);
}

template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
rational to_rational(T t) {
  return rational{t};
}

template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T> ||
                                                  std::is_integral_v<T>>>
inline ratio to_ratio(T t) {
  return to_rational(t);
}

inline std::string to_string(const boost::multiprecision::cpp_rational& i) {
  return boost::lexical_cast<std::string>(i);
}

template <typename Backend>
inline std::string to_string(const boost::multiprecision::number<Backend>& i) {
  return boost::lexical_cast<std::string>(i);
}

inline std::wstring to_wstring(const boost::multiprecision::cpp_rational& i) {
  return ::sequant::to_wstring(boost::lexical_cast<std::string>(i));
}

template <typename Backend>
inline std::wstring to_wstring(
    const boost::multiprecision::number<Backend>& i) {
  return ::sequant::to_wstring(boost::lexical_cast<std::string>(i));
}

template <typename Backend>
inline std::wstring to_latex(const boost::multiprecision::number<Backend>& t) {
  std::wstring result = L"{";
  using ::sequant::to_wstring;
  result += to_wstring(t) + L"}";
  return result;
}

inline std::wstring to_latex(const boost::multiprecision::cpp_rational& t) {
  // n.b. skip enclosing braces to make Constant::to_latex to produce same
  // output as before std::wstring result = L"{";
  std::wstring result;
  if (denominator(t) == 1)
    result += to_latex(numerator(t));
  else {
    const auto num = numerator(t);
    // n.b. extra braces around \frac and use of to_wstring instead of to_latex
    // to avoid extra braces around args to \frac
    if (num > 0) {
      result += L"{\\frac{" + to_wstring(numerator(t)) + L"}{" +
                to_wstring(denominator(t)) + L"}}";
    } else if (num < 0) {
      result += L"{-\\frac{" + to_wstring(-num) + L"}{" +
                to_wstring(denominator(t)) + L"}}";
    } else
      result += L"0";
  }
  // n.b.
  // result += L"}";
  return result;
}

template <typename Backend>
inline std::wstring to_wolfram(
    const boost::multiprecision::number<Backend>& t) {
  return ::sequant::to_wstring(t);
}

inline std::wstring to_wolfram(const boost::multiprecision::cpp_rational& t) {
  using ::sequant::to_wstring;
  if (denominator(t) == 1) {
    // n.b. use to_string to skip extra braces so that output agrees with code
    // that used scalars return to_wolfram(t.numerator());
    return to_wstring(numerator(t));
  } else
    return std::wstring(L"Rational[") + to_wolfram(numerator(t)) + L"," +
           to_wolfram(denominator(t)) + L"]";
}

}  // namespace sequant

#endif  // SEQUANT_CORE_RATIONAL_H
