#ifndef SEQUANT_CORE_RATIONAL_H
#define SEQUANT_CORE_RATIONAL_H

#include <SeQuant/core/hash.hpp>
#include <SeQuant/core/wstring.hpp>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <memory>
#include <mutex>

namespace sequant {

using rational = boost::multiprecision::cpp_rational;

/// shorter=sweeter? sometimes
using ratio = rational;

/// the integer type supporting the rational
using intmax_t = rational::value_type;

}  // namespace sequant

namespace boost {

inline auto hash_value(const sequant::rational& i) {
  auto val = sequant::hash::value(numerator(i));
  sequant::hash::combine(val, denominator(i));
  return val;
}

}  // namespace boost

namespace sequant {

/// convert a floating-point number to a rational number to a given precision

/// @param t the floating-point number to convert
/// @param eps the target precision, i.e. the largest permitted value of
/// `T(result)-t`; the default is the `sqrt` of the machine epsilon for `T`;
/// set to 0 to obtain exact representation of `t` as a rational number
/// @param max_niter the maximum number of iterations to use to construct
/// inexact representation by refining the Stern-Brocot tree
/// (see https://mathworld.wolfram.com/Stern-BrocotTree.html)
/// @return the rational number within @p eps of @p t
template <typename T, typename = std::enable_if_t<std::is_floating_point_v<T>>>
inline rational to_rational(
    T t, T eps = std::sqrt(std::numeric_limits<T>::epsilon()),
    std::size_t max_niter = 1000) {
  if (std::isnan(t) || std::isinf(t)) {
    throw std::invalid_argument(
        "sequant::to_rational: cannot make a rational out of " +
        std::to_string(t));
  }
  // e.g.
  // https://gist.github.com/mikeando/7073d62385a34a61a6f7#file-main2-cpp-L42
  auto sbtree = [max_niter](T f, T tol) {
    using fraction = rational;
    using Int = rational::value_type;
    auto base = std::floor(f);
    auto base_Int = boost::numeric_cast<Int>(base);
    f -= base;
    if (f < tol) return fraction(base, 1);
    if (1 - tol < f) return fraction(base + 1, 1);

    fraction lower(0, 1);
    fraction upper(1, 1);

    std::size_t niter = 0;
    const auto f_plus_top_exact = fraction(f + tol);
    const auto f_minus_top_exact = fraction(f - tol);
    while (niter < max_niter) {
      fraction middle(numerator(lower) + numerator(upper),
                      denominator(lower) + denominator(upper));

      if (denominator(middle) * f_plus_top_exact < numerator(middle)) {
        upper = middle;
      } else if (numerator(middle) < denominator(middle) * f_minus_top_exact) {
        lower = middle;
      } else {
        return fraction(denominator(middle) * base_Int + numerator(middle),
                        denominator(middle));
      }
      ++niter;
    }
    throw std::invalid_argument(
        "sequant::rationalize: could not rationalize " + std::to_string(f) +
        " to eps=" + std::to_string(tol) + " in " + std::to_string(max_niter) +
        " iterations");  // unreachable
  };
  if (eps == 0)  // exact conversion ... rarely what's desired
    return rational(t);
  else
    return sbtree(t, eps);
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

inline std::string to_string(const rational& i) {
  return boost::lexical_cast<std::string>(i);
}

template <typename Backend>
inline std::string to_string(const boost::multiprecision::number<Backend>& i) {
  return boost::lexical_cast<std::string>(i);
}

inline std::wstring to_wstring(const rational& i) {
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

inline std::wstring to_latex(const rational& t) {
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

inline std::wstring to_wolfram(const rational& t) {
  using ::sequant::to_wstring;
  if (denominator(t) == 1) {
    // n.b. use to_string to skip extra braces so that output agrees with code
    // that used scalars return to_wolfram(t.numerator());
    return to_wstring(numerator(t));
  } else
    return std::wstring(L"Rational[") + to_wolfram(numerator(t)) + L"," +
           to_wolfram(denominator(t)) + L"]";
}

/// @brief Returns the absolute value of a rational number
/// @param r a sequant::rational object
inline rational abs(const rational& r) { return numerator(r) < 0 ? -r : r; }

}  // namespace sequant

#endif  // SEQUANT_CORE_RATIONAL_H
