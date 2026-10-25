//
// Created by Eduard Valeyev on 3/30/18.
//

#ifndef SEQUANT_CORE_IO_LATEX_LATEX_HPP
#define SEQUANT_CORE_IO_LATEX_LATEX_HPP

#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>

namespace sequant {
class ExprPtr;
class Power;
class ResultExpr;
}  // namespace sequant

namespace sequant::io::latex {

/// @returns an empty string if @p expr is null.
std::string to_string(const ExprPtr& expr);

template <typename T>
concept has_to_latex_member = requires(const T& t) { t.to_latex(); };

template <typename T>
concept pointer_can_call_to_latex = requires(const T& t) { t->to_latex(); };

template <typename T>
  requires(has_to_latex_member<T>)
std::string to_string(const T& t) {
  return t.to_latex();
}

template <typename T>
  requires(!has_to_latex_member<T> && pointer_can_call_to_latex<T>)
std::string to_string(const T& t) {
  return t->to_latex();
}

template <typename T>
  requires(!has_to_latex_member<T> && !pointer_can_call_to_latex<T> &&
           requires(const T& t) { t._to_latex(); })
std::string to_string(const T& t) {
  return t._to_latex();
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     !std::is_floating_point_v<std::decay_t<T>>,
                 std::string>
to_string(T&& t) {
  std::string result = "{";
  using ::sequant::to_string;
  result += to_string(t) + "}";
  return result;
}

template <typename T>
std::enable_if_t<std::is_arithmetic_v<std::decay_t<T>> &&
                     std::is_floating_point_v<std::decay_t<T>> &&
                     !std::is_same_v<std::decay_t<T>, rational>,
                 std::string>
to_string(T&& t) {
  using Real = std::decay_t<T>;
  static const auto eps_sqrt = std::sqrt(std::numeric_limits<Real>::epsilon());

  std::string result = "{";
  using ::sequant::to_string;

  const long round_t = std::lround(t);
  if (std::abs(round_t - t) < eps_sqrt)  // exact integer
    result += to_string(round_t) + "}";
  else {
    const auto inv_t = Real(1) / t;
    const long round_inv_t = std::lround(inv_t);
    if (std::abs(round_inv_t - inv_t) < eps_sqrt) {  // exact inverse of an
                                                     // integer
      long denom = round_inv_t;
      using namespace std::literals;
      result += (std::signbit(t) ? "-"s : ""s) + "\\frac{1}{"s +
                std::to_string(std::abs(denom)) + "}}"s;
    } else
      result += to_string(t) + "}";
  }
  return result;
}

template <typename T>
std::string to_string(const std::complex<T>& t) {
  std::string result = "{";
  result += to_string(t.real());
  if (t.imag() > 0) {
    result += " + i " + to_string(t.imag());
  } else if (t.imag() < 0)
    result += " - i " + to_string(-t.imag());
  result += "}";
  return result;
}

std::string to_string(const rational& num);

std::string to_string(const Power& power);

std::string to_string(const ResultExpr& expr);

/// replaces certain greek characters in a string with their (math-mode) LaTeX
/// equivalents
/// @param str UTF-8 encoded input string
/// @return @p str with greek characters contained within Unicode ranges
/// `[0x3B1,0x3C9]` and `[0391,03A9]` replaced with their LaTeX equivalents
std::string greek_characters_to_string(std::string_view str);

/// replaces certain diactric marks with their (math-mode) LaTeX equivalents
/// @param str UTF-8 encoded input string
/// @return @p str with some diactrics replaced with their LaTeX equivalents
/// @note only some combined Unicode characters are currently supported
std::string diactrics_to_string(std::string_view str);

/// replaces certain Unicode characters with their (math-mode) LaTeX equivalents
/// @param str UTF-8 encoded input string
/// @return @p str with some Unicode characters replaced by their LaTeX
/// equivalents
/// @note only some combined Unicode characters are currently supported
/// @internal useful resources
/// - https://milde.users.sourceforge.net/LUCR/Math/unimathsymbols.pdf
/// - https://www.unicode.org/charts/PDF/U1D400.pdf
std::string utf_to_string(std::string_view str);

}  // namespace sequant::io::latex

#endif  // SEQUANT_CORE_IO_LATEX_LATEX_HPP
