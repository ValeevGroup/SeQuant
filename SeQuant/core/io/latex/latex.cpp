//
// Created by Eduard Valeyev on 7/18/23.
//

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <algorithm>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace sequant::io::latex {

std::string to_string(const ExprPtr& expr) { return expr.to_latex(); }

std::string to_string(const rational& t) {
  // n.b. skip enclosing braces to make Constant::to_latex to produce same
  // output as before std::string result = "{";
  std::string result;
  if (denominator(t) == 1)
    result += "{" + ::sequant::to_string(numerator(t)) + "}";
  else {
    const auto num = numerator(t);
    // n.b. extra braces around \frac and use of to_string instead of to_latex
    // to avoid extra braces around args to \frac
    if (num > 0) {
      result += "{\\frac{" + ::sequant::to_string(numerator(t)) + "}{" +
                ::sequant::to_string(denominator(t)) + "}}";
    } else if (num < 0) {
      result += "{-\\frac{" + ::sequant::to_string(decltype(num)(-num)) + "}{" +
                ::sequant::to_string(denominator(t)) + "}}";
    } else
      result += "0";
  }
  // n.b.
  // result += "}";
  return result;
}

std::string to_string(const Power& power) {
  std::string result;
  const auto& base = power.base();
  const auto& exp = power.exponent();

  // special case when base is a unit fraction 1/k with k > 1
  const auto base_nr_1 = base->is<Constant>() &&
                         base->as<Constant>().value().imag() == 0 &&
                         numerator(base->as<Constant>().value().real()) == 1 &&
                         denominator(base->as<Constant>().value().real()) > 1;

  if (exp == 1) {
    result = base->to_latex();
  } else if (base_nr_1) {
    // render (1/k)^n as \frac{1}{k^n}
    result =
        "{\\frac{1}{" +
        ::sequant::to_string(denominator(base->as<Constant>().value().real())) +
        "^" + to_string(exp) + "}}";
  } else {
    result = base->to_latex() + "^" + to_string(exp);
  }
  if (power.conjugated()) result = "{" + result + "^*}";
  return result;
}

std::string to_string(const ResultExpr& expr) {
  std::stringstream stream;
  if (expr.produces_tensor()) {
    stream << to_string(expr.result_as_tensor());
  } else {
    stream << to_string(expr.result_as_variable());
  }

  stream << " = " << to_string(expr.expression());

  return stream.str();
}

namespace {

// n.b. the implementations operate on decoded code points (UTF-32)

std::u32string greek_characters_to_string_impl(std::u32string_view str) {
  // lower-case greek characters in the order of their appearance in Unicode
  // chart https://www.unicode.org/charts/PDF/U0370.pdf
  const std::vector<std::u32string> lc = {
      U"\\alpha",   U"\\beta", U"\\gamma",    U"\\delta", U"\\epsilon",
      U"\\zeta",    U"\\eta",  U"\\theta",    U"\\iota",  U"\\kappa",
      U"\\lambda",  U"\\mu",   U"\\nu",       U"\\xi",    U"o",
      U"\\pi",      U"\\rho",  U"\\varsigma", U"\\sigma", U"\\tau",
      U"\\upsilon", U"\\phi",  U"\\chi",      U"\\psi",   U"\\omega"};
  const std::vector<std::u32string> uc = {
      U"A",         U"B",     U"\\Gamma", U"\\Delta", U"E",
      U"Z",         U"H",     U"\\Theta", U"I",       U"K",
      U"\\Lambda",  U"M",     U"N",       U"\\Xi",    U"O",
      U"\\Pi",      U"P",     U"",        U"\\Sigma", U"T",
      U"\\Upsilon", U"\\Phi", U"X",       U"\\Psi",   U"\\Omega"};
  auto is_lc = [](char32_t ch) {
    return ch >= 0x3B1     // alpha
           && ch <= 0x3C9  // omega
        ;
  };
  auto is_uc = [](char32_t ch) {
    return ch >= 0x391     // Alpha
           && ch <= 0x3A9  // Omega
           && ch != 0x3A2  // gap in the table
        ;
  };

  std::u32string result;
  result.reserve(str.size());
  /// TODO iterate over graphemes
  for (const char32_t ch : str) {
    if (is_lc(ch)) {
      const auto& lc_str = lc[ch - 0x3B1];
      SEQUANT_ASSERT(lc_str.size() > 0);
      result += lc_str;
    } else if (is_uc(ch)) {
      const auto& uc_str = uc[ch - 0x391];
      SEQUANT_ASSERT(uc_str.size() > 0);
      result += uc_str;
    } else {  // pass through all other characters
      result.push_back(ch);
    }
  }

  return result;
}

std::u32string diactrics_to_string_impl(std::u32string_view str) {
  // precomposed characters with diacritics
  static const container::map<char32_t, std::u32string_view> precomposed = {
      // lower-case characters with tilde
      {U'ã', U"\\tilde{a}"},
      {U'ẽ', U"\\tilde{e}"},
      {U'ñ', U"\\tilde{n}"},
      {U'õ', U"\\tilde{o}"},
      {U'ũ', U"\\tilde{u}"},
      {U'ṽ', U"\\tilde{v}"},
      // upper-case characters with tilde
      {U'Ã', U"\\tilde{A}"},
      {U'Ẽ', U"\\tilde{E}"},
      {U'Ñ', U"\\tilde{N}"},
      {U'Õ', U"\\tilde{O}"},
      {U'Ũ', U"\\tilde{U}"},
      {U'Ṽ', U"\\tilde{V}"},
      // lower-case characters with circumflex
      {U'â', U"\\hat{a}"},
      {U'ê', U"\\hat{e}"},
      {U'î', U"\\hat{i}"},
      {U'ô', U"\\hat{o}"},
      // upper-case characters with circumflex
      {U'Â', U"\\hat{A}"},
      {U'Ê', U"\\hat{E}"},
      {U'Î', U"\\hat{I}"},
      {U'Ô', U"\\hat{O}"},
      {U'Ŝ', U"\\hat{S}"}};

  // Combining diacritics:
  // https://www.ncbi.nlm.nih.gov/staff/beck/charents/accents.html
  auto combining_prefix = [](char32_t ch) -> std::u32string_view {
    switch (ch) {
      case 0x303:  // tilde
        return U"\\tilde{";
      case 0x301:  // acute
        return U"\\acute{";
      case 0x300:  // grave
        return U"\\grave{";
      case 0x30C:  // caron
        return U"\\check{";
      case 0x302:  // hat/circumflex
        return U"\\hat{";
      default:
        return {};
    }
  };

  std::u32string result;
  result.reserve(str.size());
  /// TODO iterate over graphemes
  for (auto it = str.begin(); it != str.end(); ++it) {
    const char32_t ch = *it;

    if (it + 1 != str.end()) {
      if (const auto prefix = combining_prefix(*(it + 1)); !prefix.empty()) {
        result += prefix;
        result.push_back(ch);
        result += U"}";
        ++it;
        continue;
      }
    }

    if (auto pc_it = precomposed.find(ch); pc_it != precomposed.end()) {
      result += pc_it->second;
    } else {  // pass through all other characters
      result.push_back(ch);
    }
  }

  return result;
}

bool is_ascii(std::string_view str) {
  return std::ranges::all_of(
      str, [](char c) { return static_cast<unsigned char>(c) <= 0x7F; });
}

}  // namespace

std::string greek_characters_to_string(std::string_view str) {
  // fast path: no non-ASCII characters to replace
  if (is_ascii(str)) return std::string(str);
  return to_utf8(greek_characters_to_string_impl(to_utf32(str)));
}

std::string diactrics_to_string(std::string_view str) {
  // fast path: no non-ASCII characters to replace
  if (is_ascii(str)) return std::string(str);
  return to_utf8(diactrics_to_string_impl(to_utf32(str)));
}

std::string utf_to_string(std::string_view str) {
  // fast path: no non-ASCII characters to replace
  if (is_ascii(str)) return std::string(str);
  // replace diacritics first since it relies on combining characters following
  // their base character
  return to_utf8(
      greek_characters_to_string_impl(diactrics_to_string_impl(to_utf32(str))));
}

}  // namespace sequant::io::latex
