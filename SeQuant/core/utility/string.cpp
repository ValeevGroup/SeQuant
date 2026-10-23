#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

#include <utf8/checked.h>
#include <utf8/unchecked.h>

namespace sequant {

bool is_valid_utf8(std::string_view str) { return utf8::is_valid(str); }

std::u32string to_utf32(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  std::u32string result;
  utf8::unchecked::utf8to32(str.begin(), str.end(), std::back_inserter(result));
  return result;
}

std::string to_utf8(std::u32string_view str) {
  std::string result;
  utf8::utf32to8(str.begin(), str.end(), std::back_inserter(result));
  return result;
}

std::string to_utf8(char32_t code_point) {
  return to_utf8(std::u32string_view(&code_point, 1));
}

std::string_view first_code_point(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  if (str.empty()) return str;
  auto it = str.begin();
  utf8::unchecked::next(it);
  return str.substr(0, std::distance(str.begin(), it));
}

std::string_view last_code_point(std::string_view str) {
  SEQUANT_ASSERT(is_valid_utf8(str));
  if (str.empty()) return str;
  auto it = str.end();
  utf8::unchecked::prior(it);
  return str.substr(std::distance(str.begin(), it));
}

std::string utf8_regex_pattern(std::string_view pattern) {
  SEQUANT_ASSERT(is_valid_utf8(pattern));

  // matches a single UTF-8 encoded code point that is not a line terminator
  // and not part of the given (ASCII) bracket expression contents
  auto any_code_point_except = [](std::string_view ascii_exclusions) {
    return "(?:[^\\n\\r\\x80-\\xFF" + std::string(ascii_exclusions) +
           "]|[\\xC0-\\xF7][\\x80-\\xBF]+)";
  };
  auto is_ascii = [](char32_t c) { return c < 0x80; };
  // non-ASCII characters are grouped such that quantifiers apply to the
  // entire (multi-byte) character
  auto literal = [](char32_t c) { return "(?:" + to_utf8(c) + ")"; };

  const std::u32string cps = to_utf32(pattern);
  std::string result;
  result.reserve(pattern.size());

  for (std::size_t i = 0; i < cps.size(); ++i) {
    const char32_t c = cps[i];

    if (c == U'\\' && i + 1 < cps.size()) {
      const char32_t escaped = cps[++i];
      if (!is_ascii(escaped)) {
        result += literal(escaped);
      } else if (escaped == U'D' || escaped == U'W' || escaped == U'S') {
        // negated character class escapes must match entire code points
        result += any_code_point_except(
            std::string("\\") + static_cast<char>(escaped + ('a' - 'A')));
      } else {
        result += '\\';
        result += static_cast<char>(escaped);
      }
    } else if (c == U'.') {
      result += any_code_point_except("");
    } else if (c == U'[') {
      std::size_t j = i + 1;
      const bool negated = j < cps.size() && cps[j] == U'^';
      if (negated) ++j;

      std::string ascii_members;
      std::vector<std::string> non_ascii_members;
      bool terminated = false;
      for (bool first = true; j < cps.size(); ++j, first = false) {
        const char32_t m = cps[j];
        if (m == U']' && !first) {
          terminated = true;
          break;
        }
        if (m == U'\\' && j + 1 < cps.size()) {
          const char32_t escaped = cps[++j];
          if (!is_ascii(escaped)) {
            non_ascii_members.push_back(to_utf8(escaped));
          } else if (escaped == U'D' || escaped == U'W' || escaped == U'S') {
            throw Exception(
                "utf8_regex_pattern: negated character class escapes inside "
                "bracket expressions are not supported");
          } else {
            ascii_members += '\\';
            ascii_members += static_cast<char>(escaped);
          }
          continue;
        }
        const bool is_range =
            j + 2 < cps.size() && cps[j + 1] == U'-' && cps[j + 2] != U']';
        if (is_range) {
          if (!is_ascii(m) || !is_ascii(cps[j + 2]))
            throw Exception(
                "utf8_regex_pattern: ranges of non-ASCII characters in "
                "bracket expressions are not supported");
          ascii_members += static_cast<char>(m);
          ascii_members += '-';
          ascii_members += static_cast<char>(cps[j + 2]);
          j += 2;
        } else if (is_ascii(m)) {
          ascii_members += static_cast<char>(m);
        } else {
          non_ascii_members.push_back(to_utf8(m));
        }
      }

      if (!terminated) {
        // malformed pattern: pass on as-is to let the regex engine complain
        result += to_utf8(std::u32string_view(cps).substr(i));
        break;
      }
      i = j;

      std::string alternatives;
      for (const auto &member : non_ascii_members) {
        if (!alternatives.empty()) alternatives += '|';
        alternatives += member;
      }

      if (negated) {
        result += alternatives.empty() ? "(?:" : "(?:(?!" + alternatives + ")";
        result += any_code_point_except(ascii_members) + ")";
      } else if (non_ascii_members.empty()) {
        result += "[" + ascii_members + "]";
      } else if (ascii_members.empty()) {
        result += "(?:" + alternatives + ")";
      } else {
        result += "(?:[" + ascii_members + "]|" + alternatives + ")";
      }
    } else if (!is_ascii(c)) {
      result += literal(c);
    } else {
      result += static_cast<char>(c);
    }
  }

  return result;
}

}  // namespace sequant
