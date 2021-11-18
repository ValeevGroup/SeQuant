//
// Created by Bimal Gaudel on 6/15/21.
//

#ifndef SEQUANT_PARSE_REGEX_SEQUANT_HPP
#define SEQUANT_PARSE_REGEX_SEQUANT_HPP

#include <string_view>

namespace sequant::parse {

class regex_patterns {
 public:
  static std::wstring_view tensor_expanded();

  static std::wstring_view tensor_terse();

  static std::wstring_view fraction();

 private:
  static std::wstring capture(std::wstring_view pat);

  static std::wstring capture_not(std::wstring_view pat);

  static std::wstring zero_or_more_non_greedy(std::wstring_view pat);

  static std::wstring this_or_that(std::wstring_view pat1,
                                   std::wstring_view pat2);

  static std::wstring look_ahead(std::wstring_view pat);

  static std::wstring_view index();

  static std::wstring_view indices();

  static std::wstring_view bra_expanded();

  static std::wstring_view ket_expanded();

  static std::wstring_view bra_expanded_capture();

  static std::wstring_view ket_expanded_capture();
};

} // namespace

#endif  // SEQUANT_PARSE_REGEX_SEQUANT_HPP
