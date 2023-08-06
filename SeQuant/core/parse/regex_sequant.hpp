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

  constexpr static std::wstring_view abs_real_num() {
    return LR"=(((?:\d+(?:\.\d*)?)|(?:\.\d+)))=";
  }

  ///
  /// \brief A label cannot be a space, punctuation, or a control character.
  ///        Also at the first character cannot be a digit.
  ///
  constexpr static std::wstring_view label() {
    return L"[^[:digit:][:punct:][:space:][:cntrl:]]"
           L"(?:[^[:punct:][:space:][:cntrl:]])*?";
  }

  /// Capture group 1: numerator
  /// If captured, group 2: denominator
  static std::wstring_view abs_real_frac();

  static std::wstring pure_index();

  static std::wstring pure_index_capture();

  static std::wstring index_capture();

 private:
  static std::wstring capture(std::wstring_view pat);

  static std::wstring capture_not(std::wstring_view pat);

  static std::wstring zero_or_more_non_greedy(std::wstring_view pat);

  static std::wstring this_or_that(std::wstring_view pat1,
                                   std::wstring_view pat2);

  static std::wstring look_ahead(std::wstring_view pat);

  static std::wstring pure_indices();

  static std::wstring proto_indices();

  static std::wstring proto_indices_capture();

  static std::wstring index();

  static std::wstring_view indices();

  static std::wstring_view bra_expanded();

  static std::wstring_view ket_expanded();

  static std::wstring_view bra_expanded_capture();

  static std::wstring_view ket_expanded_capture();
};

}  // namespace sequant::parse

#endif  // SEQUANT_PARSE_REGEX_SEQUANT_HPP
