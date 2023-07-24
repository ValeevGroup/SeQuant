//
// Created by Eduard Valeyev on 7/18/23.
//

#include "SeQuant/core/latex.ipp"

namespace sequant::detail {

#define SQ_IMPL1(CHAR)                                             \
  template std::basic_string<CHAR> greek_characters_to_latex_impl< \
      CHAR, std::char_traits<CHAR>, std::allocator<CHAR>>(         \
      std::basic_string_view<CHAR>);

SQ_IMPL1(char);
SQ_IMPL1(wchar_t);
#if __cplusplus >= 202002L
SQ_IMPL1(char8_t);
SQ_IMPL1(char16_t);
SQ_IMPL1(char32_t);
#endif

#define SQ_IMPL2(CHAR)                                      \
  template std::basic_string<CHAR> diactrics_to_latex_impl< \
      CHAR, std::char_traits<CHAR>, std::allocator<CHAR>>(  \
      std::basic_string_view<CHAR>);

SQ_IMPL2(char);
SQ_IMPL2(wchar_t);
#if __cplusplus >= 202002L
SQ_IMPL2(char8_t);
SQ_IMPL2(char16_t);
SQ_IMPL2(char32_t);
#endif

}  // namespace sequant::detail
