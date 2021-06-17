//
// Created by Bimal Gaudel on 6/15/21.
//

#include "regex_sequant.hpp"

#include <string>

namespace sequant::parse {

using namespace std::string_literals;

std::wstring regex_patterns::capture(std::wstring_view pat) {
  return L"("s + pat.data() + L")";
}

std::wstring regex_patterns::capture_not(std::wstring_view pat) {
  return L"(?:"s + pat.data() + L")";
}

std::wstring regex_patterns::zero_or_more_non_greedy(std::wstring_view pat) {
  return capture_not(pat) + L"*?";
}

std::wstring regex_patterns::this_or_that(std::wstring_view pat1,
                                          std::wstring_view pat2) {
  return capture_not(capture_not(pat1)
                         + L"|"
                         + capture_not(pat2));
}

std::wstring regex_patterns::look_ahead(std::wstring_view pat) {
  return L"(?="s + pat.data() + L")";
}

std::wstring_view regex_patterns::index() {
  static const std::wstring idx = L"[ia]_?\\d+";
  return idx;
}

std::wstring_view regex_patterns::indices() {
  static const std::wstring idxs = index().data()
                                   + zero_or_more_non_greedy(L","s
                                                             +index().data());
  return idxs;
}

std::wstring_view regex_patterns::bra_expanded() {
  static const std::wstring bra = L"\\_\\{"s
                                  + indices().data()
                                  + L"\\}";
  return bra;
}

std::wstring_view regex_patterns::ket_expanded() {
  static const std::wstring ket = L"\\^\\{"s
                                  + indices().data()
                                  + L"\\}";
  return ket;
}

std::wstring_view regex_patterns::bra_expanded_capture() {
  static const std::wstring bra = L"\\_\\{"s
                                  + capture(indices())
                                  + L"\\}";
  return bra;
}

std::wstring_view regex_patterns::ket_expanded_capture() {
  static const std::wstring ket = L"\\^\\{"s
                                  + capture(indices())
                                  + L"\\}";
  return ket;
}

std::wstring_view regex_patterns::fraction() {
  static const std::wstring frac =
      LR"=((\d+(?:\.\d*)?|\.\d+)(?:\/(\d+(?:\.\d*)?|\.\d+))?)="s;
  return frac;
}

std::wstring_view regex_patterns::tensor_expanded() {
  static const std::wstring bk = look_ahead(L"\\S*?"s
                                            + bra_expanded_capture().data())
                                 + look_ahead(L"\\S*?"s
                                              + ket_expanded_capture().data())
                                 + this_or_that(L""s
                                                    + bra_expanded().data()
                                                    + ket_expanded().data(),
                                                L""s
                                                    + ket_expanded().data()
                                                    + bra_expanded().data());

  static const std::wstring tensor = capture(L"\\w+") + bk;
  return tensor;
}

std::wstring_view regex_patterns::tensor_terse() {
  static const std::wstring tensor = capture(L"\\w+")
                                     + L"\\{"s
                                     + capture(indices())
                                     + L";"
                                     + capture(indices())
                                     + L"\\}";
  return tensor;
}

} // namespace