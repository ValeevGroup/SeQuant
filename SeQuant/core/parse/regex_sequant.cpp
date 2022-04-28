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

std::wstring regex_patterns::pure_index() {
    return L"[ia][⁺⁻]?_?\\d+";
}

std::wstring regex_patterns::pure_indices() {
  return pure_index() + zero_or_more_non_greedy(L"," + pure_index());
}

std::wstring regex_patterns::proto_indices() {
  return L"<" + pure_indices() + L">";
}

std::wstring regex_patterns::proto_indices_capture() {
  return L"<" + capture(pure_indices()) + L">";
}

std::wstring regex_patterns::index() {
  // index with optional proto-indices
  return pure_index() + capture_not(proto_indices()) + L"?";
}

std::wstring regex_patterns::index_capture() {
  return capture(pure_index()) + capture_not(proto_indices_capture()) + L"?";
}

std::wstring_view regex_patterns::indices() {
  static const std::wstring idxs = index().data()
                                   + zero_or_more_non_greedy(L","s
                                                             +index().data());
  return idxs;
}

std::wstring_view regex_patterns::bra_expanded() {
  static const std::wstring bra = L"_\\{"s
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
  static const std::wstring bra = L"_\\{"s
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

std::wstring_view regex_patterns::abs_real_frac() {
  static std::wstring frac =
      std::wstring{} + abs_real_num().data()
      + capture_not(std::wstring{}
                    + LR"(\/)"
                    + capture(abs_real_num()).data())
      + L"?";

  return frac; // guranteed numerator and optional denominator
}

std::wstring_view regex_patterns::tensor_expanded() {
  static const std::wstring tensor = capture(L"\\w+")
                                     + look_ahead(L"\\S*?"s
                                     + bra_expanded_capture().data())
                                 + look_ahead(L"\\S*?"s
                                              + ket_expanded_capture().data())
                                 + this_or_that(L""s
                                                    + bra_expanded().data()
                                                    + ket_expanded().data(),
                                                L""s
                                                    + ket_expanded().data()
                                                    + bra_expanded().data())
                                 + capture_not(L":" + capture(L"A|S|N"))
                                 + L"?";

  return tensor;
}

std::wstring_view regex_patterns::tensor_terse() {
  static const std::wstring tensor = capture(L"\\w+")
                                     + L"\\{"s
                                     + capture(indices())
                                     + L";"
                                     + capture(indices())
                                     + L"\\}"
                                     + capture_not(L":" + capture(L"A|S|N"))
                                     + L"?";
  return tensor;
}

} // namespace