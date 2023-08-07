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
  return capture_not(capture_not(pat1) + L"|" + capture_not(pat2));
}

std::wstring regex_patterns::look_ahead(std::wstring_view pat) {
  return L"(?="s + pat.data() + L")";
}

std::wstring regex_patterns::pure_index() {
    return label().data() + L"_?\\d+"s;
}

std::wstring regex_patterns::pure_index_capture() {
    return capture(label()) + L"_?" + capture(L"\\d+");
}

std::wstring regex_patterns::pure_indices() {
  return pure_index() + zero_or_more_non_greedy(L"\\s*?,\\s*?" + pure_index());
}

std::wstring regex_patterns::proto_indices() {
  return L"<\\s*?" + pure_indices() + L"\\s*?>";
}

std::wstring regex_patterns::proto_indices_capture() {
  return L"<\\s*?" + capture(pure_indices()) + L"\\s*?>";
}

std::wstring regex_patterns::index() {
  // index with optional proto-indices
  return pure_index() + L"\\s*?" + capture_not(proto_indices()) + L"?";
}

std::wstring regex_patterns::index_capture() {
  return capture(pure_index()) + L"\\s*?" + capture_not(proto_indices_capture()) + L"?";
}

std::wstring_view regex_patterns::indices() {
  static const std::wstring idxs =
      index() + zero_or_more_non_greedy(L"\\s*?,\\s*?"s + index());
  return idxs;
}

std::wstring_view regex_patterns::bra_expanded() {
  static const std::wstring bra = L"_\\s*?\\{\\s*?"s + indices().data() + L"\\s*?\\}";
  return bra;
}

std::wstring_view regex_patterns::ket_expanded() {
  static const std::wstring ket = L"\\^\\s*?\\{\\s*?"s + indices().data() + L"\\s*?\\}";
  return ket;
}

std::wstring_view regex_patterns::bra_expanded_capture() {
  static const std::wstring bra = L"_\\s*?\\{\\s*?"s + capture(indices()) + L"\\s*?\\}";
  return bra;
}

std::wstring_view regex_patterns::ket_expanded_capture() {
  static const std::wstring ket = L"\\^\\s*?\\{\\s*?"s + capture(indices()) + L"\\s*?\\}";
  return ket;
}

std::wstring_view regex_patterns::abs_real_frac() {
  static std::wstring frac =
      abs_real_num().data() +
      capture_not(LR"(\s*?\/\s*?)"s + abs_real_num().data() ) +
      L"?";

  return frac;  // guaranteed numerator and optional denominator
}

std::wstring_view regex_patterns::tensor_expanded() {
  static const std::wstring tensor =
      capture(label()) +
      look_ahead(L".*?"s + bra_expanded_capture().data()) +
      look_ahead(L".*?"s + ket_expanded_capture().data()) +
      this_or_that(L"\\s*?"s + bra_expanded().data() + ket_expanded().data(),
                   L"\\s*?"s + ket_expanded().data() + bra_expanded().data()) +
      capture_not(L"\\s*?:\\s*?" + capture(L"A|S|N")) + L"?";

  return tensor;
}

std::wstring_view regex_patterns::tensor_terse() {
  static const std::wstring tensor =
      capture(label()) + L"\\s*?\\{\\s*?"s + capture(indices()) + L"\\s*?;\\s*?" +
      capture(indices()) + L"\\s*?\\}" + capture_not(L"\\s*?:\\s*?" + capture(L"A|S|N")) +
      L"?";
  return tensor;
}

}  // namespace sequant::parse
