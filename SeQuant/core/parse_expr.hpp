#ifndef SEQUANT_PARSE_EXPR_HPP
#define SEQUANT_PARSE_EXPR_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <regex>
#include <string>

/**
 * Create SeQuant expression from string input.
 *
 * @author: Bimal Gaudel
 * @version: 04 Dec, 2020
 */

namespace sequant {

namespace detail {

/**
 * Squash all whitespace to nil.
 */
std::wstring prune_space(std::wstring_view);

/**
 * Convert string to decimal and return as double.
 * Returns 0 if empty string is passed or conversion fails.
 */
double to_decimal(std::wstring_view);

/**
 * Return numerator / denomrator.
 */
double as_fraction(std::wstring_view num, std::wstring_view denom);

/**
 * Check if a string can be parsed to a term with success.
 */
bool validate_term(std::wstring_view expr);

/**
 * Check if a string can be parsed to a sum with success.
 */
bool validate_sum(std::wstring_view expr);

/**
 * Make a tensor out of a regex match object that matched for tensor with
 * captured sub-groups: Group 1: tensor label. Group 2: tensor bra. Group 3:
 * tensor ket.
 */
ExprPtr make_tensor(const std::wsmatch& match_obj, Symmetry tensor_sym);

ExprPtr parse_tensor_term(const std::wstring&, Symmetry tensor_sym);

ExprPtr parse_product_term(const std::wstring&, Symmetry tensor_sym);

ExprPtr parse_term(const std::wstring&, Symmetry tensor_sym);

ExprPtr parse_sum(const std::wstring&, Symmetry tensor_sym);

const auto expr_rgx_pat = []() {  // a map from expr components to regex pattern
                                  // to match them
  using namespace std::string_literals;
  auto rgx_capture = [](std::wstring_view x) {
    return L"("s + x.data() + L")";
  };

  auto rgx_capture_not = [](std::wstring_view x) {
    return L"(?:"s + x.data() + L")";
  };

  auto rgx_or_nc = [&](std::wstring_view x, std::wstring_view y) {
    return rgx_capture_not(x.data() + L"|"s + y.data());
  };

  auto rgx_zero_or_more_nc = [&](std::wstring_view x) {
    return rgx_capture_not(x.data()) + L"*"s;
  };

  auto rgx_zero_or_one_nc = [&](std::wstring_view x) {
    return rgx_capture_not(x.data()) + L"?"s;
  };

  auto rgx_one_or_more_nc = [&](std::wstring_view x) {
    return rgx_capture_not(x.data()) + L"+"s;
  };

  auto rgx_look_ahead = [](std::wstring_view x) {
    return L"(?="s + x.data() + L")";
  };

  auto rgx_look_ahead_not = [](std::wstring_view x) {
    return L"(?!"s + x.data() + L")";
  };

  const std::wstring pi = L"[ia]_?[[:d:]]+";  // index pattern

  const std::wstring pii =  // (i1, i2, ..)
      LR"(\{)"              // open brace
      + rgx_capture_not(pi + rgx_zero_or_more_nc(L"," + pi)) +
      LR"(\})";  // close brace

  const std::wstring pb = L"_" + pii;      // bra
  const std::wstring pk = LR"(\^)" + pii;  // ket

  const std::wstring pt =         // tensor pattern
      rgx_capture(L"[[:w:]]+") +  // group 1 tensor label
      rgx_look_ahead(LR"([^[:s:]]*?)" +
                     rgx_capture(pb)) +  // group 2 tensor bra
      rgx_look_ahead(LR"([^[:s:]]*?)" +
                     rgx_capture(pk)) +  // group 3 tensor ket
      rgx_or_nc(pb + pk, pk + pb);

  const std::wstring ptt =  // tensor term: ie. just a tensor. no scaling
      rgx_or_nc(L"^" + pt + rgx_look_ahead_not(LR"(\*)" + pt),
                LR"(\+)" + pt + rgx_look_ahead_not(LR"(\*)" + pt));

  const std::wstring pd =  // decimal number
      LR"((?:\+|-)?)"
      L"[[:d:]]+"
      LR"((?:\.(?:[[:d:]]+)?)?)"
      L"|"
      LR"((?:\+|-)?\.[[:d:]]+)";

  const std::wstring pf =  // decimal fraction
      rgx_capture(pd) + rgx_zero_or_one_nc(LR"(\/)" + rgx_capture(pd));

  const std::wstring pmt =  // tensor following a multiplication sign
      LR"(\*)" + pt;

  const std::wstring pp =  // product
      L"(?:"               // product pattern begin
      //
      L"(?:"  // alternative 1
      L"^" +
      pt + rgx_one_or_more_nc(pmt) +
      L")"  // alternative 1 end
      //
      L"|"
      L"(?:"  // alternative 2
      L"^" +
      pf + rgx_one_or_more_nc(pmt) +
      L")"  // alternative 2 end
      //
      L"|"
      L"(?:"  // alternative 3
      L"-" +
      pt + rgx_zero_or_more_nc(pmt) +
      L")"  // alternative 3 end
      //
      L"|"
      //
      L"(?:"  // alternative 4
      LR"((?:\+|-))" +
      pt + rgx_one_or_more_nc(pmt) +
      L")"  // alternative 4 end
      //
      L"|"
      //
      L"(?:"  // alternative 5
      LR"((?:\+|-))" +
      pf + rgx_one_or_more_nc(pmt) +
      L")"  // alternative 5 end
      //
      L")";  // product pattern end

  const std::wstring pterm =  // a term: product or tensor
                              // if group 1 captured: it's a tensor term
                              // if group 2 captured: it's a product term
      rgx_or_nc(rgx_capture(ptt), rgx_capture(pp));

  const std::wstring psum =  // a sum: two or more terms
      L"^" + pterm + L"{2,}$";

  container::map<std::string, std::wstring> pat_mat{};

  pat_mat.emplace("index", pi);  // nothing captured
  pat_mat.emplace("indices",
                  pii);  // one or more index in paren. nothing captured
  pat_mat.emplace("tensor", pt);    // group 1 label. group 2 bra. group 3 ket
  pat_mat.emplace("bra", pb);       // nothing captured
  pat_mat.emplace("ket", pk);       // nothing captured
  pat_mat.emplace("decimal", pd);   // nothing captured
  pat_mat.emplace("fraction", pf);  // captured group 1 num. group 2 denom.
  pat_mat.emplace("tensor_term",
                  ptt);  // tensor pattern captured implicitly
  pat_mat.emplace("product_term",
                  pp);             // tensor patterns captured implicitly
  pat_mat.emplace("term", pterm);  // group 1 or group 2 captured. 1 implies
  // tensor term. 2 implies product term.
  pat_mat.emplace("sum", psum);  // tensor patterns captured implicitly

  return pat_mat;
}();

}  // namespace detail

ExprPtr parse_expr(std::wstring_view raw, Symmetry tensor_sym);

ExprPtr parse_expr_asymm(std::wstring_view raw);

}  // namespace sequant::utils

#endif  // SEQUANT_PARSE_EXPR_HPP
