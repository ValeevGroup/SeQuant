#include "parse_expr.hpp"

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/tensor.hpp>
#include <range/v3/view.hpp>
#include <regex>
#include <string>

namespace sequant {

namespace detail {

std::wstring prune_space(std::wstring_view raw) {
  return std::regex_replace(raw.data(), std::wregex(L"[[:space:]]+"), L"");
}

double to_decimal(std::wstring_view numstr) {
  double num;
  std::wistringstream wiss{numstr.data()};
  wiss >> num;
  return wiss ? num : 0;
}

double as_fraction(std::wstring_view num, std::wstring_view denom) {
  return (num.empty() ? 1.0 : to_decimal(num)) /
         (denom.empty() ? 1.0 : to_decimal(denom));
}

bool validate_term(std::wstring_view expr) {
  return std::regex_match(
      expr.data(),
      std::wregex(expr_rgx_pat.at("term"), std::regex_constants::nosubs));
}

bool validate_sum(std::wstring_view expr) {
  return std::regex_match(
      expr.data(),
      std::wregex(expr_rgx_pat.at("sum"), std::regex_constants::nosubs));
}

ExprPtr make_tensor(const std::wsmatch& match_obj, Symmetry tensor_sym) {
  // HELPER LAMBDAS
  /*
   * i1   : i_1
   * a10  : a_10
   * i_1  : i_1
   * a_10 : a_10
   */
  auto add_underscore = [](const std::wstring& lbl) -> std::wstring {
    auto noUscore =
        lbl |
        ranges::views::filter([](const auto& x) { return x != char{'_'}; }) |
        ranges::to<std::wstring>;
    return noUscore.substr(0, 1) + L"_" + noUscore.substr(1);
  };

  auto indices = [&](const std::wstring& bk) {
    return bk | ranges::views::tokenize(std::wregex(expr_rgx_pat.at("index"))) |
           ranges::views::transform(add_underscore) |
           ranges::to<container::vector<Index>>;
  };
  // HELPER LAMBDAS END

  const auto& lbl = match_obj[1].str();
  const auto& bra = match_obj[2].str();
  const auto& ket = match_obj[3].str();

  return ex<Tensor>(Tensor{lbl, indices(bra), indices(ket), tensor_sym});
}

ExprPtr parse_tensor_term(const std::wstring& expr, Symmetry tensor_sym) {
  std::wsmatch match_obj;
  std::regex_search(expr, match_obj, std::wregex(expr_rgx_pat.at("tensor")));
  return make_tensor(match_obj, tensor_sym);
}

ExprPtr parse_product_term(const std::wstring& expr, Symmetry tensor_sym) {
  using namespace std::string_literals;
  // find the scalar first
  std::wsmatch match_obj;
  std::regex_search(
      expr, match_obj,
      std::wregex(LR"(^(\+|-)?)"s +                       // optional +/- sign
                  L"(?:" + expr_rgx_pat.at("fraction") +  // optional fraction
                  L")?"));
  double scalar =
      (!match_obj[1].matched || match_obj[1].str() == L"+") ? 1.0 : -1.0;

  scalar *= as_fraction(match_obj[2].matched ? match_obj[2].str() : L"1.0",
                        match_obj[3].matched ? match_obj[3].str() : L"1.0");

  // find the tensors
  const auto& tregex = std::wregex(expr_rgx_pat.at("tensor"));
  std::wsregex_iterator beg(expr.cbegin(), expr.cend(), tregex);
  std::wsregex_iterator end{};

  auto tensors = ranges::make_subrange(beg, end) |
                 ranges::views::transform([&tensor_sym](const auto& x) {
                   return make_tensor(x, tensor_sym);
                 }) |
                 ranges::to<container::svector<ExprPtr>>;

  if (scalar == 1.0 && tensors.size() == 1) return tensors[0];

  return ex<Product>(Product(scalar, tensors.begin(), tensors.end()));
}

ExprPtr parse_term(const std::wstring& expr, Symmetry tensor_sym) {
  if (std::regex_match(expr, std::wregex(expr_rgx_pat.at("tensor_term"),
                                         std::regex_constants::nosubs))) {
    return parse_tensor_term(expr, tensor_sym);
  } else if (std::regex_match(expr,
                              std::wregex(expr_rgx_pat.at("product_term"),
                                          std::regex_constants::nosubs))) {
    return parse_product_term(expr, tensor_sym);
  } else {
    return nullptr;  // not reachable
  }
}

ExprPtr parse_sum(const std::wstring& expr, Symmetry tensor_sym) {
  auto summands = expr |
                  ranges::views::tokenize(std::wregex{
                      expr_rgx_pat.at("term"), std::regex_constants::nosubs}) |
                  ranges::views::transform([&tensor_sym](const auto& x) {
                    return parse_term(x, tensor_sym);
                  });
  return ex<Sum>(Sum{summands.begin(), summands.end()});
}

}  // namespace detail

ExprPtr parse_expr(std::wstring_view raw, Symmetry tensor_sym) {
  const auto& expr = detail::prune_space(raw);
  using detail::expr_rgx_pat;
  if (std::regex_match(expr, std::wregex{expr_rgx_pat.at("term"),
                                         std::regex_constants::nosubs}))
    return detail::parse_term(expr, tensor_sym);
  else if (std::regex_match(expr, std::wregex(expr_rgx_pat.at("sum"),
                                              std::regex_constants::nosubs))) {
    return detail::parse_sum(expr, tensor_sym);
  } else
    throw std::logic_error("Expression parse error. Invalid expression spec.");
}

ExprPtr parse_expr_asymm(std::wstring_view raw){
return parse_expr(raw, Symmetry::antisymm);
}

}  // namespace sequant
