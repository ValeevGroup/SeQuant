#include "SeQuant/core/parse_expr.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/parse/regex_sequant.hpp>
#include <SeQuant/core/tensor.hpp>

#include <boost/regex.hpp>

namespace sequant {

namespace details {
std::wstring deparse_expr(Tensor const& tnsr, bool annot_sym) {
  std::wstring idx_str = tnsr.label().data();

  auto add_idx = [&idx_str](Index const& idx) {
    std::wstring str = idx.label().data();
    if (idx.has_proto_indices()) {
      str += L"<";  // begin proto-indices token
      for (auto&& pi : idx.proto_indices()) {
        str += pi.label().data();
        str += L",";
      }
      // remove that trailing comma
      str.pop_back();
      str += L">";  // end proto-indices token
    }
    idx_str += str;
  };

  idx_str += L"{";  // open tensor bra-ket and begin bra
  for (auto const& idx : tnsr.bra()) {
    add_idx(idx);
    idx_str += L",";
  }
  *idx_str.rbegin() = L';';  // begin ket
  for (auto const& idx : tnsr.ket()) {
    add_idx(idx);
    idx_str += L",";
  }

  *idx_str.rbegin() = L'}';  // close tensor bra-ket
  if (annot_sym) {
    if (tnsr.symmetry() == Symmetry::antisymm)
      idx_str += L":A";
    else if (tnsr.symmetry() == Symmetry::symm)
      idx_str += L":S";
    else if (tnsr.symmetry() == Symmetry::nonsymm)
      idx_str += L":N";
    else
      idx_str += L":INVALID";
  }
  return idx_str;
}

std::wstring deparse_expr(Product const& prod, bool annot_sym) {
  auto str = std::wstring{};

  auto scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    std::wstring scalar_latex = Constant{scal}.to_latex();

    static auto const decimal_rgx =
        boost::wregex{std::wstring{} + L"(-)?" +
                      parse::regex_patterns::abs_real_num().data()};

    container::svector<std::wstring, 2> num_denom{};
    auto const end = boost::wsregex_iterator{};
    for (auto iter = boost::wsregex_iterator{scalar_latex.begin(),
                                             scalar_latex.end(), decimal_rgx};
         iter != end; ++iter) {
      std::wstring n{};
      if ((*iter)[1].matched) n += (*iter)[1].str() + L" ";
      n += (*iter)[2].str();
      num_denom.emplace_back(n);
    }

    assert(num_denom.size() == 1 || num_denom.size() == 2);

    std::wstring scalar_text =
        num_denom[0] + (num_denom.size() == 2 ? L"/" + num_denom[1] : L"");

    if (scalar_text == L"- 1")
      scalar_text.pop_back();
    else if (scalar_text == L"1")
      scalar_text.pop_back();

    if (*scalar_text.rbegin() != L' ' && !scalar_text.empty())
      scalar_text += L" * ";
    str += scalar_text;
  }

  str += deparse_expr(prod.factor(0), annot_sym);
  for (auto&& xpr : ranges::views::tail(prod.factors()))
    str += L" * " + deparse_expr(xpr, annot_sym);

  return str;
}

std::wstring deparse_expr(Sum const& sum, bool annot_sym) {
  auto str = deparse_expr(sum.summand(0), annot_sym);
  for (auto&& xpr : ranges::views::tail(sum.summands())) {
    auto to_app = deparse_expr(xpr, annot_sym);
    str += (to_app.front() == L'-' ? L" " : L" + ");
    str += to_app;
  }
  return str;
}

}  // namespace details

std::wstring deparse_expr(ExprPtr expr, bool annot_sym) {
  using namespace details;
  if (expr->is<Tensor>())
    return deparse_expr(expr->as<Tensor>(), annot_sym);
  else if (expr->is<Sum>())
    return deparse_expr(expr->as<Sum>(), annot_sym);
  else if (expr->is<Product>())
    return deparse_expr(expr->as<Product>(), annot_sym);
  else
    throw std::runtime_error("Unsupported expr type for deparse!");
}

}  // namespace sequant
