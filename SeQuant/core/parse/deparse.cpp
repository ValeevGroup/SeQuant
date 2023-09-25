#include "SeQuant/core/parse_expr.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

#include <boost/multiprecision/cpp_int.hpp>

#include <cassert>
#include <codecvt>
#include <locale>
#include <string>

namespace sequant {

namespace details {

std::wstring deparse_index(const Index& index) {
  std::wstring deparsed(index.label());

  if (index.has_proto_indices()) {
    deparsed += L"<";
    const auto& protos = index.proto_indices();
    for (std::size_t i = 0; i < protos.size(); ++i) {
      deparsed += protos[i].label();

      if (i + 1 < protos.size()) {
        deparsed += L", ";
      }
    }
    deparsed += L">";
  }

  return deparsed;
}

template <typename Range>
std::wstring deparse_indices(const Range& indices) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < indices.size(); ++i) {
    deparsed += deparse_index(indices[i]);

    if (i + 1 < indices.size()) {
      deparsed += L", ";
    }
  }

  return deparsed;
}

std::wstring deparse_sym(Symmetry sym) {
  switch (sym) {
    case Symmetry::symm:
      return L"S";
    case Symmetry::antisymm:
      return L"A";
    case Symmetry::nonsymm:
      return L"N";
    case Symmetry::invalid:
      return L"INVALID";
  }

  assert(false);
  return L"INVALIDANDUNREACHABLE";
}

std::wstring deparse_expr(Tensor const& tensor, bool annot_sym) {
  return std::wstring(tensor.label()) + L"^{" + deparse_indices(tensor.ket()) +
         L"}_{" + deparse_indices(tensor.bra()) + L"}" +
         (annot_sym ? std::wstring(L":") + deparse_sym(tensor.symmetry())
                    : std::wstring{});
}

std::wstring deparse_scalar(const Constant::scalar_type& scalar) {
  const auto& real = scalar.real();
  const auto& realNumerator = boost::multiprecision::numerator(real);
  const auto& realDenominator = boost::multiprecision::denominator(real);
  const auto& imag = scalar.imag();
  const auto& imagNumerator = boost::multiprecision::numerator(imag);
  const auto& imagDenominator = boost::multiprecision::denominator(imag);

  std::string deparsed;
  if (realNumerator != 0) {
    deparsed += realNumerator.str();

    if (realDenominator != 1) {
      deparsed += "/" + realDenominator.str();
    }
  }
  if (imagNumerator != 0) {
    if (!deparsed.empty()) {
      deparsed += imagNumerator < 0 ? " + i " : " - i ";
    }

    deparsed += imagNumerator.str();

    if (imagDenominator != 1) {
      deparsed += "/" + imagDenominator.str();
    }
  }

  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.from_bytes(deparsed);
}

std::wstring deparse_expr(const Constant& constant) {
  return deparse_scalar(constant.value());
}

std::wstring deparse_expr(const Variable& variable) {
  return std::wstring(variable.label());
}

std::wstring deparse_expr(Product const& prod, bool annot_sym) {
  std::wstring deparsed;

  const auto& scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    deparsed += deparse_scalar(scal) + L" ";
  }

  deparsed += deparse_expr(prod.factor(0), annot_sym);
  for (auto&& xpr : ranges::views::tail(prod.factors()))
    deparsed += L" * " + deparse_expr(xpr, annot_sym);

  return deparsed;
}

std::wstring deparse_expr(Sum const& sum, bool annot_sym) {
  std::wstring deparsed = deparse_expr(sum.summand(0), annot_sym);
  for (auto&& xpr : ranges::views::tail(sum.summands())) {
    std::wstring current = deparse_expr(xpr, annot_sym);

    if (current.front() == L'-') {
      deparsed += L" - " + current.substr(1);
    } else {
      deparsed += L" + " + current;
    }
  }
  return deparsed;
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
  else if (expr->is<Constant>())
    return deparse_expr(expr->as<Constant>());
  else if (expr->is<Variable>())
    return deparse_expr(expr->as<Variable>());
  else
    throw std::runtime_error("Unsupported expr type for deparse!");
}

}  // namespace sequant
