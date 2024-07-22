#include <SeQuant/core/parse.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor.hpp>

#include <range/v3/all.hpp>

#include <cassert>
#include <codecvt>
#include <cstddef>
#include <locale>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sequant {

namespace details {

template <typename Range>
std::wstring deparse_indices(const Range& indices) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < indices.size(); ++i) {
    deparsed += deparse(indices[i]);

    if (i + 1 < indices.size()) {
      deparsed += L",";
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

  SEQUANT_PRAGMA_CLANG(diagnostic push)
  SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
  SEQUANT_PRAGMA_GCC(diagnostic push)
  SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")

  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.from_bytes(deparsed);

  SEQUANT_PRAGMA_CLANG(diagnostic pop)
  SEQUANT_PRAGMA_GCC(diagnostic pop)
}

}  // namespace details

std::wstring deparse(const ExprPtr& expr, bool annot_sym) {
  using namespace details;
  if (!expr) return {};
  if (expr->is<Tensor>())
    return deparse(expr->as<Tensor>(), annot_sym);
  else if (expr->is<Sum>())
    return deparse(expr->as<Sum>(), annot_sym);
  else if (expr->is<Product>())
    return deparse(expr->as<Product>(), annot_sym);
  else if (expr->is<Constant>())
    return deparse(expr->as<Constant>());
  else if (expr->is<Variable>())
    return deparse(expr->as<Variable>());
  else
    throw std::runtime_error("Unsupported expr type for deparse!");
}

std::wstring deparse(const Index& index) {
  std::wstring deparsed(index.label());

  if (index.has_proto_indices()) {
    deparsed += L"<";
    const auto& protos = index.proto_indices();
    for (std::size_t i = 0; i < protos.size(); ++i) {
      deparsed += protos[i].label();

      if (i + 1 < protos.size()) {
        deparsed += L",";
      }
    }
    deparsed += L">";
  }

  return deparsed;
}

std::wstring deparse(Tensor const& tensor, bool annot_sym) {
  std::wstring deparsed(tensor.label());
  deparsed += L"{" + details::deparse_indices(tensor.bra());
  if (tensor.ket_rank() > 0) {
    deparsed += L";" + details::deparse_indices(tensor.ket());
  }
  deparsed += L"}";

  if (annot_sym) {
    deparsed += L":" + details::deparse_sym(tensor.symmetry());
  }

  return deparsed;
}

std::wstring deparse(const Constant& constant) {
  return details::deparse_scalar(constant.value());
}

std::wstring deparse(const Variable& variable) {
  return std::wstring(variable.label());
}

std::wstring deparse(Product const& prod, bool annot_sym) {
  std::wstring deparsed;

  const auto& scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    deparsed += details::deparse_scalar(scal) + L" ";
  }

  for (std::size_t i = 0; i < prod.size(); ++i) {
    const ExprPtr& current = prod[i];
    bool parenthesize = false;
    if (current->is<Product>() || current->is<Sum>()) {
      parenthesize = true;
      deparsed += L"(";
    }

    deparsed += deparse(current, annot_sym);

    if (parenthesize) {
      deparsed += L")";
    }

    if (i + 1 < prod.size()) {
      deparsed += L" * ";
    }
  }

  return deparsed;
}

std::wstring deparse(Sum const& sum, bool annot_sym) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < sum.size(); ++i) {
    ExprPtr& current = sum[i];

    const bool parenthesize = current->is<Sum>();

    std::wstring current_deparsed = deparse(current, annot_sym);

    bool is_negative = false;
    if (parenthesize) {
      current_deparsed += L"(" + current_deparsed + L")";
    } else {
      is_negative = current_deparsed.front() == L'-';
    }

    if (i > 0) {
      if (is_negative) {
        deparsed += L" - " + current_deparsed.substr(1);
      } else {
        deparsed += L" + " + current_deparsed;
      }
    } else {
      deparsed = std::move(current_deparsed);
    }
  }
  return deparsed;
}

}  // namespace sequant
