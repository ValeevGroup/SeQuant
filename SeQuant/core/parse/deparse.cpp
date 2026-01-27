#include <SeQuant/core/parse.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/all.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sequant {

namespace details {

template <typename Range>
std::wstring deparse_indices(Range&& indices, const DeparseOptions& options) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < indices.size(); ++i) {
    deparsed += deparse(indices[i], options);

    if (i + 1 < indices.size()) {
      deparsed += L",";
    }
  }

  return deparsed;
}

template <typename Range>
std::wstring deparse_ops(const Range& ops, const DeparseOptions& options) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < ops.size(); ++i) {
    deparsed += deparse(ops[i].index(), options);

    if (i + 1 < ops.size()) {
      deparsed += L",";
    }
  }

  return deparsed;
}

std::wstring deparse_symm(Symmetry symm, const DeparseOptions&) {
  switch (symm) {
    case Symmetry::Symm:
      return L"S";
    case Symmetry::Antisymm:
      return L"A";
    case Symmetry::Nonsymm:
      return L"N";
  }

  SEQUANT_UNREACHABLE;
}

std::wstring deparse_symm(BraKetSymmetry symm, const DeparseOptions&) {
  switch (symm) {
    case BraKetSymmetry::Conjugate:
      return L"C";
    case BraKetSymmetry::Symm:
      return L"S";
    case BraKetSymmetry::Nonsymm:
      return L"N";
  }

  SEQUANT_UNREACHABLE;
}

std::wstring deparse_symm(ColumnSymmetry symm, const DeparseOptions&) {
  switch (symm) {
    case ColumnSymmetry::Symm:
      return L"S";
    case ColumnSymmetry::Nonsymm:
      return L"N";
  }

  SEQUANT_UNREACHABLE;
}

std::wstring deparse_scalar(const Constant::scalar_type& scalar,
                            const DeparseOptions&) {
  if (scalar == 0) {
    return L"0";
  }

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

  SEQUANT_ASSERT(!deparsed.empty());

  return to_wstring(deparsed);
}

std::wstring deparse(Tensor const& tensor, const DeparseOptions& options) {
  return deparse(static_cast<const AbstractTensor&>(tensor), options);
}

std::wstring deparse(const Constant& constant, const DeparseOptions& options) {
  return details::deparse_scalar(constant.value(), options);
}

std::wstring deparse(const Variable& variable, const DeparseOptions&) {
  return std::wstring(variable.label()) + (variable.conjugated() ? L"^*" : L"");
}

std::wstring deparse(Product const& prod, const DeparseOptions& options) {
  std::wstring deparsed;

  const auto& scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    deparsed += details::deparse_scalar(scal, options) + L" ";
  }

  for (std::size_t i = 0; i < prod.size(); ++i) {
    const ExprPtr& current = prod[i];
    bool parenthesize = false;
    if (current->is<Product>() || current->is<Sum>()) {
      parenthesize = true;
      deparsed += L"(";
    }

    deparsed += deparse(current, options);

    if (parenthesize) {
      deparsed += L")";
    }

    if (i + 1 < prod.size()) {
      deparsed += L" * ";
    }
  }

  return deparsed;
}

std::wstring deparse(Sum const& sum, const DeparseOptions& options) {
  std::wstring deparsed;

  for (std::size_t i = 0; i < sum.size(); ++i) {
    ExprPtr& current = sum[i];

    const bool parenthesize = current->is<Sum>();

    std::wstring current_deparsed = deparse(current, options);

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

}  // namespace details

std::wstring deparse(const ExprPtr& expr, const DeparseOptions& options) {
  if (!expr) return {};

  return deparse(*expr, options);
}

std::wstring deparse(const Expr& expr, const DeparseOptions& options) {
  using namespace details;
  if (expr.is<Tensor>())
    return details::deparse(expr.as<Tensor>(), options);
  else if (expr.is<FNOperator>())
    return deparse(expr.as<FNOperator>(), options);
  else if (expr.is<BNOperator>())
    return deparse(expr.as<BNOperator>(), options);
  else if (expr.is<Sum>())
    return details::deparse(expr.as<Sum>(), options);
  else if (expr.is<Product>())
    return details::deparse(expr.as<Product>(), options);
  else if (expr.is<Constant>())
    return details::deparse(expr.as<Constant>(), options);
  else if (expr.is<Variable>())
    return details::deparse(expr.as<Variable>(), options);
  else
    throw std::runtime_error("Unsupported expr type for deparse!");
}

std::wstring deparse(const ResultExpr& result, const DeparseOptions& options) {
  std::wstring deparsed;
  if (result.produces_tensor()) {
    deparsed = details::deparse(result.result_as_tensor(L"?"), options);
  } else {
    deparsed = details::deparse(result.result_as_variable(L"?"), options);
  }

  return deparsed + L" = " + deparse(result.expression(), options);
}

std::wstring deparse(const Index& index, const DeparseOptions&) {
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

std::wstring deparse(AbstractTensor const& tensor,
                     const DeparseOptions& options) {
  std::wstring deparsed(tensor._label());
  deparsed += L"{" + details::deparse_indices(tensor._bra(), options);
  if (tensor._ket_rank() > 0) {
    deparsed += L";" + details::deparse_indices(tensor._ket(), options);
  }
  if (tensor._aux_rank() > 0) {
    if (tensor._ket_rank() == 0) {
      deparsed += L";";
    }
    deparsed += L";" + details::deparse_indices(tensor._aux(), options);
  }
  deparsed += L"}";

  if (options.annot_symm) {
    deparsed += L":" + details::deparse_symm(tensor._symmetry(), options);
    deparsed +=
        L"-" + details::deparse_symm(tensor._braket_symmetry(), options);
    deparsed +=
        L"-" + details::deparse_symm(tensor._column_symmetry(), options);
  }

  return deparsed;
}

template <Statistics S>
std::wstring deparse(NormalOperator<S> const& nop,
                     const DeparseOptions& options) {
  std::wstring deparsed(nop.label());
  deparsed += L"{" + details::deparse_ops(nop.annihilators(), options);
  if (nop.ncreators() > 0) {
    deparsed += L";" + details::deparse_ops(nop.creators(), options);
  }
  deparsed += L"}";

  return deparsed;
}

template std::wstring deparse<Statistics::FermiDirac>(
    NormalOperator<Statistics::FermiDirac> const& nop,
    const DeparseOptions& options);
template std::wstring deparse<Statistics::BoseEinstein>(
    NormalOperator<Statistics::BoseEinstein> const& nop,
    const DeparseOptions& options);

}  // namespace sequant
