#include <SeQuant/core/io/serialization/serialization.hpp>

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

namespace sequant::io::serialization::v1 {

namespace details {

template <typename Range>
std::wstring serialize_indices(Range&& indices,
                               const SerializationOptions& options) {
  std::wstring serialized;

  for (std::size_t i = 0; i < indices.size(); ++i) {
    serialized += v1::to_string(indices[i], options);

    if (i + 1 < indices.size()) {
      serialized += L",";
    }
  }

  return serialized;
}

template <typename Range>
std::wstring serialize_ops(const Range& ops,
                           const SerializationOptions& options) {
  std::wstring serialized;

  for (std::size_t i = 0; i < ops.size(); ++i) {
    serialized += v1::to_string(ops[i].index(), options);

    if (i + 1 < ops.size()) {
      serialized += L",";
    }
  }

  return serialized;
}

std::wstring serialize_symm(Symmetry symm, const SerializationOptions&) {
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

std::wstring serialize_symm(BraKetSymmetry symm, const SerializationOptions&) {
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

std::wstring serialize_symm(ColumnSymmetry symm, const SerializationOptions&) {
  switch (symm) {
    case ColumnSymmetry::Symm:
      return L"S";
    case ColumnSymmetry::Nonsymm:
      return L"N";
  }

  SEQUANT_UNREACHABLE;
}

std::wstring serialize_scalar(const Constant::scalar_type& scalar,
                              const SerializationOptions&) {
  if (scalar == 0) {
    return L"0";
  }

  const auto& real = scalar.real();
  const auto& realNumerator = boost::multiprecision::numerator(real);
  const auto& realDenominator = boost::multiprecision::denominator(real);
  const auto& imag = scalar.imag();
  const auto& imagNumerator = boost::multiprecision::numerator(imag);
  const auto& imagDenominator = boost::multiprecision::denominator(imag);

  std::string serialized;
  if (realNumerator != 0) {
    serialized += realNumerator.str();

    if (realDenominator != 1) {
      serialized += "/" + realDenominator.str();
    }
  }
  if (imagNumerator != 0) {
    if (!serialized.empty()) {
      serialized += imagNumerator < 0 ? " + i " : " - i ";
    }

    serialized += imagNumerator.str();

    if (imagDenominator != 1) {
      serialized += "/" + imagDenominator.str();
    }
  }

  SEQUANT_ASSERT(!serialized.empty());

  return toUtf16(serialized);
}

std::wstring to_string(Tensor const& tensor,
                       const SerializationOptions& options) {
  return to_string(static_cast<const AbstractTensor&>(tensor), options);
}

std::wstring to_string(const Constant& constant,
                       const SerializationOptions& options) {
  return details::serialize_scalar(constant.value(), options);
}

std::wstring to_string(const Variable& variable, const SerializationOptions&) {
  return std::wstring(variable.label()) + (variable.conjugated() ? L"^*" : L"");
}

std::wstring to_string(Product const& prod,
                       const SerializationOptions& options) {
  std::wstring serialized;

  const auto& scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    serialized += details::serialize_scalar(scal, options) + L" ";
  }

  for (std::size_t i = 0; i < prod.size(); ++i) {
    const ExprPtr& current = prod[i];
    bool parenthesize = false;
    if (current->is<Product>() || current->is<Sum>()) {
      parenthesize = true;
      serialized += L"(";
    }

    serialized += to_string(current, options);

    if (parenthesize) {
      serialized += L")";
    }

    if (i + 1 < prod.size()) {
      serialized += L" * ";
    }
  }

  return serialized;
}

std::wstring to_string(Sum const& sum, const SerializationOptions& options) {
  std::wstring serialized;

  for (std::size_t i = 0; i < sum.size(); ++i) {
    ExprPtr& current = sum[i];

    const bool parenthesize = current->is<Sum>();

    std::wstring current_serialized = to_string(current, options);

    bool is_negative = false;
    if (parenthesize) {
      current_serialized += L"(" + current_serialized + L")";
    } else {
      is_negative = current_serialized.front() == L'-';
    }

    if (i > 0) {
      if (is_negative) {
        serialized += L" - " + current_serialized.substr(1);
      } else {
        serialized += L" + " + current_serialized;
      }
    } else {
      serialized = std::move(current_serialized);
    }
  }
  return serialized;
}

}  // namespace details

std::wstring to_string(const ExprPtr& expr,
                       const SerializationOptions& options) {
  if (!expr) return {};

  return v1::to_string(*expr, options);
}

std::wstring to_string(const Expr& expr, const SerializationOptions& options) {
  using namespace details;
  if (expr.is<Tensor>())
    return details::to_string(expr.as<Tensor>(), options);
  else if (expr.is<FNOperator>())
    return v1::to_string(expr.as<FNOperator>(), options);
  else if (expr.is<BNOperator>())
    return v1::to_string(expr.as<BNOperator>(), options);
  else if (expr.is<Sum>())
    return details::to_string(expr.as<Sum>(), options);
  else if (expr.is<Product>())
    return details::to_string(expr.as<Product>(), options);
  else if (expr.is<Constant>())
    return details::to_string(expr.as<Constant>(), options);
  else if (expr.is<Variable>())
    return details::to_string(expr.as<Variable>(), options);
  else
    throw std::runtime_error("Unsupported expr type for serialize!");
}

std::wstring to_string(const ResultExpr& result,
                       const SerializationOptions& options) {
  std::wstring serialized;
  if (result.produces_tensor()) {
    serialized = details::to_string(result.result_as_tensor(L"?"), options);
  } else {
    serialized = details::to_string(result.result_as_variable(L"?"), options);
  }

  return serialized + L" = " + v1::to_string(result.expression(), options);
}

std::wstring to_string(const Index& index, const SerializationOptions&) {
  std::wstring serialized(index.label());

  if (index.has_proto_indices()) {
    serialized += L"<";
    const auto& protos = index.proto_indices();
    for (std::size_t i = 0; i < protos.size(); ++i) {
      serialized += protos[i].label();

      if (i + 1 < protos.size()) {
        serialized += L",";
      }
    }
    serialized += L">";
  }

  return serialized;
}

std::wstring to_string(AbstractTensor const& tensor,
                       const SerializationOptions& options) {
  std::wstring serialized(tensor._label());
  serialized += L"{" + details::serialize_indices(tensor._bra(), options);
  if (tensor._ket_rank() > 0) {
    serialized += L";" + details::serialize_indices(tensor._ket(), options);
  }
  if (tensor._aux_rank() > 0) {
    if (tensor._ket_rank() == 0) {
      serialized += L";";
    }
    serialized += L";" + details::serialize_indices(tensor._aux(), options);
  }
  serialized += L"}";

  if (options.annot_symm) {
    serialized += L":" + details::serialize_symm(tensor._symmetry(), options);
    serialized +=
        L"-" + details::serialize_symm(tensor._braket_symmetry(), options);
    serialized +=
        L"-" + details::serialize_symm(tensor._column_symmetry(), options);
  }

  return serialized;
}

template <Statistics S>
std::wstring to_string(NormalOperator<S> const& nop,
                       const SerializationOptions& options) {
  std::wstring serialized(nop.label());
  serialized += L"{" + details::serialize_ops(nop.annihilators(), options);
  if (nop.ncreators() > 0) {
    serialized += L";" + details::serialize_ops(nop.creators(), options);
  }
  serialized += L"}";

  return serialized;
}

template std::wstring to_string<Statistics::FermiDirac>(
    NormalOperator<Statistics::FermiDirac> const& nop,
    const SerializationOptions& options);
template std::wstring to_string<Statistics::BoseEinstein>(
    NormalOperator<Statistics::BoseEinstein> const& nop,
    const SerializationOptions& options);

}  // namespace sequant::io::serialization::v1
