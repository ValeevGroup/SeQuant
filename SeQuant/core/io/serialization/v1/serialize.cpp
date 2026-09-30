#include <SeQuant/core/io/serialization/serialization.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::io::serialization::v1 {

namespace details {

template <typename Range>
std::wstring serialize_indices(Range&& indices,
                               const SerializationOptions& options) {
  return join_strings<std::wstring>(indices, L",", [&](const Index& idx) {
    return v1::to_string(idx, options);
  });
}

template <typename Range>
std::wstring serialize_ops(const Range& ops,
                           const SerializationOptions& options) {
  return join_strings<std::wstring>(ops, L",", [&](const auto& op) {
    return v1::to_string(op.index(), options);
  });
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

std::wstring serialize_symm(BraKetSymmetry symm, Hermiticity hermiticity,
                            const SerializationOptions&) {
  // a definite hermiticity is the trait the observable exchange symmetry is
  // derived from, so it is spelled with its own letter and resolved back
  // against the parity and the indices' field when the tensor is rebuilt (see
  // to_string(const AbstractTensor&)); this keeps the spelling of a real-basis
  // tensor readable over a complex one too
  switch (hermiticity) {
    case Hermiticity::Hermitian:
      return L"H";
    case Hermiticity::AntiHermitian:
      return L"A";
    case Hermiticity::NonHermitian:
      break;
  }

  // an indefinite hermiticity derives Nonsymm in every basis, so `N` is the
  // only letter a tensor reaches here with; the pin letters stay for a tensor
  // kind that spells an exchange symmetry its traits do not derive
  switch (symm) {
    case BraKetSymmetry::Conjugate:
      return L"C";
    case BraKetSymmetry::Symm:
      return L"S";
    case BraKetSymmetry::Nonsymm:
      return L"N";
    case BraKetSymmetry::Antisymm:
    case BraKetSymmetry::AntiConjugate:
      // only a definite hermiticity derives a signed exchange
      throw Exception(
          "io::serialization::v1: a signed bra/ket exchange symmetry with an "
          "indefinite hermiticity has no spelling");
  }

  SEQUANT_UNREACHABLE;
}

std::wstring serialize_symm(ConjugationParity parity,
                            const SerializationOptions&) {
  switch (parity) {
    case ConjugationParity::Even:
      return L"E";
    case ConjugationParity::Odd:
      return L"O";
    case ConjugationParity::None:
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

/// @return whether @p scalar spells as a real term plus an imaginary one, so
///         that it needs parentheses wherever juxtaposition binds tighter
///         than `+` -- as a factor of a Product does
bool is_composite_scalar(const Constant::scalar_type& scalar) {
  return numerator(scalar.real()) != 0 && numerator(scalar.imag()) != 0;
}

std::wstring serialize_scalar(const Constant::scalar_type& scalar,
                              const SerializationOptions&) {
  if (scalar == 0) {
    return L"0";
  }

  // a rational spells as `n` or `n/d`; the imaginary part is the same spelling
  // with an `i` abutting it (`2i`, `1/2i`), which is the grammar's imaginary
  // literal
  const auto spell = [](const auto& num, const auto& den) {
    std::string serialized = num.str();
    if (den != 1) {
      serialized += "/" + den.str();
    }
    return serialized;
  };

  const auto& real = scalar.real();
  const auto& imag = scalar.imag();
  auto imagNumerator = numerator(imag);

  std::string serialized;
  if (numerator(real) != 0) {
    serialized += spell(numerator(real), denominator(real));
  }
  if (imagNumerator != 0) {
    if (!serialized.empty()) {
      if (imagNumerator < 0) {
        serialized += " - ";
        imagNumerator *= -1;
      } else {
        serialized += " + ";
      }
    }

    serialized += spell(imagNumerator, denominator(imag)) + "i";
  }

  SEQUANT_ASSERT(!serialized.empty());

  return toUtf16(serialized);
}

/// @return serialize_scalar(), parenthesized where the spelling is composite
std::wstring serialize_scalar_atom(const Constant::scalar_type& scalar,
                                   const SerializationOptions& options) {
  auto serialized = serialize_scalar(scalar, options);
  if (is_composite_scalar(scalar)) {
    return L"(" + std::move(serialized) + L")";
  }
  return serialized;
}

std::wstring to_string(Tensor const& tensor,
                       const SerializationOptions& options) {
  auto serialized =
      to_string(static_cast<const AbstractTensor&>(tensor), options);
  // the states are spelled as marks trailing the label, the form the
  // deserializer grammar accepts, so the round-trip is lossless
  serialized.replace(0, tensor.label().size(), tensor.decorated_label());
  return serialized;
}

std::wstring to_string(const Constant& constant,
                       const SerializationOptions& options) {
  return details::serialize_scalar(constant.value(), options);
}

std::wstring to_string(const Variable& variable, const SerializationOptions&) {
  // the state is spelled as a mark trailing the label, the form the
  // deserializer grammar accepts, so the round-trip is lossless
  return variable.decorated_label();
}

std::wstring to_string(const Power& power,
                       const SerializationOptions& options) {
  std::wstring core;
  const auto& base = power.base();
  // parenthesize bases for Constants and conjugated Variables
  const bool parenthesize_base =
      base->is<Constant>() ||
      (base->is<Variable>() && base->as<Variable>().conjugated());
  if (parenthesize_base) core += L"(";
  core += v1::to_string(*base, options);
  if (parenthesize_base) core += L")";
  core += L"^(";
  core += serialize_scalar(Constant::scalar_type{power.exponent()}, options);
  core += L")";
  if (power.conjugated()) return L"(" + std::move(core) + L")^*";
  return core;
}

/// Spells `Re[...]` / `Im[...]` with a real scalar of the wrapped expression
/// hoisted in front of the brackets, as `Re(c X) = c Re(X)` allows for a real
/// c. That is the form the real_part()/imaginary_part() builders and the
/// wrappers' canonicalization produce, so a wrapper that reached the
/// serializer by neither route still spells as its canonical self and the
/// round trip is a fixed point. A complex scalar stays inside: neither
/// wrapper is linear over it.
std::wstring serialize_wrapper(std::wstring_view projection,
                               const ExprPtr& inner,
                               const SerializationOptions& options) {
  std::wstring serialized;
  ExprPtr wrapped = inner;

  if (inner->is<Product>()) {
    const auto& prod = inner->as<Product>();
    const auto scalar = prod.scalar();
    if (scalar.imag() == 0 && scalar != Product::scalar_type{1}) {
      serialized = scalar == Product::scalar_type{-1}
                       ? std::wstring(L"-")
                       : serialize_scalar_atom(scalar, options) + L" ";
      wrapped = sequant::detail::strip_scalar(prod);
    }
  }

  serialized += projection;
  serialized += L"[" + v1::to_string(*wrapped, options) + L"]";
  return serialized;
}

std::wstring to_string(const RealPart& part,
                       const SerializationOptions& options) {
  return serialize_wrapper(L"Re", part.inner(), options);
}

std::wstring to_string(const ImagPart& part,
                       const SerializationOptions& options) {
  return serialize_wrapper(L"Im", part.inner(), options);
}

std::wstring to_string(Product const& prod,
                       const SerializationOptions& options) {
  std::wstring serialized;

  const auto& scal = prod.scalar();
  if (scal == Product::scalar_type{-1}) {
    // a negated product spells its sign rather than a `-1` factor, so that a
    // subtracted summand reads `- b` and not `- 1 b`
    serialized += L"-";
  } else if (scal != Product::scalar_type{1}) {
    serialized += details::serialize_scalar_atom(scal, options) + L" ";
  }

  for (std::size_t i = 0; i < prod.size(); ++i) {
    const ExprPtr& current = prod[i];
    bool parenthesize = false;
    // a composite scalar reads as a sum where a factor is expected, so it is
    // parenthesized here as it is in the scalar position above
    if (current->is<Product>() || current->is<Sum>() ||
        (current->is<Constant>() &&
         is_composite_scalar(current->as<Constant>().value()))) {
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
    const ExprPtr& current = sum[i];

    const bool parenthesize = current->is<Sum>();

    std::wstring current_serialized = to_string(current, options);

    bool is_negative = false;
    if (parenthesize) {
      current_serialized = L"(" + current_serialized + L")";
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
  else if (expr.is<Power>())
    return details::to_string(expr.as<Power>(), options);
  else if (expr.is<RealPart>())
    return details::to_string(expr.as<RealPart>(), options);
  else if (expr.is<ImagPart>())
    return details::to_string(expr.as<ImagPart>(), options);
  else
    throw Exception("Unsupported expr type for serialize!");
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
    serialized +=
        join_strings<std::wstring>(index.proto_indices(), L",", &Index::label);
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
    const auto braket_symmetry = tensor._braket_symmetry();
    const auto hermiticity = tensor._hermiticity();
    serialized +=
        L"-" + details::serialize_symm(braket_symmetry, hermiticity, options);
    serialized +=
        L"-" + details::serialize_symm(tensor._column_symmetry(), options);
    // the fourth letter is optional: emitted only when the parity the parser
    // back-fills from the bra/ket letter is not the one this tensor carries.
    // A definite hermiticity is spelled with its trait letter (H/A), which
    // says nothing about the parity, so the parser leaves it at the Tensor
    // default; with an indefinite hermiticity the exchange letter is what the
    // parser reads the parity from, against the same base field.
    const ConjugationParity implied_parity =
        hermiticity == Hermiticity::NonHermitian
            ? to_conjugation_parity(braket_symmetry, tensor._base_field())
            : Tensor::Defaults::conjugation_parity;
    if (tensor._conjugation_parity() != implied_parity) {
      serialized +=
          L"-" + details::serialize_symm(tensor._conjugation_parity(), options);
    }
  }

  return serialized;
}

std::wstring to_string(Tensor const& tensor,
                       const SerializationOptions& options) {
  return v1::to_string(static_cast<const Expr&>(tensor), options);
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
