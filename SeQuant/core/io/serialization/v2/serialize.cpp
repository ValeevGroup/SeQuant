#include <SeQuant/core/io/serialization/serialization.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/function_registry.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::io::serialization::v2 {

namespace {

std::wstring serialize_index(const Index &index) {
  if (!index) return {};  // a blank (null) index occupies its position silently

  std::wstring result(index.label());
  if (auto instance = index.basis().unnamed_instance())
    result += L"#" + std::to_wstring(*instance);
  if (index.has_proto_indices()) {
    result += L"<";
    result += join_strings<std::wstring>(index.proto_indices(), L",",
                                         serialize_index);
    result += L">";
  }
  return result;
}

template <typename Range>
std::wstring serialize_indices(Range &&indices) {
  return join_strings<std::wstring>(
      indices, L",", [](const Index &idx) { return serialize_index(idx); });
}

template <typename Range>
std::wstring serialize_ops(const Range &ops) {
  return join_strings<std::wstring>(
      ops, L",", [](const auto &op) { return serialize_index(op.index()); });
}

std::wstring serialize_symm(Symmetry symm) {
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

std::wstring serialize_symm(BraKetSymmetry symm) {
  switch (symm) {
    case BraKetSymmetry::Conjugate:
      return L"bkC";
    case BraKetSymmetry::Symm:
      return L"bkS";
    case BraKetSymmetry::Nonsymm:
      return L"bkN";
  }
  SEQUANT_UNREACHABLE;
}

std::wstring serialize_symm(ColumnSymmetry symm) {
  switch (symm) {
    case ColumnSymmetry::Symm:
      return L"pS";
    case ColumnSymmetry::Nonsymm:
      return L"pN";
  }
  SEQUANT_UNREACHABLE;
}

std::wstring serialize_scalar(const Constant::scalar_type &scalar) {
  if (scalar == 0) return L"0";

  // the imaginary unit goes onto the numerator ("1i/2"): "1/2i" would read as
  // 1/(2i)
  auto part = [](const auto &rational, std::string_view unit = "") {
    std::string result = numerator(rational).str() + std::string(unit);
    if (denominator(rational) != 1) result += "/" + denominator(rational).str();
    return result;
  };

  const auto &real = scalar.real();
  const auto &imag = scalar.imag();

  if (imag == 0) return toUtf16(part(real));
  if (real == 0) return toUtf16(part(imag, "i"));

  // mixed complex constants have no single-token form; emit a parenthesized sum
  std::string result = "(" + part(real);
  result += numerator(imag) < 0 ? " - " : " + ";
  auto imag_mag = imag < 0 ? -imag : imag;
  result += part(imag_mag, "i") + ")";
  return toUtf16(result);
}

std::wstring serialize_tensor(const AbstractTensor &tensor,
                              const SerializationOptions &options) {
  std::wstring result(tensor._label());
  result += L"[" + serialize_indices(tensor._bra());
  if (tensor._ket_rank() > 0 || tensor._aux_rank() > 0)
    result += L";" + serialize_indices(tensor._ket());
  if (tensor._aux_rank() > 0) result += L";" + serialize_indices(tensor._aux());
  result += L"]";

  if (options.annot_symm) {
    result += L":" + serialize_symm(tensor._symmetry());
    result += L"-" + serialize_symm(tensor._braket_symmetry());
    result += L"-" + serialize_symm(tensor._column_symmetry());
  }
  return result;
}

template <Statistics S>
std::wstring serialize_normal_operator(const NormalOperator<S> &nop) {
  std::wstring result(nop.label());
  result += L"{" + serialize_ops(nop.annihilators());
  if (nop.ncreators() > 0) result += L";" + serialize_ops(nop.creators());
  result += L"}";
  return result;
}

template <Statistics S>
std::wstring serialize_raw_operator(const Operator<S> &op) {
  std::wstring result = L"~op(";
  result += join_strings<std::wstring>(op, L", ", [](const Op<S> &o) {
    return serialize_index(o.index()) +
           (o.action() == Action::Create ? L"+" : L"-");
  });
  result += S == Statistics::FermiDirac ? L"):F" : L"):B";
  return result;
}

template <Statistics S>
std::wstring serialize_operator_sequence(const NormalOperatorSequence<S> &seq) {
  std::wstring result = L"~nopseq(";
  result +=
      join_strings<std::wstring>(seq, L", ", [](const NormalOperator<S> &nop) {
        return serialize_normal_operator(nop);
      });
  result += L")";
  return result;
}

}  // namespace

std::wstring to_string(const Index &index, const SerializationOptions &) {
  return serialize_index(index);
}

std::wstring to_string(const Tensor &tensor,
                       const SerializationOptions &options) {
  return serialize_tensor(tensor, options);
}

std::wstring to_string(const AbstractTensor &tensor,
                       const SerializationOptions &options) {
  return serialize_tensor(tensor, options);
}

template <Statistics S>
std::wstring to_string(const NormalOperator<S> &nop,
                       const SerializationOptions &) {
  return serialize_normal_operator(nop);
}
template std::wstring to_string(const NormalOperator<Statistics::FermiDirac> &,
                                const SerializationOptions &);
template std::wstring to_string(
    const NormalOperator<Statistics::BoseEinstein> &,
    const SerializationOptions &);

std::wstring to_string(const Expr &expr, const SerializationOptions &options) {
  if (expr.is<Tensor>()) return serialize_tensor(expr.as<Tensor>(), options);
  if (expr.is<NormalOperator<Statistics::FermiDirac>>())
    return serialize_normal_operator(
        expr.as<NormalOperator<Statistics::FermiDirac>>());
  if (expr.is<NormalOperator<Statistics::BoseEinstein>>())
    return serialize_normal_operator(
        expr.as<NormalOperator<Statistics::BoseEinstein>>());
  if (expr.is<Operator<Statistics::FermiDirac>>())
    return serialize_raw_operator(expr.as<Operator<Statistics::FermiDirac>>());
  if (expr.is<Operator<Statistics::BoseEinstein>>())
    return serialize_raw_operator(
        expr.as<Operator<Statistics::BoseEinstein>>());
  if (expr.is<NormalOperatorSequence<Statistics::FermiDirac>>())
    return serialize_operator_sequence(
        expr.as<NormalOperatorSequence<Statistics::FermiDirac>>());
  if (expr.is<NormalOperatorSequence<Statistics::BoseEinstein>>())
    return serialize_operator_sequence(
        expr.as<NormalOperatorSequence<Statistics::BoseEinstein>>());
  if (expr.is<Constant>()) return serialize_scalar(expr.as<Constant>().value());
  if (expr.is<Variable>()) {
    const Variable &var = expr.as<Variable>();
    return std::wstring(var.label()) + (var.conjugated() ? L"^*" : L"");
  }
  if (expr.is<Power>()) {
    const Power &power = expr.as<Power>();
    std::wstring base = v2::to_string(*power.base(), options);
    if (power.base()->is<Constant>() ||
        (power.base()->is<Variable>() &&
         power.base()->as<Variable>().conjugated()))
      base = L"(" + base + L")";
    std::wstring result =
        base + L"^(" +
        serialize_scalar(Constant::scalar_type{power.exponent()}) + L")";
    return power.conjugated() ? L"(" + result + L")^*" : result;
  }
  if (expr.is<Product>()) {
    const Product &product = expr.as<Product>();
    std::wstring result;
    // a single-factor Product needs its prefactor to not read back as the bare
    // factor
    if (product.scalar() != Constant::scalar_type(1) ||
        product.factors().size() == 1)
      result += serialize_scalar(product.scalar()) + L" ";
    result += join_strings<std::wstring>(
        product.factors(), L" ", [&](const ExprPtr &factor) {
          std::wstring f = v2::to_string(*factor, options);
          return (factor->is<Sum>() || factor->is<Product>()) ? L"(" + f + L")"
                                                              : f;
        });
    return result;
  }
  if (expr.is<Sum>()) {
    const Sum &sum = expr.as<Sum>();
    std::wstring result;
    bool first = true;
    for (const ExprPtr &summand : sum.summands()) {
      std::wstring current = v2::to_string(*summand, options);
      if (first) {
        result = std::move(current);
        first = false;
      } else if (!current.empty() && current.front() == L'-') {
        result += L" - " + current.substr(1);
      } else {
        result += L" + " + current;
      }
    }
    return result;
  }

  // last resort: a type registered with the ~name(...) extension registry
  if (auto serialized = FunctionRegistry::instance().serialize(expr, options))
    return *serialized;

  throw SerializationError("unsupported expression type for v2 serialization");
}

std::wstring to_string(const ExprPtr &expr,
                       const SerializationOptions &options) {
  return expr ? v2::to_string(*expr, options) : std::wstring{};
}

std::wstring to_string(const ResultExpr &result,
                       const SerializationOptions &options) {
  std::wstring lhs =
      result.produces_tensor()
          ? serialize_tensor(result.result_as_tensor(L"?"), options)
          : [&] {
              const Variable &var = result.result_as_variable(L"?");
              return std::wstring(var.label()) +
                     (var.conjugated() ? L"^*" : L"");
            }();
  return lhs + L" = " + v2::to_string(result.expression(), options);
}

}  // namespace sequant::io::serialization::v2
