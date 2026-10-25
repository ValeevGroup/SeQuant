#include <SeQuant/core/io/serialization/serialization.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace sequant::io::serialization::v1 {

namespace detail {

template <typename Range>
std::string serialize_indices(Range&& indices,
                              const SerializationOptions& options) {
  return join_strings(indices, ",", [&](const Index& idx) {
    return v1::to_string(idx, options);
  });
}

template <typename Range>
std::string serialize_ops(const Range& ops,
                          const SerializationOptions& options) {
  return join_strings(ops, ",", [&](const auto& op) {
    return v1::to_string(op.index(), options);
  });
}

std::string serialize_symm(Symmetry symm, const SerializationOptions&) {
  switch (symm) {
    case Symmetry::Symm:
      return "S";
    case Symmetry::Antisymm:
      return "A";
    case Symmetry::Nonsymm:
      return "N";
  }

  SEQUANT_UNREACHABLE;
}

std::string serialize_symm(BraKetSymmetry symm, const SerializationOptions&) {
  switch (symm) {
    case BraKetSymmetry::Conjugate:
      return "C";
    case BraKetSymmetry::Symm:
      return "S";
    case BraKetSymmetry::Nonsymm:
      return "N";
  }

  SEQUANT_UNREACHABLE;
}

std::string serialize_symm(ColumnSymmetry symm, const SerializationOptions&) {
  switch (symm) {
    case ColumnSymmetry::Symm:
      return "S";
    case ColumnSymmetry::Nonsymm:
      return "N";
  }

  SEQUANT_UNREACHABLE;
}

std::string serialize_scalar(const Constant::scalar_type& scalar,
                             const SerializationOptions&) {
  if (scalar == 0) {
    return "0";
  }

  const auto& real = scalar.real();
  const auto& realNumerator = numerator(real);
  const auto& realDenominator = denominator(real);
  const auto& imag = scalar.imag();
  auto imagNumerator = numerator(imag);
  const auto& imagDenominator = denominator(imag);

  std::string serialized;
  if (realNumerator != 0) {
    serialized += realNumerator.str();

    if (realDenominator != 1) {
      serialized += "/" + realDenominator.str();
    }
  }
  if (imagNumerator != 0) {
    if (!serialized.empty()) {
      if (imagNumerator < 0) {
        serialized += " - i ";
        imagNumerator *= -1;
      } else {
        serialized += " + i ";
      }
    }

    serialized += imagNumerator.str();

    if (imagDenominator != 1) {
      serialized += "/" + imagDenominator.str();
    }
  }

  SEQUANT_ASSERT(!serialized.empty());

  return serialized;
}

std::string to_string(Tensor const& tensor,
                      const SerializationOptions& options) {
  return to_string(static_cast<const AbstractTensor&>(tensor), options);
}

std::string to_string(const Constant& constant,
                      const SerializationOptions& options) {
  return detail::serialize_scalar(constant.value(), options);
}

std::string to_string(const Variable& variable, const SerializationOptions&) {
  return std::string(variable.label()) + (variable.conjugated() ? "^*" : "");
}

std::string to_string(const Power& power, const SerializationOptions& options) {
  std::string core;
  const auto& base = power.base();
  // parenthesize bases for Constants and conjugated Variables
  const bool parenthesize_base =
      base->is<Constant>() ||
      (base->is<Variable>() && base->as<Variable>().conjugated());
  if (parenthesize_base) core += "(";
  core += v1::to_string(*base, options);
  if (parenthesize_base) core += ")";
  core += "^(";
  core += serialize_scalar(Constant::scalar_type{power.exponent()}, options);
  core += ")";
  if (power.conjugated()) return "(" + std::move(core) + ")^*";
  return core;
}

std::string to_string(Product const& prod,
                      const SerializationOptions& options) {
  std::string serialized;

  const auto& scal = prod.scalar();
  if (scal != Product::scalar_type{1}) {
    serialized += detail::serialize_scalar(scal, options) + " ";
  }

  for (std::size_t i = 0; i < prod.size(); ++i) {
    const ExprPtr& current = prod[i];
    bool parenthesize = false;
    if (current->is<Product>() || current->is<Sum>()) {
      parenthesize = true;
      serialized += "(";
    }

    serialized += to_string(current, options);

    if (parenthesize) {
      serialized += ")";
    }

    if (i + 1 < prod.size()) {
      serialized += " * ";
    }
  }

  return serialized;
}

std::string to_string(Sum const& sum, const SerializationOptions& options) {
  std::string serialized;

  for (std::size_t i = 0; i < sum.size(); ++i) {
    const ExprPtr& current = sum[i];

    const bool parenthesize = current->is<Sum>();

    std::string current_serialized = to_string(current, options);

    bool is_negative = false;
    if (parenthesize) {
      current_serialized = "(" + current_serialized + ")";
    } else {
      is_negative = current_serialized.front() == '-';
    }

    if (i > 0) {
      if (is_negative) {
        serialized += " - " + current_serialized.substr(1);
      } else {
        serialized += " + " + current_serialized;
      }
    } else {
      serialized = std::move(current_serialized);
    }
  }
  return serialized;
}

}  // namespace detail

std::string to_string(const ExprPtr& expr,
                      const SerializationOptions& options) {
  if (!expr) return {};

  return v1::to_string(*expr, options);
}

std::string to_string(const Expr& expr, const SerializationOptions& options) {
  using namespace detail;
  if (expr.is<Tensor>())
    return detail::to_string(expr.as<Tensor>(), options);
  else if (expr.is<FNOperator>())
    return v1::to_string(expr.as<FNOperator>(), options);
  else if (expr.is<BNOperator>())
    return v1::to_string(expr.as<BNOperator>(), options);
  else if (expr.is<Sum>())
    return detail::to_string(expr.as<Sum>(), options);
  else if (expr.is<Product>())
    return detail::to_string(expr.as<Product>(), options);
  else if (expr.is<Constant>())
    return detail::to_string(expr.as<Constant>(), options);
  else if (expr.is<Variable>())
    return detail::to_string(expr.as<Variable>(), options);
  else if (expr.is<Power>())
    return detail::to_string(expr.as<Power>(), options);
  else
    throw Exception(
        expr.type_name() +
        " cannot be serialized; use LaTeX (to_latex or {:l}) if the "
        "type provides it");
}

std::string to_string(const ResultExpr& result,
                      const SerializationOptions& options) {
  std::string serialized;
  if (result.produces_tensor()) {
    serialized = detail::to_string(result.result_as_tensor("?"), options);
  } else {
    serialized = detail::to_string(result.result_as_variable("?"), options);
  }

  return serialized + " = " + v1::to_string(result.expression(), options);
}

std::string to_string(const Index& index, const SerializationOptions&) {
  std::string serialized(index.label());

  const bool unnamed_instance = index.basis().unnamed_instance().has_value();
  if (index.has_proto_indices() || unnamed_instance) {
    serialized += "<";
    serialized += join_strings(
        index.proto_indices(), ",", [](const Index& proto) -> std::string {
          std::string label(proto.label());
          if (proto.basis().unnamed_instance())
            label += "<" + proto.basis().instance_suffix() + ">";
          return label;
        });
    if (unnamed_instance) serialized += index.basis().instance_suffix();
    serialized += ">";
  }

  return serialized;
}

std::string to_string(AbstractTensor const& tensor,
                      const SerializationOptions& options) {
  std::string serialized(tensor._label());
  serialized += "{" + detail::serialize_indices(tensor._bra(), options);
  if (tensor._ket_rank() > 0) {
    serialized += ";" + detail::serialize_indices(tensor._ket(), options);
  }
  if (tensor._aux_rank() > 0) {
    if (tensor._ket_rank() == 0) {
      serialized += ";";
    }
    serialized += ";" + detail::serialize_indices(tensor._aux(), options);
  }
  serialized += "}";

  if (options.annot_symm) {
    serialized += ":" + detail::serialize_symm(tensor._symmetry(), options);
    serialized +=
        "-" + detail::serialize_symm(tensor._braket_symmetry(), options);
    serialized +=
        "-" + detail::serialize_symm(tensor._column_symmetry(), options);
  }

  return serialized;
}

std::string to_string(Tensor const& tensor,
                      const SerializationOptions& options) {
  return v1::to_string(static_cast<const Expr&>(tensor), options);
}

template <Statistics S>
std::string to_string(NormalOperator<S> const& nop,
                      const SerializationOptions& options) {
  std::string serialized(nop.label());
  serialized += "{" + detail::serialize_ops(nop.annihilators(), options);
  if (nop.ncreators() > 0) {
    serialized += ";" + detail::serialize_ops(nop.creators(), options);
  }
  serialized += "}";

  return serialized;
}

template std::string to_string<Statistics::FermiDirac>(
    NormalOperator<Statistics::FermiDirac> const& nop,
    const SerializationOptions& options);
template std::string to_string<Statistics::BoseEinstein>(
    NormalOperator<Statistics::BoseEinstein> const& nop,
    const SerializationOptions& options);

}  // namespace sequant::io::serialization::v1
