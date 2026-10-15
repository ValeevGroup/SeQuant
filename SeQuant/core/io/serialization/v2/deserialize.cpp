#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/io/serialization/detail/construct.hpp>
#include <SeQuant/core/io/serialization/function_registry.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/serialization/v2/error.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/algorithm/find.hpp>

#include <peglib.h>

#include <any>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant::io::serialization::v2 {

// Defines a variable peg_serialization_grammar that contains the definition of
// the grammar (as a string_view)
#include "peg_grammar.ipp"

namespace {

/// a parsed index group list; a std::nullopt entry is a blank index
using IndexGroupsValue = std::vector<std::vector<std::optional<Index>>>;

/// @return the first value of type T stored in @p vs, or nullptr
template <typename T>
const T *find_value(const peg::SemanticValues &vs) {
  for (const std::any &value : vs) {
    if (value.type() == typeid(T)) return std::any_cast<T>(&value);
  }
  return nullptr;
}

template <typename T>
bool has_value(const peg::SemanticValues &vs) {
  return find_value<T>(vs) != nullptr;
}

const DeserializationOptions &options_of(std::any &dt) {
  return *std::any_cast<const DeserializationOptions *>(dt);
}

container::vector<Index> to_index_vector(
    const std::vector<std::optional<Index>> &group) {
  container::vector<Index> result;
  result.reserve(group.size());
  for (const std::optional<Index> &slot : group) {
    result.push_back(slot ? *slot : Index{});
  }
  return result;
}

ExprPtr negate(ExprPtr expr) {
  if (expr->is<Constant>()) {
    return ex<Constant>(Constant(-expr->as<Constant>().value()));
  }
  // fold the sign into the prefactor, which is where the serializer puts it
  if (expr->is<Product>()) {
    ExprPtr product = expr->clone();
    product->as<Product>().scale(-1);
    return product;
  }
  return ex<Product>(Constant::scalar_type(-1),
                     std::vector<ExprPtr>{std::move(expr)},
                     Product::Flatten::No);
}

ExprPtr build_normal_operator(const std::wstring &name,
                              const IndexGroupsValue &groups) {
  container::vector<Index> annihilators = groups.size() > 0
                                              ? to_index_vector(groups[0])
                                              : container::vector<Index>{};
  container::vector<Index> creators = groups.size() > 1
                                          ? to_index_vector(groups[1])
                                          : container::vector<Index>{};
  if (groups.size() > 2) {
    throw Exception(
        "a normal-ordered operator takes at most two index groups "
        "(annihilators;creators)");
  }

  const auto &flabels = NormalOperator<Statistics::FermiDirac>::labels();
  if (auto it = ranges::find(flabels, name); it != flabels.end()) {
    const Vacuum vac = it == flabels.begin()
                           ? Vacuum::Physical
                           : detail::tilde_vacuum(Statistics::FermiDirac);
    return ex<NormalOperator<Statistics::FermiDirac>>(
        cre(std::move(creators)), ann(std::move(annihilators)), vac);
  }
  const auto &blabels = NormalOperator<Statistics::BoseEinstein>::labels();
  if (auto it = ranges::find(blabels, name); it != blabels.end()) {
    const Vacuum vac = it == blabels.begin()
                           ? Vacuum::Physical
                           : detail::tilde_vacuum(Statistics::BoseEinstein);
    return ex<NormalOperator<Statistics::BoseEinstein>>(
        cre(std::move(creators)), ann(std::move(annihilators)), vac);
  }
  throw Exception("'" + toUtf8(name) +
                  "' is not a recognized normal-ordered operator label");
}

ExprPtr build_tensor(const std::wstring &name, const IndexGroupsValue &groups,
                     const std::vector<std::string> &symmetries,
                     const DeserializationOptions &options) {
  container::vector<Index> bra = groups.size() > 0 ? to_index_vector(groups[0])
                                                   : container::vector<Index>{};
  container::vector<Index> ket = groups.size() > 1 ? to_index_vector(groups[1])
                                                   : container::vector<Index>{};
  container::vector<Index> aux = groups.size() > 2 ? to_index_vector(groups[2])
                                                   : container::vector<Index>{};
  if (groups.size() > 3) {
    throw Exception("a tensor takes at most three index groups (bra;ket;aux)");
  }

  const auto [def_perm, def_braket, def_column] =
      detail::resolve_default_symmetries(options);
  Symmetry perm = def_perm;
  std::variant<BraKetSymmetry, Hermiticity> braket = def_braket;
  ColumnSymmetry column = def_column;
  bool perm_spec = false, braket_spec = false, column_spec = false;

  // records that a field is given, which it may only be once
  auto specify = [](bool &specified, const char *field) {
    if (specified)
      throw Exception(std::string("the ") + field +
                      " symmetry is specified more than once");
    specified = true;
  };
  for (const std::string &token : symmetries) {
    if (token == "A") {
      specify(perm_spec, "permutational"), perm = Symmetry::Antisymm;
    } else if (token == "S") {
      specify(perm_spec, "permutational"), perm = Symmetry::Symm;
    } else if (token == "N") {
      specify(perm_spec, "permutational"), perm = Symmetry::Nonsymm;
    } else if (token == "bkS") {
      specify(braket_spec, "bra-ket"), braket = BraKetSymmetry::Symm;
    } else if (token == "bkC") {
      specify(braket_spec, "bra-ket"), braket = BraKetSymmetry::Conjugate;
    } else if (token == "bkN") {
      specify(braket_spec, "bra-ket"), braket = BraKetSymmetry::Nonsymm;
    } else if (token == "pS") {
      specify(column_spec, "particle"), column = ColumnSymmetry::Symm;
    } else if (token == "pN") {
      specify(column_spec, "particle"), column = ColumnSymmetry::Nonsymm;
    } else {
      throw Exception("invalid symmetry specifier '" + token + "'");
    }
  }

  return detail::make_tensor(name, std::move(bra), std::move(ket),
                             std::move(aux), perm, braket, column, perm_spec,
                             braket_spec, column_spec);
}

/// Installs \p action for \p rule such that an Exception it throws (e.g. for
/// an unknown index space) is reported as an Error at the rule's match
template <typename F>
void on(peg::parser &parser, const char *rule, F action) {
  parser[rule] = [rule, action = std::move(action)](
                     const peg::SemanticValues &vs, std::any &dt) {
    try {
      if constexpr (std::is_invocable_v<const F &, const peg::SemanticValues &,
                                        std::any &>)
        return action(vs, dt);
      else
        return action(vs);
    } catch (const Error &) {
      throw;
    } catch (const Exception &e) {
      const auto [line, column] = vs.line_info();
      throw Error(line, column, rule, e.what());
    }
  };
}

void register_actions(peg::parser &parser) {
  using peg::SemanticValues;

  on(parser, "BareIdentifier",
     [](const SemanticValues &vs) { return std::string(vs.token()); });
  on(parser, "QuotedIdentifier",
     [](const SemanticValues &vs) { return toUtf8(parse_label(vs.token())); });
  on(parser, "Identifier", [](const SemanticValues &vs) {
    return std::any_cast<std::string>(vs[0]);
  });
  on(parser, "Real",
     [](const SemanticValues &vs) { return string_to<rational>(vs.token()); });
  on(parser, "Imaginary", [](const SemanticValues &vs) {
    return Constant::scalar_type(0, std::any_cast<rational>(vs[0]));
  });
  on(parser, "Complex", [](const SemanticValues &vs) {
    return vs.choice() == 0
               ? std::any_cast<Constant::scalar_type>(vs[0])
               : Constant::scalar_type(std::any_cast<rational>(vs[0]), 0);
  });
  on(parser, "Constant", [](const SemanticValues &vs) {
    return ExprPtr(ex<Constant>(std::any_cast<Constant::scalar_type>(vs[0])));
  });

  on(parser, "Conjugation", [](const SemanticValues &) { return true; });
  on(parser, "Sign", [](const SemanticValues &vs) { return vs.token()[0]; });
  on(parser, "MulOp", [](const SemanticValues &vs) { return vs.token()[0]; });
  on(parser, "AddOp", [](const SemanticValues &vs) { return vs.token()[0]; });
  on(parser, "Symmetry",
     [](const SemanticValues &vs) { return std::string(vs.token()); });
  on(parser, "SymmetrySpec", [](const SemanticValues &vs) {
    std::vector<std::string> result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::string>(vs[i]));
    return result;
  });
  on(parser, "BasisInstance", [](const SemanticValues &vs) {
    return string_to<IndexBasis::instance_type>(vs.token());
  });
  on(parser, "ProtoIndex", [](const SemanticValues &vs) {
    Index proto(toUtf16(*find_value<std::string>(vs)));
    if (auto instance = find_value<IndexBasis::instance_type>(vs))
      proto = proto.replace_basis_instance(*instance);
    return proto;
  });
  on(parser, "ProtoList", [](const SemanticValues &vs) {
    container::vector<Index> result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<Index>(vs[i]));
    return result;
  });
  on(parser, "IndexDomain", [](const SemanticValues &vs) {
    container::vector<Index> protos;
    if (auto p = find_value<container::vector<Index>>(vs)) protos = *p;
    return protos;
  });
  on(parser, "Index", [](const SemanticValues &vs) {
    const std::wstring label = toUtf16(*find_value<std::string>(vs));
    const auto *protos = find_value<container::vector<Index>>(vs);
    Index index =
        protos && !protos->empty() ? Index(label, *protos) : Index(label);
    if (auto instance = find_value<IndexBasis::instance_type>(vs))
      index = index.replace_basis_instance(*instance);
    return index;
  });
  on(parser, "IndexItem", [](const SemanticValues &vs) {
    return vs.choice() == 1 ? std::optional<Index>(std::any_cast<Index>(vs[0]))
                            : std::optional<Index>{};
  });
  on(parser, "Indices", [](const SemanticValues &vs) {
    std::vector<std::optional<Index>> result;
    // a group without any text (t[] / t[;a1]) is empty rather than one blank
    if (vs.sv().find_first_not_of(" \t\r\n") == std::string_view::npos)
      return result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::optional<Index>>(vs[i]));
    return result;
  });
  on(parser, "IndexGroups", [](const SemanticValues &vs) {
    IndexGroupsValue result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::vector<std::optional<Index>>>(vs[i]));
    return result;
  });

  on(parser, "Variable", [](const SemanticValues &vs) {
    Variable var(toUtf16(*find_value<std::string>(vs)));
    if (has_value<bool>(vs)) var.conjugate();
    return ExprPtr(ex<Variable>(std::move(var)));
  });

  on(parser, "Tensor", [](const SemanticValues &vs, std::any &dt) -> ExprPtr {
    const std::wstring name = toUtf16(*find_value<std::string>(vs));
    const IndexGroupsValue &groups = *find_value<IndexGroupsValue>(vs);
    std::vector<std::string> symmetries;
    if (auto s = find_value<std::vector<std::string>>(vs)) symmetries = *s;
    return build_tensor(name, groups, symmetries, options_of(dt));
  });

  on(parser, "NormalOrderedOperator", [](const SemanticValues &vs) {
    const std::wstring name = toUtf16(*find_value<std::string>(vs));
    const IndexGroupsValue &groups = *find_value<IndexGroupsValue>(vs);
    return build_normal_operator(name, groups);
  });

  on(parser, "PowerBase",
     [](const SemanticValues &vs) { return std::any_cast<ExprPtr>(vs[0]); });
  on(parser, "Exponent",
     [](const SemanticValues &vs) { return string_to<rational>(vs.token()); });
  on(parser, "Power", [](const SemanticValues &vs) {
    ExprPtr power = vs.choice() == 0 ? ex<Power>(*find_value<ExprPtr>(vs),
                                                 *find_value<rational>(vs))
                                     : find_value<ExprPtr>(vs)->clone();
    if (has_value<bool>(vs)) power->as<Power>().conjugate();
    return power;
  });

  on(parser, "RawArguments",
     [](const SemanticValues &vs) { return std::string(vs.token()); });
  // the handler's errors (including those of nested parses of its arguments,
  // whose positions are relative to the argument) are reported at the call
  parser["Function"] = [](const SemanticValues &vs, std::any &dt) {
    FunctionCall call;
    call.name = toUtf16(*find_value<std::string>(vs));
    // find_value would also match the name, which comes first
    call.arguments = std::any_cast<std::string>(vs[1]);
    try {
      return FunctionRegistry::instance().deserialize(call, options_of(dt));
    } catch (const Exception &e) {
      const auto [line, column] = vs.line_info();
      throw Error(line, column, "Function",
                  "~" + toUtf8(call.name) + "(...): " + e.what());
    }
  };

  on(parser, "Nullary",
     [](const SemanticValues &vs) { return std::any_cast<ExprPtr>(vs[0]); });
  on(parser, "Factor", [](const SemanticValues &vs) {
    ExprPtr expr = std::any_cast<ExprPtr>(vs.back());
    std::size_t minus = 0;
    for (std::size_t i = 0; i + 1 < vs.size(); ++i)
      if (std::any_cast<char>(vs[i]) == '-') ++minus;
    return (minus % 2 == 1) ? negate(std::move(expr)) : expr;
  });
  on(parser, "ProductTail", [](const SemanticValues &vs) {
    return vs.choice() == 0
               ? std::pair<char, ExprPtr>(std::any_cast<char>(vs[0]),
                                          std::any_cast<ExprPtr>(vs[1]))
               : std::pair<char, ExprPtr>('*', std::any_cast<ExprPtr>(vs[0]));
  });
  on(parser, "Product", [](const SemanticValues &vs) -> ExprPtr {
    Constant::scalar_type prefactor(1);
    bool explicit_prefactor = false;
    std::vector<ExprPtr> factors;

    auto apply = [&](char op, ExprPtr operand) {
      if (operand->is<Constant>()) {
        explicit_prefactor = true;
        if (op == '/')
          prefactor = prefactor / operand->as<Constant>().value();
        else
          prefactor *= operand->as<Constant>().value();
      } else if (op == '/') {
        factors.push_back(ex<Power>(std::move(operand), rational(-1)));
      } else {
        factors.push_back(std::move(operand));
      }
    };

    apply('*', std::any_cast<ExprPtr>(vs[0]));
    for (std::size_t i = 1; i < vs.size(); ++i) {
      auto [op, operand] = std::any_cast<std::pair<char, ExprPtr>>(vs[i]);
      apply(op, std::move(operand));
    }

    if (factors.empty()) return ex<Constant>(Constant(prefactor));
    // a spelled-out prefactor (even 1) makes the factor a Product, as written
    if (factors.size() == 1 && !explicit_prefactor) return factors.front();
    return ex<Product>(prefactor, std::move(factors), Product::Flatten::No);
  });
  on(parser, "Sum", [](const SemanticValues &vs) -> ExprPtr {
    if (vs.size() == 1) return std::any_cast<ExprPtr>(vs[0]);
    std::vector<ExprPtr> summands;
    summands.push_back(std::any_cast<ExprPtr>(vs[0]));
    for (std::size_t i = 1; i + 1 < vs.size(); i += 2) {
      const char op = std::any_cast<char>(vs[i]);
      ExprPtr operand = std::any_cast<ExprPtr>(vs[i + 1]);
      summands.push_back(op == '-' ? negate(std::move(operand))
                                   : std::move(operand));
    }
    ExprPtr sum = ex<Sum>(std::move(summands));
    // Sum adds up constant summands; one that folded entirely is a Constant,
    // e.g. the mixed complex "3 + 2i", which has no single-token form
    if (sum->size() == 1 && sum->as<Sum>().summand(0)->is<Constant>())
      return sum->as<Sum>().summand(0);
    return sum;
  });

  on(parser, "Result",
     [](const SemanticValues &vs) { return std::any_cast<ExprPtr>(vs[0]); });
  on(parser, "ResultExpression", [](const SemanticValues &vs) {
    ExprPtr lhs = std::any_cast<ExprPtr>(vs[0]);
    ExprPtr rhs = std::any_cast<ExprPtr>(vs[1]);
    if (lhs->is<Tensor>())
      return std::make_shared<ResultExpr>(lhs->as<Tensor>(), rhs);
    return std::make_shared<ResultExpr>(lhs->as<Variable>(), rhs);
  });
}

std::string describe_error(const peg::ErrorReport &report) {
  if (!report.message.empty()) return report.message;

  std::string msg = "syntax error";
  if (!report.unexpected_token.empty())
    msg += " at '" + report.unexpected_token + "'";

  std::vector<std::string> expected = report.expected_literals;
  expected.insert(expected.end(), report.expected_rules.begin(),
                  report.expected_rules.end());
  if (!expected.empty()) {
    msg += ", expected ";
    for (std::size_t i = 0; i < expected.size(); ++i) {
      if (i != 0) msg += i + 1 == expected.size() ? " or " : ", ";
      msg += expected[i];
    }
  }
  return msg;
}

peg::parser &get_parser(const char *start_rule) {
  static thread_local std::map<std::string, std::shared_ptr<peg::parser>>
      parsers;

  if (auto it = parsers.find(start_rule); it != parsers.end())
    return *it->second;

  auto parser = std::make_shared<peg::parser>();
  // grammar-definition errors surface through the logger during load; the
  // "rule is not referenced" notices for the unused top-level rules (e.g.
  // ResultExpression when a leaf start rule is used) are expected, not errors
  parser->set_logger([](std::size_t line, std::size_t column,
                        const std::string &msg, const std::string &rule) {
    if (msg.find("is not referenced") != std::string::npos) return;
    throw Error(line, column, rule,
                "internal grammar error: " + msg + " (" + rule + ")");
  });
  const bool ok = parser->load_grammar(peg_serialization_grammar, start_rule);
  SEQUANT_ASSERT(ok);
  (void)ok;

  parser->set_error_reporter([](const peg::ErrorReport &report) {
    throw Error(report.line, report.col, report.label, describe_error(report));
  });
  register_actions(*parser);

  auto [it, _] = parsers.emplace(start_rule, std::move(parser));
  return *it->second;
}

// Parses \p input beginning at \p start_rule, which must yield a value of the
// default-constructible type \p T.
template <typename T>
T parse_value(const char *start_rule, std::string_view input,
              const DeserializationOptions &options) {
  peg::parser &parser = get_parser(start_rule);
  std::any dt = &options;
  T value{};
  const bool ok = parser.parse(input, dt, value);
  if (!ok) {
    // the error reporter normally throws first; this is the fallback
    throw Error(0, 0, start_rule, "failed to parse input");
  }
  return value;
}

ExprPtr parse_expr(const char *start_rule, std::string_view input,
                   const DeserializationOptions &options) {
  return parse_value<ExprPtr>(start_rule, input, options);
}

// Parses \p input beginning at \p start_rule and requires the result to be a
// \p T, called \p name in the error otherwise
template <typename T>
T parse_as(const char *start_rule, std::string_view input,
           const DeserializationOptions &options, const char *name) {
  ExprPtr expr = parse_expr(start_rule, input, options);
  if (!expr->is<T>())
    throw Error(0, 0, start_rule,
                std::string("input does not denote a ") + name);
  return expr->as<T>();
}

}  // namespace

// All grammar rules that back a leaf type (Constant, Variable, Tensor) produce
// an ExprPtr; the typed entry points extract the concrete object from it.
#define SEQUANT_DESERIALIZE_EXPR_AS(Type, Rule)                   \
  template <>                                                     \
  Type from_string<Type>(std::string_view input,                  \
                         const DeserializationOptions &options) { \
    return parse_as<Type>(Rule, input, options, #Type);           \
  }                                                               \
  template <>                                                     \
  Type from_string<Type>(std::wstring_view input,                 \
                         const DeserializationOptions &options) { \
    return v2::from_string<Type>(toUtf8(input), options);         \
  }

SEQUANT_DESERIALIZE_EXPR_AS(Variable, "Variable")
SEQUANT_DESERIALIZE_EXPR_AS(Tensor, "Tensor")
SEQUANT_DESERIALIZE_EXPR_AS(Product, "Expression")
SEQUANT_DESERIALIZE_EXPR_AS(Sum, "Expression")

// A Constant is parsed as a full expression, so that an arithmetic combination
// of numbers -- e.g. the fraction "1/2" or the mixed complex "3 + 2i", neither
// of which is a single Constant token -- still yields one.
template <>
Constant from_string<Constant>(std::string_view input,
                               const DeserializationOptions &options) {
  return parse_as<Constant>("Expression", input, options, "Constant");
}
template <>
Constant from_string<Constant>(std::wstring_view input,
                               const DeserializationOptions &options) {
  return v2::from_string<Constant>(toUtf8(input), options);
}

template <>
Index from_string<Index>(std::string_view input,
                         const DeserializationOptions &options) {
  return parse_value<Index>("Index", input, options);
}
template <>
Index from_string<Index>(std::wstring_view input,
                         const DeserializationOptions &options) {
  return v2::from_string<Index>(toUtf8(input), options);
}

template <>
ExprPtr from_string<ExprPtr>(std::string_view input,
                             const DeserializationOptions &options) {
  return parse_expr("Expression", input, options);
}
template <>
ExprPtr from_string<ExprPtr>(std::wstring_view input,
                             const DeserializationOptions &options) {
  return v2::from_string<ExprPtr>(toUtf8(input), options);
}

template <>
ResultExpr from_string<ResultExpr>(std::string_view input,
                                   const DeserializationOptions &options) {
  return *parse_value<std::shared_ptr<ResultExpr>>("ResultExpression", input,
                                                   options);
}
template <>
ResultExpr from_string<ResultExpr>(std::wstring_view input,
                                   const DeserializationOptions &options) {
  return v2::from_string<ResultExpr>(toUtf8(input), options);
}

}  // namespace sequant::io::serialization::v2
