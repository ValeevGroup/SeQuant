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

/// a parsed index group list; a std::nullopt entry is a positionally blank
/// index, and a group that is a single blank is an empty group
using IndexGroupsValue = std::vector<std::vector<std::optional<Index>>>;

/// an index as matched, before it is validated against the index registry; the
/// deferred construction lets a speculative ActionIndex match fail (and the
/// parser backtrack) without a registry lookup throwing
struct RawIndex {
  std::wstring label;
  std::optional<IndexBasis::instance_type> instance;
  container::vector<Index> protos;
};

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
  // a single blank slot is how an empty group is spelled (t[] / t[;])
  if (group.size() == 1 && !group.front().has_value()) return {};

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

Index make_index(const RawIndex &raw) {
  Index index =
      raw.protos.empty() ? Index(raw.label) : Index(raw.label, raw.protos);
  if (raw.instance) index = index.replace_basis_instance(*raw.instance);
  return index;
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
    throw Error(0, 0, "NormalOrderedOperator",
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
  throw Error(0, 0, "NormalOrderedOperator",
              "'" + toUtf8(name) +
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
    throw Error(0, 0, "Tensor",
                "a tensor takes at most three index groups (bra;ket;aux)");
  }

  const auto [def_perm, def_braket, def_column] =
      detail::resolve_default_symmetries(options);
  Symmetry perm = def_perm;
  std::variant<BraKetSymmetry, Hermiticity> braket = def_braket;
  ColumnSymmetry column = def_column;
  bool perm_spec = false, braket_spec = false, column_spec = false;

  for (const std::string &token : symmetries) {
    if (token == "A") {
      perm = Symmetry::Antisymm, perm_spec = true;
    } else if (token == "S") {
      perm = Symmetry::Symm, perm_spec = true;
    } else if (token == "N") {
      perm = Symmetry::Nonsymm, perm_spec = true;
    } else if (token == "bkS") {
      braket = BraKetSymmetry::Symm, braket_spec = true;
    } else if (token == "bkC") {
      braket = BraKetSymmetry::Conjugate, braket_spec = true;
    } else if (token == "bkN") {
      braket = BraKetSymmetry::Nonsymm, braket_spec = true;
    } else if (token == "pS") {
      column = ColumnSymmetry::Symm, column_spec = true;
    } else if (token == "pN") {
      column = ColumnSymmetry::Nonsymm, column_spec = true;
    } else {
      throw Error(0, 0, "Symmetry",
                  "invalid symmetry specifier '" + token + "'");
    }
  }

  return detail::make_tensor(name, std::move(bra), std::move(ket),
                             std::move(aux), perm, braket, column, perm_spec,
                             braket_spec, column_spec);
}

void register_actions(peg::parser &parser) {
  using peg::SemanticValues;

  parser["Identifier"] = [](const SemanticValues &vs) {
    return std::string(vs.token());
  };
  parser["Real"] = [](const SemanticValues &vs) {
    return string_to<rational>(vs.token());
  };
  parser["Imaginary"] = [](const SemanticValues &vs) {
    return Constant::scalar_type(0, std::any_cast<rational>(vs[0]));
  };
  parser["Complex"] = [](const SemanticValues &vs) {
    return vs.choice() == 0
               ? std::any_cast<Constant::scalar_type>(vs[0])
               : Constant::scalar_type(std::any_cast<rational>(vs[0]), 0);
  };
  parser["Constant"] = [](const SemanticValues &vs) {
    return ExprPtr(ex<Constant>(std::any_cast<Constant::scalar_type>(vs[0])));
  };

  parser["Conjugation"] = [](const SemanticValues &) { return true; };
  parser["Sign"] = [](const SemanticValues &vs) { return vs.token()[0]; };
  parser["MulOp"] = [](const SemanticValues &vs) { return vs.token()[0]; };
  parser["AddOp"] = [](const SemanticValues &vs) { return vs.token()[0]; };
  parser["Symmetry"] = [](const SemanticValues &vs) {
    return std::string(vs.token());
  };
  parser["SymmetrySpec"] = [](const SemanticValues &vs) {
    std::vector<std::string> result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::string>(vs[i]));
    return result;
  };
  parser["Statistic"] = [](const SemanticValues &vs) {
    return vs.choice() == 0 ? Statistics::FermiDirac : Statistics::BoseEinstein;
  };
  parser["CreatorMark"] = [](const SemanticValues &) { return Action::Create; };
  parser["AnnihilatorMark"] = [](const SemanticValues &) {
    return Action::Annihilate;
  };

  parser["BasisInstance"] = [](const SemanticValues &vs) {
    return string_to<IndexBasis::instance_type>(vs.token());
  };
  parser["ProtoIndex"] = [](const SemanticValues &vs) {
    Index proto(toUtf16(*find_value<std::string>(vs)));
    if (auto instance = find_value<IndexBasis::instance_type>(vs))
      proto = proto.replace_basis_instance(*instance);
    return proto;
  };
  parser["ProtoList"] = [](const SemanticValues &vs) {
    container::vector<Index> result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<Index>(vs[i]));
    return result;
  };
  parser["IndexDomain"] = [](const SemanticValues &vs) {
    container::vector<Index> protos;
    if (auto p = find_value<container::vector<Index>>(vs)) protos = *p;
    return protos;
  };
  parser["Index"] = [](const SemanticValues &vs) {
    RawIndex raw;
    raw.label = toUtf16(*find_value<std::string>(vs));
    if (auto instance = find_value<IndexBasis::instance_type>(vs))
      raw.instance = *instance;
    if (auto protos = find_value<container::vector<Index>>(vs))
      raw.protos = *protos;
    return raw;
  };
  parser["IndexItem"] = [](const SemanticValues &vs) {
    return vs.choice() == 0 ? std::optional<Index>(
                                  make_index(std::any_cast<RawIndex>(vs[0])))
                            : std::optional<Index>{};
  };
  parser["Indices"] = [](const SemanticValues &vs) {
    std::vector<std::optional<Index>> result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::optional<Index>>(vs[i]));
    return result;
  };
  parser["IndexGroups"] = [](const SemanticValues &vs) {
    IndexGroupsValue result;
    for (std::size_t i = 0; i < vs.size(); ++i)
      result.push_back(std::any_cast<std::vector<std::optional<Index>>>(vs[i]));
    return result;
  };

  parser["Variable"] = [](const SemanticValues &vs) {
    Variable var(toUtf16(*find_value<std::string>(vs)));
    if (has_value<bool>(vs)) var.conjugate();
    return ExprPtr(ex<Variable>(std::move(var)));
  };

  parser["Tensor"] = [](const SemanticValues &vs, std::any &dt) -> ExprPtr {
    const std::wstring name = toUtf16(*find_value<std::string>(vs));
    const IndexGroupsValue &groups = *find_value<IndexGroupsValue>(vs);
    std::vector<std::string> symmetries;
    if (auto s = find_value<std::vector<std::string>>(vs)) symmetries = *s;
    return build_tensor(name, groups, symmetries, options_of(dt));
  };

  parser["NormalOrderedOperator"] = [](const SemanticValues &vs) {
    const std::wstring name = toUtf16(*find_value<std::string>(vs));
    const IndexGroupsValue &groups = *find_value<IndexGroupsValue>(vs);
    return build_normal_operator(name, groups);
  };

  parser["PowerBase"] = [](const SemanticValues &vs) {
    return std::any_cast<ExprPtr>(vs[0]);
  };
  parser["Exponent"] = [](const SemanticValues &vs) {
    return string_to<rational>(vs.token());
  };
  parser["Power"] = [](const SemanticValues &vs) {
    ExprPtr power = vs.choice() == 0 ? ex<Power>(*find_value<ExprPtr>(vs),
                                                 *find_value<rational>(vs))
                                     : find_value<ExprPtr>(vs)->clone();
    if (has_value<bool>(vs)) power->as<Power>().conjugate();
    return power;
  };

  parser["ActionIndex"] = [](const SemanticValues &vs) {
    return std::pair<Index, Action>(make_index(std::any_cast<RawIndex>(vs[0])),
                                    std::any_cast<Action>(vs[1]));
  };
  parser["FunctionArg"] = [](const SemanticValues &vs) {
    FunctionArgument arg;
    if (vs.choice() == 0)
      arg.action_index = std::any_cast<std::pair<Index, Action>>(vs[0]);
    else
      arg.expr = std::any_cast<ExprPtr>(vs[0]);
    return arg;
  };
  parser["ArgList"] = [](const SemanticValues &vs) {
    std::vector<FunctionArgument> args;
    for (std::size_t i = 0; i < vs.size(); ++i)
      args.push_back(std::any_cast<FunctionArgument>(vs[i]));
    return args;
  };
  parser["Function"] = [](const SemanticValues &vs, std::any &dt) {
    FunctionCall call;
    call.name = toUtf16(*find_value<std::string>(vs));
    if (auto args = find_value<std::vector<FunctionArgument>>(vs))
      call.arguments = *args;
    if (auto stat = find_value<Statistics>(vs)) call.statistic = *stat;
    return FunctionRegistry::instance().deserialize(call, options_of(dt));
  };

  parser["Nullary"] = [](const SemanticValues &vs) {
    return std::any_cast<ExprPtr>(vs[0]);
  };
  parser["Factor"] = [](const SemanticValues &vs) {
    ExprPtr expr = std::any_cast<ExprPtr>(vs.back());
    std::size_t minus = 0;
    for (std::size_t i = 0; i + 1 < vs.size(); ++i)
      if (std::any_cast<char>(vs[i]) == '-') ++minus;
    return (minus % 2 == 1) ? negate(std::move(expr)) : expr;
  };
  parser["ProductTail"] = [](const SemanticValues &vs) {
    return vs.choice() == 0
               ? std::pair<char, ExprPtr>(std::any_cast<char>(vs[0]),
                                          std::any_cast<ExprPtr>(vs[1]))
               : std::pair<char, ExprPtr>('*', std::any_cast<ExprPtr>(vs[0]));
  };
  parser["Product"] = [](const SemanticValues &vs) -> ExprPtr {
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
  };
  parser["Sum"] = [](const SemanticValues &vs) -> ExprPtr {
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
  };

  parser["Result"] = [](const SemanticValues &vs) {
    return std::any_cast<ExprPtr>(vs[0]);
  };
  parser["ResultExpression"] = [](const SemanticValues &vs) {
    ExprPtr lhs = std::any_cast<ExprPtr>(vs[0]);
    ExprPtr rhs = std::any_cast<ExprPtr>(vs[1]);
    if (lhs->is<Tensor>())
      return std::make_shared<ResultExpr>(lhs->as<Tensor>(), rhs);
    return std::make_shared<ResultExpr>(lhs->as<Variable>(), rhs);
  };
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

}  // namespace

// All grammar rules that back a leaf type (Constant, Variable, Tensor) produce
// an ExprPtr; the typed entry points extract the concrete object from it.
#define SEQUANT_DESERIALIZE_EXPR_AS(Type, Rule)                   \
  template <>                                                     \
  Type from_string<Type>(std::string_view input,                  \
                         const DeserializationOptions &options) { \
    return parse_expr(Rule, input, options)->as<Type>();          \
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
  ExprPtr expr = parse_expr("Expression", input, options);
  if (!expr->is<Constant>())
    throw Error(0, 0, "Expression", "input does not denote a constant");
  return expr->as<Constant>();
}
template <>
Constant from_string<Constant>(std::wstring_view input,
                               const DeserializationOptions &options) {
  return v2::from_string<Constant>(toUtf8(input), options);
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
