#include <SeQuant/core/irreps/config.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/irreps/propagation.hpp>
#include <SeQuant/core/irreps/table_algebra.hpp>
#include <SeQuant/core/irreps/z2_algebra.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr_filter.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <boost/algorithm/string/trim.hpp>

#include <algorithm>
#include <cctype>
#include <memory>
#include <optional>
#include <regex>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace sequant {

namespace {

/// Converts the given IRREP names into a set, using @p parse to convert a
/// single name
template <typename Irrep, typename Parser>
IrrepSet<Irrep> to_irrep_set(const IrrepContextSpec::IrrepNames &names,
                             const Parser &parse) {
  IrrepSet<Irrep> set;
  for (const std::string &name : names) {
    if (name == "*") {
      set = IrrepSet<Irrep>::any();
    } else {
      set.insert(parse(name));
    }
  }

  return set;
}

/// Creates an IRREP context for the given algebra with the leaf tensor IRREPs
/// as specified in the 'tensors' and default entries of the given
/// specification
template <IrrepAlgebra A, typename Parser>
AnyIrrepContext make_context(A algebra, const IrrepContextSpec &spec,
                             const Parser &parse) {
  using Set = IrrepSet<typename A::irrep_type>;

  struct TensorRule {
    std::string key;
    std::shared_ptr<const ExpressionFilter::Rule> rule;
    Set irreps;
  };
  std::vector<TensorRule> blocks;
  std::vector<TensorRule> labels;
  container::map<std::wstring, Set> literal_labels;
  for (const auto &[key, names] : spec.tensors) {
    Set set = to_irrep_set<typename A::irrep_type>(names, parse);

    if (key.find('{') == std::string::npos) {
      literal_labels.emplace(toUtf16(key), set);

      std::vector<std::wregex> patterns;
      try {
        patterns.emplace_back(toUtf16(key));
      } catch (const std::regex_error &e) {
        throw Exception("Invalid label regex '" + key + "': " + e.what());
      }
      labels.push_back(
          TensorRule{.key = key,
                     .rule = std::make_shared<LabelRule>(std::move(patterns)),
                     .irreps = std::move(set)});
      continue;
    }

    const io::serialization::DeserializationOptions options{
        .def_perm_symm = Symmetry::Nonsymm,
        .def_braket_symm = BraKetSymmetry::Nonsymm,
        .def_col_symm = ColumnSymmetry::Nonsymm};
    ExprPtr tensor = io::serialization::from_string<ExprPtr>(key, options);
    if (!tensor->is<Tensor>()) {
      throw Exception("Expected '" + key + "' to be a tensor");
    }
    blocks.push_back(TensorRule{
        .key = key,
        .rule = std::make_shared<ContainsRule>(ExprMatcher(
            std::move(*tensor), {.tensor_cmp = TensorComparison::Block})),
        .irreps = std::move(set)});
  }

  // An empty optional means that unspecified tensors are an error
  std::optional<Set> fallback;
  if (!spec.error_on_unspecified) {
    fallback = to_irrep_set<typename A::irrep_type>(spec.default_irreps, parse);
  }

  auto resolver = [blocks = std::move(blocks), labels = std::move(labels),
                   literal_labels = std::move(literal_labels),
                   fallback = std::move(fallback)](const Tensor &tensor) {
    // @returns The rule matching the tensor (if any)
    // @throws Exception if several rules with different IRREPs match
    auto unique_match =
        [&tensor](const std::vector<TensorRule> &rules) -> const TensorRule * {
      const TensorRule *match = nullptr;
      for (const TensorRule &current : rules) {
        if (!current.rule->matches(tensor)) {
          continue;
        }

        if (match && match->irreps != current.irreps) {
          throw Exception("Ambiguous IRREPs for tensor " +
                          toUtf8(tensor.to_latex()) + ": matched by '" +
                          match->key + "' and '" + current.key + "'");
        }
        match = match ? match : &current;
      }

      return match;
    };

    if (const TensorRule *match = unique_match(blocks)) {
      return match->irreps;
    }

    // A key spelling out the label literally takes precedence over patterns
    if (auto it = literal_labels.find(std::wstring(tensor.label()));
        it != literal_labels.end()) {
      return it->second;
    }

    if (const TensorRule *match = unique_match(labels)) {
      return match->irreps;
    }

    if (!fallback.has_value()) {
      throw Exception("No IRREPs specified for tensor " +
                      toUtf8(tensor.to_latex()));
    }

    return fallback.value();
  };

  return IrrepContext<A>(std::move(algebra), std::move(resolver));
}

/// Creates an IRREP context for a Z2 group (D2h or one of its subgroups, or
/// a group consisting of symbols only)
/// @param symbols_must_be_declared Whether all symbolic IRREPs must be listed
///        in the specification's symbols
template <Z2Group Group>
AnyIrrepContext make_z2_context(const IrrepContextSpec &spec,
                                bool symbols_must_be_declared) {
  std::optional<std::set<std::string, std::less<>>> declared;
  if (symbols_must_be_declared || spec.symbols.has_value()) {
    const std::vector<std::string> symbols =
        spec.symbols.value_or(std::vector<std::string>{});
    declared.emplace(symbols.begin(), symbols.end());
  }

  auto parse = [declared](std::string_view name) {
    Z2Irrep<Group> irrep = make_irrep<Group>(name);
    if (declared.has_value()) {
      for (const std::string &symbol : irrep.symbols()) {
        if (!declared->contains(symbol)) {
          throw Exception("Unknown IRREP '" + symbol +
                          "': neither an IRREP of the point group nor a "
                          "declared symbol");
        }
      }
    }

    return irrep;
  };

  std::vector<std::pair<Z2Irrep<Group>, Z2Irrep<Group>>> relations;
  for (const std::string &relation :
       spec.relations.value_or(std::vector<std::string>{})) {
    const auto pos = relation.find('=');
    if (pos == std::string::npos ||
        relation.find('=', pos + 1) != std::string::npos) {
      throw Exception(
          "Expected IRREP relation of the form 'lhs = rhs' but got '" +
          relation + "'");
    }

    std::string lhs = relation.substr(0, pos);
    std::string rhs = relation.substr(pos + 1);
    boost::trim(lhs);
    boost::trim(rhs);
    if (lhs.empty() || rhs.empty()) {
      throw Exception("Incomplete IRREP relation '" + relation + "'");
    }

    relations.emplace_back(parse(lhs), parse(rhs));
  }

  return make_context(Z2Algebra<Group>(relations), spec, parse);
}

/// Creates an IRREP context for a group given by an explicit multiplication
/// table
AnyIrrepContext make_table_context(const IrrepContextSpec &spec) {
  if (spec.relations.has_value() || spec.symbols.has_value()) {
    throw Exception(
        "'relations' and 'symbols' are only supported for the shipped point "
        "groups and 'symbolic'");
  }

  const IrrepContextSpec::Table &group = spec.table.value();

  std::set<std::string, std::less<>> known{group.totally_symmetric};
  auto add_known = [&known](std::string_view name) {
    return std::string(*known.emplace(name).first);
  };

  TableAlgebra::table_type table;
  for (const auto &[key, names] : group.products) {
    const auto pos = key.find('*');
    if (pos == std::string::npos ||
        key.find('*', pos + 1) != std::string::npos) {
      throw Exception("Expected product of the form 'X*Y' but got '" + key +
                      "'");
    }

    std::string lhs = add_known(key.substr(0, pos));
    std::string rhs = add_known(key.substr(pos + 1));
    table.emplace(std::make_pair(std::move(lhs), std::move(rhs)),
                  to_irrep_set<std::string>(names, add_known));
  }

  TableAlgebra::conjugation_map conjugates;
  for (const auto &[irrep, conjugate] : group.conjugates) {
    conjugates.emplace(add_known(irrep), add_known(conjugate));
  }

  auto parse = [known = std::move(known)](std::string_view name) {
    if (!known.contains(name)) {
      throw Exception("Unknown IRREP '" + std::string(name) + "'");
    }
    return std::string(name);
  };

  return make_context(TableAlgebra(group.totally_symmetric, std::move(table),
                                   std::move(conjugates)),
                      spec, parse);
}

}  // namespace

std::string normalize_irrep_spec(std::string_view spec) {
  std::string result;
  std::size_t start = 0;
  while (start <= spec.size()) {
    std::size_t end = spec.find_first_of("*=", start);
    const bool last = end == std::string_view::npos;
    std::string name(spec.substr(start, last ? spec.npos : end - start));
    boost::trim(name);
    if (std::ranges::any_of(name,
                            [](unsigned char c) { return std::isspace(c); })) {
      throw Exception("IRREP names must not contain whitespace but got '" +
                      name + "'");
    }

    result += name;
    if (last) {
      break;
    }
    result += spec[end];
    start = end + 1;
  }

  return result;
}

AnyIrrepContext make_irrep_context(const IrrepContextSpec &input) {
  IrrepContextSpec spec = input;
  auto normalize = [](std::string &name) { name = normalize_irrep_spec(name); };
  auto normalize_all = [&](auto &names) {
    std::ranges::for_each(names, normalize);
  };
  normalize(spec.group);
  if (spec.table.has_value()) {
    normalize(spec.table->totally_symmetric);
    for (auto &[product, names] : spec.table->products) {
      normalize(product);
      normalize_all(names);
    }
    for (auto &[irrep, conjugate] : spec.table->conjugates) {
      normalize(irrep);
      normalize(conjugate);
    }
  }
  if (spec.symbols.has_value()) {
    normalize_all(*spec.symbols);
  }
  if (spec.relations.has_value()) {
    normalize_all(*spec.relations);
  }
  for (auto &[tensor, names] : spec.tensors) {
    normalize_all(names);
  }
  normalize_all(spec.default_irreps);

  if (spec.table.has_value()) {
    return make_table_context(spec);
  }

  const std::string &name = spec.group;
  if (name == "symbolic") {
    return make_z2_context<NoConcreteGroup>(spec, false);
  } else if (name == "C1") {
    return make_z2_context<point_groups::C1>(spec, true);
  } else if (name == "Ci") {
    return make_z2_context<point_groups::Ci>(spec, true);
  } else if (name == "C2") {
    return make_z2_context<point_groups::C2>(spec, true);
  } else if (name == "Cs") {
    return make_z2_context<point_groups::Cs>(spec, true);
  } else if (name == "C2v") {
    return make_z2_context<point_groups::C2v>(spec, true);
  } else if (name == "C2h") {
    return make_z2_context<point_groups::C2h>(spec, true);
  } else if (name == "D2") {
    return make_z2_context<point_groups::D2>(spec, true);
  } else if (name == "D2h") {
    return make_z2_context<point_groups::D2h>(spec, true);
  }

  throw Exception("Unknown point group '" + name + "'");
}

}  // namespace sequant
