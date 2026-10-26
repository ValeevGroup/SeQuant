#include "export_step.hpp"
#include "filter_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"
#include "utils.hpp"

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/export/context.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/expression_group.hpp>
#include <SeQuant/core/export/generation_optimizer.hpp>
#include <SeQuant/core/export/generator.hpp>
#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/export/name_recycling.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/irreps/any.hpp>
#include <SeQuant/core/irreps/propagation.hpp>
#include <SeQuant/core/irreps/table_algebra.hpp>
#include <SeQuant/core/irreps/z2_algebra.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>
#include <SeQuant/core/utility/topological.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <range/v3/view/concat.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

class MetaAwareIftContext : public ItfContext {
 public:
  MetaAwareIftContext(const ExportStep::ItfMeta &meta) : meta_(&meta) {
    set_index_id_offset(meta.idx_id_offset);
  }

  std::string get_tag(const IndexSpace &space) const override {
    SEQUANT_ASSERT(meta_);
    auto it = meta_->index_spaces.find(space);
    if (it == meta_->index_spaces.end()) {
      return ItfContext::get_name(space);
    }

    return it->second.tag;
  }

  std::string get_name(const IndexSpace &space) const override {
    SEQUANT_ASSERT(meta_);
    auto it = meta_->index_spaces.find(space);
    if (it == meta_->index_spaces.end()) {
      return ItfContext::get_tag(space);
    }

    return it->second.name;
  }

  std::string irrep_name(const AnyIrrep &irrep) const override {
    SEQUANT_ASSERT(meta_);
    std::string name = irrep.to_string();
    auto it = meta_->irrep_names.find(name);

    return it == meta_->irrep_names.end() ? name : it->second;
  }

  std::string any_irreps_name() const override {
    SEQUANT_ASSERT(meta_);
    auto it = meta_->irrep_names.find("*");

    return it == meta_->irrep_names.end() ? ItfContext::any_irreps_name()
                                          : it->second;
  }

 private:
  const ExportStep::ItfMeta *meta_ = nullptr;
};

namespace {

using json = nlohmann::json;

/// Parses a set of IRREPs given either as a single IRREP name or as an array
/// of IRREP names. The name "*" refers to any IRREP.
template <typename Irrep, typename Parser>
IrrepSet<Irrep> parse_irrep_set(const json &value, const Parser &parse,
                                std::string_view what) {
  auto add = [&](const json &entry, IrrepSet<Irrep> &set) {
    if (!entry.is_string()) {
      throw Exception("Expected IRREPs of " + std::string(what) +
                      " to be given as string(s)");
    }

    const std::string name = entry.get<std::string>();
    if (name == "*") {
      set = IrrepSet<Irrep>::any();
    } else {
      set.insert(parse(name));
    }
  };

  IrrepSet<Irrep> set;
  if (value.is_array() && !value.empty()) {
    for (const json &entry : value) {
      add(entry, set);
    }
  } else if (value.is_string()) {
    add(value, set);
  } else {
    throw Exception("Expected IRREPs of " + std::string(what) +
                    " to be given as a string or a non-empty array of strings");
  }

  return set;
}

/// Creates an IRREP context for the given algebra with the leaf tensor IRREPs
/// as specified in the 'tensors' and 'default' entries of the given
/// configuration
template <IrrepAlgebra A, typename Parser>
AnyIrrepContext make_irrep_context(A algebra, const json &config,
                                   const Parser &parse) {
  using Set = IrrepSet<typename A::irrep_type>;

  // Tensors are matched consistently with the other steps: keys containing
  // a '{' are serialized tensors matched by block (cmp. imports), all other
  // keys are full-match regular expressions for the tensor label (cmp. label
  // rules of the filter step)
  struct TensorRule {
    std::string key;
    std::shared_ptr<const ExpressionFilter::Rule> rule;
    Set irreps;
  };
  std::vector<TensorRule> blocks;
  std::vector<TensorRule> labels;
  container::map<std::wstring, Set> literal_labels;
  if (config.contains("tensors")) {
    const json &tensors = config.at("tensors");
    if (!tensors.is_object()) {
      throw Exception("Expected 'tensors' to be an object");
    }

    for (const auto &[key, value] : tensors.items()) {
      Set set = parse_irrep_set<typename A::irrep_type>(value, parse,
                                                        "tensor '" + key + "'");

      if (key.find('{') == std::string::npos) {
        literal_labels.emplace(toUtf16(key), set);
        labels.push_back(
            TensorRule{.key = key,
                       .rule = parse_contains_label(json{{"label", key}}),
                       .irreps = std::move(set)});
        continue;
      }

      ExprPtr tensor = io::serialization::from_string<ExprPtr>(key);
      if (!tensor->is<Tensor>()) {
        throw Exception("Expected '" + key + "' to be a tensor");
      }
      blocks.push_back(TensorRule{
          .key = key,
          .rule =
              parse_contains_expr(json{{"expr", key}}, TensorComparison::Block),
          .irreps = std::move(set)});
    }
  }

  // An empty optional means that unspecified tensors are an error
  std::optional<Set> fallback = Set::any();
  if (config.contains("default")) {
    const json &value = config.at("default");
    if (value.is_string() && value.get<std::string>() == "error") {
      fallback.reset();
    } else {
      fallback =
          parse_irrep_set<typename A::irrep_type>(value, parse, "the default");
    }
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

/// @returns The list of strings in the given (optional) array entry
std::vector<std::string> string_list(const json &config, const char *key) {
  std::vector<std::string> list;
  if (!config.contains(key)) {
    return list;
  }

  const json &value = config.at(key);
  if (!value.is_array()) {
    throw Exception("Expected '" + std::string(key) +
                    "' to be an array of strings");
  }
  for (const json &entry : value) {
    if (!entry.is_string()) {
      throw Exception("Expected '" + std::string(key) +
                      "' to be an array of strings");
    }
    list.push_back(entry.get<std::string>());
  }

  return list;
}

/// Creates an IRREP context for a Z2 group (D2h or one of its subgroups, or
/// a group consisting of symbols only)
/// @param symbols_must_be_declared Whether all symbolic IRREPs must be listed
///        in the configuration's 'symbols' entry
template <Z2Group Group>
AnyIrrepContext make_z2_context(const json &config,
                                bool symbols_must_be_declared) {
  std::optional<std::set<std::string, std::less<>>> declared;
  if (symbols_must_be_declared || config.contains("symbols")) {
    std::vector<std::string> symbols = string_list(config, "symbols");
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
  for (const std::string &relation : string_list(config, "relations")) {
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

  return make_irrep_context(Z2Algebra<Group>(relations), config, parse);
}

/// Creates an IRREP context for a group given by an explicit multiplication
/// table
AnyIrrepContext make_table_context(const json &group, const json &config) {
  if (config.contains("relations") || config.contains("symbols")) {
    throw Exception(
        "'relations' and 'symbols' are only supported for the shipped point "
        "groups and 'symbolic'");
  }

  for (const auto &[key, value] : group.items()) {
    if (key != "totally_symmetric" && key != "products" &&
        key != "conjugates") {
      throw Exception("Unknown key '" + key + "' in multiplication table");
    }
  }

  if (!group.contains("totally_symmetric") ||
      !group.at("totally_symmetric").is_string()) {
    throw Exception(
        "Multiplication table requires 'totally_symmetric' to be a string");
  }
  const std::string totally_symmetric =
      group.at("totally_symmetric").get<std::string>();

  std::set<std::string, std::less<>> known{totally_symmetric};
  auto add_known = [&known](std::string_view name) {
    return std::string(*known.emplace(name).first);
  };

  TableAlgebra::table_type table;
  if (group.contains("products")) {
    const json &products = group.at("products");
    if (!products.is_object()) {
      throw Exception("Expected 'products' to be an object");
    }

    for (const auto &[key, value] : products.items()) {
      const auto pos = key.find('*');
      if (pos == std::string::npos ||
          key.find('*', pos + 1) != std::string::npos) {
        throw Exception("Expected product of the form 'X*Y' but got '" + key +
                        "'");
      }

      std::string lhs = add_known(key.substr(0, pos));
      std::string rhs = add_known(key.substr(pos + 1));
      table.emplace(std::make_pair(std::move(lhs), std::move(rhs)),
                    parse_irrep_set<std::string>(value, add_known,
                                                 "product '" + key + "'"));
    }
  }

  TableAlgebra::conjugation_map conjugates;
  if (group.contains("conjugates")) {
    const json &value = group.at("conjugates");
    if (!value.is_object()) {
      throw Exception("Expected 'conjugates' to be an object");
    }

    for (const auto &[irrep, conjugate] : value.items()) {
      if (!conjugate.is_string()) {
        throw Exception("Expected conjugate of '" + irrep + "' to be a string");
      }
      conjugates.emplace(add_known(irrep),
                         add_known(conjugate.get<std::string>()));
    }
  }

  auto parse = [known = std::move(known)](std::string_view name) {
    if (!known.contains(name)) {
      throw Exception("Unknown IRREP '" + std::string(name) + "'");
    }
    return std::string(name);
  };

  return make_irrep_context(
      TableAlgebra(totally_symmetric, std::move(table), std::move(conjugates)),
      config, parse);
}

}  // namespace

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ExportStep);

std::string ExportStep::kind() const { return "export"; }

bool ExportStep::accepts_options() const { return true; }

bool ExportStep::requires_options() const { return true; }

void ExportStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "language") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      language_ = value.get<std::string>();
    } else if (key == "optimize") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      optimize_ = value.get<bool>();
    } else if (key == "output") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      filepath_ = value.get<std::string>();
    } else if (key == "grouping") {
      if (!value.is_object()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be an object");
      }

      for (const auto &[name, ids] : value.items()) {
        if (!ids.is_string()) {
          throw Exception("Expected IDs for group '" + name +
                          "' to be given as a string");
        }

        group_assoc_.emplace(ids.get<std::string>(), name);
      }
    } else if (key == "irreps") {
      irrep_context_ = parse_irreps(value);
    } else if (key == "meta") {
      // Processing happens later as we need to ensure we know the target
      // language beforehand
    } else if (key == "relative_order") {
      if (!value.is_array()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be an array");
      }

      for (const nlohmann::json &current : value) {
        if (!current.is_string()) {
          throw Exception("Entries in " + kind() + "." + key +
                          " are expected to be strings");
        }

        relative_order_.emplace_back(current.get<std::string>());
      }
    } else if (key == "imports") {
      if (!value.is_object()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be an object");
      }

      for (const auto &[tensor, name] : value.items()) {
        if (!name.is_string()) {
          throw Exception("Expected import name for '" + tensor +
                          "' to be a string");
        }

        import_names_.emplace(
            ExprMatcher(*io::serialization::from_string<ExprPtr>(tensor),
                        {.tensor_cmp = TensorComparison::Block}),
            name.get<std::string>());
      }
    } else if (key == "recycle_intermediates") {
      recycle_intermediates_ = parse_recycling(value);
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  if (language_.empty()) {
    throw Exception("The 'language' option for " + kind() + " is mandatory");
  }

  if (filepath_.empty()) {
    throw Exception("The 'output' option for " + kind() + " is mandatory");
  }

  if (options.contains("meta")) {
    meta_ = parse_meta(language_, options.at("meta"));
  }
}

struct ExportTreeDataCompare {
  std::size_t min_order(const ExecutionContext::Data<ProcessingData> entry) {
    std::size_t order = std::numeric_limits<std::size_t>::max();

    for (std::string_view current : ranges::views::concat(
             entry.associated_ids, entry.associated_group_ids)) {
      auto pos = current.find('.');
      if (pos != std::string_view::npos) {
        // Strip step name from ID to allow referencing e.g. 'res' instead of
        // 'step5.res'
        current = current.substr(pos + 1);
      }

      auto it = std::ranges::find(rel_order, current);
      order = std::min(order, static_cast<std::size_t>(std::ranges::distance(
                                  rel_order.begin(), it)));
    }

    return order;
  }

  bool operator()(const ExecutionContext::Data<ProcessingData> &lhs,
                  const ExecutionContext::Data<ProcessingData> &rhs) {
    const std::size_t lhs_order = min_order(lhs);
    const std::size_t rhs_order = min_order(rhs);

    if (lhs_order != rhs_order) {
      return lhs_order < rhs_order;
    }

    return (*this)(convert_data<ExportTreeData>(lhs.data.get()),
                   convert_data<ExportTreeData>(rhs.data.get()));
  }

  bool operator()(const ExportTreeData &lhs, const ExportTreeData &rhs) const {
    if (lhs.entries.size() != rhs.entries.size()) {
      return lhs.entries.size() < rhs.entries.size();
    }

    for (std::size_t i = 0; i < lhs.entries.size(); ++i) {
      if (!equal(lhs.entries.at(i), rhs.entries.at(i))) {
        return (*this)(lhs.entries.at(i), rhs.entries.at(i));
      }
    }

    return false;
  }

  bool equal(const ExportTreeData::Entry &lhs,
             const ExportTreeData::Entry &rhs) const {
    if (lhs.symm_contribution_target != rhs.symm_contribution_target) {
      return false;
    }

    if (lhs.tree->is_tensor() != rhs.tree->is_tensor()) {
      return false;
    }

    return lhs.tree->is_tensor()
               ? lhs.tree->as_tensor() == rhs.tree->as_tensor()
               : lhs.tree->as_variable() == rhs.tree->as_variable();
  }

  bool operator()(const ExportTreeData::Entry &lhs,
                  const ExportTreeData::Entry &rhs) const {
    if (lhs.symm_contribution_target != rhs.symm_contribution_target) {
      return lhs.symm_contribution_target < rhs.symm_contribution_target;
    }

    if (lhs.tree->is_tensor() != rhs.tree->is_tensor()) {
      return rhs.tree->is_tensor();
    }

    return lhs.tree->is_tensor()
               ? lhs.tree->as_tensor() < rhs.tree->as_tensor()
               : lhs.tree->as_variable() < rhs.tree->as_variable();
  }

  const std::vector<std::string> &rel_order;
};

std::size_t ExportStep::run(std::string_view, ExecutionContext &exctx,
                            const std::vector<std::string_view> &inputs) {
  auto perform_export = [&](auto &&generator, auto &&genctx)
    requires(
        std::derived_from<std::remove_cvref_t<decltype(genctx)>,
                          ExportContext> &&
        std::derived_from<std::remove_cvref_t<decltype(generator)>,
                          Generator<std::remove_cvref_t<decltype(genctx)>>>)
  {
    // Setup load/store/create/import strategies
    std::set<std::size_t> always_create_ids;
    std::vector<ExpressionGroup<>> groups =
        prepare_expressions(exctx, inputs, genctx, always_create_ids);

    if (recycle_intermediates_.has_value()) {
      // Entries flagged to overwrite their previous result (CSE
      // intermediates) are the ones whose names may be recycled
      auto overwrites = [&always_create_ids](const ExportNode<> &tree) {
        return always_create_ids.contains(tree->id());
      };

      recycle_intermediate_names(groups, genctx, overwrites,
                                 recycle_intermediates_.value());

      for (const ExpressionGroup<> &group : groups) {
        // Renamed results need the same zero strategy as the one set for
        // their original names in prepare_expressions
        for (const ExportNode<> &tree : group) {
          if (!overwrites(tree)) {
            continue;
          }

          if (tree->is_tensor()) {
            Tensor result = tree->as_tensor();
            genctx.rewrite(result);
            genctx.setZeroStrategy(result, ZeroStrategy::ZeroOnCreate);
          } else {
            Variable result = tree->as_variable();
            genctx.rewrite(result);
            genctx.setZeroStrategy(result, ZeroStrategy::ZeroOnCreate);
          }
        }
      }
    }

    for (const ExpressionGroup<> &current_group : groups) {
      std::set<Tensor> created_tensors;
      std::set<Variable> created_variables;

      for (const ExportNode<> &current_node : current_group) {
        std::string import_name;
        if (auto it = import_names_.find(current_node->expr());
            it != import_names_.end()) {
          import_name = it->second;
        }

        // Entries flagged to overwrite their previous result (e.g. CSE
        // intermediates) are always created fresh and intentionally not
        // recorded in the shared-result tracking, so that two unrelated
        // subexpressions reusing the same result name+block each get their own
        // allocation instead of accumulating onto one another.
        const bool force_create = always_create_ids.find(current_node->id()) !=
                                  always_create_ids.end();

        if (current_node->is_tensor()) {
          Tensor result = current_node->as_tensor();
          genctx.rewrite(result);

          if (!import_name.empty()) {
            genctx.set_import_name(result, import_name);
          }

          if (force_create ||
              created_tensors.find(result) == created_tensors.end()) {
            genctx.setLoadStrategy(result, LoadStrategy::Create,
                                   current_node->id());
            if (!force_create) created_tensors.emplace(std::move(result));
          } else {
            genctx.setLoadStrategy(result, LoadStrategy::Load,
                                   current_node->id());
          }
        } else {
          SEQUANT_ASSERT(current_node->is_variable());
          Variable result = current_node->as_variable();
          genctx.rewrite(result);

          if (!import_name.empty()) {
            genctx.set_import_name(result, import_name);
          }

          if (force_create ||
              created_variables.find(result) == created_variables.end()) {
            genctx.setLoadStrategy(result, LoadStrategy::Create,
                                   current_node->id());
            if (!force_create) created_variables.emplace(std::move(result));
          } else {
            genctx.setLoadStrategy(result, LoadStrategy::Load,
                                   current_node->id());
          }
        }
      }
    }

    // TODO: handle index batching

    export_groups(groups, generator, genctx);

    return generator.get_generated_code();
  };

  std::string generated;

  if (language_ == "itf") {
    ItfGenerator<MetaAwareIftContext> generator;
    MetaAwareIftContext genctx(std::get<ItfMeta>(meta_));
    if (irrep_context_.has_value()) {
      genctx.set_irrep_context(irrep_context_.value());
    }

    if (optimize_) {
      generated = perform_export(
          GenerationOptimizer<decltype(generator)>(std::move(generator)),
          genctx);
    } else {
      generated = perform_export(generator, genctx);
    }
  } else {
    throw Exception("Unsupported export language '" + language_ + "'");
  }

  SEQUANT_ASSERT(!generated.empty());

  std::ofstream stream(filepath_);
  stream << generated;

  return 0;
}

std::optional<opt::ShapeCompatibility> ExportStep::parse_recycling(
    const nlohmann::json &value) {
  if (value.is_boolean()) {
    if (!value.get<bool>()) {
      return std::nullopt;
    }
    return opt::ShapeCompatibility::IndexMultiset;
  }

  if (!value.is_object()) {
    throw Exception(
        "Value for export option 'recycle_intermediates' must be a boolean or "
        "an object");
  }

  opt::ShapeCompatibility compatibility =
      opt::ShapeCompatibility::IndexMultiset;
  for (const auto &[key, entry] : value.items()) {
    if (key != "shape") {
      throw Exception(
          "Unknown key in export option 'recycle_intermediates': '" + key +
          "'");
    }

    const std::string shape = entry.is_string() ? entry.get<std::string>() : "";
    if (shape == "multiset") {
      compatibility = opt::ShapeCompatibility::IndexMultiset;
    } else if (shape == "block") {
      compatibility = opt::ShapeCompatibility::ExactBlock;
    } else {
      throw Exception(
          "Value for 'recycle_intermediates.shape' must be either 'multiset' "
          "or 'block'");
    }
  }

  return compatibility;
}

ExportStep::meta_type ExportStep::parse_meta(std::string_view language,
                                             const nlohmann::json &meta) {
  if (language == "itf") {
    ItfMeta data;

    for (const auto &[key, value] : meta.items()) {
      if (key == "index_spaces") {
        if (!value.is_object()) {
          throw Exception("Expected '" + key + "' to be an object");
        }

        for (const auto &[space_label, info] : value.items()) {
          IndexSpace space(space_label);

          data.index_spaces.emplace(
              std::move(space),
              ItfMeta::IndexSpaceMeta{.name = info.value("name", std::string{}),
                                      .tag = info.value("tag", std::string{})});
        }
      } else if (key == "min_index_id") {
        if (!value.is_number_integer()) {
          throw Exception("Expected '" + key + "' to be an integer");
        }
        if (!value.is_number_unsigned()) {
          throw Exception("Expected '" + key + "' to be non-negative");
        }

        data.idx_id_offset = value.get<std::size_t>();
      } else if (key == "irrep_names") {
        if (!value.is_object()) {
          throw Exception("Expected '" + key + "' to be an object");
        }

        for (const auto &[irrep, name] : value.items()) {
          if (!name.is_string()) {
            throw Exception("Expected the name of IRREP '" + irrep +
                            "' to be a string");
          }
          data.irrep_names.emplace(irrep, name.get<std::string>());
        }
      } else {
        throw Exception("Unknown meta data type '" + key + "' for language '" +
                        std::string(language) + "'");
      }
    }

    return data;
  } else {
    if (!meta.empty()) {
      throw Exception("Exporting to '" + std::string(language) +
                      "' does not accept meta data");
    }

    return {};
  }
}

AnyIrrepContext ExportStep::parse_irreps(const nlohmann::json &irreps) {
  try {
    if (!irreps.is_object()) {
      throw Exception("Expected an object");
    }

    for (const auto &[key, value] : irreps.items()) {
      if (key != "group" && key != "relations" && key != "symbols" &&
          key != "tensors" && key != "default") {
        throw Exception("Unknown key '" + key + "'");
      }
    }

    if (!irreps.contains("group")) {
      throw Exception("The 'group' entry is mandatory");
    }

    const json &group = irreps.at("group");
    if (group.is_object()) {
      return make_table_context(group, irreps);
    }
    if (!group.is_string()) {
      throw Exception(
          "Expected 'group' to be the name of a point group or a "
          "multiplication table");
    }

    const std::string name = group.get<std::string>();
    if (name == "symbolic") {
      return make_z2_context<NoConcreteGroup>(irreps, false);
    } else if (name == "C1") {
      return make_z2_context<point_groups::C1>(irreps, true);
    } else if (name == "Ci") {
      return make_z2_context<point_groups::Ci>(irreps, true);
    } else if (name == "C2") {
      return make_z2_context<point_groups::C2>(irreps, true);
    } else if (name == "Cs") {
      return make_z2_context<point_groups::Cs>(irreps, true);
    } else if (name == "C2v") {
      return make_z2_context<point_groups::C2v>(irreps, true);
    } else if (name == "C2h") {
      return make_z2_context<point_groups::C2h>(irreps, true);
    } else if (name == "D2") {
      return make_z2_context<point_groups::D2>(irreps, true);
    } else if (name == "D2h") {
      return make_z2_context<point_groups::D2h>(irreps, true);
    }

    throw Exception("Unknown point group '" + name + "'");
  } catch (const Exception &e) {
    throw Exception("Invalid 'irreps' option for export: " +
                    std::string(e.what()));
  }
}

std::vector<ExpressionGroup<>> ExportStep::prepare_expressions(
    ExecutionContext &exctx, const std::vector<std::string_view> &inputs,
    ExportContext &genctx, std::set<std::size_t> &always_create_ids) const {
  using CtxData = ExecutionContext::Data<ProcessingData>;

  std::vector<CtxData> data;
  for (std::string_view current_input : inputs) {
    for (CtxData current : exctx.get_data(current_input)) {
      data.push_back(std::move(current));
    }
  }

  auto tree_data = data | std::views::transform([](const CtxData &dat) {
                     return convert_data<ExportTreeData>(dat.data.get());
                   });

  // Pre-compute CtxData inter-dependencies
  std::map<const CtxData *, container::svector<const CtxData *>> dependencies;
  for (const CtxData &base : data) {
    const ExportTreeData &current =
        convert_data<ExportTreeData>(base.data.get());

    container::svector<const CtxData *> deps;

    auto all_dep_ids =
        current.entries |
        std::views::transform(&ExportTreeData::Entry::dependencies) |
        std::views::join;

    for (std::size_t id : all_dep_ids) {
      auto it = std::ranges::find_if(
          tree_data, [id](const ExportTreeData &tree_data) {
            for (const auto &entry : tree_data.entries) {
              if (entry.tree->id() == id) {
                return true;
              }
            }
            return false;
          });

      SEQUANT_ASSERT(it != end(tree_data));

      const CtxData *producer = &(*it.base());

      // A dependency resolving to an entry in the same group (e.g. after a
      // merge_inputs cse pulled dependent results into a single group) is
      // satisfied by the intra-group ordering; turning it into an inter-group
      // edge would be a spurious self-edge and break the topological sort.
      if (producer == &base) {
        continue;
      }

      deps.emplace_back(producer);
    }

    dependencies.emplace(&base, std::move(deps));
  }

  auto to_deps = [&](const CtxData &ctxdata) {
    auto it = dependencies.find(&ctxdata);
    SEQUANT_ASSERT(it != end(dependencies));
    return it->second | std::views::transform(
                            [](auto ptr) -> decltype(auto) { return *ptr; });
  };

  // Sort inputs in partitions of contributions to the same result
  // while respecting hard dependencies as well as requested relative order
  std::vector<std::size_t> access_order =
      topological_order(data, to_deps, ExportTreeDataCompare{relative_order_});

  std::vector<ExpressionGroup<>> groups;
  std::vector<std::vector<std::pair<Tensor, Tensor>>> required_symmetrizations;
  for (std::size_t idx : access_order) {
    const ExportTreeData &current_data = tree_data[idx];

    std::string name = [&]() -> std::string {
      for (std::string_view current :
           ranges::views::concat(data.at(idx).associated_ids,
                                 data.at(idx).associated_group_ids)) {
        auto it = group_assoc_.find(detail::strip_step_prefix(current));
        if (it != group_assoc_.end()) {
          return it->second;
        }
      }

      return "Default";
    }();

    auto it =
        std::ranges::find_if(groups, [&name](const ExpressionGroup<> &group) {
          return group.is_named() && group.name() == name;
        });
    if (it == groups.end()) {
      groups.emplace_back(std::move(name));
      it = groups.end() - 1;

      required_symmetrizations.emplace_back();
    }

    ExpressionGroup<> &group = *it;
    std::vector<std::pair<Tensor, Tensor>> &symms =
        required_symmetrizations.at(std::ranges::distance(groups.begin(), it));

    for (const ExportTreeData::Entry &current : current_data.entries) {
      group.add(current.tree);

      if (current.overwrite_previous.has_value()) {
        // Never ZeroOnReuse: the summands of a sum all accumulate into the same
        // result, so zeroing on every reappearance would discard what the
        // earlier summands have already added.
        constexpr ZeroStrategy strategy = ZeroStrategy::ZeroOnCreate;
        // Key the strategy by the rewritten result, as every other strategy
        // set-site does (and as the code generator's lookup expects):
        // genctx.rewrite() may reorder a tensor's slots, so setting it on the
        // unrewritten result would store it under a key the lookup never hits.
        // rewrite() is the generic per-generator hook, so this holds for every
        // backend.
        if (current.tree->is_tensor()) {
          Tensor result = current.tree->as_tensor();
          genctx.rewrite(result);
          genctx.setZeroStrategy(result, strategy);
        } else {
          SEQUANT_ASSERT(current.tree->is_variable());
          Variable result = current.tree->as_variable();
          genctx.rewrite(result);
          genctx.setZeroStrategy(result, strategy);
        }

        // An entry flagged to overwrite its previous result (e.g. a CSE
        // intermediate) must be allocated freshly rather than accumulated onto.
        // Different equations can reuse the same intermediate name+block for
        // unrelated subexpressions; without this, the shared-result tracking
        // below would treat the second definition as a continuation of the
        // first and emit a load-and-accumulate, corrupting the value.
        if (current.overwrite_previous.value()) {
          always_create_ids.insert(current.tree->id());
        }
      }

      if (!current.batch_indices.empty()) {
        // Note: this throws on generators that don't support batching
        // we leave choosing to batch in the user's responsibility
        genctx.set_batch_indices(current.batch_indices, current.tree->id());
      }

      if (current.symm_contribution_target.has_value()) {
        auto symm_it = std::ranges::find(
            symms, current.tree->as_tensor(),
            &std::remove_cvref_t<decltype(symms)>::value_type::first);

        if (symm_it == symms.end()) {
          symms.emplace_back(current.tree->as_tensor(),
                             current.symm_contribution_target.value());
        }
      }
    }
  }

  SEQUANT_ASSERT(groups.size() == required_symmetrizations.size());
  for (std::size_t i = 0; i < groups.size(); ++i) {
    ExpressionGroup<> &group = groups.at(i);

    for (const auto &[unsymm_res, symm_res] : required_symmetrizations.at(i)) {
      auto it = std::find_if(group.rbegin(), group.rend(),
                             [&unsymm_res](const ExportNode<> &node) {
                               return node->is_tensor() &&
                                      node->as_tensor() == unsymm_res;
                             });

      SEQUANT_ASSERT(it != group.rend());

      ResultExpr symmetrization(
          symm_res, generateResultSymmetrization(symm_res, unsymm_res.label()));

      const auto batch_indices = genctx.batch_indices((*it)->id());

      ExportNode<> symm_tree = to_export_tree(std::move(symmetrization));
      std::size_t symm_id = symm_tree->id();

      group.insert(it.base(), std::move(symm_tree));

      if (!batch_indices.empty()) {
        genctx.set_batch_indices(batch_indices, symm_id);
      }
    }
  }

  return groups;
}

}  // namespace sequant::util::extint
