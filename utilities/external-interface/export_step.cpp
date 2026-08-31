#include "export_step.hpp"
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
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/topological.hpp>

#include <range/v3/view/concat.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
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

 private:
  const ExportStep::ItfMeta *meta_ = nullptr;
};

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ExportStep, "export");

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
    std::vector<ExpressionGroup<>> groups =
        prepare_expressions(exctx, inputs, genctx);

    for (const ExpressionGroup<> &current_group : groups) {
      std::set<Tensor> created_tensors;
      std::set<Variable> created_variables;

      for (const ExportNode<> &current_node : current_group) {
        std::string import_name;
        if (auto it = import_names_.find(current_node->expr());
            it != import_names_.end()) {
          import_name = it->second;
        }

        if (current_node->is_tensor()) {
          Tensor result = current_node->as_tensor();
          genctx.rewrite(result);

          if (!import_name.empty()) {
            genctx.set_import_name(result, import_name);
          }

          if (created_tensors.find(result) == created_tensors.end()) {
            genctx.setLoadStrategy(result, LoadStrategy::Create,
                                   current_node->id());

            created_tensors.emplace(std::move(result));
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

          if (created_variables.find(result) == created_variables.end()) {
            genctx.setLoadStrategy(result, LoadStrategy::Create,
                                   current_node->id());

            created_variables.emplace(std::move(result));
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

std::vector<ExpressionGroup<>> ExportStep::prepare_expressions(
    ExecutionContext &exctx, const std::vector<std::string_view> &inputs,
    ExportContext &genctx) const {
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

      deps.emplace_back(&(*it.base()));
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
        auto pos = current.find('.');
        if (pos != std::string_view::npos) {
          // Strip step name from ID to allow referencing e.g. 'res' instead of
          // 'step5.res'
          current = current.substr(pos + 1);
        }

        auto it = group_assoc_.find(current);
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
        ZeroStrategy strategy =
            current.overwrite_previous.value()
                ? ZeroStrategy::ZeroOnReuse | ZeroStrategy::ZeroOnCreate
                : ZeroStrategy::ZeroOnCreate;
        if (current.tree->is_tensor()) {
          genctx.setZeroStrategy(current.tree->as_tensor(), strategy);
        } else {
          SEQUANT_ASSERT(current.tree->is_variable());
          genctx.setZeroStrategy(current.tree->as_variable(), strategy);
        }
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

      group.insert(it.base(), to_export_tree(symmetrization));
    }
  }

  return groups;
}

}  // namespace sequant::util::extint
