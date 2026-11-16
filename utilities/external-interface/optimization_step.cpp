#include "optimization_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"
#include "utils.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <utility>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(OptimizationStep);

std::string OptimizationStep::kind() const { return "optimize"; }

bool OptimizationStep::accepts_options() const { return true; }

bool OptimizationStep::requires_options() const { return false; }

void OptimizationStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  bool batch_given = false;

  for (const auto &[key, value] : options.items()) {
    if (key == "objective") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      if (value == "DenseFLOPs") {
        options_.objective_function = ObjectiveFunction::DenseFLOPs;
      } else if (value == "DenseTimeSpace") {
        options_.objective_function = ObjectiveFunction::DenseTimeSpace;
      } else if (value == "DenseTimeSpaceBatched") {
        options_.objective_function = ObjectiveFunction::DenseTimeSpaceBatched;
      } else {
        throw Exception("Invalid value for " + kind() + " option '" + key +
                        "': '" + value.get<std::string>() + "'");
      }
    } else if (key == "reorder_sums") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      options_.reorder =
          value.get<bool>() ? ReorderSum::Reorder : ReorderSum::NoReorder;
    } else if (key == "cse") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      if (value == "none") {
        options_.CSE.subnet = false;
      } else if (value == "subnet") {
        options_.CSE.subnet = true;
      } else {
        throw Exception("Invalid value for " + kind() + " option '" + key +
                        "': '" + value.get<std::string>() + "'");
      }
    } else if (key == "intermediate_size_penalty") {
      if (!value.is_number()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a number");
      }

      options_.footprint_weight = value.get<double>();
    } else if (key == "replay") {
      set_replay_options(value);
    } else if (key == "batch") {
      set_batch_options(value);
      batch_given = true;
    } else if (key == "prune_outer_products") {
      if (!value.is_boolean()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a boolean");
      }

      options_.prune_outer_products = value.get<bool>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }

  const ObjectiveFunction objective = options_.objective_function;
  const bool batched = objective == ObjectiveFunction::DenseTimeSpaceBatched;

  if (options_.footprint_weight != 0.0 &&
      objective != ObjectiveFunction::DenseFLOPs) {
    throw Exception(kind() +
                    " option 'intermediate_size_penalty' requires the "
                    "'DenseFLOPs' objective");
  }
  if (options_.CSE.subnet && objective != ObjectiveFunction::DenseFLOPs) {
    throw Exception(kind() +
                    " option 'cse' = 'subnet' requires the 'DenseFLOPs' "
                    "objective");
  }
  if (options_.batch_policy.persistent_only && !replay_.has_value()) {
    throw Exception(kind() +
                    " option 'batch.persistent_only' requires the 'replay' "
                    "option");
  }
  if (batch_given && !batched) {
    throw Exception(kind() +
                    " option 'batch' requires the 'DenseTimeSpaceBatched' "
                    "objective");
  }
  if (batched && !batch_given) {
    throw Exception(kind() +
                    " objective 'DenseTimeSpaceBatched' requires the 'batch' "
                    "option");
  }
}

OptimizationStep::OptimizationStep() {
  options_.objective_function = ObjectiveFunction::DenseTimeSpace;
}

void OptimizationStep::set_batch_options(const nlohmann::json &options) {
  const std::string name = kind() + ".batch";

  if (!options.is_object()) {
    throw Exception(name + " expects a JSON object");
  }

  auto get_labels = [&](const std::string &key, const nlohmann::json &value) {
    std::set<std::string> labels;
    const std::string where = "Value for " + name + " option '" + key + "'";
    if (!value.is_array()) {
      throw Exception(where + " must be an array of index space labels");
    }
    for (const nlohmann::json &entry : value) {
      if (!entry.is_string() || entry.get<std::string>().empty()) {
        throw Exception(where + " must only contain non-empty strings");
      }
      labels.insert(entry.get<std::string>());
    }
    return labels;
  };
  auto get_size = [&](const std::string &key, const nlohmann::json &value) {
    if (!value.is_number_integer() || value.get<long long>() <= 0) {
      throw Exception("Value for " + name + " option '" + key +
                      "' must be a positive integer");
    }
    return value.get<std::size_t>();
  };
  auto get_bool = [&](const std::string &key, const nlohmann::json &value) {
    if (!value.is_boolean()) {
      throw Exception("Value for " + name + " option '" + key +
                      "' must be a boolean");
    }
    return value.get<bool>();
  };

  BatchPolicy &policy = options_.batch_policy;
  std::set<std::string> contracted;
  std::set<std::string> external;
  std::map<std::string, std::size_t> target_sizes;
  std::optional<std::size_t> uniform_target_size;
  bool threshold_given = false;

  for (const auto &[key, value] : options.items()) {
    if (key == "peak_threshold") {
      if (!value.is_number() || !(value.get<double>() > 0.0)) {
        throw Exception("Value for " + name + " option '" + key +
                        "' must be a positive number (bytes)");
      }
      policy.peak_threshold = value.get<double>();
      threshold_given = true;
    } else if (key == "contracted") {
      contracted = get_labels(key, value);
    } else if (key == "external") {
      external = get_labels(key, value);
    } else if (key == "target_size") {
      if (value.is_object()) {
        for (const auto &[label, size] : value.items()) {
          target_sizes[label] = get_size(key + "." + label, size);
        }
      } else {
        uniform_target_size = get_size(key, value);
      }
    } else if (key == "spectator_indices") {
      policy.batch_spectator_indices = get_bool(key, value);
    } else if (key == "persistent_only") {
      policy.persistent_only = get_bool(key, value);
    } else if (key == "accumulation_factor") {
      if (!value.is_number() || value.get<double>() < 0.0) {
        throw Exception("Value for " + name + " option '" + key +
                        "' must be a non-negative number");
      }
      policy.accumulation_factor = value.get<double>();
    } else {
      throw Exception("Unknown option key for " + name + ": '" + key + "'");
    }
  }

  if (!threshold_given) {
    throw Exception("The 'peak_threshold' option for " + name +
                    " is mandatory");
  }
  if (contracted.empty() && external.empty()) {
    throw Exception(name +
                    " requires at least one of 'contracted' and 'external'");
  }
  if (!uniform_target_size && target_sizes.empty()) {
    throw Exception("The 'target_size' option for " + name + " is mandatory");
  }
  for (const std::set<std::string> *labels : {&contracted, &external}) {
    for (const std::string &label : *labels) {
      if (!uniform_target_size && !target_sizes.contains(label)) {
        throw Exception("No 'target_size' given for batchable space '" + label +
                        "' in " + name);
      }
    }
  }

  if (!external.empty() && !policy.batch_spectator_indices) {
    throw Exception(name +
                    " option 'external' has no effect without "
                    "'spectator_indices'");
  }
  if (policy.batch_spectator_indices && external.empty()) {
    throw Exception(name +
                    " option 'spectator_indices' requires batchable "
                    "'external' spaces");
  }

  batch_labels_.clear();
  batch_labels_.insert(contracted.begin(), contracted.end());
  batch_labels_.insert(external.begin(), external.end());
  for (const auto &entry : target_sizes) {
    batch_labels_.insert(entry.first);
  }

  using WLabelSet = std::set<std::wstring, std::less<>>;
  auto to_wide = [](const std::set<std::string> &labels) {
    WLabelSet res;
    for (const std::string &label : labels) {
      res.insert(toUtf16(label));
    }
    return res;
  };
  policy.is_batchable_contracted_index =
      [labels = to_wide(contracted)](const Index &ix) {
        return labels.contains(ix.space().base_key());
      };
  policy.is_batchable_external_index =
      [labels = to_wide(external)](const Index &ix) {
        return labels.contains(ix.space().base_key());
      };
  std::map<std::wstring, std::size_t, std::less<>> wide_sizes;
  for (const auto &[label, size] : target_sizes) {
    wide_sizes.emplace(toUtf16(label), size);
  }
  policy.batch_target_size =
      [wide_sizes = std::move(wide_sizes),
       uniform_target_size](const Index &ix) -> std::size_t {
    if (uniform_target_size) {
      return *uniform_target_size;
    }
    auto it = wide_sizes.find(ix.space().base_key());
    return it != wide_sizes.end() ? it->second : ix.space().dimension();
  };
}

void OptimizationStep::check_batch_spaces() const {
  if (batch_labels_.empty()) {
    return;
  }

  std::set<std::wstring> known;
  for (const IndexSpace &space :
       get_default_context().index_basis_registry()->spaces()) {
    known.insert(std::wstring(space.base_key()));
  }
  for (const std::string &label : batch_labels_) {
    if (!known.contains(toUtf16(label))) {
      throw Exception("Unknown index space '" + label + "' in " + kind() +
                      ".batch");
    }
  }
}

void OptimizationStep::set_replay_options(const nlohmann::json &options) {
  const std::string name = kind() + ".replay";

  if (!options.is_object()) {
    throw Exception(name + " expects a JSON object");
  }

  auto is_volatile = std::make_shared<ExpressionFilter>();
  is_volatile->set_require_all(false);
  bool volatile_given = false;

  ReplayOptions replay;

  auto get_string = [&](const std::string &key, const nlohmann::json &value) {
    if (!value.is_string() || value.get<std::string>().empty()) {
      throw Exception("Value for " + name + " option '" + key +
                      "' must be a non-empty string");
    }

    return value.get<std::string>();
  };

  for (const auto &[key, value] : options.items()) {
    if (key == "volatile") {
      auto add_rule = [&](const nlohmann::json &entry) {
        if (entry.is_string()) {
          is_volatile->add_rule(
              parse_contains_label(nlohmann::json{{"label", entry}}));
        } else if (entry.is_object() && entry.contains("expr")) {
          is_volatile->add_rule(
              parse_contains_expr(entry, TensorComparison::Block));
        } else if (entry.is_object() && entry.contains("label")) {
          is_volatile->add_rule(parse_contains_label(entry));
        } else {
          throw Exception("Entries of " + name +
                          ".volatile must be strings or objects containing "
                          "either 'expr' or 'label'");
        }
      };

      if (value.is_array()) {
        if (value.empty()) {
          throw Exception(name + ".volatile must not be empty");
        }

        for (const nlohmann::json &entry : value) {
          add_rule(entry);
        }
      } else {
        add_rule(value);
      }

      volatile_given = true;
    } else if (key == "weight") {
      if (!value.is_number()) {
        throw Exception("Value for " + name + " option '" + key +
                        "' must be a number");
      }

      options_.volatile_weight = value.get<double>();
    } else if (key == "intermediate_label") {
      replay.intermediate_label = get_string(key, value);
    } else if (key == "persistent_id") {
      replay.persistent_id = get_string(key, value);
    } else if (key == "volatile_id") {
      replay.volatile_id = get_string(key, value);
    } else {
      throw Exception("Unknown option key for " + name + ": '" + key + "'");
    }
  }

  if (!volatile_given) {
    throw Exception("The 'volatile' option for " + name + " is mandatory");
  }

  for (const std::string &id : {replay.persistent_id, replay.volatile_id}) {
    if (!ExecutionContext::is_valid_id(id, false)) {
      throw Exception("Invalid output ID '" + id + "' in " + name);
    }
  }
  if (replay.persistent_id == replay.volatile_id) {
    throw Exception("persistent_id and volatile_id in " + name +
                    " must differ");
  }

  options_.batch_policy.is_volatile_leaf = [is_volatile](const Tensor &tensor) {
    return is_volatile->matches(tensor);
  };
  replay.is_volatile = std::move(is_volatile);

  replay_ = std::move(replay);
}

/// Optimizes @p expr, returning the result symmetrizer that has been split
/// off beforehand (see push_symmetrizer)
std::optional<ExprPtr> optimize_unsymmetrized(ResultExpr &expr,
                                              const OptimizeOptions &options) {
  std::optional<ExprPtr> symmetrizer = pop_symmetrizer(expr);

  optimize(expr, options);

  return symmetrizer;
}

std::size_t OptimizationStep::run(std::string_view step_id,
                                  ExecutionContext &ctx,
                                  const std::vector<std::string_view> &inputs) {
  check_batch_spaces();

  if (!replay_.has_value()) {
    return OneByOneProcessingStep<ExpressionData, true,
                                  InputAccess::WorkingCopy>::run(step_id, ctx,
                                                                 inputs);
  }

  using Node = ExportNode<>;

  auto to_tree = [](const auto &expr) {
    return to_export_tree(expr, /*retain_braket=*/true);
  };
  auto to_result = [](const Node &tree) {
    if (tree->is_tensor()) {
      return ResultExpr(tree->as_tensor(), to_expr(tree));
    }
    return ResultExpr(tree->as_variable(), to_expr(tree));
  };

  // All inputs are optimized up front, so that persistent intermediates can
  // be shared across all of them
  std::vector<Node> trees;
  std::vector<std::size_t> tree_to_input;
  std::vector<std::optional<ExprPtr>> symmetrizers;
  // The inputs, which may have been consumed by now (their IDs are still valid)
  std::vector<ExecutionContext::Data<ProcessingData>> input_data;
  for (std::string_view input : inputs) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         ctx.get_data(input)) {
      SEQUANT_ASSERT(!current.associated_ids.empty());

      try {
        WorkingCopy<ExpressionData> working =
            ctx.working_copy<ExpressionData>(current);

        for (ResultExpr &expr : working.data.expressions) {
          symmetrizers.push_back(optimize_unsymmetrized(expr, options_));
          trees.push_back(to_tree(expr));
          tree_to_input.push_back(input_data.size());
        }
      } catch (const std::exception &e) {
        throw Exception("Error in " + kind() + " on input " +
                        std::string(current.associated_ids.front()) + ": " +
                        e.what());
      }

      input_data.push_back(current);
    }
  }
  const std::size_t n_inputs = input_data.size();

  const std::vector<std::size_t> definitions =
      opt::extract_persistent_subexpressions(
          trees, std::function<bool(const Node &)>([this](const Node &leaf) {
            return replay_->is_volatile->matches(*leaf->expr());
          }),
          to_tree,
          std::function<std::string(const Node &, std::size_t)>(
              [this](const Node &, std::size_t counter) {
                return replay_->intermediate_label + std::to_string(counter);
              }));

  // Definitions are inserted right before their first use and are assigned
  // to the persistent part of the input that use belongs to
  std::vector<std::pair<ExpressionData, ExpressionData>> parts(n_inputs);
  std::vector<ResultExpr> unassigned;
  std::size_t next_definition = 0;
  std::size_t expr_idx = 0;
  for (std::size_t i = 0; i < trees.size(); ++i) {
    if (next_definition < definitions.size() &&
        definitions[next_definition] == i) {
      unassigned.push_back(to_result(trees[i]));
      ++next_definition;
      continue;
    }

    auto &[persistent_part, volatile_part] =
        parts.at(tree_to_input.at(expr_idx));

    std::ranges::move(unassigned,
                      std::back_inserter(persistent_part.expressions));
    unassigned.clear();

    ResultExpr result = to_result(trees[i]);
    push_symmetrizer(result, std::move(symmetrizers.at(expr_idx)));
    volatile_part.expressions.push_back(std::move(result));

    ++expr_idx;
  }
  SEQUANT_ASSERT(unassigned.empty());
  SEQUANT_ASSERT(expr_idx == tree_to_input.size());

  detail::AliasPropagator propagator(ctx, step_id);
  for (std::size_t i = 0; i < n_inputs; ++i) {
    ctx.set_data(step_id, 2 * i, std::move(parts[i].first));
    ctx.set_data(step_id, 2 * i + 1, std::move(parts[i].second));

    propagator.record(input_data[i].associated_ids,
                      input_data[i].associated_group_ids, 2 * i, 2,
                      /* alias_source = */ false);
  }
  propagator.flush();

  const std::size_t n_outputs = 2 * n_inputs;

  // Expose the two parts of every alias created above as sub-IDs
  std::map<std::string, std::vector<std::size_t>> alias_inputs;
  for (std::size_t input_idx = 0; input_idx < n_inputs; ++input_idx) {
    const ExecutionContext::Data<ProcessingData> &current =
        input_data[input_idx];

    std::vector<std::string> aliases = {std::string(step_id)};
    for (std::string_view id : current.associated_ids) {
      if (std::optional<std::string_view> suffix =
              detail::strip_autogenerated_id(id)) {
        aliases.push_back(std::string(step_id) + std::string(*suffix));
      }
    }
    for (std::string_view id : current.associated_group_ids) {
      aliases.push_back(detail::make_group_alias(step_id, id));
    }

    for (std::string &alias : aliases) {
      if (alias == step_id || ctx.has_data(alias)) {
        alias_inputs[std::move(alias)].push_back(input_idx);
      }
    }
  }

  for (const auto &[alias, input_indices] : alias_inputs) {
    for (const auto &[sub_id, offset] :
         {std::pair{replay_->persistent_id, 0}, {replay_->volatile_id, 1}}) {
      std::string sub_alias = alias + "." + sub_id;
      if (ctx.has_data(sub_alias)) {
        throw Exception("Output ID '" + sub_alias + "' of " + kind() +
                        " is ambiguous");
      }

      std::vector<std::string> ids;
      for (std::size_t idx : input_indices) {
        ids.push_back(std::string(step_id) + "." +
                      std::to_string(2 * idx + offset));
      }

      ctx.add_data_alias(ids, std::move(sub_alias));
    }
  }

  return n_outputs;
}

std::size_t OptimizationStep::process(std::string_view id_prefix,
                                      std::size_t id_start,
                                      ExecutionContext &ctx,
                                      WorkingCopy<ExpressionData> &input) {
  SEQUANT_ASSERT(!replay_.has_value());

  for (ResultExpr &expr : input.data.expressions) {
    std::optional<ExprPtr> symmetrizer = optimize_unsymmetrized(expr, options_);
    push_symmetrizer(expr, std::move(symmetrizer));
  }

  ctx.set_data(id_prefix, id_start, std::move(input.data));

  return 1;
}

}  // namespace sequant::util::extint
