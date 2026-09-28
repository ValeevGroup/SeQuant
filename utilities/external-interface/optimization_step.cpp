#include "optimization_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"
#include "utils.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iterator>
#include <map>
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

  for (const auto &[key, value] : options.items()) {
    if (key == "objective") {
      if (!value.is_string()) {
        throw Exception("Value for " + kind() + " option '" + key +
                        "' must be a string");
      }

      if (value == "DenseFLOPs") {
        options_.objective_function = ObjectiveFunction::DenseFLOPs;
      } else if (value == "DenseSize") {
        options_.objective_function = ObjectiveFunction::DenseSize;
      } else if (value == "DensePeakSize") {
        options_.objective_function = ObjectiveFunction::DensePeakSize;
      } else if (value == "DensePeakSizeBatched") {
        options_.objective_function = ObjectiveFunction::DensePeakSizeBatched;
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
/// off beforehand (see reattach_symmetrizer)
std::optional<ExprPtr> optimize_unsymmetrized(ResultExpr &expr,
                                              const OptimizeOptions &options) {
  std::optional<ExprPtr> symmetrizer = pop_symmetrizer(expr);

  optimize(expr, options);

  return symmetrizer;
}

void reattach_symmetrizer(ResultExpr &expr,
                          std::optional<ExprPtr> symmetrizer) {
  if (symmetrizer.has_value()) {
    expr.expression() = ex<Product>(ExprPtrList{std::move(symmetrizer.value()),
                                                std::move(expr.expression())},
                                    Product::Flatten::No);
  }
}

/// Adds the number of tensor leaves of @p node to @p num_tensors and returns
/// whether none of its leaves is volatile.
// Note: visit_leaf can't be used as it doesn't stop at non-root nodes
bool count_tensors_if_persistent(const ExportNode<> &node,
                                 const ExpressionFilter &is_volatile,
                                 std::size_t &num_tensors) {
  if (node.leaf()) {
    num_tensors += node->is_tensor();
    return !is_volatile.matches(*node->expr());
  }

  return count_tensors_if_persistent(node.left(), is_volatile, num_tensors) &&
         count_tensors_if_persistent(node.right(), is_volatile, num_tensors);
}

std::size_t OptimizationStep::run(std::string_view step_id,
                                  ExecutionContext &ctx,
                                  const std::vector<std::string_view> &inputs) {
  if (!replay_.has_value()) {
    return OneByOneProcessingStep<ExpressionData>::run(step_id, ctx, inputs);
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
  std::size_t n_inputs = 0;
  for (std::string_view input : inputs) {
    for (const ExecutionContext::Data<const ProcessingData> &current :
         std::as_const(ctx).get_data(input)) {
      for (const ResultExpr &expr :
           convert_data<ExpressionData>(current.data.get()).expressions) {
        ResultExpr copy = expr.clone();
        symmetrizers.push_back(optimize_unsymmetrized(copy, options_));
        trees.push_back(to_tree(copy));
        tree_to_input.push_back(n_inputs);
      }

      ++n_inputs;
    }
  }

  auto is_persistent = [this](const Node &node) {
    std::size_t num_tensors = 0;
    if (!count_tensors_if_persistent(node, *replay_->is_volatile,
                                     num_tensors)) {
      return false;
    }

    // Same threshold as in the cse step: extracting a single tensor is not
    // worth it
    return num_tensors >= 2;
  };

  // Maximal persistent subtrees (their subtrees are not visited)
  SubexpressionUsageCounts<Node> persistent;
  for (const Node &tree : trees) {
    tree.visit_internal([&](const Node &node) {
      if (!is_persistent(node)) {
        return true;
      }

      ++persistent[node];
      return false;
    });
  }

  std::vector<std::size_t> definitions;
  if (!persistent.empty()) {
    auto label_gen = [this](const Node &, std::size_t counter) {
      return replay_->intermediate_label + std::to_string(counter);
    };

    opt::cse::SubexpressionReplacer<std::vector<Node>, Node, decltype(to_tree),
                                    decltype(label_gen)>
        replacer(trees, persistent, to_tree, label_gen);
    replacer.perform_replacements(0, trees.size());

    definitions = replacer.cse_indices();
  }

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
    reattach_symmetrizer(result, std::move(symmetrizers.at(expr_idx)));
    volatile_part.expressions.push_back(std::move(result));

    ++expr_idx;
  }
  SEQUANT_ASSERT(unassigned.empty());
  SEQUANT_ASSERT(expr_idx == tree_to_input.size());

  pending_.assign(std::make_move_iterator(parts.begin()),
                  std::make_move_iterator(parts.end()));

  const std::size_t n_outputs =
      OneByOneProcessingStep<ExpressionData>::run(step_id, ctx, inputs);
  SEQUANT_ASSERT(pending_.empty());
  SEQUANT_ASSERT(n_outputs == 2 * n_inputs);

  // Expose the two parts of every alias created above as sub-IDs
  std::map<std::string, std::vector<std::size_t>> alias_inputs;
  std::size_t input_idx = 0;
  for (std::string_view input : inputs) {
    for (const ExecutionContext::Data<const ProcessingData> &current :
         std::as_const(ctx).get_data(input)) {
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

      ++input_idx;
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
                                      const ExpressionData &data) {
  if (replay_.has_value()) {
    SEQUANT_ASSERT(!pending_.empty());
    auto [persistent_part, volatile_part] = std::move(pending_.front());
    pending_.pop_front();

    ctx.set_data(id_prefix, id_start, std::move(persistent_part));
    ctx.set_data(id_prefix, id_start + 1, std::move(volatile_part));

    return 2;
  }

  ExpressionData result;
  result.expressions.reserve(data.expressions.size());
  for (const ResultExpr &input : data.expressions) {
    result.expressions.emplace_back(input.clone());

    std::optional<ExprPtr> symmetrizer =
        optimize_unsymmetrized(result.expressions.back(), options_);
    reattach_symmetrizer(result.expressions.back(), std::move(symmetrizer));
  }

  ctx.set_data(id_prefix, id_start, std::move(result));

  return 1;
}

}  // namespace sequant::util::extint
