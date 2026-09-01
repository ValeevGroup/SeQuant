#include "cse_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <functional>
#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CSEStep, "cse");

std::string CSEStep::kind() const { return "cse"; }

bool CSEStep::accepts_options() const { return true; }

bool CSEStep::requires_options() const { return false; }

void CSEStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "min_usage") {
      if (!value.is_number_integer()) {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " expects an unsigned integer");
      }

      min_usage_ = value.get<std::size_t>();
    } else if (key == "merge_inputs") {
      if (!value.is_boolean()) {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " has to be a boolean");
      }

      merge_inputs_ = value.get<bool>();
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

template <std::ranges::range ExprRange>
std::vector<std::size_t> perform_cse(ExprRange &&exprs, std::size_t min_usage) {
  opt::CSEOptions<ExportNode<>> cse_opts;
  cse_opts.filter_predicate = [min_usage](const ExportNode<> &tree,
                                          std::size_t usage_count) {
    if (usage_count < min_usage) {
      return false;
    }

    std::size_t num_tensors = 0;
    tree.visit_leaf([&num_tensors](const ExportNode<> &node) {
      num_tensors += node->is_tensor();
    });

    if (num_tensors < 2) {
      // A subexpression that contains less than two tensors is not
      // worth the hassle of creating, storing and reusing it.
      // Specifically, this avoids CSE in symmetrization expressions
      // such as
      // 1/2 * R2u:eecc[abij] + 1/2 * R2u:eecc[baji]
      // where 1/2 * R2u:eecc would be the kind of subexpression we
      // don't want)
      return false;
    }

    return true;
  };

  return opt::eliminate_common_subexpressions(
      exprs,
      [](const auto &expr) {
        // Note: the lambda is needed to make the callable usable for
        // ExprPtr as well as ResultExpr objects
        return to_export_tree(expr);
      },
      cse_opts);
}

std::size_t CSEStep::run(std::string_view step_id, ExecutionContext &ctx,
                         const std::vector<std::string_view> &input_ids) {
  if (!merge_inputs_) {
    return OneByOneProcessingStep<ExportTreeData>::run(step_id, ctx, input_ids);
  }

  std::vector<ExecutionContext::Data<ProcessingData>> inputs;
  std::vector<ExportNode<>> expressions;
  std::vector<std::size_t> expr_to_input;

  for (std::string_view current_input : input_ids) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         ctx.get_data(current_input)) {
      const ExportTreeData &data =
          convert_data<ExportTreeData>(current.data.get());

      for (std::size_t i = 0; i < data.entries.size(); ++i) {
        expr_to_input.emplace_back(inputs.size());
      }

      std::ranges::transform(data.entries, std::back_inserter(expressions),
                             &ExportTreeData::Entry::tree);

      inputs.emplace_back(std::move(current));
    }
  }

  SEQUANT_ASSERT(expressions.size() == expr_to_input.size());

  std::vector<std::size_t> cse_positions = perform_cse(expressions, min_usage_);

  if (cse_positions.empty()) {
    ctx.add_data_alias(input_ids, std::string(step_id) + ".0");
    return 1;  // we "produced" an alias
  }

  std::ranges::sort(cse_positions, std::greater<>{});

  ExportTreeData output;
  output.entries.reserve(expressions.size());

  std::size_t expr_offset = 0;
  std::size_t last_input_idx = 0;
  std::size_t entry_offset = 0;

  for (std::size_t i = 0; i < expressions.size(); ++i) {
    std::optional<Tensor> symm_target;

    std::optional<bool> overwrite;
    if (cse_positions.empty() || cse_positions.back() != i) {
      SEQUANT_ASSERT(i >= expr_offset);
      const std::size_t expr_idx = i - expr_offset;
      const std::size_t input_idx = expr_to_input.at(expr_idx);
      const ExecutionContext::Data<ProcessingData> &assoc_data =
          inputs.at(input_idx);

      const ExportTreeData &assoc_tree_data =
          convert_data<ExportTreeData>(assoc_data.data.get());

      if (last_input_idx != input_idx) {
        // Reset offset to ensure it is a per-input offset into entries
        entry_offset = 0;
        last_input_idx = input_idx;
      }

      SEQUANT_ASSERT(entry_offset < assoc_tree_data.entries.size());
      symm_target =
          assoc_tree_data.entries.at(entry_offset).symm_contribution_target;
      ++entry_offset;
    } else if (!cse_positions.empty()) {
      cse_positions.pop_back();
      ++expr_offset;
      overwrite = true;
    }

    ExportTreeData::Entry current{
        .tree = std::move(expressions.at(i)),
        .symm_contribution_target = std::move(symm_target),
        .overwrite_previous = std::move(overwrite)};

    output.entries.emplace_back(std::move(current));
  }

  ctx.set_data(step_id, 0, std::move(output));

  return 1;
}

std::size_t CSEStep::process(std::string_view id_prefix, std::size_t id_start,
                             ExecutionContext &ctx,
                             const ExportTreeData &input) {
  std::vector<ExportNode<>> expressions;
  expressions.reserve(std::ranges::size(input.entries));
  std::ranges::transform(input.entries, std::back_inserter(expressions),
                         &ExportTreeData::Entry::tree);

  std::vector<std::size_t> cse_positions = perform_cse(expressions, min_usage_);

  if (cse_positions.empty()) {
    return 0;
  }

  std::ranges::sort(cse_positions, std::greater<>{});

  ExportTreeData output;
  output.entries.reserve(expressions.size());

  std::size_t offset = 0;

  for (std::size_t i = 0; i < expressions.size(); ++i) {
    std::optional<Tensor> symm_target;

    std::optional<bool> overwrite;
    if (cse_positions.empty() || cse_positions.back() != i) {
      SEQUANT_ASSERT(i >= offset);
      symm_target = input.entries.at(i - offset).symm_contribution_target;
    } else if (!cse_positions.empty()) {
      cse_positions.pop_back();
      ++offset;
      overwrite = true;
    }

    ExportTreeData::Entry current{
        .tree = std::move(expressions.at(i)),
        .symm_contribution_target = std::move(symm_target),
        .overwrite_previous = std::move(overwrite)};

    output.entries.emplace_back(std::move(current));
  }

  ctx.set_data(id_prefix, id_start, std::move(output));

  return 1;
}

bool CSEStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
