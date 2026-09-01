#include "cse_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <range/v3/view/repeat_n.hpp>

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

template <std::ranges::random_access_range ExprRange,
          std::ranges::random_access_range InputRange,
          std::ranges::random_access_range Expr2Inp>
std::optional<ExportTreeData> perform_cse(ExprRange &&exprs,
                                          InputRange &&inputs,
                                          Expr2Inp &&expr_to_input,
                                          std::size_t min_usage) {
  ExportTreeData output;
  output.entries.reserve(exprs.size());

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

  std::vector<std::size_t> cse_positions;
  std::vector<std::vector<Index>> used_batch_indices;
  std::vector<std::pair<std::size_t, std::size_t>> batches;

  // Determine CSEs in batches of expressions that share the same
  // batching indices
  std::pair<std::size_t, std::size_t> partition = {0, 0};
  while (partition.second < exprs.size()) {
    partition.first = partition.second;

    SEQUANT_ASSERT(partition.first >= cse_positions.size());
    auto batch_indices =
        inputs[expr_to_input[partition.first - cse_positions.size()]]
            .entries.front()
            .batch_indices;

    // Note: we explicitly also check the first ==  second case
    // as the first set of entries could already have different indices
    for (; partition.second < exprs.size(); ++partition.second) {
      SEQUANT_ASSERT(partition.second >= cse_positions.size());
      const auto &entries =
          inputs[expr_to_input[partition.second - cse_positions.size()]]
              .entries;
      if (std::ranges::any_of(
              entries,
              [&](const auto &indices) { return indices != batch_indices; },
              &ExportTreeData::Entry::batch_indices)) {
        break;
      }
    }

    if (partition.first == partition.second) {
      ++partition.second;
      batch_indices.clear();
    }

    cse_opts.batch_indices = std::move(batch_indices);

    std::vector<std::size_t> new_cses = opt::eliminate_common_subexpressions(
        exprs, partition.first, partition.second,
        [](const auto &expr) {
          // Note: the lambda is needed to make the callable usable for
          // ExprPtr as well as ResultExpr objects
          return to_export_tree(expr);
        },
        cse_opts);

    cse_positions.insert(cse_positions.end(), new_cses.begin(), new_cses.end());

    partition.second += new_cses.size();

    used_batch_indices.emplace_back(std::move(cse_opts.batch_indices));
    batches.emplace_back(partition);
  }

  SEQUANT_ASSERT(batches.back().second == exprs.size());
  SEQUANT_ASSERT(batches.size() == used_batch_indices.size());

  if (cse_positions.empty()) {
    return {};
  }

  std::ranges::sort(cse_positions, std::greater<>{});

  std::size_t expr_offset = 0;
  std::size_t last_input_idx = 0;
  std::size_t entry_offset = 0;

  for (std::size_t batch = 0; batch < batches.size(); ++batch) {
    const auto &current_batch_indices = used_batch_indices.at(batch);
    const auto &current_batch = batches.at(batch);

    for (std::size_t i = current_batch.first; i < current_batch.second; ++i) {
      std::optional<Tensor> symm_target;

      std::optional<bool> overwrite;
      if (cse_positions.empty() || cse_positions.back() != i) {
        SEQUANT_ASSERT(i >= expr_offset);
        const std::size_t expr_idx = i - expr_offset;
        const std::size_t input_idx = expr_to_input[expr_idx];

        const ExportTreeData &assoc_tree_data = inputs[input_idx];

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

      std::vector<Index> indices;
      if (exprs[i]->is_tensor()) {
        // Only batch over indices that are actually present in the result
        std::ranges::copy_if(
            exprs[i]->as_tensor().const_indices(), std::back_inserter(indices),
            [&current_batch_indices](const Index &idx) {
              return std::ranges::find(current_batch_indices, idx) !=
                     current_batch_indices.end();
            });
      } else {
        SEQUANT_ASSERT(current_batch_indices.empty());
      }

      ExportTreeData::Entry current{
          .tree = std::move(exprs[i]),
          .symm_contribution_target = std::move(symm_target),
          .overwrite_previous = std::move(overwrite),
          .batch_indices = std::move(indices)};

      output.entries.emplace_back(std::move(current));
    }

    // TODO: stable_sort output entry to group by batching index
  }

  return output;
}

decltype(auto) to_tree_data(
    const ExecutionContext::Data<ProcessingData> &data_obj) {
  return convert_data<ExportTreeData>(data_obj.data.get());
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

  std::optional<ExportTreeData> result =
      perform_cse(expressions, inputs | std::views::transform(&to_tree_data),
                  expr_to_input, min_usage_);

  if (!result.has_value()) {
    ctx.add_data_alias(input_ids, std::string(step_id) + ".0");
    return 1;  // we "produced" an alias
  }

  ctx.set_data(step_id, 0, std::move(result.value()));

  return 1;
}

std::size_t CSEStep::process(std::string_view id_prefix, std::size_t id_start,
                             ExecutionContext &ctx,
                             const ExportTreeData &input) {
  std::vector<ExportNode<>> expressions;
  expressions.reserve(std::ranges::size(input.entries));
  std::ranges::transform(input.entries, std::back_inserter(expressions),
                         &ExportTreeData::Entry::tree);

  std::optional<ExportTreeData> result =
      perform_cse(expressions, std::ranges::subrange(&input, &input + 1),
                  ranges::views::repeat_n(0, input.entries.size()), min_usage_);

  if (!result.has_value()) {
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(result.value()));

  return 1;
}

bool CSEStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
