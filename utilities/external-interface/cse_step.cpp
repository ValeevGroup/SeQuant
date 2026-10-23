#include "cse_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <range/v3/view/concat.hpp>
#include <range/v3/view/repeat_n.hpp>

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <functional>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CSEStep);

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
      merge_all_ = false;
      merge_groups_.clear();

      if (value.is_boolean()) {
        merge_all_ = value.get<bool>();
      } else if (value.is_array() &&
                 std::ranges::all_of(value, [](const nlohmann::json &id) {
                   return id.is_string() && !id.get<std::string>().empty();
                 })) {
        merge_groups_ = value.get<std::vector<std::string>>();
      } else {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " has to be a boolean or an array of non-empty IDs");
      }
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
  if (std::ranges::empty(exprs)) {
    return {};
  }

  ExportTreeData output;
  output.entries.reserve(exprs.size());

  opt::CSEOptions<ExportNode<>> cse_opts;
  cse_opts.filter_predicate = opt::min_usage_filter<ExportNode<>>(min_usage);

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
      std::vector<std::size_t> dependencies;

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
        const ExportTreeData::Entry &assoc_entry =
            assoc_tree_data.entries.at(entry_offset);
        symm_target = assoc_entry.symm_contribution_target;
        // Preserve dependencies set on the original expression (e.g. by a
        // preceding dependencies step). CSE only introduces new intermediate
        // trees and leaves the root IDs these refer to intact.
        dependencies = assoc_entry.dependencies;
        ++entry_offset;
      } else if (!cse_positions.empty()) {
        cse_positions.pop_back();
        ++expr_offset;
        overwrite = true;
      }

      std::vector<Index> indices;
      if (exprs[i]->is_tensor()) {
        auto tensor_indices = exprs[i]->as_tensor().const_indices();
        // Only batch over indices that are actually present in the result
        std::ranges::copy_if(current_batch_indices, std::back_inserter(indices),
                             [&tensor_indices](const Index &idx) {
                               return std::ranges::find(tensor_indices, idx) !=
                                      tensor_indices.end();
                             });
      } else {
        SEQUANT_ASSERT(current_batch_indices.empty());
      }

      ExportTreeData::Entry current{
          .tree = std::move(exprs[i]),
          .symm_contribution_target = std::move(symm_target),
          .overwrite_previous = std::move(overwrite),
          .dependencies = std::move(dependencies),
          .batch_indices = std::move(indices)};

      output.entries.emplace_back(std::move(current));
    }

    SEQUANT_ASSERT(output.entries.size() == current_batch.second);

    auto batch_entries =
        std::ranges::subrange(output.entries.begin() + current_batch.first,
                              output.entries.begin() + current_batch.second);

    // Sort by used batch indices per expression. Batch indices are treated as
    // having a "priority". We start at the highest priority (first) batch index
    // and compare on the first batch index that only one of the expressions
    // has. Expressions with the lower priority batching indices (or even none)
    // come before the ones with the higher priority (or more) batching indices.
    // Since all non-CSEs are chose to have the same batching indices,
    // their relative order remains unchanged. And since CSEs can't depend
    // on each other (due to the way we determine them), their order is
    // irrelevant as long as they are computed before they are used, which
    // is guaranteed because the CSE can at most have as many batched
    // indices as the associated main expression and hence they compare
    // equivalent with the used comparator (and stable_sort guarantees
    // retained relative order of equivalent elements).
    std::ranges::stable_sort(
        batch_entries,
        [&current_batch_indices](const ExportNode<> &lhs,
                                 const ExportNode<> &rhs) {
          if (!rhs->is_tensor()) {
            return true;
          }
          if (!lhs->is_tensor()) {
            return std::ranges::any_of(
                rhs->as_tensor().const_indices(), [&](const Index &idx) {
                  return std::ranges::find(current_batch_indices, idx) !=
                         current_batch_indices.end();
                });
          }

          auto lhs_indices = lhs->as_tensor().const_indices();
          auto rhs_indices = rhs->as_tensor().const_indices();

          for (const Index &idx : current_batch_indices) {
            const bool lhs_contains =
                std::ranges::find(lhs_indices, idx) != lhs_indices.end();
            const bool rhs_contains =
                std::ranges::find(rhs_indices, idx) != rhs_indices.end();

            if (lhs_contains != rhs_contains) {
              return rhs_contains;
            }
          }

          return false;
        },
        &ExportTreeData::Entry::tree);
  }

  return output;
}

decltype(auto) to_tree_data(
    const ExecutionContext::Data<ProcessingData> &data_obj) {
  return convert_data<ExportTreeData>(data_obj.data.get());
}

std::optional<ExportTreeData> perform_merged_cse(
    const std::vector<ExecutionContext::Data<ProcessingData>> &datasets,
    const std::vector<std::size_t> &members, std::size_t min_usage) {
  std::vector<ExecutionContext::Data<ProcessingData>> inputs;
  std::vector<ExportNode<>> expressions;
  std::vector<std::size_t> expr_to_input;

  for (std::size_t member : members) {
    const ExportTreeData &data = to_tree_data(datasets[member]);

    for (std::size_t i = 0; i < data.entries.size(); ++i) {
      expr_to_input.emplace_back(inputs.size());
    }

    std::ranges::transform(data.entries, std::back_inserter(expressions),
                           &ExportTreeData::Entry::tree);

    inputs.push_back(datasets[member]);
  }

  SEQUANT_ASSERT(expressions.size() == expr_to_input.size());

  return perform_cse(expressions, inputs | std::views::transform(&to_tree_data),
                     expr_to_input, min_usage);
}

std::size_t CSEStep::run(std::string_view step_id, ExecutionContext &ctx,
                         const std::vector<std::string_view> &input_ids) {
  if (!merge_all_ && merge_groups_.empty()) {
    return OneByOneProcessingStep<ExportTreeData>::run(step_id, ctx, input_ids);
  }

  std::vector<ExecutionContext::Data<ProcessingData>> datasets;
  for (std::string_view current_input : input_ids) {
    for (const ExecutionContext::Data<ProcessingData> &current :
         ctx.get_data(current_input)) {
      datasets.push_back(current);
    }
  }

  // Index into merge_groups_ for each dataset (merge_groups_.size() if the
  // dataset isn't part of any group and thus processed individually). With
  // merge_all_, everything belongs to the single group 0.
  const std::size_t ungrouped = merge_all_ ? 1 : merge_groups_.size();
  std::vector<std::size_t> group_of(datasets.size(), ungrouped);
  std::vector<std::vector<std::size_t>> members(ungrouped);

  for (std::size_t i = 0; i < datasets.size(); ++i) {
    if (merge_all_) {
      group_of[i] = 0;
    } else {
      for (std::string_view id : ranges::views::concat(
               datasets[i].associated_ids, datasets[i].associated_group_ids)) {
        auto it =
            std::ranges::find(merge_groups_, detail::strip_step_prefix(id));
        if (it == merge_groups_.end()) {
          continue;
        }

        std::size_t group = it - merge_groups_.begin();
        if (group_of[i] != ungrouped && group_of[i] != group) {
          throw Exception("Input '" +
                          std::string(datasets[i].associated_ids.front()) +
                          "' matches multiple merge_inputs groups of " +
                          kind() + ": '" + merge_groups_[group_of[i]] +
                          "' and '" + merge_groups_[group] + "'");
        }
        group_of[i] = group;
      }
    }

    if (group_of[i] != ungrouped) {
      members[group_of[i]].push_back(i);
    }
  }

  detail::AliasPropagator propagator(ctx, step_id);
  std::vector<bool> group_done(members.size(), false);
  std::size_t total_outputs = 0;

  for (std::size_t i = 0; i < datasets.size(); ++i) {
    const ExecutionContext::Data<ProcessingData> &current = datasets[i];
    const std::size_t group = group_of[i];

    if (group == ungrouped) {
      std::size_t new_outputs =
          process(step_id, total_outputs, ctx, to_tree_data(current));

      const bool alias_source = new_outputs == 0;
      if (alias_source) {
        new_outputs = 1;
      }

      propagator.record(current.associated_ids, current.associated_group_ids,
                        total_outputs, new_outputs, alias_source);
      total_outputs += new_outputs;
      continue;
    }

    if (group_done[group]) {
      continue;
    }
    group_done[group] = true;

    std::optional<ExportTreeData> result =
        perform_merged_cse(datasets, members[group], min_usage_);

    const std::string output_id =
        std::string(step_id) + "." + std::to_string(total_outputs);

    if (result.has_value()) {
      ctx.set_data(step_id, total_outputs, std::move(result.value()));
    } else if (merge_all_) {
      ctx.add_data_alias(input_ids, output_id);
    } else {
      std::vector<std::string_view> member_ids;
      for (std::size_t member : members[group]) {
        member_ids.push_back(datasets[member].associated_ids.front());
      }
      ctx.add_data_alias(member_ids, output_id);
    }

    if (!merge_all_) {
      const std::string group_alias =
          std::string(step_id) + "." + merge_groups_[group];
      if (ctx.has_data(group_alias)) {
        throw Exception("Duplicate group ID name '" + merge_groups_[group] +
                        "' leading to ambiguity for '" + group_alias + "'");
      }
      ctx.add_data_alias(output_id, group_alias);
    }

    ++total_outputs;
  }

  propagator.flush();

  return total_outputs;
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
