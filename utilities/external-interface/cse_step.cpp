#include "cse_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <range/v3/view/concat.hpp>
#include <range/v3/view/repeat_n.hpp>

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/common_subexpression_elimination.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/topological.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CSEStep);

std::string CSEStep::kind() const { return "cse"; }

bool CSEStep::accepts_options() const { return true; }

bool CSEStep::requires_options() const { return false; }

namespace {

bool is_valid_id(const nlohmann::json &id) {
  return id.is_string() && !id.get<std::string>().empty();
}

bool is_valid_merge_group(const nlohmann::json &group) {
  return is_valid_id(group) || (group.is_array() && !group.empty() &&
                                std::ranges::all_of(group, is_valid_id));
}

std::string describe_merge_group(const std::vector<std::string> &ids) {
  return ids.front() + (ids.size() > 1 ? ", ..." : "");
}

}  // namespace

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
                 std::ranges::all_of(value, is_valid_merge_group)) {
        for (const nlohmann::json &group : value) {
          merge_groups_.push_back(group.is_string()
                                      ? std::vector{group.get<std::string>()}
                                      : group.get<std::vector<std::string>>());
        }

        std::set<std::string_view> seen;
        for (const std::vector<std::string> &ids : merge_groups_) {
          for (const std::string &id : ids) {
            if (!seen.insert(id).second) {
              throw Exception("ID '" + id + "' in option '" + key + "' for " +
                              kind() + " is listed more than once");
            }
          }
        }
      } else {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " has to be a boolean or an array of non-empty IDs or "
                        "non-empty arrays of non-empty IDs");
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

namespace {
/// Moves the trees out of the entries of @p data (which keep all other
/// information) and appends them to @p trees. Their IDs are appended to
/// @p root_ids.
void take_trees(ExportTreeData &data, std::vector<ExportNode<>> &trees,
                std::vector<std::size_t> &root_ids) {
  for (ExportTreeData::Entry &entry : data.entries) {
    root_ids.push_back(entry.tree->id());
    trees.push_back(std::move(entry.tree));
  }
}

/// Counterpart to take_trees. Moves the trees from the range starting at
/// @p trees back into the entries of @p data, restoring their original IDs.
void restore_trees(ExportTreeData &data, std::ranges::range auto &&trees,
                   std::ranges::range auto &&root_ids) {
  auto tree_it = std::ranges::begin(trees);
  auto id_it = std::ranges::begin(root_ids);

  for (ExportTreeData::Entry &entry : data.entries) {
    SEQUANT_ASSERT(tree_it != std::ranges::end(trees));
    (*tree_it)->set_id(*id_it);
    entry.tree = std::move(*tree_it);
    ++tree_it;
    ++id_it;
  }
}
}  // namespace

template <std::ranges::random_access_range ExprRange,
          std::ranges::random_access_range InputRange,
          std::ranges::random_access_range Expr2Inp>
std::optional<ExportTreeData> perform_cse(
    ExprRange &&exprs, InputRange &&inputs, Expr2Inp &&expr_to_input,
    const std::vector<std::size_t> &root_ids, std::size_t min_usage) {
  if (std::ranges::empty(exprs)) {
    return {};
  }

  ExportTreeData output;
  output.entries.reserve(exprs.size());

  opt::CSEOptions<ExportNode<>> cse_opts;
  cse_opts.filter_predicate = opt::min_usage_filter<ExportNode<>>(min_usage);

  // Expressions can only share subexpressions if they are batched identically
  std::vector<std::optional<std::vector<Index>>> expr_batch_indices;
  expr_batch_indices.reserve(exprs.size());
  for (std::size_t input_idx : expr_to_input) {
    const auto &entries = inputs[input_idx].entries;
    const auto &batch_indices = entries.front().batch_indices;

    if (std::ranges::all_of(
            entries,
            [&](const auto &indices) { return indices == batch_indices; },
            &ExportTreeData::Entry::batch_indices)) {
      expr_batch_indices.emplace_back(batch_indices);
    } else {
      expr_batch_indices.emplace_back(std::nullopt);
    }
  }

  opt::BatchedCSEResult cse_result =
      opt::eliminate_common_subexpressions_batched(
          exprs, expr_batch_indices,
          [](const auto &expr) {
            // Note: the lambda is needed to make the callable usable for
            // ExprPtr as well as ResultExpr objects
            return to_export_tree(expr);
          },
          cse_opts);

  std::vector<std::size_t> &cse_positions = cse_result.cse_positions;
  const auto &batches = cse_result.partitions;
  const auto &used_batch_indices = cse_result.batch_indices;

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
        // Rewriting a root (e.g. replacing it by a CSE) creates a new node.
        // Restore the original ID, as other entries may depend on it.
        exprs[i]->set_id(root_ids.at(expr_idx));
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

    // Since all non-CSEs are chosen to have the same batching indices, their
    // relative order remains unchanged. And since CSEs can't depend on each
    // other (due to the way we determine them), their order is irrelevant as
    // long as they are computed before they are used, which is guaranteed
    // because the CSE can at most have as many batched indices as the
    // associated main expression and hence they compare equivalent with the
    // used comparator (and stable_sort guarantees retained relative order of
    // equivalent elements).
    std::ranges::stable_sort(
        batch_entries,
        [&current_batch_indices](const ExportNode<> &lhs,
                                 const ExportNode<> &rhs) {
          return opt::batch_priority_less(lhs, rhs, current_batch_indices);
        },
        &ExportTreeData::Entry::tree);
  }

  return output;
}

std::optional<ExportTreeData> perform_merged_cse(
    std::vector<WorkingCopy<ExportTreeData>> &inputs, std::size_t min_usage) {
  std::vector<ExportNode<>> expressions;
  std::vector<std::size_t> root_ids;
  std::vector<std::size_t> expr_to_input;

  for (std::size_t i = 0; i < inputs.size(); ++i) {
    expr_to_input.insert(expr_to_input.end(), inputs[i].data.entries.size(), i);
    take_trees(inputs[i].data, expressions, root_ids);
  }

  SEQUANT_ASSERT(expressions.size() == expr_to_input.size());

  std::optional<ExportTreeData> result = perform_cse(
      expressions,
      inputs | std::views::transform(&WorkingCopy<ExportTreeData>::data),
      expr_to_input, root_ids, min_usage);

  if (!result.has_value()) {
    // The inputs remain unchanged
    std::size_t offset = 0;
    for (WorkingCopy<ExportTreeData> &input : inputs) {
      const std::size_t n = input.data.entries.size();
      restore_trees(input.data,
                    std::ranges::subrange(expressions.begin() + offset,
                                          expressions.begin() + offset + n),
                    std::ranges::subrange(root_ids.begin() + offset,
                                          root_ids.begin() + offset + n));
      offset += n;
    }
  }

  return result;
}

std::size_t CSEStep::run(std::string_view step_id, ExecutionContext &ctx,
                         const std::vector<std::string_view> &input_ids) {
  if (!merge_all_ && merge_groups_.empty()) {
    return OneByOneProcessingStep<ExportTreeData, true,
                                  InputAccess::WorkingCopy>::run(step_id, ctx,
                                                                 input_ids);
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
        const std::string_view stripped = detail::strip_step_prefix(id);
        auto it = std::ranges::find_if(
            merge_groups_, [stripped](const std::vector<std::string> &ids) {
              return std::ranges::find(ids, stripped) != ids.end();
            });
        if (it == merge_groups_.end()) {
          continue;
        }

        std::size_t group = it - merge_groups_.begin();
        if (group_of[i] != ungrouped && group_of[i] != group) {
          throw Exception(
              "Input '" + std::string(datasets[i].associated_ids.front()) +
              "' matches multiple merge_inputs groups of " + kind() + ": '" +
              describe_merge_group(merge_groups_[group_of[i]]) + "' and '" +
              describe_merge_group(merge_groups_[group]) + "'");
        }
        group_of[i] = group;
      }
    }

    if (group_of[i] != ungrouped) {
      members[group_of[i]].push_back(i);
    }
  }

  // Merging puts the members into a single dataset, for which the export only
  // honors the intra-dataset order. Dependencies between members (e.g. set by a
  // dependencies step) would otherwise be lost, so order the members such that
  // producers precede their consumers.
  for (std::vector<std::size_t> &group_members : members) {
    std::map<std::size_t, std::size_t> producer_of;
    for (std::size_t pos = 0; pos < group_members.size(); ++pos) {
      for (const ExportTreeData::Entry &entry :
           convert_data<ExportTreeData>(datasets[group_members[pos]].data.get())
               .entries) {
        producer_of.emplace(entry.tree->id(), pos);
      }
    }

    std::vector<std::size_t> order =
        topological_order_indexed(group_members.size(), [&](std::size_t pos) {
          std::set<std::size_t> deps;
          for (const ExportTreeData::Entry &entry :
               convert_data<ExportTreeData>(
                   datasets[group_members[pos]].data.get())
                   .entries) {
            for (std::size_t id : entry.dependencies) {
              auto it = producer_of.find(id);
              if (it != producer_of.end() && it->second != pos) {
                deps.insert(it->second);
              }
            }
          }
          return deps;
        });

    std::vector<std::size_t> sorted;
    sorted.reserve(order.size());
    for (std::size_t pos : order) sorted.push_back(group_members[pos]);
    group_members = std::move(sorted);
  }

  detail::AliasPropagator propagator(ctx, step_id);
  std::vector<bool> group_done(members.size(), false);
  std::size_t total_outputs = 0;

  for (std::size_t i = 0; i < datasets.size(); ++i) {
    const ExecutionContext::Data<ProcessingData> &current = datasets[i];
    const std::size_t group = group_of[i];

    if (group == ungrouped) {
      total_outputs +=
          process_input(step_id, ctx, current, total_outputs, propagator);
      continue;
    }

    if (group_done[group]) {
      continue;
    }
    group_done[group] = true;

    std::vector<WorkingCopy<ExportTreeData>> group_inputs;
    for (std::size_t member : members[group]) {
      group_inputs.push_back(
          ctx.working_copy<ExportTreeData>(datasets[member]));
    }

    std::optional<ExportTreeData> result =
        perform_merged_cse(group_inputs, min_usage_);

    const std::string output_id =
        std::string(step_id) + "." + std::to_string(total_outputs);

    if (result.has_value()) {
      ctx.set_data(step_id, total_outputs, std::move(result.value()));
    } else if (!std::ranges::all_of(group_inputs,
                                    &WorkingCopy<ExportTreeData>::is_copy)) {
      // At least one input has been consumed and thus can't be aliased. All
      // inputs are passed on instead, as a group.
      std::vector<ProcessingData> passed_on;
      for (WorkingCopy<ExportTreeData> &input : group_inputs) {
        passed_on.emplace_back(std::move(input.data));
      }
      ctx.set_data(step_id, std::vector<std::size_t>{total_outputs},
                   std::move(passed_on));
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
      for (const std::string &group_id : merge_groups_[group]) {
        const std::string group_alias = std::string(step_id) + "." + group_id;
        if (ctx.has_data(group_alias)) {
          throw Exception("Duplicate group ID name '" + group_id +
                          "' leading to ambiguity for '" + group_alias + "'");
        }
        ctx.add_data_alias(output_id, group_alias);
      }
    }

    ++total_outputs;
  }

  propagator.flush();

  return total_outputs;
}

std::size_t CSEStep::process(std::string_view id_prefix, std::size_t id_start,
                             ExecutionContext &ctx,
                             WorkingCopy<ExportTreeData> &input) {
  std::vector<ExportNode<>> expressions;
  std::vector<std::size_t> root_ids;
  expressions.reserve(input.data.entries.size());
  root_ids.reserve(input.data.entries.size());
  take_trees(input.data, expressions, root_ids);

  std::optional<ExportTreeData> result = perform_cse(
      expressions, std::ranges::subrange(&input.data, &input.data + 1),
      ranges::views::repeat_n(0, input.data.entries.size()), root_ids,
      min_usage_);

  if (result.has_value()) {
    ctx.set_data(id_prefix, id_start, std::move(result.value()));
    return 1;
  }

  restore_trees(input.data, expressions, root_ids);

  return 0;
}

bool CSEStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
