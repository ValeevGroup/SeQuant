#include "index_batching_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(IndexBatchingStep, "batch_indices");

std::string IndexBatchingStep::kind() const { return "batch_indices"; }

bool IndexBatchingStep::accepts_options() const { return true; }

bool IndexBatchingStep::requires_options() const { return false; }

void IndexBatchingStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception(kind() + " expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "min_unbatched") {
      if (!value.is_number_integer()) {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " expects an unsigned integer");
      }

      min_unbatched_ = value.get<std::size_t>();
    } else if (key == "max_batched") {
      if (!value.is_number_integer()) {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " expects an unsigned integer");
      }

      max_batched_ = value.get<std::size_t>();
    } else if (key == "select_strategy") {
      if (!value.is_string()) {
        throw Exception("Value for option '" + key + "' for " + kind() +
                        " expects a string");
      }

      if (value == "largest") {
        start_with_largest_ = true;
      } else if (value == "smallest") {
        start_with_largest_ = false;
      } else {
        throw Exception("Unknown index selection strategy '" +
                        value.get<std::string>() + "'");
      }
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t IndexBatchingStep::process(std::string_view id_prefix,
                                       std::size_t id_start,
                                       ExecutionContext &ctx,
                                       const ExportTreeData &input) {
  ExportTreeData output;

  bool batched_some = false;
  for (const ExportTreeData::Entry &current : input.entries) {
    if (!current.tree->is_tensor()) {
      output.entries.emplace_back(current);
      continue;
    }

    const Tensor &result = current.tree->as_tensor();

    std::vector<Index> indices;
    indices.insert(indices.end(), result.const_indices().begin(),
                   result.const_indices().end());

    std::ranges::sort(indices, std::less<>{}, [](const Index &idx) {
      return idx.space().approximate_size();
    });

    std::vector<Index> batch_indices;
    if (start_with_largest_) {
      auto selected = indices | std::ranges::views::reverse |
                      std::ranges::views::drop(min_unbatched_);
      batch_indices.insert(batch_indices.end(), selected.begin(),
                           selected.end());
    } else {
      auto selected = indices | std::ranges::views::drop(min_unbatched_);
      batch_indices.insert(batch_indices.end(), selected.begin(),
                           selected.end());
    }

    if (batch_indices.size() > max_batched_) {
      batch_indices.resize(max_batched_);
    }

    batched_some |= !batch_indices.empty();

    output.entries.emplace_back(current);
    output.entries.back().batch_indices = std::move(batch_indices);
  }

  if (!batched_some) {
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(output));

  return 1;
}

bool IndexBatchingStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
