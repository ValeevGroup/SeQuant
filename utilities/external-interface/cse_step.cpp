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
    } else {
      throw Exception("Unknown option key for " + kind() + ": '" + key + "'");
    }
  }
}

std::size_t CSEStep::process(std::string_view id_prefix, std::size_t id_start,
                             ExecutionContext &ctx,
                             const ExportTreeData &input) {
  std::vector<ExportNode<>> expressions;
  expressions.reserve(input.entries.size());
  std::ranges::transform(input.entries, std::back_inserter(expressions),
                         &ExportTreeData::Entry::tree);

  opt::CSEOptions<ExportNode<>> cse_opts;
  cse_opts.filter_predicate = [this](const ExportNode<> &tree,
                                     std::size_t usage_count) {
    if (usage_count < min_usage_) {
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

  std::vector<std::size_t> cse_positions = opt::eliminate_common_subexpressions(
      expressions,
      [](const auto &expr) {
        // Note: the lambda is needed to make the callable usable for
        // ExprPtr as well as ResultExpr objects
        return to_export_tree(expr);
      },
      cse_opts);

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

}  // namespace sequant::util::extint
