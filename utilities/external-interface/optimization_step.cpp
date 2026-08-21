#include "optimization_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize/optimize.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(OptimizationStep, "optimize");

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

std::size_t OptimizationStep::process(std::string_view id_prefix,
                                      std::size_t id_start,
                                      ExecutionContext &ctx,
                                      const ExpressionData &data) {
  ExpressionData result;
  result.expressions.reserve(data.expressions.size());
  for (const ResultExpr &input : data.expressions) {
    result.expressions.emplace_back(input.clone());

    optimize(result.expressions.back(), options_);
  }

  ctx.set_data(id_prefix, id_start, std::move(result));

  return 1;
}

}  // namespace sequant::util::extint
