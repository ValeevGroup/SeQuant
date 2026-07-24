#include "output_step.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <iostream>
#include <string>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(OutputStep, "output");

std::string OutputStep::kind() const { return "Output"; }

bool OutputStep::accepts_options() const { return true; }

bool OutputStep::requires_options() const { return false; }

void OutputStep::set_options(const nlohmann::json &options) {
  if (!options.is_object()) {
    throw Exception("OutputStep expects a JSON object for its options!");
  }

  for (const auto &[key, value] : options.items()) {
    if (key == "format") {
      if (!value.is_string()) {
        throw Exception(
            "Option 'format' for OutputStep is expected to be a string");
      }

      if (value == "latex") {
        latex_ = true;
      } else if (value == "serialize") {
        latex_ = false;
      } else {
        throw Exception("Invalid output format for OutputStep: '" +
                        value.get<std::string>() + "'");
      }
    } else if (key == "annotate_symmetry") {
      if (!value.is_boolean()) {
        throw Exception(
            "Option 'annotate_symmetry' for OutputStep is expected to be a "
            "boolean");
      }

      annot_symm_ = value.get<bool>();
    } else {
      throw Exception("Unknown option key for OutputStep: '" + key + "'");
    }
  }
}

std::size_t OutputStep::run(std::string_view, ExecutionContext &ctx,
                            const std::vector<std::string_view> &inputs) {
  for (std::string_view current_input : inputs) {
    for (const ProcessingData &current_data : ctx.get_data(current_input)) {
      const ExpressionData &data = convert_data<ExpressionData>(current_data);

      for (const auto &expr : data.expressions) {
        if (latex_) {
          std::wcout << io::latex::to_string(expr) << std::endl;
        } else {
          std::wcout << io::serialization::to_string(
                            expr, {.annot_symm = annot_symm_})
                     << std::endl;
        }
      }
    }
  }

  return 0;
}

}  // namespace sequant::util::extint
