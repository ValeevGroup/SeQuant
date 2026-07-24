#include "validate_step.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>

#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ValidateStep, "validate");

std::string ValidateStep::kind() const { return "validate"; }

bool ValidateStep::accepts_options() const { return false; }

bool ValidateStep::requires_options() const { return false; }

void ValidateStep::set_options(const nlohmann::json &) {
  throw Exception("validate doesn't take any options");
}

std::size_t ValidateStep::run(std::string_view, ExecutionContext &ctx,
                              const std::vector<std::string_view> &inputs) {
  for (std::string_view current_input : inputs) {
    std::size_t outer_counter = 1;
    for (const ProcessingData &current_data : ctx.get_data(current_input)) {
      const ExpressionData &data = convert_data<ExpressionData>(current_data);

      std::size_t inner_counter = 1;
      for (const ResultExpr &expr : data.expressions) {
        std::string msg;
        if (!is_valid(expr, &msg)) {
          throw Exception("Expression " + std::string(current_input) + " #" +
                          std::to_string(outer_counter) + "." +
                          std::to_string(inner_counter) +
                          " is invalid: " + msg);
        }
        ++inner_counter;
      }
      ++outer_counter;
    }
  }

  return 0;
}

}  // namespace sequant::util::extint
