#include "validate_step.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ValidateStep, "validate");

std::string ValidateStep::kind() const { return "validate"; }

bool ValidateStep::accepts_options() const { return false; }

bool ValidateStep::requires_options() const { return false; }

void ValidateStep::set_options(const nlohmann::json &) {
  throw Exception("validate doesn't take any options");
}

std::size_t ValidateStep::process(std::string_view, std::size_t,
                                  ExecutionContext &,
                                  const ExpressionData &data) {
  std::size_t expr_counter = 1;
  for (const ResultExpr &expr : data.expressions) {
    std::string msg;
    if (!is_valid(expr, &msg)) {
      throw Exception("Expression #" + std::to_string(expr_counter) +
                      " is invalid: " + msg);
    }
    ++expr_counter;
  }

  return 0;
}

}  // namespace sequant::util::extint
