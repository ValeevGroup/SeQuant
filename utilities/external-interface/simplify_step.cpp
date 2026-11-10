#include "simplify_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(SimplifyStep);

std::string SimplifyStep::kind() const { return "simplify"; }

bool SimplifyStep::accepts_options() const { return false; }

bool SimplifyStep::requires_options() const { return false; }

void SimplifyStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t SimplifyStep::process(std::string_view id_prefix,
                                  std::size_t id_start, ExecutionContext &ctx,
                                  WorkingCopy<ExpressionData> input) {
  for (ResultExpr &expr : input.data.expressions) {
    simplify(expr);
  }

  if (input.unchanged()) {
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(input.data));

  return 1;
}

bool SimplifyStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
