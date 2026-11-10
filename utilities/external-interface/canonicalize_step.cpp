#include "canonicalize_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(CanonicalizeStep);

std::string CanonicalizeStep::kind() const { return "canonicalize"; }

bool CanonicalizeStep::accepts_options() const { return false; }

bool CanonicalizeStep::requires_options() const { return false; }

void CanonicalizeStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t CanonicalizeStep::process(std::string_view id_prefix,
                                      std::size_t id_start,
                                      ExecutionContext &ctx,
                                      WorkingCopy<ExpressionData> input) {
  for (ResultExpr &expr : input.data.expressions) {
    canonicalize(expr);
  }

  if (input.unchanged()) {
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(input.data));

  return 1;
}

bool CanonicalizeStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
