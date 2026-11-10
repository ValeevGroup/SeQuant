#include "sort_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"

#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(SortStep);

std::string SortStep::kind() const { return "sort"; }

bool SortStep::accepts_options() const { return false; }

bool SortStep::requires_options() const { return false; }

void SortStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t SortStep::process(std::string_view id_prefix, std::size_t id_start,
                              ExecutionContext &ctx,
                              WorkingCopy<ExpressionData> &input) {
  const LexicographicExprCompare cmp;

  for (ResultExpr &expr : input.data.expressions) {
    sort_expr(expr, cmp);
  }

  if (input.unchanged()) {
    return 0;
  }

  ctx.set_data(id_prefix, id_start, std::move(input.data));

  return 1;
}

bool SortStep::alias_unchanged_inputs() const { return true; }

}  // namespace sequant::util::extint
