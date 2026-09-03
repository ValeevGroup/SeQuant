#include "to_export_tree_step.hpp"
#include "processing_data.hpp"
#include "processing_step_factory.hpp"
#include "utils.hpp"

#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/export/export_expr.hpp>
#include <SeQuant/core/export/export_node.hpp>
#include <SeQuant/core/expr.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

namespace sequant::util::extint {

SEQUANT_EXTINT_REGISTER_STEP_TYPE(ToExportTreeStep, "to_export_tree");

std::string ToExportTreeStep::kind() const { return "to_export_tree"; }

bool ToExportTreeStep::accepts_options() const { return false; }

bool ToExportTreeStep::requires_options() const { return false; }

void ToExportTreeStep::set_options(const nlohmann::json &) {
  throw Exception(kind() + " doesn't take any options");
}

std::size_t ToExportTreeStep::process(std::string_view id_prefix,
                                      std::size_t id_start,
                                      ExecutionContext &ctx,
                                      const ExpressionData &data) {
  ExportTreeData output;

  for (const ResultExpr &current : data.expressions) {
    std::optional<Tensor> symmetrized_result;

    ExportNode<> tree = [&]() {
      if (!needsSymmetrization(current.expression())) {
        return to_export_tree(current);
      }

      ResultExpr copy = current.clone();
      [[maybe_unused]] std::optional<ExprPtr> symmetrizer =
          pop_symmetrizer(copy);

      SEQUANT_ASSERT(copy.produces_tensor());
      SEQUANT_ASSERT(copy.has_label());

      SEQUANT_ASSERT(symmetrizer.has_value());
      SEQUANT_ASSERT(symmetrizer->is<Tensor>());
      SEQUANT_ASSERT(symmetrizer->as<Tensor>().bra() == copy.ket());
      SEQUANT_ASSERT(symmetrizer->as<Tensor>().ket() == copy.bra());

      symmetrized_result = copy.result_as_tensor();

      copy.set_label(std::wstring(copy.label()) + L"u");

      return to_export_tree(copy);
    }();

    output.entries.push_back(
        {.tree = std::move(tree),
         .symm_contribution_target = std::move(symmetrized_result)});
  }

  ctx.set_data(id_prefix, id_start, std::move(output));

  return 1;
}

}  // namespace sequant::util::extint
