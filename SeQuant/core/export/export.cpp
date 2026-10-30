#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string>

namespace sequant::detail {

bool rename(Tensor &tensor, PreprocessResult &result,
            const std::function<bool(const Tensor &)> &acceptable) {
  SEQUANT_ASSERT(acceptable || result.tensorReferences.at(tensor) > 0);

  auto allocation = result.names.lowest_available(
      opt::detail::label_basename(tensor.label()), tensor, acceptable);
  tensor.set_label(std::move(allocation.label));

  return allocation.recycled;
}

bool rename(Variable &variable, PreprocessResult &result,
            const std::function<bool(const Variable &)> &acceptable) {
  SEQUANT_ASSERT(acceptable || result.variableReferences.at(variable) > 0);

  auto allocation = result.names.lowest_available(
      opt::detail::label_basename(variable.label()), variable, acceptable);
  variable.set_label(std::move(allocation.label));

  return allocation.recycled;
}

RegisteredIrrepContext::RegisteredIrrepContext(const ExportContext &ctx)
    : m_ctx(ctx) {}

AnyIrrepSet RegisteredIrrepContext::irreps_of(const Tensor &tensor) const {
  // The registry is keyed by the rewritten tensors
  Tensor rewritten(tensor);
  m_ctx.rewrite(rewritten);
  if (const AnyIrrepSet *computed = m_ctx.computed_irreps_of(rewritten)) {
    return *computed;
  }

  return m_ctx.irrep_context().irreps_of(tensor);
}

AnyIrrepSet RegisteredIrrepContext::irreps_of(const Variable &variable) const {
  return m_ctx.irrep_context().irreps_of(variable);
}

AnyIrrepSet RegisteredIrrepContext::totally_symmetric() const {
  return m_ctx.irrep_context().totally_symmetric();
}

AnyIrrepSet RegisteredIrrepContext::empty_set() const {
  return m_ctx.irrep_context().empty_set();
}

AnyIrrepSet RegisteredIrrepContext::any_set() const {
  return m_ctx.irrep_context().any_set();
}

AnyIrrepSet RegisteredIrrepContext::product(const AnyIrrepSet &lhs,
                                            const AnyIrrepSet &rhs) const {
  return m_ctx.irrep_context().product(lhs, rhs);
}

AnyIrrepSet RegisteredIrrepContext::conjugate(const AnyIrrepSet &set) const {
  return m_ctx.irrep_context().conjugate(set);
}

}  // namespace sequant::detail
