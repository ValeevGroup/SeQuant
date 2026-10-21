#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/export/export.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <cwctype>
#include <string>

namespace sequant::detail {

auto get_with_default(const auto &map, const auto &key, auto def_value) {
  auto it = map.find(key);

  return it == map.end() ? def_value : it->second;
}

bool rename(Tensor &tensor, PreprocessResult &result,
            const std::function<bool(const Tensor &)> &acceptable) {
  SEQUANT_ASSERT(acceptable || result.tensorReferences.at(tensor) > 0);

  std::wstring basename = [&tensor]() {
    std::size_t size = tensor.label().size();

    while (size > 0 && std::iswdigit(tensor.label().at(size - 1))) {
      size--;
    }
    SEQUANT_ASSERT(size > 0);

    return std::wstring(tensor.label().substr(0, size));
  }();

  tensor.set_label(basename);

  std::size_t counter = 2;
  while (get_with_default(result.tensorReferences, tensor, 0) > 0 ||
         (acceptable && !acceptable(tensor))) {
    tensor.set_label(basename + std::to_wstring(counter));
    counter++;
  }

  return result.tensorReferences.find(tensor) != result.tensorReferences.end();
}

bool rename(Variable &variable, PreprocessResult &result,
            const std::function<bool(const Variable &)> &acceptable) {
  SEQUANT_ASSERT(acceptable || result.variableReferences.at(variable) > 0);

  std::wstring basename = [&variable]() {
    std::size_t size = variable.label().size();

    while (size > 0 && std::iswdigit(variable.label().at(size - 1))) {
      size--;
    }
    SEQUANT_ASSERT(size > 0);

    return std::wstring(variable.label().substr(0, size));
  }();

  variable.set_label(basename);

  std::size_t counter = 2;
  while (get_with_default(result.variableReferences, variable, 0) > 0 ||
         (acceptable && !acceptable(variable))) {
    variable.set_label(basename + std::to_wstring(counter));
    counter++;
  }

  return result.variableReferences.find(variable) !=
         result.variableReferences.end();
}

RegisteredIrrepContext::RegisteredIrrepContext(const ExportContext &ctx)
    : m_ctx(ctx) {}

AnyIrrepSet RegisteredIrrepContext::irreps_of(const Tensor &tensor) const {
  // The registry is keyed by the rewritten tensors
  Tensor rewritten(tensor);
  m_ctx.rewrite(rewritten);
  if (const AnyIrrepSet *known = m_ctx.irreps(rewritten)) {
    return *known;
  }

  return m_ctx.irrep_context().irreps_of(tensor);
}

AnyIrrepSet RegisteredIrrepContext::totally_symmetric() const {
  return m_ctx.irrep_context().totally_symmetric();
}

AnyIrrepSet RegisteredIrrepContext::empty_set() const {
  return m_ctx.irrep_context().empty_set();
}

AnyIrrepSet RegisteredIrrepContext::product(const AnyIrrepSet &lhs,
                                            const AnyIrrepSet &rhs) const {
  return m_ctx.irrep_context().product(lhs, rhs);
}

AnyIrrepSet RegisteredIrrepContext::conjugate(const AnyIrrepSet &set) const {
  return m_ctx.irrep_context().conjugate(set);
}

}  // namespace sequant::detail
