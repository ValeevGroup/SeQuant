#include "contexts.hpp"

#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

namespace sequant::deriv::test {

ContextScope real_spaces(std::vector<std::wstring> keys) {
  auto registry = mbpt::make_sr_spaces();
  if (keys.empty())
    for (const auto& s : *registry)
      if (!s.base_key().empty()) keys.push_back(s.base_key());
  for (const auto& k : keys) registry->field(k, Field::Real);
  auto context = get_default_context();
  context.set(registry);
  return set_scoped_default_context(context);
}

ContextScope with_batching_spaces() {
  auto registry = mbpt::make_sr_spaces();
  mbpt::add_batching_spaces(registry);
  auto context = get_default_context();
  context.set(registry);
  return set_scoped_default_context(context);
}

}  // namespace sequant::deriv::test
