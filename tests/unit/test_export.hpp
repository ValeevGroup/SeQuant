#ifndef SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP
#define SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP

#include <SeQuant/core/index_space_registry.hpp>

#include "SeQuant/core/context.hpp"

namespace {
[[nodiscard]] [[maybe_unused]] inline auto to_export_context() {
  using namespace sequant;
  auto reg = std::make_shared<IndexSpaceRegistry>();
  reg->add(L"i", 0b001, is_particle, 10);
  reg->add(L"a", 0b010, is_vacuum_occupied, is_reference_occupied, is_hole,
           100);
  reg->add(L"u", 0b100, is_reference_occupied, is_hole, is_particle, 5);

  return set_scoped_default_context(
      {.index_space_registry_shared_ptr = std::move(reg)});
}

}  // anonymous namespace

#endif  // SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP
