#ifndef SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP
#define SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP

#include <SeQuant/core/index_basis_registry.hpp>

#include "SeQuant/core/context.hpp"

namespace {
[[nodiscard]] [[maybe_unused]] inline auto to_export_context() {
  using namespace sequant;
  auto reg = std::make_shared<IndexBasisRegistry>();
  reg->add(L"i", 0b001, {.dimension = 10, .particle = true});
  reg->add(L"a", 0b010,
           {.dimension = 100,
            .vacuum_occupied = true,
            .reference_occupied = true,
            .hole = true});
  reg->add(L"u", 0b100,
           {.dimension = 5,
            .reference_occupied = true,
            .hole = true,
            .particle = true});
  reg->add(L"x", 0b1000, {.dimension = 50, .particle = true});

  return set_scoped_default_context(
      {.index_basis_registry_shared_ptr = std::move(reg)});
}

}  // anonymous namespace

#endif  // SEQUANT_UNIT_TESTS_TEST_EXPORT_HPP
