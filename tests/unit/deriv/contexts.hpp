#ifndef SEQUANT_TESTS_UNIT_DERIV_CONTEXTS_HPP
#define SEQUANT_TESTS_UNIT_DERIV_CONTEXTS_HPP

#include <SeQuant/core/context.hpp>

#include <string>
#include <vector>

namespace sequant::deriv::test {

/// restores the default context it replaced when destroyed
using ContextScope = sequant::detail::ImplicitContextResetter<
    container::map<Statistics, Context>>;

/// the default context, while the result lives, with the index spaces of
/// @p keys real, or every space if @p keys is empty
[[nodiscard]] ContextScope real_spaces(std::vector<std::wstring> keys = {});

/// the default context with the batching space z, while the result lives
[[nodiscard]] ContextScope with_batching_spaces();

}  // namespace sequant::deriv::test

#endif  // SEQUANT_TESTS_UNIT_DERIV_CONTEXTS_HPP
