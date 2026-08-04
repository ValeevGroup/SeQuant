//
// Created by Ajay Melekamburath on 8/4/26.
//

#ifndef SEQUANT_CORE_UTILITY_AGGREGATE_HPP
#define SEQUANT_CORE_UTILITY_AGGREGATE_HPP

namespace sequant::detail {

/// Tag type for a hidden first member of an aggregate `Options`-style struct,
/// to force designated initialization at every call site, e.g.
/// `Options{.foo = true}` rather than `Options{true, false}`. Positional
/// aggregate init is otherwise still legal C++, so renaming, reordering, or
/// adding a field silently reinterprets any positional-init caller instead of
/// failing to compile.
/// @note declare it `[[no_unique_address]] designated_only designated_only_ =
/// {};` as the struct's first member; designated initializers may skip a
/// defaulted member, so this does not need to be named at any call site.
struct designated_only {};

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_AGGREGATE_HPP
