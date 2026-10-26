//
// Created by Eduard Valeyev on 10/21/26.
//

#ifndef SEQUANT_CORE_UTILITY_DEPRECATED_MEMBER_HPP
#define SEQUANT_CORE_UTILITY_DEPRECATED_MEMBER_HPP

#include <concepts>
#include <type_traits>
#include <utility>

namespace sequant {

namespace detail {
/// the tag the member type SEQUANT_DEPRECATED_MEMBER declares passes to
/// DeprecatedMember's value constructor, so that only it can call that one
struct deprecated_member_init_t {
  explicit deprecated_member_init_t() = default;
};
}  // namespace detail

/// @brief the base of the type of a deprecated data member of an aggregate,
/// declared with SEQUANT_DEPRECATED_MEMBER: holds a T, and can be given one
/// only by the derived type, whose converting constructor carries the
/// deprecation
///
/// Marking an aggregate's data member `[[deprecated]]` directly is not
/// portable: gcc 15 counts the default member initializer of a deprecated
/// member as a use of it, in every designated initialization of the aggregate
/// that omits the member (that is, in every caller's expression, which no
/// pragma in the header that declares the member reaches). Carrying the
/// deprecation on the converting constructor of the member's type instead
/// diagnoses the member where it is given a value, `{.old_name = ...}`, and
/// nowhere else, on every compiler.
/// @tparam T the type of the value the member holds
template <typename T>
struct DeprecatedMember {
  static_assert(std::is_object_v<T> && !std::is_const_v<T> &&
                    std::is_default_constructible_v<T>,
                "DeprecatedMember<T>: T is a non-const, default-constructible "
                "object type (the member holds a T, value-initialized when the "
                "alias is not spelled)");

  /// the member's type must stay default-constructible for the aggregate's
  /// default member initializer, which is not a use of the deprecated spelling
  DeprecatedMember() = default;

  /// what the member was given, or a value-initialized T if it was not
  T value{};

 protected:
  template <typename U>
    requires std::constructible_from<T, U&&>
  constexpr DeprecatedMember(detail::deprecated_member_init_t, U&& u)
      : value(std::forward<U>(u)) {}
};

// clang-format off
/// Declares a deprecated data member of an aggregate: `Name`, of a type
/// (named `Name##_deprecated_t`, derived from DeprecatedMember<Type>) that
/// holds a `Type` in its `value` and whose one converting constructor is
/// `[[deprecated(Msg)]]`, so that `{.Name = x}` is diagnosed with @p Msg where
/// it is written and the aggregate is otherwise untouched by the deprecation:
/// @code
///   struct Options {
///     SEQUANT_DESIGNATED_INIT_ONLY;
///     std::shared_ptr<const Registry> basis_registry = nullptr;
///     /// @deprecated use basis_registry
///     SEQUANT_DEPRECATED_MEMBER(space_registry, "use basis_registry",
///                               std::shared_ptr<const Registry>);
///   };
///   // reading the alias: options.space_registry.value
/// @endcode
/// @param Name the member's name
/// @param Msg the deprecation message, a string literal
/// @param ... the member's type: a non-const, default-constructible object
/// type (variadic, so that it may contain commas)
/// @note the member's default member initializer is `{}`; a designated
/// initializer that gives it a value, and nothing else, is what warns, so
/// `{.Name = nullptr}` warns even where `nullptr` is the default: the
/// spelling is what is deprecated
/// @note the member is the declared type, not a `Type`: code that read the
/// member as a `Type` reads `Name.value` instead; the converting constructor
/// never takes the declared type itself, so copying or moving the member is
/// not a use of the deprecated spelling
// clang-format on
#ifdef SEQUANT_DOXYGEN
// the documentation's view of the member: a plain member of the value's type,
// so that the @deprecated comment above an invocation attaches to it (the
// Doxyfile predefines SEQUANT_DOXYGEN and expands this macro by name)
#define SEQUANT_DEPRECATED_MEMBER(Name, Msg, ...) __VA_ARGS__ Name
#else
#define SEQUANT_DEPRECATED_MEMBER(Name, Msg, ...)                          \
  struct Name##_deprecated_t final                                         \
      : ::sequant::DeprecatedMember<__VA_ARGS__> {                         \
    Name##_deprecated_t() = default;                                       \
    template <typename U>                                                  \
      requires(std::constructible_from<__VA_ARGS__, U &&> &&               \
               !std::same_as<std::remove_cvref_t<U>, Name##_deprecated_t>) \
    [[deprecated(Msg)]] Name##_deprecated_t(U&& u)                         \
        : ::sequant::DeprecatedMember<__VA_ARGS__>(                        \
              ::sequant::detail::deprecated_member_init_t{},               \
              std::forward<U>(u)) {}                                       \
  };                                                                       \
  Name##_deprecated_t Name = {}
#endif  // SEQUANT_DOXYGEN

namespace detail {
/// self-test of SEQUANT_DEPRECATED_MEMBER: the aggregate stays an aggregate,
/// its default construction and designated initialization of the other
/// members do not touch the deprecated constructor, and the alias holds what
/// it is given
struct deprecated_member_probe {
  int current = 0;
  SEQUANT_DEPRECATED_MEMBER(old, "use current", int);
};
static_assert(std::is_aggregate_v<deprecated_member_probe>);
static_assert(std::is_default_constructible_v<deprecated_member_probe>);
static_assert(deprecated_member_probe{.current = 1}.old.value == 0);
static_assert(
    !std::is_constructible_v<deprecated_member_probe::old_deprecated_t,
                             detail::deprecated_member_init_t, int>,
    "the value constructor is for the declared type only");

/// a value type constructible from anything: were the member's converting
/// constructor to take the member's own type, copying the member from a
/// non-const lvalue would route through it (and warn) instead of the copy
/// constructor, and the copy would hold the member instead of its value
struct deprecated_member_greedy {
  int from_alias = 0;
  constexpr deprecated_member_greedy() = default;
  template <typename U>
  constexpr deprecated_member_greedy(U&&) : from_alias(1) {}
};
struct deprecated_member_greedy_probe {
  SEQUANT_DEPRECATED_MEMBER(old, "use current", deprecated_member_greedy);
};
static_assert(
    [] {
      deprecated_member_greedy_probe p{};
      auto copy = p.old;
      return copy.value.from_alias == 0;
    }(),
    "copying the member is a copy, not a use of the deprecated constructor");
}  // namespace detail

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_DEPRECATED_MEMBER_HPP
