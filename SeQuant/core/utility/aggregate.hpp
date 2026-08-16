//
// Created by Ajay Melekamburath on 8/4/26.
//

#ifndef SEQUANT_CORE_UTILITY_AGGREGATE_HPP
#define SEQUANT_CORE_UTILITY_AGGREGATE_HPP

#include <compare>
#include <type_traits>

namespace sequant::detail {

/// Tag type used by SEQUANT_DESIGNATED_INIT_ONLY; not meant to be named
/// directly. Its default constructor is private, so the only way to produce one
/// is `make()`, which is what the macro's default member initializer calls.
/// That is what makes both `Opts{a, b}` (`bool` does not convert to this type)
/// and `Opts{{}, a, b}` (no accessible default constructor) ill-formed.
class designated_init_only {
 public:
  constexpr designated_init_only(const designated_init_only&) noexcept =
      default;
  /// declaring the copy constructor deprecates the *implicit* copy assignment
  /// (-Wdeprecated-copy, an error under this project's -Werror), so declare it
  /// too -- a guarded aggregate must stay copy-assignable
  constexpr designated_init_only& operator=(
      const designated_init_only&) noexcept = default;
  static constexpr designated_init_only make() noexcept { return {}; }

  /// all tags compare equal, so that a guarded aggregate can still default its
  /// comparison operators
  friend constexpr bool operator==(designated_init_only,
                                   designated_init_only) noexcept {
    return true;
  }
  friend constexpr std::strong_ordering operator<=>(
      designated_init_only, designated_init_only) noexcept {
    return std::strong_ordering::equal;
  }

 private:
  constexpr designated_init_only() noexcept = default;
};

}  // namespace sequant::detail

// clang-format off
/// Declares a hidden first member of an `Options`-style aggregate that makes
/// positional (non-designated) aggregate initialization ill-formed, so that
/// every call site must name the fields it sets:
/// @code
///   struct LSTOptions {
///     SEQUANT_DESIGNATED_INIT_ONLY;
///     bool unitary = false;
///     bool use_connected_form = false;
///   };
/// @endcode
/// Positional aggregate init is otherwise legal C++, and it is not merely
/// unreadable: reordering the fields, or changing the meaning of one while
/// keeping its position, silently reinterprets any positional caller instead of
/// failing to compile. Designated initializers do not have that failure mode --
/// renaming a field is a hard error at every call site that sets it.
///
/// Rejected: `Opts{a, b}` and `Opts{{}, a, b}`.
/// Still accepted, and unaffected: `Opts{}`, `Opts{.field = a}`, copy and move,
/// `Opts` as a defaulted function parameter (`f(Opts opts = {})`), and use in
/// constant expressions.
///
/// @note Must be the *first* member, and requires C++20 designated
/// initializers.
/// @note The member is declared `[[no_unique_address]]`, so the aggregate's
/// size is unchanged, and it stays an aggregate and trivially copyable.
// clang-format on
#define SEQUANT_DESIGNATED_INIT_ONLY                            \
  [[no_unique_address]] ::sequant::detail::designated_init_only \
      sequant_designated_init_only_ =                           \
          ::sequant::detail::designated_init_only::make()

namespace sequant::detail {

/// self-test of SEQUANT_DESIGNATED_INIT_ONLY: the guard must not change what
/// the aggregate *is*, only how it may be written
struct designated_init_only_probe {
  SEQUANT_DESIGNATED_INIT_ONLY;
  bool a = false;
  bool b = false;
};
static_assert(std::is_aggregate_v<designated_init_only_probe>);
static_assert(std::is_trivially_copyable_v<designated_init_only_probe>);
// [[no_unique_address]] elision is permitted, not required, and MSVC ignores
// the standard spelling outright (it has [[msvc::no_unique_address]], since
// honoring this one would change its ABI). A toolchain that does not elide
// still gets a working guard -- the tag merely costs a byte -- so check the
// size only where elision is guaranteed rather than failing the build over it.
#if !defined(_MSC_VER) || defined(__clang__)
static_assert(sizeof(designated_init_only_probe) == 2 * sizeof(bool),
              "[[no_unique_address]] must elide the tag member");
#endif
/// odr-uses the copy assignment, which is what diagnoses -Wdeprecated-copy;
/// a trait check alone would not, since it never defines the operator
static_assert([] {
  designated_init_only_probe p{};
  p = designated_init_only_probe{.a = true};
  return p.a;
}());

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_AGGREGATE_HPP
