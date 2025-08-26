//
// Created by Eduard Valeyev on 7/24/24.
//

#ifndef SEQUANT_CORE_UTILITY_STRONG_HPP
#define SEQUANT_CORE_UTILITY_STRONG_HPP

#include <SeQuant/core/meta.hpp>

#include <utility>

namespace sequant::detail {

/// \brief Base class for strong types

/// This class is modeled after the `strong_type` class of the excellent
/// [type_safe](https://type_safe.foonathan.net) library, but with
/// a stripped-down design covering our intended use cases. Namely, this class:
/// - provides implicit cast to the underlying type to make many uses less
/// verbose
///   once the API surface has been breached
/// - provides comparability and max/min functions if the underlying type
/// supports them
/// - provides range semantics if the underlying type supports it
///
/// Normally this is used as the base for the strong type, with the
/// CRPT-like use of the strong type as \p Tag:
/// \code
/// struct StrongInt : strong_type_base<int, StrongInt> {
///   using strong_type_base::strong_type_base;
/// };
/// \endcode
/// \tparam T the underlying type
/// \tparam Tag tag type
template <typename T, typename Tag>
class strong_type_base {
  template <typename U, typename TagU>
  friend class strong_type_base;

 public:
  /// \name standard constructors and assignment operators
  /// @{
  constexpr strong_type_base() noexcept(
      std::is_nothrow_default_constructible_v<T>) = default;
  constexpr strong_type_base(const strong_type_base&) noexcept(
      std::is_nothrow_copy_constructible_v<T>) = default;
  constexpr strong_type_base(strong_type_base&&) noexcept(
      std::is_nothrow_move_constructible_v<T>) = default;
  constexpr strong_type_base& operator=(const strong_type_base&) noexcept(
      std::is_nothrow_copy_assignable_v<T>) = default;
  constexpr strong_type_base& operator=(strong_type_base&&) noexcept(
      std::is_nothrow_move_assignable_v<T>) = default;
  /// @}

  /// \name extension of standard constructors and assignment operators to
  /// handle different tags as well as payload conversions
  /// @{
  template <
      typename U, typename TagU,
      typename = std::enable_if_t<
          (!std::is_same_v<T, U> ||
           !std::is_same_v<Tag, TagU>)&&std::is_constructible_v<T, const U&>>>
  explicit constexpr strong_type_base(
      const strong_type_base<U, TagU>&
          other) noexcept(std::is_nothrow_constructible_v<T, const U&>)
      : value_(other.value_) {}
  template <typename U, typename TagU,
            typename = std::enable_if_t<
                (!std::is_same_v<T, U> ||
                 !std::is_same_v<Tag, TagU>)&&std::is_constructible_v<T, U&&>>>
  explicit constexpr strong_type_base(
      strong_type_base<U, TagU>&&
          other) noexcept(std::is_nothrow_constructible_v<T, U&&>)
      : value_(std::move(other.value_)) {}
  template <
      typename U, typename TagU,
      typename = std::enable_if_t<
          (!std::is_same_v<T, U> ||
           !std::is_same_v<Tag, TagU>)&&std::is_assignable_v<T&, const U&>>>
  constexpr strong_type_base&
  operator=(const strong_type_base<U, TagU>& other) noexcept(
      std::is_nothrow_assignable_v<T, const U&>) {
    value_ = other.value_;
    return *this;
  }
  template <typename U, typename TagU,
            typename = std::enable_if_t<
                (!std::is_same_v<T, U> ||
                 !std::is_same_v<Tag, TagU>)&&std::is_assignable_v<T&, U&&>>>
  constexpr strong_type_base&
  operator=(strong_type_base<U, TagU>&& other) noexcept(
      std::is_nothrow_assignable_v<T, U&&>) {
    value_ = std::move(other.value_);
    return *this;
  }
  /// @}

  /// \name constructors that take the value of underlying type
  /// @{

  /// \param value the value to this will hold
  /// \note this is disabled if \p T is the castable_to_any placeholder to allow
  /// GCC resolve string_type({}) construction ambiguity
  template <typename T_ = T,
            typename = std::enable_if_t<!std::is_same_v<
                meta::remove_cvref_t<T_>, meta::castable_to_any>>>
  explicit constexpr strong_type_base(const T& value) : value_(value) {}

  /// \param value the value to this will hold
  /// \note this is disabled if \p T is the castable_to_any placeholder to allow
  /// GCC resolve string_type({}) construction ambiguity
  template <typename T_ = T,
            typename = std::enable_if_t<!std::is_same_v<
                meta::remove_cvref_t<T_>, meta::castable_to_any>>>
  explicit constexpr strong_type_base(T&& value) noexcept(
      std::is_nothrow_move_constructible<T>::value)
      : value_(std::move(value)) {}

  /// if \p T is a range, this will initialize the contained value using the
  /// provided initializer_list \param elements the elements to initialize the
  /// value with
  template <typename U, typename T_ = T,
            typename = std::enable_if_t<
                meta::is_range_v<T_> &&
                meta::is_statically_castable_v<U, meta::range_value_t<T_>>>>
  explicit constexpr strong_type_base(std::initializer_list<U> elements) {
    if constexpr (meta::has_memfn_push_back_v<T, U>)
      std::copy(elements.begin(), elements.end(), std::back_inserter(value_));
    else {
      assert(elements.size() == value_.size());
      std::copy(elements.begin(), elements.end(), value_.begin());
    }
  }

  /// if \p T is a std::array, this will initialize the contained value using
  /// the provided parameter pack \param elements the elements to initialize the
  /// value with
  template <typename... Elements,
            typename = std::enable_if_t<meta::is_std_array_v<T> &&
                                        meta::std_array_size_v<T> ==
                                            sizeof...(Elements)>>
  explicit constexpr strong_type_base(Elements&&... elements)
      : value_{std::forward<Elements>(elements)...} {}
  /// @}

  /// \name implicit conversion operators provide direct access the value
  /// @{
  constexpr operator T&() & noexcept { return value_; }
  constexpr operator const T&() const& noexcept { return value_; }
  constexpr operator T&&() && noexcept { return std::move(value_); }
  constexpr operator const T&&() const&& noexcept { return std::move(value_); }
  /// @}

  /// \name explicit access to the value
  /// @{
  constexpr T& value() & noexcept { return value_; }
  constexpr const T& value() const& noexcept { return value_; }
  constexpr T&& value() && noexcept { return std::move(value_); }
  constexpr const T&& value() const&& noexcept { return std::move(value_); }
  /// @}

  /// \name pass-through the range interface of \p T , if available
  /// @{
  template <typename T_ = T,
            typename = std::enable_if_t<meta::has_memfn_size_v<const T_>>>
  decltype(auto) size() const {
    return value_.size();
  }

  template <typename T_ = T,
            typename = std::enable_if_t<meta::has_memfn_empty_v<const T_>>>
  decltype(auto) empty() const {
    return value_.empty();
  }

  template <typename T_ = T,
            typename = std::enable_if_t<meta::is_range_v<const T_>>>
  decltype(auto) begin() const {
    using ranges::begin;
    return begin(value_);
  }
  template <typename T_ = T, typename = std::enable_if_t<meta::is_range_v<T_>>>
  decltype(auto) begin() {
    using ranges::begin;
    return begin(value_);
  }
  template <typename T_ = T,
            typename = std::enable_if_t<meta::is_range_v<const T_>>>
  decltype(auto) end() const {
    using ranges::end;
    return end(value_);
  }
  template <typename T_ = T, typename = std::enable_if_t<meta::is_range_v<T_>>>
  decltype(auto) end() {
    using ranges::end;
    return end(value_);
  }

  template <typename T_ = T, typename = std::enable_if_t<
                                 meta::has_operator_subscript_v<const T_>>>
  decltype(auto) operator[](std::size_t i) const {
    return value_[i];
  }

  template <typename T_ = T,
            typename = std::enable_if_t<meta::has_operator_subscript_v<T_>>>
  decltype(auto) operator[](std::size_t i) {
    return value_[i];
  }

  template <typename T_ = T,
            typename = std::enable_if_t<meta::has_memfn_at_v<const T_>>>
  decltype(auto) at(std::size_t i) const {
    return value_.at(i);
  }

  template <typename T_ = T,
            typename = std::enable_if_t<meta::has_memfn_at_v<T_>>>
  decltype(auto) at(std::size_t i) {
    return value_.at(i);
  }

  template <typename T_ = T, typename U,
            typename = std::enable_if_t<meta::has_memfn_push_back_v<T_, U&&>>>
  decltype(auto) push_back(U&& u) const {
    return value_.push_back(std::forward<U>(u));
  }

  /// @}

  friend void swap(strong_type_base& a, strong_type_base& b) noexcept {
    using std::swap;
    swap(static_cast<T&>(a), static_cast<T&>(b));
  }

  /// \name pass-through the comparison interface of \p T , if available
  /// @{
  template <typename U, typename TagU,
            typename = std::enable_if_t<meta::are_equality_comparable_v<T, U>>>
  friend bool operator==(const strong_type_base& a,
                         const strong_type_base<U, TagU>& b) noexcept {
    return static_cast<const T&>(a) == static_cast<const U&>(b);
  }

  template <typename U, typename TagU,
            typename = std::enable_if_t<meta::are_less_than_comparable_v<T, U>>>
  friend bool operator<(const strong_type_base& a,
                        const strong_type_base<U, TagU>& b) noexcept {
    return static_cast<const T&>(a) < static_cast<const U&>(b);
  }

  template <
      typename U, typename TagU,
      typename = std::enable_if_t<meta::are_greater_than_comparable_v<T, U> ||
                                  (meta::are_equality_comparable_v<T, U> &&
                                   meta::are_less_than_comparable_v<T, U>)>>
  friend bool operator>(const strong_type_base& a,
                        const strong_type_base<U, TagU>& b) noexcept {
    if constexpr (meta::are_greater_than_comparable_v<T, U>)
      return static_cast<const T&>(a) > static_cast<const U&>(b);
    else
      return !(static_cast<const T&>(a) < static_cast<const U&>(b) ||
               static_cast<const T&>(a) == static_cast<const U&>(b));
  }

  template <typename U, typename TagU,
            typename = meta::are_less_than_comparable<T, U>>
  friend std::common_type_t<T, U> max(
      const strong_type_base& a, const strong_type_base<U, TagU>& b) noexcept {
    return a < b ? b.value() : a.value();
  }

  template <typename U, typename TagU,
            typename = meta::are_less_than_comparable<T, U>>
  friend std::common_type_t<T, U> min(
      const strong_type_base& a, const strong_type_base<U, TagU>& b) noexcept {
    return a < b ? a.value() : b.value();
  }
  /// @}

  /// \name pass-through the arithmetic interface of \p T , if available
  /// @{
  template <typename T_ = T,
            typename = std::enable_if_t<std::is_integral_v<T_> &&
                                        std::is_signed_v<T_>>>
  strong_type_base operator-() const noexcept {
    return strong_type_base(-value());
  }
  /// @}

 private:
  T value_{};
};  // class strong_type_base

}  // namespace sequant::detail

#ifndef DEFINE_STRONG_TYPE_FOR_RANGE
// N.B. default deduction guide with older gcc does not defer to the
// initializer_list guide, so end up with non-range instantiations clang and
// recent gcc (14) handle this fine ... this means we have to explicitly define
// the list of types that the default deduction guide will use as is (i.e.
// interpret as T)
#define DEFINE_STRONG_TYPE_FOR_RANGE(ID)                                       \
  template <typename T>                                                        \
  struct ID : detail::strong_type_base<T, ID<T>> {                             \
    using value_type = T;                                                      \
    using base_type = detail::strong_type_base<T, ID<T>>;                      \
    using base_type::base_type;                                                \
    using base_type::operator=;                                                \
  };                                                                           \
                                                                               \
  template <typename T>                                                        \
  ID(T&& t) -> ID<std::conditional_t<                                          \
      ((meta::is_range_v<meta::remove_cvref_t<T>> &&                           \
        !meta::is_char_range_v<meta::remove_cvref_t<T>>) ||                    \
       std::is_arithmetic_v<meta::remove_cvref_t<T>> ||                        \
       std::is_same_v<meta::remove_cvref_t<T>, IndexSpace>),                   \
      meta::remove_cvref_t<T>,                                                 \
      container::svector<                                                      \
          meta::literal_to_string_t<meta::remove_cvref_t<T>>>>>;               \
  template <typename T = meta::castable_to_any>                                \
  ID(std::initializer_list<T> t) -> ID<                                        \
      container::svector<meta::literal_to_string_t<meta::remove_cvref_t<T>>>>; \
  template <                                                                   \
      typename T, typename... U,                                               \
      typename = std::enable_if_t<                                             \
          sizeof...(U) != 0 &&                                                 \
          (std::is_same_v<meta::remove_cvref_t<T>, meta::remove_cvref_t<U>> && \
           ...)>>                                                              \
  ID(T&& t, U&&... rest)                                                       \
      -> ID<std::array<meta::literal_to_string_t<meta::remove_cvref_t<T>>,     \
                       1 + sizeof...(U)>>;                                     \
  template <typename T = meta::castable_to_any>                                \
  ID() -> ID<std::array<T, 0>>;
#endif  // DEFINE_STRONG_TYPE_FOR_RANGE

#ifndef DEFINE_STRONG_TYPE_FOR_INTEGER
#define DEFINE_STRONG_TYPE_FOR_INTEGER(ID, IntType)          \
  struct ID : detail::strong_type_base<IntType, ID> {        \
    using value_type = IntType;                              \
    using base_type = detail::strong_type_base<IntType, ID>; \
    using base_type::base_type;                              \
    using base_type::operator=;                              \
    using base_type::operator-;                              \
    ID(const base_type& b) noexcept : base_type(b) {}        \
    ID(base_type&& b) noexcept : base_type(std::move(b)) {}  \
  };
#endif  // DEFINE_STRONG_TYPE_FOR_INTEGER
#ifndef DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE
#define DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE(ID)               \
  DEFINE_STRONG_TYPE_FOR_INTEGER(SEQUANT_CONCAT(n, ID), std::size_t) \
  DEFINE_STRONG_TYPE_FOR_RANGE(ID)
#endif  // DEFINE_STRONG_TYPE_FOR_RANGE_AND_RANGESIZE

#endif  // SEQUANT_CORE_UTILITY_STRONG_HPP
