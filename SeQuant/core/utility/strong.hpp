//
// Created by Eduard Valeyev on 7/24/24.
//

#ifndef SEQUANT_CORE_UTILITY_STRONG_HPP
#define SEQUANT_CORE_UTILITY_STRONG_HPP

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

  /// \name constructors that take the value of underlying type
  /// @{

  /// \param value the value to this will hold
  explicit constexpr strong_type_base(const T& value) : value_(value) {}

  /// \param value the value to this will hold
  explicit constexpr strong_type_base(T&& value) noexcept(
      std::is_nothrow_move_constructible<T>::value)
      : value_(std::move(value)) {}

  /// if \p T is a range, this will initialize the contained value using the
  /// provided initializer_list \param elements the elements to initialize the
  /// value with
  template <typename U, typename T_ = T,
            typename = std::enable_if_t<meta::is_range_v<T_>>>
  explicit constexpr strong_type_base(std::initializer_list<U> elements)
      : value_{elements.begin(), elements.end()} {}

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

 private:
  T value_{};
};  // class strong_type_base

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_STRONG_HPP
