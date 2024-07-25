//
// Created by Eduard Valeyev on 7/24/24.
//

#ifndef SEQUANT_CORE_UTILITY_STRONG_HPP
#define SEQUANT_CORE_UTILITY_STRONG_HPP

#include <utility>

namespace sequant::detail {

template <typename T, typename Tag>
class strong_type_base {
 public:
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

  explicit constexpr strong_type_base(const T& value) : value_(value) {}

  explicit constexpr strong_type_base(T&& value) noexcept(
      std::is_nothrow_move_constructible<T>::value)
      : value_(std::move(value)) {}

  constexpr operator T&() & noexcept { return value_; }
  constexpr operator const T&() const& noexcept { return value_; }
  constexpr operator T&&() && noexcept { return std::move(value_); }
  constexpr operator const T&&() const&& noexcept { return std::move(value_); }

  constexpr T& value() & noexcept { return value_; }
  constexpr const T& value() const& noexcept { return value_; }
  constexpr T&& value() && noexcept { return std::move(value_); }
  constexpr const T&& value() const&& noexcept { return std::move(value_); }

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

  friend void swap(strong_type_base& a, strong_type_base& b) noexcept {
    using std::swap;
    swap(static_cast<T&>(a), static_cast<T&>(b));
  }

  template <typename U, typename TagU,
            typename = meta::are_less_than_comparable<T, U>>
  friend std::common_type_t<T, U> const max(
      const strong_type_base& a, const strong_type_base<U, TagU>& b) noexcept {
    return static_cast<const T&>(a) < static_cast<const U&>(b)
               ? static_cast<const U&>(b)
               : static_cast<const T&>(a);
  }

  template <typename U, typename TagU,
            typename = meta::are_less_than_comparable<T, U>>
  friend std::common_type_t<T, U> const min(
      const strong_type_base& a, const strong_type_base<U, TagU>& b) noexcept {
    return static_cast<const T&>(a) < static_cast<const U&>(b)
               ? static_cast<const T&>(a)
               : static_cast<const U&>(b);
  }

 private:
  T value_{};
};  // class strong_type_base

}  // namespace sequant::detail

#endif  // SEQUANT_CORE_UTILITY_STRONG_HPP
