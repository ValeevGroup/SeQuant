#ifndef SEQUANT_EXPRESSIONS_EXPR_ITERATOR_HPP
#define SEQUANT_EXPRESSIONS_EXPR_ITERATOR_HPP

#include <SeQuant/core/utility/macros.hpp>

#include <compare>
#include <concepts>
#include <iterator>
#include <type_traits>

namespace sequant {

class ExprPtr;

namespace detail {

template <bool is_const>
class ExprIteratorImpl {
 public:
  using value_type = ExprPtr;
  using reference = std::add_lvalue_reference_t<
      std::conditional_t<is_const, std::add_const_t<value_type>, value_type>>;
  using const_reference =
      std::add_lvalue_reference_t<std::add_const_t<value_type>>;
  using pointer = std::add_pointer_t<
      std::conditional_t<is_const, std::add_const_t<value_type>, value_type>>;
  using difference_type = std::ptrdiff_t;

  explicit ExprIteratorImpl(pointer ptr = nullptr) : ptr_(ptr) {}

  /// converting constructor: a mutable iterator converts to a const iterator
  /// (but not the other way around)
  /// @note this is a constructor *template* on purpose: that way it is never
  ///       considered a copy constructor and thus never suppresses the
  ///       implicitly-declared one
  template <bool other_is_const>
    requires(is_const && !other_is_const)
  ExprIteratorImpl(const ExprIteratorImpl<other_is_const> &other)
      : ptr_(other.ptr_) {}

  ExprIteratorImpl &operator+=(difference_type val) {
    ptr_ += val;
    return *this;
  }

  friend ExprIteratorImpl operator+(const ExprIteratorImpl &it,
                                    difference_type val) {
    return ExprIteratorImpl(it.ptr_ + val);
  }

  friend ExprIteratorImpl operator+(difference_type val,
                                    const ExprIteratorImpl &it) {
    return ExprIteratorImpl(it.ptr_ + val);
  }

  ExprIteratorImpl &operator++() {
    ++ptr_;
    return *this;
  }

  ExprIteratorImpl operator++(int) {
    ExprIteratorImpl copy = *this;

    ++ptr_;

    return copy;
  }

  ExprIteratorImpl &operator-=(difference_type val) {
    ptr_ -= val;
    return *this;
  }

  friend ExprIteratorImpl operator-(const ExprIteratorImpl &it,
                                    difference_type val) {
    return ExprIteratorImpl(it.ptr_ - val);
  }

  ExprIteratorImpl &operator--() {
    --ptr_;
    return *this;
  }

  ExprIteratorImpl operator--(int) {
    ExprIteratorImpl copy = *this;

    --ptr_;

    return copy;
  }

  reference operator*() const {
    SEQUANT_ASSERT(ptr_);
    return *ptr_;
  }

  pointer operator->() const {
    SEQUANT_ASSERT(ptr_);
    return ptr_;
  }

  template <bool other_is_const>
  difference_type operator-(
      const ExprIteratorImpl<other_is_const> &other) const {
    return ptr_ - other.ptr_;
  }

  template <bool other_is_const>
  bool operator==(const ExprIteratorImpl<other_is_const> &other) const {
    return ptr_ == other.ptr_;
  }

  reference operator[](difference_type offset) const {
    SEQUANT_ASSERT(ptr_);
    return *(ptr_ + offset);
  }

  template <bool other_is_const>
  std::strong_ordering operator<=>(
      const ExprIteratorImpl<other_is_const> &other) const {
    return ptr_ <=> other.ptr_;
  }

 private:
  // needed so that the const and the non-const specializations can access each
  // other's ptr_ (see the converting constructor and the heterogeneous
  // comparison/difference operators above)
  template <bool>
  friend class ExprIteratorImpl;

  pointer ptr_ = nullptr;
};

}  // namespace detail

using ExprIterator = detail::ExprIteratorImpl<false>;
using ConstExprIterator = detail::ExprIteratorImpl<true>;

static_assert(std::bidirectional_iterator<ExprIterator>);
static_assert(std::random_access_iterator<ExprIterator>);
static_assert(std::bidirectional_iterator<ConstExprIterator>);
static_assert(std::random_access_iterator<ConstExprIterator>);

// mutable iterators must interoperate with (and convert to) const iterators,
// but not vice versa
static_assert(std::convertible_to<ExprIterator, ConstExprIterator>);
static_assert(!std::convertible_to<ConstExprIterator, ExprIterator>);
static_assert(std::equality_comparable_with<ExprIterator, ConstExprIterator>);
static_assert(std::totally_ordered_with<ExprIterator, ConstExprIterator>);

namespace detail {
// `it - n` is a valid random-access-iterator expression, `n - it` is not
// (unlike `n + it`, which is)
template <typename It>
concept subtractable_from_difference =
    requires(It it, std::ptrdiff_t n) { n - it; };
}  // namespace detail
static_assert(!detail::subtractable_from_difference<ExprIterator>);
static_assert(!detail::subtractable_from_difference<ConstExprIterator>);

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_EXPR_ITERATOR_HPP
