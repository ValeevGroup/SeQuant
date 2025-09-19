//
// Created by Eduard Valeyev on 2019-01-30.
//

#ifndef SEQUANT_TAG_HPP
#define SEQUANT_TAG_HPP

#include <SeQuant/core/meta.hpp>

#include <any>
#include <cassert>

namespace sequant {

class bad_any_comparable_cast : public std::bad_any_cast {
 public:
  bad_any_comparable_cast() = default;
  virtual ~bad_any_comparable_cast() {}
  virtual const char *what() const noexcept {
    return "Bad any_comparable_cast";
  }
};

namespace detail {

class any_comparable {
 public:
  // this is constexpr in the standard
  any_comparable() : impl_(nullptr) {}
  any_comparable(const any_comparable &other)
      : impl_(other.impl_ ? other.impl_->clone() : nullptr) {}
  any_comparable(any_comparable &&other) = default;
  template <
      typename ValueType,
      typename = std::enable_if_t<
          !std::is_base_of<any_comparable, std::decay_t<ValueType> >::value &&
          meta::is_less_than_comparable_v<std::decay_t<ValueType> > > >
  any_comparable(ValueType &&value)
      : impl_(new impl<typename std::decay<ValueType>::type>(
            std::forward<ValueType>(value))) {}
  ~any_comparable() = default;

  any_comparable &operator=(const any_comparable &rhs) {
    impl_ = rhs.impl_ ? decltype(impl_)(rhs.impl_->clone()) : nullptr;
    return *this;
  }
  any_comparable &operator=(any_comparable &&rhs) {
    impl_ = std::move(rhs.impl_);
    return *this;
  }
  template <
      typename ValueType,
      typename = std::enable_if_t<
          !std::is_base_of<any_comparable, std::decay_t<ValueType> >::value &&
          meta::is_less_than_comparable_v<std::decay_t<ValueType> > > >
  any_comparable &operator=(ValueType &&rhs) {
    impl_ = decltype(impl_)(new impl<typename std::decay<ValueType>::type>(
        std::forward<ValueType>(rhs)));
    return *this;
  }

  template <class ValueType, class... Args>
  std::enable_if_t<meta::is_less_than_comparable_v<std::decay_t<ValueType> >,
                   std::decay_t<ValueType> &>
  emplace(Args &&...args) {
    reset();
    impl_ = new impl<typename std::decay<ValueType>::type>(
        std::forward<Args>(args)...);
    return (impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }
  template <class ValueType, class U, class... Args>
  std::enable_if_t<meta::is_less_than_comparable_v<std::decay_t<ValueType> >,
                   std::decay_t<ValueType> &>
  emplace(std::initializer_list<U> il, Args &&...args) {
    reset();
    impl_ = new impl<typename std::decay<ValueType>::type>(
        il, std::forward<Args>(args)...);
    return (impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  void reset() { impl_.reset(); }

  void swap(any_comparable &other) { std::swap(impl_, other.impl_); }

  bool has_value() const { return static_cast<bool>(impl_); }

  const std::type_info &type() const {
    if (has_value())
      return impl_->type();
    else
      return typeid(void);
  }

  /// @return true if @c a is less than @c b
  /// @param a first object to compare
  /// @param b second object to compare
  /// @note will check that the types match if NDEBUG is not defined.
  friend bool operator<(const any_comparable &a, const any_comparable &b) {
    return *a.impl_ < *b.impl_;
  }

  /// @return true if @c a is equivalent to @c b
  /// @param a first object to compare
  /// @param b second object to compare
  /// @note will check that the types match if NDEBUG is not defined.
  friend bool operator==(const any_comparable &a, const any_comparable &b) {
    return *a.impl_ == *b.impl_;
  }

  /// @return true if @p a is not equivalent to @p b
  /// @param a first object to compare
  /// @param b second object to compare
  /// @note will check that the types match if NDEBUG is not defined.
  friend bool operator!=(const any_comparable &a, const any_comparable &b) {
    return !(a == b);
  }

 private:
  template <typename T,
            typename = std::enable_if_t<meta::is_less_than_comparable_v<T> > >
  struct impl;

  struct impl_base {
    virtual ~impl_base() {}
    virtual impl_base *clone() const = 0;

    virtual const std::type_info &type() const = 0;

    // static if NDEBUG is defined, dynamic otherwise
    template <typename T>
    impl<T> *cast() {
#ifndef NDEBUG
      return this->cast_static<T>();
#else
      return dynamic_cast<impl<T> *>(this);
#endif
    }
    // static if NDEBUG is defined, dynamic otherwise
    template <typename T>
    const impl<T> *cast() const {
#ifndef NDEBUG
      return this->cast_static<T>();
#else
      return dynamic_cast<const impl<T> *>(this);
#endif
    }
    // static always
    template <typename T>
    impl<T> *cast_static() {
      return static_cast<impl<T> *>(this);
    }
    // static always
    template <typename T>
    const impl<T> *cast_static() const {
      return static_cast<const impl<T> *>(this);
    }

    virtual bool operator<(const impl_base &other) const = 0;
    virtual bool operator==(const impl_base &other) const = 0;
  };
  template <typename T, typename>
  struct impl : public impl_base {
    template <typename U>
    explicit impl(U &&v) : value(std::forward<U>(v)) {}
    impl_base *clone() const override { return new impl{value}; }

    const std::type_info &type() const override { return typeid(T); }

    bool operator<(const impl_base &other) const override {
      assert(type() == other.type());
      return value < other.cast_static<T>()->value;
    }
    bool operator==(const impl_base &other) const override {
      assert(type() == other.type());
      return value == other.cast_static<T>()->value;
    }

    T value;
  };

  template <typename ValueType>
  friend typename std::decay<ValueType>::type *any_comparable_cast(
      any_comparable *operand);
  template <typename ValueType>
  friend const typename std::decay<ValueType>::type *any_comparable_cast(
      const any_comparable *operand);

  template <typename ValueType>
  typename std::decay<ValueType>::type *value_ptr() {
    return &(impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  template <typename ValueType>
  const typename std::decay<ValueType>::type *value_ptr() const {
    return &(impl_->cast_static<typename std::decay<ValueType>::type>()->value);
  }

  std::unique_ptr<impl_base> impl_;
};

template <typename ValueType>
typename std::decay<ValueType>::type *any_comparable_cast(
    any_comparable *operand) {
  if (operand->type() == typeid(typename std::decay<ValueType>::type))
    return operand->value_ptr<typename std::decay<ValueType>::type>();
  else
    return nullptr;
}

template <typename ValueType>
const typename std::decay<ValueType>::type *any_comparable_cast(
    const any_comparable *operand) {
  if (operand->type() == typeid(typename std::decay<ValueType>::type))
    return operand->value_ptr<typename std::decay<ValueType>::type>();
  else
    return nullptr;
}

template <typename ValueType>
ValueType any_comparable_cast(const any_comparable &operand) {
  const auto *cast_ptr =
      any_comparable_cast<typename std::decay<ValueType>::type>(&operand);
  if (cast_ptr != nullptr) return *cast_ptr;
  throw bad_any_comparable_cast();
}

template <typename ValueType>
ValueType any_comparable_cast(any_comparable &operand) {
  auto *cast_ptr =
      any_comparable_cast<typename std::decay<ValueType>::type>(&operand);
  if (cast_ptr != nullptr) return *cast_ptr;
  throw bad_any_comparable_cast();
}

}  // namespace detail

/// @brief type-erasing holder for Comparable objects
/// @note tag is not part of the state
/// (since tag are temporary and do not affect objects identity except for
/// temporarily labeling them), hence all methods are const
class Taggable {
 public:
  using any_comparable = ::sequant::detail::any_comparable;

  Taggable() noexcept : tag_{} { assert(!has_value()); }

  /// tags this object with tag @c t
  /// @param t tag to assign
  /// @return reference to this (for chaining)
  /// @pre `!this->has_value()`
  /// @post `this->value() == t`
  template <typename T>
  const Taggable &assign(const T &t) const {
    assert(!tag_.has_value());
    tag_ = t;
    assert(tag_.has_value());
    return *this;
  }

  /// @return this tag's value
  /// @throw bad_any_comparable_cast if the contained value is not convertible
  /// to T
  template <typename T>
  const T &value() const {
    assert(tag_.has_value());
    using detail::any_comparable_cast;
    return *any_comparable_cast<T>(&tag_);
  }

  /// @return true if tag has been assigned
  bool has_value() const { return tag_.has_value(); }

  /// resets this tag
  /// @return reference to this (for chaining)
  const Taggable &reset() const {
    tag_.reset();
    return *this;
  }

 private:
  mutable any_comparable tag_;

  /// @param a first object to compare
  /// @param b second object to compare
  /// @return true if @p a is less than @p b
  friend bool operator<(const Taggable &a, const Taggable &b) {
    return a.tag_ < b.tag_;
  }

  /// @param a first object to compare
  /// @param b second object to compare
  /// @return true if @p a is equal to @p b
  friend bool operator==(const Taggable &a, const Taggable &b) {
    return a.tag_ == b.tag_;
  }

  /// @param a first object to compare
  /// @param b second object to compare
  /// @return true if @p a is not equal to @p b
  friend bool operator!=(const Taggable &a, const Taggable &b) {
    return !(a == b);
  }
};

}  // namespace sequant

#endif  // SEQUANT_TAG_HPP
