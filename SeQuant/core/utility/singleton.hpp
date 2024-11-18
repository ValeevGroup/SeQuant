//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_CORE_UTILITY_SINGLETON_HPP
#define SEQUANT_CORE_UTILITY_SINGLETON_HPP

#include <cassert>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>

namespace sequant {

// clang-format off
/// @brief Singleton base class
/// Singleton is a CRTP base class that provides a thread-safe way to create an immutable singleton object.
///
/// To create a singleton class @c A do:
/// \code
/// class A : Singleton<A> {
///   private:
///     friend class Singleton<A>;
///     A(...);  // define (private) constructors
/// };
/// \endcode
/// Here's how to use it:
/// \code
/// // optional: create the instance of A
/// A::set_instance(Args...);  // this call-once method creates the instance of A; if A is default-constructible, can skip this as the instance will be created by first call to instance()
/// // access the instance of A
/// A& the_instance = A::instance();      // throws if A is not default-constructible or A::set_instance() had not been called to create the instance
/// // ... or ...
/// A* the_instance_ptr = A::instance_ptr();  // alternative to A::instance(), instead of throwing, returns nullptr
/// // the instance of A will be destroyed with other static-linkage objects
/// \endcode
/// Singleton is thread-safe, in that multiple threads can call `A::instance()`, even if `A::set_instance()` has not been called yet.
// clang-format on
template <typename Derived>
class Singleton {
  // can't use std::is_default_constructible since Derived's ctors should be
  // private
  template <typename T, typename Enabler = void>
  struct is_default_constructible_helper : public std::false_type {};
  template <typename T>
  struct is_default_constructible_helper<T, std::void_t<decltype(T{})>>
      : public std::true_type {};
  constexpr static bool derived_is_default_constructible =
      is_default_constructible_helper<Derived>::value;

 public:
  /// @return reference to the instance
  /// @throw std::logic_error if the instance has not been constructed (because
  /// Derived is not default-constructible and set_instance() had not been
  /// called)
  static Derived& instance() {
    const auto& result_ptr = instance_accessor();
    if (result_ptr != nullptr) return *result_ptr;
    if constexpr (derived_is_default_constructible) {
      return set_default_instance();
    } else
      throw std::logic_error(
          "sequant::Singleton: is not default-constructible and set_instance() "
          "has not been called");
  }

  /// same as instance(), but returns pointer to the instance, or nullptr if
  /// the instance has not been constructed
  /// @return pointer to the instance, or nullptr if it has not yet been
  /// constructed
  static Derived* instance_ptr() {
    const auto& result_ptr = instance_accessor();
    if (result_ptr != nullptr) return result_ptr.get();
    if constexpr (derived_is_default_constructible) {
      return set_default_instance();
    } else
      return nullptr;
  }

  /// Constructs the instance. This must be called if Derived is not
  /// default-constructible.
  /// @tparam Args a parameter pack type
  /// @param args a parameter pack
  /// @return reference to the newly-created instance
  /// @throw std::logic_error if the instance has already been constructed
  template <typename... Args>
  static Derived& set_instance(Args&&... args) {
    //    WARNING: can't check constructibility since the ctor may be private
    //    static_assert(std::is_constructible_v<Derived, Args...>,
    //                  "sequant::Singleton::set_instance: Derived is not
    //                  constructible with Args");
    std::scoped_lock lock(instance_mutex());
    if (instance_accessor() != nullptr)
      throw std::logic_error(
          "sequant::Singleton::set_instance: instance has already been "
          "constructed");
    instance_accessor() = std::move(
        std::unique_ptr<Derived>(new Derived(std::forward<Args>(args)...)));
    return *instance_accessor();
  }

 protected:
  template <typename... Args>
  Singleton(Args&&... args) {}  // all constructors are private

  static auto& instance_accessor() {
    static std::unique_ptr<Derived> instance;
    return instance;
  }
  // provides mutex that controls access to instance_accessor()'s object
  static auto& instance_mutex() {
    static std::mutex mtx;
    return mtx;
  }

 private:
  /// Constructs a default-constructed instance. If called from multiple threads
  /// only the first call will construct the instance.
  static Derived& set_default_instance() {
    std::scoped_lock lock(instance_mutex());
    if (!instance_accessor()) {
      instance_accessor() = std::move(std::unique_ptr<Derived>(new Derived));
    }
    return *instance_accessor();
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_SINGLETON_HPP
