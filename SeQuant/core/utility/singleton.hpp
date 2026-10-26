//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_CORE_UTILITY_SINGLETON_HPP
#define SEQUANT_CORE_UTILITY_SINGLETON_HPP

#include <SeQuant/core/utility/exception.hpp>
#include <atomic>
#include <cassert>
#include <memory>
#include <mutex>

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
  template <typename T>
  static constexpr bool has_default_ctor = requires { T{}; };

 public:
  /// @return reference to the instance
  /// @throw Exception if the instance has not been constructed (because
  /// Derived is not default-constructible and set_instance() had not been
  /// called)
  static Derived& instance() {
    if (auto* result_ptr = published_instance().load(std::memory_order_acquire))
      return *result_ptr;
    if constexpr (has_default_ctor<Derived>) {
      return set_default_instance();
    } else
      throw Exception(
          "sequant::Singleton: is not default-constructible and set_instance() "
          "has not been called");
  }

  /// same as instance(), but returns pointer to the instance, or nullptr if
  /// the instance has not been constructed
  /// @return pointer to the instance, or nullptr if it has not yet been
  /// constructed
  static Derived* instance_ptr() {
    if (auto* result_ptr = published_instance().load(std::memory_order_acquire))
      return result_ptr;
    if constexpr (has_default_ctor<Derived>) {
      return set_default_instance();
    } else
      return nullptr;
  }

  /// Constructs the instance. This must be called if Derived is not
  /// default-constructible.
  /// @tparam Args a parameter pack type
  /// @param args a parameter pack
  /// @return reference to the newly-created instance
  /// @throw Exception if the instance has already been constructed
  template <typename... Args>
  static Derived& set_instance(Args&&... args) {
    //    WARNING: can't check constructibility since the ctor may be private
    //    static_assert(std::is_constructible_v<Derived, Args...>,
    //                  "sequant::Singleton::set_instance: Derived is not
    //                  constructible with Args");
    std::scoped_lock lock(instance_mutex());
    if (instance_accessor() != nullptr)
      throw Exception(
          "sequant::Singleton::set_instance: instance has already been "
          "constructed");
    // Note: std::make_unique can't deal with private friend ctors
    return publish(
        std::unique_ptr<Derived>(new Derived(std::forward<Args>(args)...)));
  }

 protected:
  template <typename... Args>
  Singleton(Args&&...) {}  // all constructors are private

  /// provides mutex that serializes construction of the instance
  static auto& instance_mutex() {
    static std::mutex mtx;
    return mtx;
  }

 private:
  /// owns the instance; read and written only under instance_mutex(), and
  /// written only by publish()
  static auto& instance_accessor() {
    static std::unique_ptr<Derived> instance;
    return instance;
  }
  /// the constructed instance, published for lock-free reads; null until
  /// publish() has been called
  static std::atomic<Derived*>& published_instance() {
    static std::atomic<Derived*> instance{nullptr};
    return instance;
  }

  /// Takes ownership of @p instance and makes it visible to instance() and
  /// instance_ptr(). Must be called under instance_mutex(), at most once.
  static Derived& publish(std::unique_ptr<Derived> instance) {
    auto& owner = instance_accessor();
    owner = std::move(instance);
    published_instance().store(owner.get(), std::memory_order_release);
    return *owner;
  }

  /// Constructs a default-constructed instance. If called from multiple threads
  /// only the first call will construct the instance.
  static Derived& set_default_instance()
    requires(has_default_ctor<Derived>)
  {
    std::scoped_lock lock(instance_mutex());
    if (instance_accessor()) return *instance_accessor();
    // Note: std::make_unique can't deal with private friend ctors
    return publish(std::unique_ptr<Derived>(new Derived()));
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_SINGLETON_HPP
