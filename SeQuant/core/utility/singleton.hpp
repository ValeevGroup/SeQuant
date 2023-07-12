//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_CORE_UTILITY_SINGLETON_HPP
#define SEQUANT_CORE_UTILITY_SINGLETON_HPP

#include <memory>
#include <type_traits>

namespace sequant {

// clang-format off
/// @brief Singleton base class
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
/// A::set_instance(Args...);  // this creates the first instance of A, if A is not default-constructible, otherwise can skip this
/// A& the_instance_ref = A::get_instance();      // throws if the instance of A had not been created
/// A* the_instance_ptr = A::get_instance_ptr();  // returns nullptr if the instance of A had not been created
/// // the instance of A will be destroyed with other static-linkage objects
/// \endcode
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
  /// @throw std::logic_error if the reference has not been contructed (because
  /// Derived is not default-constructible and set_instance() had not been
  /// called)
  static Derived& get_instance() {
    const auto& result_ptr = instance_accessor();
    if (result_ptr != nullptr) return *result_ptr;
    if constexpr (derived_is_default_constructible) {
      set_instance();
      return *instance_accessor();
    } else
      throw std::logic_error(
          "sequant::Singleton: is not default-constructible and set_instance() "
          "has not been called");
  }

  /// @return pointer to the instance, or nullptr if it has not yet been
  /// constructed
  static Derived* get_instance_ptr() {
    const auto result_ptr = instance_accessor();
    if (result_ptr != nullptr) return result_ptr.get();
    if constexpr (derived_is_default_constructible) {
      set_instance();
      return instance_accessor();
    } else
      return nullptr;
  }

  /// Constructs the instance. This must be called if Derived is not
  /// default-constructible.
  /// @tparam Args a parameter pack type
  /// @param args a parameter pack
  template <typename... Args>
  static void set_instance(Args&&... args) {
    assert(instance_accessor() == nullptr);
    instance_accessor() = std::move(
        std::unique_ptr<Derived>(new Derived(std::forward<Args>(args)...)));
  }

 protected:
  template <typename... Args>
  Singleton(Args&&... args) {}  // all constructors are private

  static auto& instance_accessor() {
    static std::unique_ptr<Derived> instance(nullptr);
    return instance;
  }
};

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_SINGLETON_HPP
