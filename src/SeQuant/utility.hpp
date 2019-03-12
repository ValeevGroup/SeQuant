//
// Created by Eduard Valeyev on 2019-03-11.
//

#ifndef SEQUANT_UTILITY_HPP
#define SEQUANT_UTILITY_HPP

namespace sequant {

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
template <typename Derived>
class Singleton {
 public:
  static Derived& get_instance() {
    const auto& result_ptr = instance_accessor();
    if (result_ptr != nullptr)
      return *result_ptr;
    if constexpr (std::is_default_constructible_v<Derived>) {
      set_instance();
      return *instance_accessor();
    }
    else
      throw std::logic_error("sequant::Singleton: is not default-constructible and set_instance() has not been called");
  }

  static Derived* get_instance_ptr() {
    const auto result_ptr = instance_accessor();
    if (result_ptr != nullptr)
      return result_ptr.get();
    if constexpr (std::is_default_constructible_v<Derived>) {
      set_instance();
      return instance_accessor();
    }
    else
      return nullptr;
  }

  template <typename ... Args> static void set_instance(Args&&...args) {
    assert(instance_accessor() == nullptr);
    instance_accessor() = std::move(std::unique_ptr<Derived>(new Derived(std::forward<Args>(args)...)));
  }

 protected:
  template <typename ... Args> Singleton(Args&& ... args) {}  // all constructors are private

  static auto& instance_accessor() {
    static std::unique_ptr<Derived> instance(nullptr);
    return instance;
  }
};

struct Logger : public Singleton<Logger> {
  bool wick_contract = false;
  bool wick_reduce = false;
  bool canonicalize = false;
 private:
  friend class Singleton<Logger>;
  Logger(int log_level = 0) {
    if (log_level > 0) {
      wick = true;
      canonicalize = true;
    }
  }
};

}

#endif //SEQUANT_UTILITY_HPP
