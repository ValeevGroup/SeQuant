#ifndef SEQUANT_CORE_UTILITY_OVERLOADS_HPP
#define SEQUANT_CORE_UTILITY_OVERLOADS_HPP

#include <initializer_list>

namespace sequant {

/// A class that can't be instantiated due to the private constructor.
/// Useful for enforcing initializer_lists can only ever be empty
class Nonconstructible {
 private:
  Nonconstructible() = default;
};

/// Tag for providing unambiguous overloads for catching function(bla, {})
/// where the initializer list could otherwise be matched to multiple other
/// overloads.
using EmptyInitializerList = std::initializer_list<Nonconstructible>;

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_OVERLOADS_HPP
