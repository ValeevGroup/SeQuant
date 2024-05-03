//
// Created by Nakul Teke on 2019-09-03.
//

#ifndef SEQUANT_ANY_HPP
#define SEQUANT_ANY_HPP

// Include C++17 any header, if available AND functional
#if __cplusplus >= 201703L  // c++17

// macos < 10.14 do not have any_cast in their libc++
#include <ciso646>  // see https://stackoverflow.com/a/31658120
#if defined(_LIBCPP_VERSION) && defined(__APPLE__)
#include <Availability.h>
#ifdef __MAC_OS_X_VERSION_MIN_ALLOWED
#if __MAC_OS_X_VERSION_MIN_ALLOWED >= 10140
#define SEQUANT_HAS_CXX17_ANY
#endif  //  10.14 or later
#endif  // have macos version
// #else   // libc++ on macos
// #define SEQUANT_HAS_CXX17_ANY
#endif  // libc++ on macos
#endif  // c++17

#ifdef SEQUANT_HAS_CXX17_ANY
#include <any>
#endif

namespace sequant {

#ifdef SEQUANT_HAS_CXX17_ANY
using std::any;
using std::any_cast;
using std::bad_any_cast;
#else

//  #error "37 any.hpp"
class bad_any_comparable_cast;

namespace detail {

#ifdef SEQUANT_HAS_CXX17_ANY
class bad_any_comparable_cast : public std::bad_any_cast {
 public:
  bad_any_comparable_cast() = default;
  virtual ~bad_any_comparable_cast() {}
  virtual const char *what() const noexcept {
    return "Bad any_comparable_cast";
  }
};
#else
class bad_any_comparable_cast : public std::bad_cast {
 public:
  bad_any_comparable_cast() = default;
  virtual ~bad_any_comparable_cast() {}
  virtual const char *what() const noexcept {
    return "Bad any_comparable_cast";
  }
};
#endif

}  // namespace detail
#endif  // SEQUANT_HAS_CXX17_ANY

}  // namespace sequant

#endif  // SEQUANT_ANY_HPP
