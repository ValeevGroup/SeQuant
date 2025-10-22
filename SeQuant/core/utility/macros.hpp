//
// Created by Eduard Valeyev on 7/31/23
//

#ifndef SEQUANT_CORE_UTILITY_MACROS_H
#define SEQUANT_CORE_UTILITY_MACROS_H

#include <SeQuant/core/utility/exception.hpp>

#include <cstdlib>
#include <iostream>
#include <source_location>
#include <sstream>
#include <utility>

/* detect C++ compiler id:
- ids taken from CMake
- macros are discussed at https://sourceforge.net/p/predef/wiki/Compilers/
*/
#define SEQUANT_CXX_COMPILER_ID_GNU 0
#define SEQUANT_CXX_COMPILER_ID_Clang 1
#define SEQUANT_CXX_COMPILER_ID_AppleClang 2
#define SEQUANT_CXX_COMPILER_ID_XLClang 3
#define SEQUANT_CXX_COMPILER_ID_Intel 4
#if defined(__INTEL_COMPILER_BUILD_DATE) /* macros like __ICC and even       \
                                            __INTEL_COMPILER can be affected \
                                            by command options like -no-icc */
#define SEQUANT_CXX_COMPILER_ID SEQUANT_CXX_COMPILER_ID_Intel
#define SEQUANT_CXX_COMPILER_IS_ICC 1
#endif
#if defined(__clang__) && !defined(SEQUANT_CXX_COMPILER_IS_ICC)
#define SEQUANT_CXX_COMPILER_IS_CLANG 1
#if defined(__apple_build_version__)
#define SEQUANT_CXX_COMPILER_ID SEQUANT_CXX_COMPILER_ID_AppleClang
#elif defined(__ibmxl__)
#define SEQUANT_CXX_COMPILER_ID SEQUANT_CXX_COMPILER_ID_XLClang
#else
#define SEQUANT_CXX_COMPILER_ID SEQUANT_CXX_COMPILER_ID_Clang
#endif
#endif
#if defined(__GNUG__) && !defined(SEQUANT_CXX_COMPILER_IS_ICC) && \
    !defined(SEQUANT_CXX_COMPILER_IS_CLANG)
#define SEQUANT_CXX_COMPILER_ID SEQUANT_CXX_COMPILER_ID_GNU
#define SEQUANT_CXX_COMPILER_IS_GCC 1
#endif

/* ----------- pragma helpers ---------------*/
#define SEQUANT_PRAGMA(x) _Pragma(#x)
/* same as SEQUANT_PRAGMA(x), but expands x */
#define SEQUANT_XPRAGMA(x) SEQUANT_PRAGMA(x)

#define SEQUANT_CONCAT_IMPL(x, y) x##y
/* "concats" a and b without a space in between */
#define SEQUANT_CONCAT(x, y) SEQUANT_CONCAT_IMPL(x, y)
/* "concats" a and b with a space in between */
#define SEQUANT_CONCAT_W_SPACE(a, b) a b
#if defined(SEQUANT_CXX_COMPILER_IS_CLANG)
#define SEQUANT_PRAGMA_CLANG(x) \
  SEQUANT_XPRAGMA(SEQUANT_CONCAT_W_SPACE(clang, x))
#else
#define SEQUANT_PRAGMA_CLANG(x)
#endif
#if defined(SEQUANT_CXX_COMPILER_IS_GCC)
#define SEQUANT_PRAGMA_GCC(x) SEQUANT_XPRAGMA(SEQUANT_CONCAT_W_SPACE(GCC, x))
#else
#define SEQUANT_PRAGMA_GCC(x)
#endif

/* Defines the default error checking behavior */
#define SEQUANT_ASSERT_THROW 2
#define SEQUANT_ASSERT_ABORT 3
#define SEQUANT_ASSERT_IGNORE 4
#define SEQUANT_STRINGIFY(x) #x
#define SEQUANT_ASSERT_BEHAVIOR \
  SEQUANT_CONCAT(SEQUANT_ASSERT_, SEQUANT_ASSERT_BEHAVIOR_)
#if SEQUANT_ASSERT_BEHAVIOR != SEQUANT_ASSERT_IGNORE
#define SEQUANT_ASSERT_ENABLED
#endif

namespace sequant {

#ifdef SEQUANT_ASSERT_ENABLED
[[noreturn]]
#endif
inline void
assert_failed([[maybe_unused]] const std::string &errmsg,
              [[maybe_unused]] const std::source_location location =
                  std::source_location::current()) {
#ifdef SEQUANT_ASSERT_ENABLED
#if SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_THROW
  std::ostringstream oss;
  oss
#elif SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_ABORT
  std::cerr
#endif  // SEQUANT_ASSERT_BEHAVIOR
      << errmsg << " at " << location.file_name() << ":" << location.line()
      << " in function '" << location.function_name() << "'";
#if SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_THROW
  throw sequant::Exception(oss.str());
#elif SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_ABORT
  std::abort();
#endif  // SEQUANT_ASSERT_BEHAVIOR
#endif  // SEQUANT_ASSERT_ENABLED
}
}  // namespace sequant

#ifdef SEQUANT_ASSERT_ENABLED
#define SEQUANT_ASSERT_MESSAGE(EXPR, ...)                          \
  "SEQUANT_ASSERT(" SEQUANT_STRINGIFY(EXPR) ") failed" __VA_OPT__( \
      " with message '" __VA_ARGS__ "'")

#define SEQUANT_ASSERT(EXPR, ...)                                        \
  do {                                                                   \
    if (!(EXPR)) {                                                       \
      sequant::assert_failed(SEQUANT_ASSERT_MESSAGE(EXPR, __VA_ARGS__)); \
    }                                                                    \
  } while (0)
#else
#define SEQUANT_ASSERT(...) \
  do {                      \
  } while (0)
#endif

namespace sequant {
[[noreturn]] inline void abort_msg(
    const std::string &errmsg,
    const std::source_location location = std::source_location::current()) {
  std::cerr << errmsg << " at " << location.file_name() << ":"
            << location.line() << " in function '" << location.function_name()
            << "'";
  std::abort();
}
}  // namespace sequant

#define SEQUANT_ABORT(msg) sequant::abort_msg(msg)

#if defined(__cpp_lib_unreachable)

#define SEQUANT_UNREACHABLE_TOKEN std::unreachable()
#elif defined(SEQUANT_CXX_COMPILER_IS_GCC) || \
    defined(SEQUANT_CXX_COMPILER_IS_CLANG)
#define SEQUANT_UNREACHABLE_TOKEN __builtin_unreachable()
#else
// Compiler-independent fallback
#define SEQUANT_UNREACHABLE_TOKEN std::abort()
#endif

#define SEQUANT_UNREACHABLE                              \
  do {                                                   \
    SEQUANT_ASSERT(false && "reached unreachable code"); \
    SEQUANT_UNREACHABLE_TOKEN;                           \
  } while (0)

#endif  // SEQUANT_CORE_UTILITY_MACROS_H
