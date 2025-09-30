//
// Created by Eduard Valeyev on 7/31/23
//

#ifndef SEQUANT_CORE_UTILITY_MACROS_H
#define SEQUANT_CORE_UTILITY_MACROS_H

#include <cassert>
#include <cstdlib>
#include <iostream>

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

#define SEQUANT_ABORT(msg)       \
  assert(false && msg);          \
  std::cerr << msg << std::endl; \
  std::abort();

#if defined(__cpp_lib_unreachable)
#include <utility>

#define SEQUANT_UNREACHABLE_TOKEN std::unreachable()
#elif defined(SEQUANT_CXX_COMPILER_IS_GCC) || \
    defined(SEQUANT_CXX_COMPILER_IS_CLANG)
#define SEQUANT_UNREACHABLE_TOKEN __builtin_unreachable()
#else
// Compiler-independent fallback
#define SEQUANT_UNREACHABLE_TOKEN std::abort()
#endif

#define SEQUANT_UNREACHABLE                      \
  do {                                           \
    assert(false && "reached unreachable code"); \
    SEQUANT_UNREACHABLE_TOKEN;                   \
  } while (0)

#endif  // SEQUANT_CORE_UTILITY_MACROS_H
