#ifndef SEQUANT_TESTS_UTILS_HPP_
#define SEQUANT_TESTS_UTILS_HPP_

#include <SeQuant/core/math.hpp>

#include <algorithm>
#include <string>

/// Checks all permutations of summands in order to match the LaTeX
/// representation to the given string. We do this as the order of summands is
/// dependent on the used hash function (which messes with tests) but ultimately
/// the order of summands does not matter whatsoever. Thus, it is sufficient to
/// test whether some permutation of summands yields the desired LaTeX
/// representation.
#define REQUIRE_SUM_EQUAL(sum, str)                                           \
  REQUIRE(sum.is<Sum>());                                                     \
  if (to_latex(sum) != str) {                                                 \
    for (sequant::intmax_t i = 0; i < sequant::factorial(sum->size()); ++i) { \
      std::next_permutation(sum->begin(), sum->end());                        \
      if (to_latex(sum) == str) {                                             \
        break;                                                                \
      }                                                                       \
    }                                                                         \
  }                                                                           \
  REQUIRE(to_latex(sum) == std::wstring(str))

#endif
