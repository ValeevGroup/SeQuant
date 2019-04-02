//
// Created by Eduard Valeyev on 2019-04-02.
//

#ifndef SEQUANT_OPTIONAL_HPP
#define SEQUANT_OPTIONAL_HPP

#if __cplusplus >= 201703L
# if __has_include(<optional>)
#  include <optional>
# elif __has_include(<experimental/optional>)
#  include <experimental/optional>
namespace std {
using std::experimental::optional;
using std::experimental::nullopt;
}
# else
#  error "std::optional is not available even though C++17 is supported"
# endif
#else
# error "C++17 is required to include optional.hpp"
#endif

#endif  // SEQUANT_OPTIONAL_HPP
