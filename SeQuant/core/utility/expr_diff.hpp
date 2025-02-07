#ifndef SEQUANT_EXPR_DIFF_HPP
#define SEQUANT_EXPR_DIFF_HPP

#include <SeQuant/core/expr_fwd.hpp>

#include <string>

namespace sequant {

/// @returns A string describing (some of) the difference between the given
/// expressions. An empty diff means that they are equal. The produced diff is
/// meant to be (resonably) human-readable.
std::string diff(const Expr &lhs, const Expr &rhs);

}  // namespace sequant

#endif  // SEQUANT_EXPR_DIFF_HPP
