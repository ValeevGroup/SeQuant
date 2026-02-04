#ifndef SEQUANT_CORE_IO_SHORTHANDS_HPP
#define SEQUANT_CORE_IO_SHORTHANDS_HPP

#include <SeQuant/core/io/concepts.hpp>
#include <SeQuant/core/io/latex/latex.hpp>

/// @file This header includes a couple of shorthands that essentially consist
/// of simple wrapper functions directly in the sequant namespace which delegate
/// to the respective functions in the sequant::io namespace.

namespace sequant {

/// Shorthand for io::latex::to_string
template <typename T>
  requires(io::convertible_to_latex<T>)
decltype(auto) to_latex(T &&expr) {
  return io::latex::to_string(std::forward<T>(expr));
}

}  // namespace sequant

#endif  // SEQUANT_CORE_IO_SHORTHANDS_HPP
