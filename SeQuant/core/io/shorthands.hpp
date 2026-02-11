#ifndef SEQUANT_CORE_IO_SHORTHANDS_HPP
#define SEQUANT_CORE_IO_SHORTHANDS_HPP

/// @file This header includes a couple of shorthands that essentially consist
/// of simple wrapper functions directly in the sequant namespace which delegate
/// to the respective functions in the sequant::io namespace.

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/io/concepts.hpp>
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

#include <string_view>

namespace sequant {

/// Shorthand for io::latex::to_string
template <typename T>
  requires(io::convertible_to_latex<T>)
decltype(auto) to_latex(T &&expr) {
  return io::latex::to_string(std::forward<T>(expr));
}

/// Shorthand for io::serialization::to_string
template <typename T>
  requires(io::serializable<T>)
decltype(auto) serialize(
    T &&t, const io::serialization::SerializationOptions &options = {}) {
  return io::serialization::to_string(std::forward<T>(t), options);
}

/// Shorthand for io::serialization::from_string
template <typename T = ExprPtr>
  requires(io::deserializable<T>)
decltype(auto) deserialize(
    std::string_view input,
    const io::serialization::DeserializationOptions &options = {}) {
  return io::serialization::from_string<T>(input, options);
}

template <typename T = ExprPtr>
  requires(io::deserializable<T>)
decltype(auto) deserialize(
    std::wstring_view input,
    const io::serialization::DeserializationOptions &options = {}) {
  return io::serialization::from_string<T>(input, options);
}

}  // namespace sequant

#endif  // SEQUANT_CORE_IO_SHORTHANDS_HPP
