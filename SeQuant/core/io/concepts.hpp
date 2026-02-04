#ifndef SEQUANT_CORE_IO_CONCEPTS_HPP
#define SEQUANT_CORE_IO_CONCEPTS_HPP

#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

namespace sequant::io {

template <typename T>
concept convertible_to_latex =
    requires(const T &t) { io::latex::to_string(t); };

template <typename T>
concept serializable =
    requires(const T &t) { io::serialization::to_string(t); };

template <typename T>
concept deserializable = requires { io::serialization::from_string<T>(""); };

}  // namespace sequant::io

#endif  // SEQUANT_CORE_IO_CONCEPTS_HPP
