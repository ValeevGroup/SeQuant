#ifndef SEQUANT_CORE_IO_CONCEPTS_HPP
#define SEQUANT_CORE_IO_CONCEPTS_HPP

#include <SeQuant/core/io/latex/latex.hpp>

namespace sequant::io {

template <typename T>
concept convertible_to_latex =
    requires(const T &t) { io::latex::to_string(t); };

}  // namespace sequant::io

#endif  // SEQUANT_CORE_IO_CONCEPTS_HPP
