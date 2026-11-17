#ifndef SEQUANT_CORE_MEMORY_LAYOUT_HPP
#define SEQUANT_CORE_MEMORY_LAYOUT_HPP

#include <cstdint>

/// Different ways of laying out a multidimensional tensor/array in memory
enum class MemoryLayout : std::uint8_t {
  RowMajor,
  ColumnMajor,
  Unspecified,
};

#endif  // SEQUANT_CORE_MEMORY_LAYOUT_HPP
