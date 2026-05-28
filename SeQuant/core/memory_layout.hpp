#ifndef SEQUANT_CORE_MEMORY_LAYOUT_HPP
#define SEQUANT_CORE_MEMORY_LAYOUT_HPP

/// Different ways of laying out a multidimensional tensor/array in memory
enum class MemoryLayout {
  RowMajor,
  ColumnMajor,
  Unspecified,
};

#endif  // SEQUANT_CORE_MEMORY_LAYOUT_HPP
