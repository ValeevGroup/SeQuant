#ifndef SEQUANT_CORE_IO_SERIALIZATION_V2_ERROR_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_V2_ERROR_HPP

#include <SeQuant/core/io/serialization/serialization.hpp>

#include <cstddef>
#include <string>

namespace sequant::io::serialization::v2 {

struct Error : SerializationError {
  /// 1-based position in the input the error refers to; both are 0 if it
  /// concerns the input as a whole
  std::size_t line;
  std::size_t column;
  /// the grammar rule that failed
  std::string rule;

  Error(std::size_t line, std::size_t column, std::string rule,
        std::string msg);
};

}  // namespace sequant::io::serialization::v2

#endif  // SEQUANT_CORE_IO_SERIALIZATION_V2_ERROR_HPP
