#ifndef SEQUANT_CORE_IO_SERIALIZATION_V1_ERROR_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_V1_ERROR_HPP

#include <SeQuant/core/io/serialization/serialization.hpp>

#include <cstddef>
#include <string>

namespace sequant::io::serialization::v1 {

struct Error : SerializationError {
  std::size_t offset;
  std::size_t length;

  Error(std::size_t offset, std::size_t length, std::string msg);
};

}  // namespace sequant::io::serialization::v1

#endif  // SEQUANT_CORE_IO_SERIALIZATION_V1_ERROR_HPP
