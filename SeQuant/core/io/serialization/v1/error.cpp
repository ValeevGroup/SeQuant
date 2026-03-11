#include <SeQuant/core/io/serialization/v1/error.hpp>

namespace sequant::io::serialization::v1 {

Error::Error(std::size_t offset, std::size_t length, std::string msg)
    : SerializationError(std::move(msg)), offset(offset), length(length) {}

}  // namespace sequant::io::serialization::v1
