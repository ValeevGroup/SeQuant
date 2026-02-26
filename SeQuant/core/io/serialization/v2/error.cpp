#include <SeQuant/core/io/serialization/v2/error.hpp>

namespace sequant::io::serialization::v2 {

Error::Error(std::size_t line, std::size_t column, std::string rule,
             std::string msg)
    : SerializationError(std::move(msg)),
      line(line),
      column(column),
      rule(std::move(rule)) {}

}  // namespace sequant::io::serialization::v2
