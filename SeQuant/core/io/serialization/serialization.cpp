#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string_view>

namespace sequant::io::serialization {

SerializationError::SerializationError(std::size_t offset, std::size_t length,
                                       std::string message)
    : Exception(std::move(message)), offset(offset), length(length) {}

#define SEQUANT_RESOLVE_DESERIALIZATION_FUNC(stringType, ExprType)        \
  template <>                                                             \
  ExprType from_string<ExprType>(stringType input,                        \
                                 const DeserializationOptions &options) { \
    switch (options.syntax) {                                             \
      case SerializationSyntax::V1:                                       \
        return serialization::v1::from_string<ExprType>(input, options);  \
    }                                                                     \
                                                                          \
    SEQUANT_UNREACHABLE;                                                  \
  }

SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, ExprPtr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, ExprPtr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, ResultExpr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, ResultExpr)

#define RESOLVE_SERIALIZE_FUNC(argType)                         \
  std::wstring to_string(const argType &arg,                    \
                         const SerializationOptions &options) { \
    switch (options.syntax) {                                   \
      case SerializationSyntax::V1:                             \
        return serialization::v1::to_string(arg, options);      \
    }                                                           \
                                                                \
    SEQUANT_UNREACHABLE;                                        \
  }

RESOLVE_SERIALIZE_FUNC(ResultExpr)
RESOLVE_SERIALIZE_FUNC(ExprPtr)
RESOLVE_SERIALIZE_FUNC(Expr)
RESOLVE_SERIALIZE_FUNC(AbstractTensor)
RESOLVE_SERIALIZE_FUNC(Index)
RESOLVE_SERIALIZE_FUNC(NormalOperator<Statistics::BoseEinstein>)
RESOLVE_SERIALIZE_FUNC(NormalOperator<Statistics::FermiDirac>)

#undef SEQUANT_RESOLVE_PARSE_FUNC
#undef SEQUANT_RESOLVE_DEPARSE_FUNC

}  // namespace sequant::io::serialization
