#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <concepts>
#include <string_view>

namespace sequant::io::serialization {

SerializationError::SerializationError(std::string message)
    : Exception(std::move(message)) {}

template <typename StringType, typename ExprType>
concept v1_can_deserialize = requires(StringType input) {
  serialization::v1::from_string<ExprType>(input);
};

template <typename StringType, typename ExprType>
concept v2_can_deserialize =
    requires(StringType input, DeserializationOptions opts) {
      serialization::v2::from_string<ExprType>(input);
    };

static_assert(!v2_can_deserialize<std::string, ExprPtr>);

template <typename StringType, typename ExprType>
ExprType from_string_indirection(StringType input,
                                 const DeserializationOptions &options) {
  // Note: We need this indirection as the if constexpr construct only works (in
  // the way we need) if the condition depends on a template parameter. If it
  // didn't, we'd get a bunch of "call to deleted function" errors as the body
  // of the if is checked even though the if is false.
  switch (options.syntax) {
    case SerializationSyntax::V1:
      if constexpr (v1_can_deserialize<StringType, ExprType>) {
        return serialization::v1::from_string<ExprType>(input, options);
      } else {
        throw Exception(
            "Deserialization of this type is not supported with syntax V1");
      }
    case SerializationSyntax::V2:
      if constexpr (v2_can_deserialize<StringType, ExprType>) {
        return serialization::v2::from_string<ExprType>(input, options);
      } else {
        throw Exception(
            "Deserialization of this type is not supported with syntax V2");
      }
  }

  SEQUANT_UNREACHABLE;
}

#define SEQUANT_RESOLVE_DESERIALIZATION_FUNC(StringType, ExprType)        \
  template <>                                                             \
  ExprType from_string<ExprType>(StringType input,                        \
                                 const DeserializationOptions &options) { \
    return from_string_indirection<StringType, ExprType>(input, options); \
  }

SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, ExprPtr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, ExprPtr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, ResultExpr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, ResultExpr)
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, Constant);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, Constant);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, Variable);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, Variable);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, Tensor);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, Tensor);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, Product);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, Product);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::string_view, Sum);
SEQUANT_RESOLVE_DESERIALIZATION_FUNC(std::wstring_view, Sum);

template <typename ExprType>
concept v1_can_serialize = requires(const ExprType &expr) {
  { serialization::v1::to_string(expr) } -> std::convertible_to<std::wstring>;
};

template <typename ExprType>
concept v2_can_serialize = requires(const ExprType &expr) {
  { serialization::v2::to_string(expr) } -> std::convertible_to<std::wstring>;
};

template <typename Arg>
std::wstring to_string_indirection(const Arg &arg,
                                   const SerializationOptions &options) {
  switch (options.syntax) {
    case SerializationSyntax::V1:
      if constexpr (v1_can_serialize<Arg>) {
        return serialization::v1::to_string(arg, options);
      } else {
        throw Exception(
            "Serialization of this type is not supported with syntax V1");
      }
    case SerializationSyntax::V2:
      if constexpr (v2_can_serialize<Arg>) {
        return serialization::v2::to_string(arg, options);
      } else {
        throw Exception(
            "Serialization of this type is not supported with syntax V2");
      }
  }

  SEQUANT_UNREACHABLE;
}

#define SEQUANT_RESOLVE_SERIALIZE_FUNC(argType)                 \
  std::wstring to_string(const argType &arg,                    \
                         const SerializationOptions &options) { \
    return to_string_indirection<argType>(arg, options);        \
  }

SEQUANT_RESOLVE_SERIALIZE_FUNC(ResultExpr)
SEQUANT_RESOLVE_SERIALIZE_FUNC(ExprPtr)
SEQUANT_RESOLVE_SERIALIZE_FUNC(Expr)
SEQUANT_RESOLVE_SERIALIZE_FUNC(AbstractTensor)
SEQUANT_RESOLVE_SERIALIZE_FUNC(Index)
SEQUANT_RESOLVE_SERIALIZE_FUNC(NormalOperator<Statistics::BoseEinstein>)
SEQUANT_RESOLVE_SERIALIZE_FUNC(NormalOperator<Statistics::FermiDirac>)

#undef SEQUANT_RESOLVE_SERIALIZE_FUNC
#undef SEQUANT_RESOLVE_DESERIALIZATION_FUNC

}  // namespace sequant::io::serialization
