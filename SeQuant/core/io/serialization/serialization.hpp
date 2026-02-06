#ifndef SEQUANT_CORE_IO_SERIALIZATION_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/exception.hpp>

#include <optional>
#include <string>
#include <string_view>

namespace sequant::io::serialization {

struct SerializationError : Exception {
  std::size_t offset;
  std::size_t length;

  SerializationError(std::size_t offset, std::size_t length,
                     std::string message);
};

/// Specifies the syntax of the textual input/representation to use. All
/// potential changes made to the syntax within a given version is understood to
/// be backwards compatible in the from_string(…) sense. That is older inputs
/// will continue to work as before. However, representations generated via
/// to_string(…) may not necessarily be deserializable by older version of
/// SeQuant.
///
/// @note Everything but Latest is considered to be deprecated by default
///       and support for it may be removed in future versions.
enum class SerializationSyntax {
  V1,

  Latest = V1
};

struct DeserializationOptions {
  /// The Symmetry (within bra and ket) to use, if none is specified in the
  /// input explicitly. The Context is queried in case this is not provided
  /// explicitly.
  std::optional<Symmetry> def_perm_symm = {};
  /// The BraKetSymmetry to use, if none is specified in the input explicitly.
  /// The Context is queried in case this is not provided explicitly.
  std::optional<BraKetSymmetry> def_braket_symm = {};
  /// The ColumnSymmetry to use, if none is specified in the input explicitly.
  /// The Context is queried in case this is not provided explicitly.
  std::optional<ColumnSymmetry> def_col_symm = {};
  /// The expected syntax version of the input
  SerializationSyntax syntax = SerializationSyntax::Latest;
};

struct SerializationOptions {
  /// Whether to explicitly annotate tensor symmetries
  bool annot_symm = true;
  /// The syntax version of the produced output
  SerializationSyntax syntax = SerializationSyntax::Latest;
};

#define SEQUANT_DECLARE_DESERIALIZATION_FUNC                          \
  template <typename T>                                               \
  T from_string(std::string_view input,                               \
                const DeserializationOptions &options = {}) = delete; \
  template <typename T>                                               \
  T from_string(std::wstring_view input,                              \
                const DeserializationOptions &options = {}) = delete;

#define SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION(ExprType)    \
  template <>                                                            \
  ExprType from_string<ExprType>(std::string_view input,                 \
                                 const DeserializationOptions &options); \
  template <>                                                            \
  ExprType from_string<ExprType>(std::wstring_view input,                \
                                 const DeserializationOptions &options);

/// \brief Construct expressions from string representations
///
/// \param input The input to deserialize
/// \param options Customization options
/// \return SeQuant expression.
SEQUANT_DECLARE_DESERIALIZATION_FUNC;

SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION(ExprPtr);
SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION(ResultExpr);

#define SEQUANT_DECLARE_SERIALIZATION_FUNC                          \
  std::wstring to_string(const ResultExpr &expr,                    \
                         const SerializationOptions &options = {}); \
  std::wstring to_string(const ExprPtr &expr,                       \
                         const SerializationOptions &options = {}); \
  std::wstring to_string(const Expr &expr,                          \
                         const SerializationOptions &options = {}); \
  std::wstring to_string(const AbstractTensor &expr,                \
                         const SerializationOptions &options = {}); \
  std::wstring to_string(const Index &index,                        \
                         const SerializationOptions &options = {}); \
  template <Statistics S>                                           \
  std::wstring to_string(const NormalOperator<S> &nop,              \
                         const SerializationOptions &options = {});

///
/// Get a serialized string from an expression.
///
/// An expression that has a flat structure (ie. product does not contain
/// sum or product subexpressions) is guaranteed to be reparsable to itself so
/// that equality comparison holds: x == from_string(to_string(x)). For nested
/// expressions the equality comparison is not guaranteed, however, for such
/// an expression x, its corresponding evaluation node, eval_node(x), and
/// the deserialized evaluation node, eval_node(from_string(to_string(x))) are
/// equivalent.
///
/// \param expr Expression to serialize
/// \param options Customization options
/// \return wstring of the expression.
SEQUANT_DECLARE_SERIALIZATION_FUNC

// Versioned variants
namespace v1 {
SEQUANT_DECLARE_DESERIALIZATION_FUNC;

SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION(ExprPtr);
SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION(ResultExpr);

SEQUANT_DECLARE_SERIALIZATION_FUNC
}  // namespace v1

#undef SEQUANT_DECLARE_DESERIALIZATION_FUNC
#undef SEQUANT_DECLARE_DESERIALIZATION_FUNC_SPECIALIZATION
#undef SEQUANT_DECLARE_SERIALIZATION_FUNC

}  // namespace sequant::io::serialization

#endif  // SEQUANT_CORE_IO_SERIALIZATION_HPP
