#ifndef SEQUANT_PARSE_HPP
#define SEQUANT_PARSE_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

///
///  Create SeQuant expression from string input.
///
/// @author: Bimal Gaudel
/// version: 21 July, 2021
///

namespace sequant {

struct ParseError : std::runtime_error {
  std::size_t offset;
  std::size_t length;

  ParseError(std::size_t offset, std::size_t length, std::string message);
};

struct ParseOptions {
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
};

struct DeparseOptions {
  /// Whether to explicitly annotate tensor symmetries
  bool annot_symm = true;
};

#define SEQUANT_DECLARE_PARSE_FUNC(name, returnType)                          \
  returnType name(std::wstring_view input, const ParseOptions &options = {}); \
  returnType name(std::string_view input, const ParseOptions &options = {});

// clang-format off
/// \brief Construct expressions from string representations
///
/// \param input The input to parse
/// \param options Customization options
/// \return SeQuant expression.
SEQUANT_DECLARE_PARSE_FUNC(parse_expr, ExprPtr);

/// \sa parse_expr
SEQUANT_DECLARE_PARSE_FUNC(parse_result_expr, ResultExpr);

#undef SEQUANT_DECLARE_PARSE_FUNC


#define SEQUANT_DECLARE_DEPARSE_FUNC(name) \
	std::wstring name(const ResultExpr &expr, const DeparseOptions &options = {}); \
	std::wstring name(const ExprPtr &expr, const DeparseOptions &options = {}); \
	std::wstring name(const Expr &expr, const DeparseOptions &options = {}); \
	std::wstring name(const AbstractTensor &expr, const DeparseOptions &options = {}); \
	std::wstring name(const Index &index, const DeparseOptions &options = {}); \
	template< Statistics S> \
	std::wstring name(const NormalOperator<S> &nop, const DeparseOptions &options = {}); \

///
/// Get a parsable string from an expression.
///
/// An expression that has a flat structure (ie. product does not contain
/// sum or product subexpressions) is guaranteed to be reparsable to itself so
/// that equality comparison holds: x == parse_expr(deparse_expr(x)). For nested
/// expressions the equality comparison is not guaranteed, however, for such
/// an expression x, its corresponding evaluation node, eval_node(x), and
/// the parsed evaluation node, eval_node(parse_expr(deparse_expr(x))) are
/// equivalent.
///
/// \param expr Expression to stringify that can be re-parsed to itself.
/// \param options Customization options
/// \return wstring of the expression.
SEQUANT_DECLARE_DEPARSE_FUNC(deparse)

}  // namespace sequant

#endif  // SEQUANT_PARSE_HPP
