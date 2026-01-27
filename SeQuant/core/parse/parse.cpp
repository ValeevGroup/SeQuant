#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string_view>

namespace sequant {

ParseError::ParseError(std::size_t offset, std::size_t length,
                       std::string message)
    : std::runtime_error(std::move(message)), offset(offset), length(length) {}

#define SEQUANT_RESOLVE_PARSE_FUNC(name, stringType, returnType)   \
  returnType name(stringType input, const ParseOptions &options) { \
    switch (options.syntax) {                                      \
      case ParseSyntax::V1:                                        \
        return parse::v1::name(input, options);                    \
    }                                                              \
                                                                   \
    SEQUANT_UNREACHABLE;                                           \
  }

SEQUANT_RESOLVE_PARSE_FUNC(parse_expr, std::wstring_view, ExprPtr)
SEQUANT_RESOLVE_PARSE_FUNC(parse_expr, std::string_view, ExprPtr)
SEQUANT_RESOLVE_PARSE_FUNC(parse_result_expr, std::wstring_view, ResultExpr)
SEQUANT_RESOLVE_PARSE_FUNC(parse_result_expr, std::string_view, ResultExpr)

#define RESOLVE_DEPARSE_FUNC(argType)                                       \
  std::wstring deparse(const argType &arg, const DeparseOptions &options) { \
    switch (options.syntax) {                                               \
      case ParseSyntax::V1:                                                 \
        return parse::v1::deparse(arg, options);                            \
    }                                                                       \
                                                                            \
    SEQUANT_UNREACHABLE;                                                    \
  }

RESOLVE_DEPARSE_FUNC(ResultExpr)
RESOLVE_DEPARSE_FUNC(ExprPtr)
RESOLVE_DEPARSE_FUNC(Expr)
RESOLVE_DEPARSE_FUNC(AbstractTensor)
RESOLVE_DEPARSE_FUNC(Index)
RESOLVE_DEPARSE_FUNC(NormalOperator<Statistics::BoseEinstein>)
RESOLVE_DEPARSE_FUNC(NormalOperator<Statistics::FermiDirac>)

#undef SEQUANT_RESOLVE_PARSE_FUNC
#undef SEQUANT_RESOLVE_DEPARSE_FUNC

}  // namespace sequant
