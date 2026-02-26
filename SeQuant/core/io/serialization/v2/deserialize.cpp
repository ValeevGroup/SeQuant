#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/io/serialization/v2/ast_conversions.hpp>
#include <SeQuant/core/io/serialization/v2/error.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <peglib.h>

#include <string_view>

namespace sequant::io::serialization::v2 {

// Defines a variable peg_serialization_grammar that contains the definition of
// the grammar (as a string_view)
#include "peg_grammar.ipp"

static peg::parser &get_parser(std::string_view start_rule) {
  static thread_local bool initialized = false;
  static thread_local std::string last_start;
  static thread_local peg::parser parser;

  peg::Log logger = [](std::size_t line, std::size_t column,
                       const std::string &msg, const std::string &rule) {
    throw Error(line, column, rule,
                "Deserialization failed at line " + std::to_string(line) + ":" +
                    std::to_string(column) + ": " + msg + " (" + rule + ")");
  };

  if (!initialized) {
    initialized = true;
    // This first logger will be notified in case there are issues with parsing
    // the grammar itself
    parser.set_logger([](std::size_t line, std::size_t column,
                         const std::string &msg, const std::string &rule) {
      throw Error(line, column, rule,
                  "Input grammar invalid at line " + std::to_string(line) +
                      ":" + std::to_string(column) + ": " + msg + " (" + rule +
                      ")");
    });

    // Note: we have to use ResultExpression in order to avoid being notified
    // about this rule not being referenced
    parser.load_grammar(peg_serialization_grammar, "ResultExpression");
    last_start = "ResultExpression";

    // This is the logger that we want to use from now on (for SeQuant syntax)
    parser.set_logger(logger);

    parser.enable_ast();
  }

  if (last_start != start_rule) {
    // Ability to change the start rule on-the-fly is implemented in
    // https://github.com/yhirose/cpp-peglib/pull/332
    // Until that is merged, we have to reload the grammar every time
    // we want to use a different start rule.

    // We assume that the grammar is valid (it has been loaded before) and in
    // order to ignore 'rule xy not referenced' logs (which our logger turns
    // into exceptions) we uninstall the logger before reloading the grammar
    parser.set_logger(peg::Log{});
    parser.load_grammar(peg_serialization_grammar, start_rule);
    parser.set_logger(logger);
    parser.enable_ast();
  }

  return parser;
}

#define SEQUANT_DESERIALIZATION_SPECIALIZATION(Type, Rule)        \
  template <>                                                     \
  Type from_string<Type>(std::string_view input,                  \
                         const DeserializationOptions &options) { \
    auto &parser = get_parser(Rule);                              \
                                                                  \
    std::shared_ptr<peg::Ast> ast;                                \
    parser.parse(input, ast);                                     \
                                                                  \
    SEQUANT_ASSERT(ast);                                          \
                                                                  \
    return ast_to<Type>(*ast, options);                           \
  }                                                               \
  template <>                                                     \
  Type from_string<Type>(std::wstring_view input,                 \
                         const DeserializationOptions &options) { \
    return v2::from_string<Type>(toUtf8(input), options);         \
  }

SEQUANT_DESERIALIZATION_SPECIALIZATION(Constant, "Constant")
SEQUANT_DESERIALIZATION_SPECIALIZATION(ExprPtr, "Expression")

}  // namespace sequant::io::serialization::v2
