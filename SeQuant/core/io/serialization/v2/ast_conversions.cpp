#include <SeQuant/core/complex.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/serialization/v2/ast_conversions.hpp>
#include <SeQuant/core/utility/conversion.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <peglib.h>

#include <cstdint>

namespace sequant::io::serialization::v2 {

std::int64_t to_int(const peg::Ast &ast) {
  SEQUANT_ASSERT(ast.name == "Integer");

  return string_to<std::int64_t>(ast.token);
}

double to_double(const peg::Ast &ast) {
  SEQUANT_ASSERT(ast.name == "Float");

  return string_to<double>(ast.token);
}

double to_real(const peg::Ast &ast) {
  SEQUANT_ASSERT(ast.name == "Real");
  SEQUANT_ASSERT(ast.nodes.size() == 1);
  SEQUANT_ASSERT(ast.nodes[0]);

  if (ast.nodes[0]->name == "Integer") {
    return to_int(*ast.nodes[0]);
  } else {
    return to_double(*ast.nodes[0]);
  }
}

Constant::scalar_type to_complex(const peg::Ast &ast) {
  SEQUANT_ASSERT(ast.name == "Complex");
  SEQUANT_ASSERT(ast.nodes.size() == 1);
  SEQUANT_ASSERT(ast.nodes[0]);

  if (ast.nodes[0]->name == "Imaginary") {
    const auto &imag_node = *ast.nodes[0];
    SEQUANT_ASSERT(ast.nodes.size() == 1);
    SEQUANT_ASSERT(ast.nodes[0]);
    return Constant::scalar_type(0, to_real(*imag_node.nodes[0]));
  } else {
    return Constant::scalar_type(to_real(*ast.nodes[0]), 0);
  }
}

template <>
Constant ast_to<Constant>(const peg::Ast &ast, const DeserializationOptions &) {
  SEQUANT_ASSERT(ast.name == "Constant");
  SEQUANT_ASSERT(ast.nodes.size() == 1);
  SEQUANT_ASSERT(ast.nodes[0]);

  return Constant(to_complex(*ast.nodes[0]));
}

template <>
ExprPtr ast_to<ExprPtr>(const peg::Ast &ast,
                        const DeserializationOptions &options) {
  (void)ast;
  (void)options;
  return {};
}

}  // namespace sequant::io::serialization::v2
