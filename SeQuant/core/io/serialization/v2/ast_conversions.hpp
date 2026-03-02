#ifndef SEQUANT_CORE_IO_SERIALIZATION_V2_AST_CONVERSIONS_HPP
#define SEQUANT_CORE_IO_SERIALIZATION_V2_AST_CONVERSIONS_HPP

#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>

#include <peglib.h>

namespace sequant::io::serialization::v2 {

template <typename T>
T ast_to(const peg::Ast &ast, const DeserializationOptions &options) = delete;

template <>
Constant ast_to<Constant>(const peg::Ast &ast,
                          const DeserializationOptions &options);

template <>
ExprPtr ast_to<ExprPtr>(const peg::Ast &ast,
                        const DeserializationOptions &options);

}  // namespace sequant::io::serialization::v2

#endif  // SEQUANT_CORE_IO_SERIALIZATION_V2_AST_CONVERSIONS_HPP
