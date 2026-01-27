//
// Created by Robert Adam on 2023-09-21
//

#include <SeQuant/core/parse/v1/ast.hpp>

namespace sequant::parse::v1::ast {

Product::Product(std::vector<NullaryValue> factors)
    : factors(std::move(factors)) {}

Sum::Sum(std::vector<Product> summands) : summands(std::move(summands)) {}

}  // namespace sequant::parse::v1::ast
