//
// Created by Eduard Valeyev on 2019-07-29.
//

#ifndef SEQUANT_SRC_SEQUANT_CORE_EXPR_WFD_H
#define SEQUANT_SRC_SEQUANT_CORE_EXPR_WFD_H

#include <memory>

namespace sequant {

class Expr;
class ResultExpr;
class ExprPtr;

class Labeled;
class Constant;
using ConstantPtr = std::shared_ptr<Constant>;
class Product;
using ProductPtr = std::shared_ptr<Product>;
class CProduct;
using CProductPtr = std::shared_ptr<CProduct>;
class NCProduct;
using NCProductPtr = std::shared_ptr<NCProduct>;
class Sum;
using SumPtr = std::shared_ptr<Sum>;
class Variable;
using VariablePtr = std::shared_ptr<Variable>;
class Tensor;
using TensorPtr = std::shared_ptr<Tensor>;

}  // namespace sequant

#endif  // SEQUANT_SRC_SEQUANT_CORE_EXPR_WFD_H
