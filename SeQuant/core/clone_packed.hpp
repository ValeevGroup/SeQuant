//
// Created by Bimal Gaudel on 9/22/21.
//

#ifndef  SEQUANT_CLONE_PACKED_HPP
#define  SEQUANT_CLONE_PACKED_HPP

#include "expr.hpp"
#include "tensor.hpp"

namespace sequant {

///
/// Clone an expression by preserving nested structures.
///
/// \param expr expression to be cloned
/// \return a cloned copy of \c expr
ExprPtr clone_packed(ExprPtr expr);

ExprPtr clone_packed(Sum const&);

ExprPtr clone_packed(Product const&);

ExprPtr clone_packed(Tensor const&);

} // namespace

#endif // SEQUANT_CLONE_PACKED_HPP
