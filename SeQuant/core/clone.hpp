//
// Created by Bimal Gaudel on 9/22/21.
//

#ifndef  SEQUANT_CLONE_HPP
#define  SEQUANT_CLONE_HPP

#include "expr_fwd.hpp"

namespace sequant {

///
/// Clone an expression by preserving nested structures.
///
/// \param expr expression to be cloned
/// \return a cloned copy of \c expr
ExprPtr clone(ExprPtr expr);

} // namespace

#endif // SEQUANT_CLONE_HPP
