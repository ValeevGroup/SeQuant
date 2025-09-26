//
// Created by Ajay Melekamburath on 9/26/25.
//

#ifndef SEQUANT_DOMAIN_MBPT_UTILS_HPP
#define SEQUANT_DOMAIN_MBPT_UTILS_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>

#include <SeQuant/domain/mbpt/op.hpp>

namespace sequant::mbpt {

/// @brief Computes the similarity transformation e^(-B) A e^B using the
/// commutator expansion up to the specified order.
/// The expansion is given by: e^(-B) A e^B = A + [A,B] + (1/2!)[[A,B],B] +
/// (1/3!)[[[A,B],B],B] + ...
///
/// @param A Operator A
/// @param B Operator B
/// @param commutator_rank The order at which to truncate the expansion (number
/// of nested commutators)
/// @param unitary If true, uses unitary ansatz with B-B^+
/// @pre This function expects \p A and \p B to be composed of mbpt::Operators
inline ExprPtr sim_tr(const ExprPtr& A, const ExprPtr& B,
                      size_t commutator_rank, bool unitary = false) {
  assert(commutator_rank >= 1 && "Truncation order must be at least 1");

  auto expr = A.clone();  // work on a copy

  // takes a product or an operator and applies the similarity transformation
  auto transform = [&B, commutator_rank, unitary](const ExprPtr& e) {
    assert(e.is<op_t>() || e.is<Product>());
    auto result = e;  // start with expr
    auto op_Sk = result;

    for (size_t k = 1; k <= commutator_rank; ++k) {
      ExprPtr op_Sk_comm_w_S;
      op_Sk_comm_w_S = op_Sk * B;  // traditional ansatz: [O,B] = (O B)_c

      if (unitary)  // unitary ansatz: [O,B-B^+] = (O B)_c + (B^+ O)_c
        op_Sk_comm_w_S += adjoint(B) * op_Sk;

      op_Sk = ex<Constant>(rational{1, k}) * op_Sk_comm_w_S;
      simplify(op_Sk);
      result += op_Sk;
    }
    return result;
  };

  // expression type dependent dispatch
  if (expr.is<op_t>()) {
    return transform(expr);
  } else if (expr.is<Product>()) {
    auto& product = expr.as<Product>();
    // Expand product as sum
    if (ranges::any_of(product.factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = sequant::expand(expr);
      simplify(expr);
      return sim_tr(expr, B, commutator_rank, unitary);
    } else {
      return transform(expr);
    }
  } else if (expr.is<Sum>()) {
    auto result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [=](const auto& op_product) { return transform(op_product); });
    return result;
  } else if (expr.is<Constant>() || expr.is<Variable>())
    return expr;
  else
    throw std::invalid_argument(
        "mbpt::sim_tr(expr): Unsupported expression type");
}

}  // namespace sequant::mbpt

#endif  // SEQUANT_DOMAIN_MBPT_UTILS_HPP
