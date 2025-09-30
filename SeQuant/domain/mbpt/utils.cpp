//
// Created by Ajay Melekamburath on 9/30/25.
//

#include <SeQuant/domain/mbpt/utils.hpp>

namespace sequant::mbpt {

ExprPtr sim_tr(ExprPtr A, const ExprPtr& B, size_t commutator_rank,
               bool unitary, bool skip_clone) {
  assert(commutator_rank >= 1 && "Truncation order must be at least 1");

  // use cloned expr to avoid side effects
  if (!skip_clone) A = A->clone();

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

  // expression type dispatch
  if (A.is<op_t>()) {
    return transform(A);
  } else if (A.is<Product>()) {
    auto& product = A.as<Product>();
    // Expand product as sum
    if (ranges::any_of(product.factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      A = sequant::expand(A);
      simplify(A);
      return sim_tr(A, B, commutator_rank, unitary, /*skip_clone*/ true);
    } else {
      return transform(A);
    }
  } else if (A.is<Sum>()) {
    auto result = sequant::transform_reduce(
        *A, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [=](const auto& term) {
          return sim_tr(term, B, commutator_rank, unitary, /*skip_clone*/ true);
        });
    return result;
  } else if (A.is<Constant>() || A.is<Variable>())
    return A;
  else
    throw std::invalid_argument(
        "mbpt::sim_tr(A, B, commutator_rank, unitary): Unsupported expression "
        "type");
}

ExprPtr screen_terms(ExprPtr expr, bool skip_clone) {
  // use cloned expr to avoid side effects
  if (!skip_clone) expr = expr->clone();

  auto screen = [](const ExprPtr& term) {
    if (!(term->is<op_t>() || term->is<Product>())) {
      throw std::invalid_argument("op::screen_terms: Unsupported term type");
    }
    qns_t term_qns;
    if (term->is<op_t>()) {
      term_qns = term->as<op_t>()();
    } else if (term->is<Product>()) {
      term_qns = mbpt::detail::compute_qnc(term->as<Product>());
    }
    return is_vacuum(term_qns) ? term : ex<Constant>(0);
  };

  // expression type dispatch
  if (expr.is<op_t>()) {
    return screen(expr);
  } else if (expr.is<Product>()) {
    auto& product = expr.as<Product>();
    // Expand product as sum
    if (ranges::any_of(product.factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = sequant::expand(expr);
      simplify(expr);
      return screen_terms(expr, /*skip_clone*/ true);
    } else {
      return screen(expr);
    }
  } else if (expr.is<Variable>() || expr.is<Constant>()) {
    return expr;
  } else if (expr.is<Sum>()) {
    auto result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [=](const auto& term) {
          return screen_terms(term, /*skip_clone*/ true);
        });
    return result;
  } else
    throw std::invalid_argument(
        "mbpt::screen_terms(expr): Unsupported expression type");
}

}  // namespace sequant::mbpt
