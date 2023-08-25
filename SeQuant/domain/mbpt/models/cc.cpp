#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>
#include <iostream>

#include <SeQuant/core/math.hpp>

#include <SeQuant/core/op.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/domain/mbpt/sr.hpp>

namespace sequant::mbpt::sr {

cceqs::cceqs(size_t n, size_t p, size_t pmin)
    : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}

ExprPtr cceqs::sim_tr(ExprPtr expr, size_t r) {
  auto transform_op_op_pdt = [this, &r](const ExprPtr& expr) {
    // TODO: find the order at which the commutator expression should truncate
    // from op/op product
    assert(expr.is<op_t>() || expr.is<Product>());
    auto result = expr;
    auto op_Tk = result;
    for (int64_t k = 1; k <= r; ++k) {
      op_Tk = simplify(ex<Constant>(rational{1, k}) * op_Tk * op::T(N));
      result += op_Tk;
    }
    return result;
  };

  if (expr.is<op_t>()) {
    return transform_op_op_pdt(expr);
  } else if (expr.is<Product>()) {
    auto& product = expr.as<Product>();
    // Expand product as sum
    if (ranges::any_of(product.factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = sequant::expand(expr);
      simplify(expr);
      return sim_tr(expr, r);
    } else {
      return transform_op_op_pdt(expr);
    }
  } else if (expr.is<Sum>()) {
    auto result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [=](const auto& op_product) {
          return transform_op_op_pdt(op_product);
        });
    return result;
  } else if (expr.is<Constant>() || expr.is<Variable>())
    return expr;
  else
    throw std::invalid_argument(
        "cceqs::sim_tr(expr): Unsupported expression type");
}

std::vector<ExprPtr> cceqs::t(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // 1. construct hbar(op) in canonical form
  auto hbar = sim_tr(op::H(), 4);

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <p|
    std::shared_ptr<Sum>
        hbar_p;  // products that can produce excitations of rank p
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p
    for (auto& term : *hbar) {
      assert(term->is<Product>() || term->is<op_t>());

      if (op::raises_vacuum_up_to_rank(term, p)) {
        if (!hbar_le_p)
          hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
        else
          hbar_le_p->append(term);
        if (op::raises_vacuum_to_rank(term, p)) {
          if (!hbar_p)
            hbar_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_p->append(term);
        }
      }
    }
    hbar = hbar_le_p;

    // 2.b project onto <p|, i.e. multiply by P(p) and evaluate
    result.at(p) = op::op_evaluate(op::P(p) * hbar_p);
  }

  return result;
}

std::vector<ExprPtr> cceqs::λ(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // construct hbar
  auto hbar = sim_tr(op::H(), 3);

  const auto One = ex<Constant>(1);
  auto lhbar = simplify((One + op::Lambda(N)) * hbar);

  std::vector<std::pair<std::wstring, std::wstring>> op_connect = {
      {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"}, {L"h", L"A"}, {L"f", L"A"},
      {L"g", L"A"}, {L"h", L"S"}, {L"f", L"S"}, {L"g", L"S"}};

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <P|
    std::shared_ptr<Sum>
        hbar_p;  // products that can produce excitations of rank p
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p
    for (auto& term : *lhbar) {  // pick terms from lhbar
      assert(term->is<Product>() || term->is<op_t>());

      if (op::lowers_rank_or_lower_to_vacuum(term, p)) {
        if (!hbar_le_p)
          hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
        else
          hbar_le_p->append(term);
        if (op::lowers_rank_to_vacuum(term, p)) {
          if (!hbar_p)
            hbar_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_p->append(term);
        }
      }
    }
    lhbar = hbar_le_p;

    // 2.b multiply by adjoint of P(p) (i.e., P(-p)) on the right side and
    // evaluate
    result.at(p) = op::op_evaluate(hbar_p * op::P(-p), op_connect);
  }
  return result;
}

}  // namespace sequant::mbpt::sr
