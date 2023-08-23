#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>
#include <iostream>

#include <SeQuant/core/math.hpp>

#include <SeQuant/core/op.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/domain/mbpt/sr.hpp>

namespace sequant::mbpt::sr {

cceqs::cceqs(size_t n, size_t p, size_t pmin)
    : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}

ExprPtr cceqs::sim_transform(ExprPtr expr, size_t r) {
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
      return sim_transform(expr, r);
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
        "cceqs::sim_transform(expr): Unsupported expression type");
}

std::vector<ExprPtr> cceqs::t(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // 1. construct hbar(op) in canonical form
  auto hbar = sim_transform(op::H(), 4);

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

    // 2.b project onto <p|, i.e. multiply by P(p)
    auto P_hbar = simplify(op::P(p) * hbar_p);

    // 2.c compute
    result.at(p) = op::vac_av(P_hbar);
    simplify(result.at(p));
  }

  return result;
}

std::vector<ExprPtr> cceqs::λ(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // construct hbar
  auto hbar = sim_transform(op::H(), 3);

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

    // 2.b multiply by adjoint of P(p) (i.e., P(-p)) on the right side
    auto hbar_P = simplify(hbar_p * op::P(-p));

    // 2.c compute vacuum average
    result.at(p) = op::vac_av(hbar_P, op_connect);
    simplify(result.at(p));
  }
  return result;
}

auto make_pert_tnsr = [](const std::wstring& label, const size_t n) {
  std::wstring bra_annot, ket_annot;
  for (size_t i = 1; i <= n; ++i) {
    bra_annot += L"a_" + std::to_wstring(i) + L",";
    ket_annot += L"i_" + std::to_wstring(i) + L",";
  }
  bra_annot.pop_back();
  ket_annot.pop_back();
  auto result = label + L"{" + bra_annot + L";" + ket_annot + L"}";
  // TODO: Don't hard-code antisymm
  return sequant::parse_expr(result, Symmetry::antisymm);
};

std::vector<sequant::ExprPtr> cceqs::pert_t1() {
  using namespace sequant::mbpt;

  // construct (V * e^T)_c = V + V * T + V * T^2/2!
  auto mu_bar = op::mu(1);
  auto mu_Tk = mu_bar;
  for (int64_t k = 1; k <= 2; ++k) {
    mu_Tk = simplify(ex<Constant>(rational{1, k}) * mu_Tk * op::T(N));
    mu_bar += mu_Tk;
  }

  // construct (H * e^T * pT1)_c = H * pT1 + H * pT1 * T + H * pT1 * T^2/2! + H
  // * pT1 * T^3/3!
  auto hbar_pert = (op::H() * op::pertT1(N));
  auto H_Tk_pert = hbar_pert;
  for (int64_t k = 1; k <= 3; ++k) {
    H_Tk_pert = simplify(ex<Constant>(rational{1, k}) * H_Tk_pert * op::T(N));
    hbar_pert += H_Tk_pert;
  }

  std::vector<std::pair<std::wstring, std::wstring>> op_connect = {
      {L"h", L"t"},  {L"f", L"t"},  {L"g", L"t"}, {L"h", L"t¹"},
      {L"f", L"t¹"}, {L"g", L"t¹"}, {L"μ", L"t"}};

  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    auto eq = simplify(op::P(p) * (mu_bar + hbar_pert));

    result.at(p) = op::vac_av(eq, op_connect);
    simplify(result.at(p));
  }
  // add frequency scaled terms
  auto omega = ex<Variable>(L"ω");
  for (auto i = 1; i <= N; ++i) {
    ExprPtr t_tnsr = make_pert_tnsr(optype2label.at(OpType::t_1), i);
    result[i] -= omega * t_tnsr;
  }

  return result;
}

std::vector<ExprPtr> cceqs::pert_λ1() {
  /// Eqn: <0| Λ^{(1)}(H * e^T * \tau_{mu}) |0> + <0| (1 + Λ^{(0)})( H * e^T *
  /// e^pT1 * \tau_{mu} ) |0>

  // construct unperturbed H_bar (reusable code)
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }

  // construct V_bar (including perturbed T amplitudes also)
  auto mu_bar = op::mu(1);
  auto mu_Tk = mu_bar;
  for (int64_t k = 1; k <= 2; ++k) {
    mu_Tk = simplify(ex<Constant>(rational{1, k}) * mu_Tk * op::T(N) *
                     op::pertT1(N));
    mu_bar += mu_Tk;
  }

  // first term
  auto eq1 =
      simplify(op::pertLambda1(N) * hbar); /* times excitation operator */

  // second term
  const auto One = ex<Constant>(1);
  auto eq2 =
      simplify((One + op::Lambda(N)) * mu_bar); /* times excitation operator */

  std::vector<std::pair<std::wstring, std::wstring>> op_connect = {
      {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"}, {L"μ", L"t"},  {L"μ", L"t¹"},
      {L"h", L"A"}, {L"f", L"A"}, {L"g", L"A"}, {L"h", L"S"},  {L"f", L"S"},
      {L"g", L"S"}, {L"t", L"A"}, {L"t", L"S"}, {L"t¹", L"A"}, {L"t¹", L"S"}};

  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    auto eq = simplify((eq1 + eq2) * op::P(-p));

    result.at(p) = op::vac_av(eq, op_connect);
    simplify(result.at(p));
  }
  // add frequency scaled terms
  auto omega = ex<Variable>(L"ω");
  for (auto i = 1; i <= N; ++i) {
    ExprPtr λ_tnsr = make_pert_tnsr(optype2label.at(OpType::λ_1), i);
    result[i] += omega * λ_tnsr;
  }

  return result;
}

}  // namespace sequant::mbpt::sr
