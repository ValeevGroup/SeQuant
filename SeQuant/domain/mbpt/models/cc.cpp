#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>

#include <clocale>
#include <iostream>

#include <SeQuant/core/math.hpp>

#include <SeQuant/core/op.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/domain/mbpt/sr.hpp>

namespace sequant::mbpt::sr {

cceqs::cceqs(size_t n, size_t p, size_t pmin)
    : N(n), P(p == std::numeric_limits<size_t>::max() ? n : p), PMIN(pmin) {}

std::vector<ExprPtr> cceqs::t(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // 1. construct hbar(op) in canonical form
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }

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

    // 2.b multiply by A(P)
    auto A_hbar = simplify(op::A(p) * hbar_p);

    // 2.c compute vacuum average
    result.at(p) = op::vac_av(A_hbar);
    simplify(result.at(p));
  }

  return result;
}

std::vector<ExprPtr> cceqs::λ(bool screen, bool use_topology,
                              bool use_connectivity, bool canonical_only) {
  // construct hbar
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 3; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }

  const auto One = ex<Constant>(1);
  auto lhbar = simplify((One + op::Lambda(N)) * hbar);

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
    lhbar = hbar_le_p;  // not needed

    // 2.b multiply by adjoint of A(P) on the right side

    auto A_hbar = simplify(hbar_p * adjoint(op::A(p)));

    // temp
    std::vector<std::pair<std::wstring, std::wstring>> new_op_connect = {
        {L"h", L"t"}, {L"f", L"t"}, {L"g", L"t"},
        {L"h", L"A"}, {L"f", L"A"}, {L"g", L"A"}};

    // 2.c compute vacuum average
    result.at(p) = op::vac_av(A_hbar, new_op_connect);
    simplify(result.at(p));
  }
  return result;
}

std::vector<sequant::ExprPtr> cceqs::pert_t1() {
  using namespace sequant::mbpt;

  // construct unperturbed H_bar (reusable code)
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }

  // construct (V * e^T)_c = V + V * T + V * T^2/2!
  auto Vbar = op::V(1);
  auto V_Tk = Vbar;
  for (int64_t k = 1; k <= 2; ++k) {
    V_Tk = simplify(ex<Constant>(rational{1, k}) * V_Tk * op::T(N));
    Vbar += V_Tk;
  }

  // construct (H * e^T * pT1)_c = H * pT1 + H * pT1 * T + H * pT1 * T^2/2! + H
  // * pT1 * T^3/3!
  auto hbar_pert = (op::H() * op::pertT1(N));
  auto H_Tk_pert = hbar_pert;
  for (int64_t k = 1; k <= 3; ++k) {
    H_Tk_pert = simplify(ex<Constant>(rational{1, k}) * H_Tk_pert * op::T(N));
    hbar_pert += H_Tk_pert;
  }

  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p >= PMIN; --p) {
    auto eq = simplify(op::A(p) * (Vbar + hbar_pert));

    std::vector<std::pair<std::wstring, std::wstring>> new_op_connect = {
        {L"h", L"t"},  {L"f", L"t"},  {L"g", L"t"}, {L"h", L"t¹"},
        {L"f", L"t¹"}, {L"g", L"t¹"}, {L"V", L"t"}};
    result.at(p) = op::vac_av(eq, new_op_connect);
    simplify(result.at(p));
  }

  return result;
}

 std::vector<ExprPtr> cceqs::pert_λ1() {
  abort(); // not implemented yet
 }

}  // namespace sequant::mbpt::sr
