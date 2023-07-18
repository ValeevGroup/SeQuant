//
// Created by Ajay Melekamburath on 7/12/23.
//

#include <SeQuant/domain/mbpt/models/ccresponse.hpp>

namespace sequant::mbpt::sr {
// response/perturbation related operators

ExprPtr V() { return make_op(OpType::V, 1)(); }

ExprPtr pT1_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return make_op(OpType::t, Nbra, Nket)();
}

ExprPtr pLambda1_(std::size_t Nbra, std::size_t Nket) {
  assert(Nbra > 0);
  assert(Nket > 0);
  return make_op(OpType::lambda, Nbra, Nket)();
}

namespace op {

ExprPtr V() {
  return ex<op_t>([]() -> std::wstring_view { return L"V"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return sr::V();
                  },
                  [=](qnc_t& qns) {
                    qns = combine(qnc_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}}, qns);
                  });
}

ExprPtr pT1_(std::size_t K) {
  return ex<op_t>([]() -> std::wstring_view { return L"t1"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return sr::pT1_(K);
                  },
                  [=](qnc_t& qns) {
                    qns = combine(qnc_t{0ul, K, K, 0ul}, qns);
                  });
}

ExprPtr pT1(std::size_t K) {
  assert(K > 0);

  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + pT1_(k) : pT1_(k);
  }
  return result;
}

ExprPtr pLambda1_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"λ1"; },
                  [=]() -> ExprPtr {
                    using namespace sequant::mbpt::sr;
                    return sr::pLambda1_(K);
                  },
                  [=](qnc_t& qns) {
                    qns = combine(qnc_t{K, 0ul, 0ul, K}, qns);
                  });
}

ExprPtr pLambda1(std::size_t K) {
  assert(K > 0);

  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + pLambda1_(k) : pLambda1_(k);
  }
  return result;
}

}  // namespace op

ccresponse::ccresponse(size_t n, size_t p, size_t pmin, size_t r)
    : N(n),
      P(p == std::numeric_limits<size_t>::max() ? n : p),
      PMIN(pmin),
      R(r) {}

std::vector<sequant::ExprPtr> ccresponse::t() {
  using namespace sequant::mbpt;
  // Equation: <0|A (e^T V)_c|0> + <0|A (H_bar^(0) e^T^(1) V)_c|0>

  // construct unperturbed H_bar (reusable code)
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }

  auto V = op::V();
  std::vector<ExprPtr> result (P + 1);
  for (auto p = P; p <= PMIN; p++){
    // try without screening
    auto eq = simplify((op::A(p) * op::T(N) * V) +
                       (op::A(p) * hbar * op::T(N) * op::pT1(N)));

    // mention connectivity here
    result.at(p) = op::vac_av(eq);
  }
  return result;
}

std::vector<ExprPtr> ccresponse::lambda() {
  using namespace sequant::mbpt;
  // Equation: <0| (1 + Lambda^(0) (H_doublebar(1) *  e^T)_c * A|0> +
  // <0| (Lambda^(1) (H_bar(0) *  e^T)_c * A|0>

  // construct unperturbed H_bar (reusable code)
  auto hbar = op::H();
  auto H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T(N));
    hbar += H_Tk;
  }
  auto V = op::V();
  const auto One = ex<Constant>(1);
  std::vector<ExprPtr> result(P + 1);
  for (auto p = P; p <= PMIN; p++) {
    auto eq =
        simplify((One + op::Lambda(N)) * (op::pT1_(N) * (hbar + V) * A(p)) +
                 (op::pLambda1_(N) * hbar * op::T(N)) * A(p));

    // mention connectivity here
    result.at(p) = op::vac_av(eq);
  }
  return result;
}
}  // namespace sequant::mbpt::sr