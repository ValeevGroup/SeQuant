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

ccresponse::ccresponse(size_t n, size_t r) : N(n), R(r){};

}  // namespace sequant::mbpt::sr