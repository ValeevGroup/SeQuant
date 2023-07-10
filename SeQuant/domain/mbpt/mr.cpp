//
// Created by Eduard Valeyev on 2019-02-19.
//

#include "SeQuant/domain/mbpt/mr.hpp"

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/math.hpp"
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/wick.hpp"
#include "SeQuant/domain/mbpt/formalism.hpp"

namespace sequant {
namespace mbpt {
namespace mr {

qninterval_t ncre(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied || s == IndexSpace::active ||
         s == IndexSpace::active_unoccupied);
  if (s == IndexSpace::active_occupied)
    return qns[0];
  else if (s == IndexSpace::active)
    return qns[2];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[4];
}

qninterval_t ncre(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  if (s == IndexSpace::active_occupied)
    return qns[0];
  else if (s == IndexSpace::active)
    return qns[2];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[4];
}

qninterval_t ncre_occ(qns_t qns) {
  return ncre(qns, IndexSpace::active_occupied);
}

qninterval_t ncre_act(qns_t qns) { return ncre(qns, IndexSpace::active); }

qninterval_t ncre_uocc(qns_t qns) {
  return ncre(qns, IndexSpace::active_unoccupied);
}

qninterval_t ncre(qns_t qns) { return qns[0] + qns[2] + qns[4]; }

qninterval_t nann(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied || s == IndexSpace::active ||
         s == IndexSpace::active_unoccupied);
  if (s == IndexSpace::active_occupied)
    return qns[1];
  else if (s == IndexSpace::active)
    return qns[3];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[5];
}

qninterval_t nann(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  if (s == IndexSpace::active_occupied)
    return qns[1];
  else if (s == IndexSpace::active)
    return qns[3];
  else  // if (s == IndexSpace::active_unoccupied)
    return qns[5];
}

qninterval_t nann_occ(qns_t qns) {
  return nann(qns, IndexSpace::active_occupied);
}

qninterval_t nann_act(qns_t qns) { return nann(qns, IndexSpace::active); }

qninterval_t nann_uocc(qns_t qns) {
  return nann(qns, IndexSpace::active_unoccupied);
}

qninterval_t nann(qns_t qns) { return qns[1] + qns[3] + qns[5]; }

qns_t combine(qns_t a, qns_t b) {
  const auto ncontr_uocc =
      qninterval_t{0, std::min(ncre(b, IndexSpace::active_unoccupied).upper(),
                               nann(a, IndexSpace::active_unoccupied).upper())};
  const auto ncontr_occ =
      qninterval_t{0, std::min(nann(b, IndexSpace::active_occupied).upper(),
                               ncre(a, IndexSpace::active_occupied).upper())};
  const auto nc_occ =
      nonnegative(ncre(a, IndexSpace::active_occupied) +
                  ncre(b, IndexSpace::active_occupied) - ncontr_occ);
  const auto nc_uocc =
      nonnegative(ncre(a, IndexSpace::active_unoccupied) +
                  ncre(b, IndexSpace::active_unoccupied) - ncontr_uocc);
  const auto na_occ =
      nonnegative(nann(a, IndexSpace::active_occupied) +
                  nann(b, IndexSpace::active_occupied) - ncontr_occ);
  const auto na_uocc =
      nonnegative(nann(a, IndexSpace::active_unoccupied) +
                  nann(b, IndexSpace::active_unoccupied) - ncontr_uocc);
  return qns_t{nc_occ, na_occ, nc_uocc, na_uocc};
}

}  // namespace mr
}  // namespace mbpt

mbpt::mr::qns_t adjoint(mbpt::mr::qns_t qns) {
  return mbpt::mr::qns_t{nann(qns, IndexSpace::active_occupied),
                         ncre(qns, IndexSpace::active_occupied),
                         nann(qns, IndexSpace::active),
                         ncre(qns, IndexSpace::active),
                         nann(qns, IndexSpace::active_unoccupied),
                         ncre(qns, IndexSpace::active_unoccupied)};
}

namespace mbpt {
namespace mr {

ExprPtr H1() {
  return get_default_context().vacuum() == Vacuum::Physical
             ? OpMaker<Statistics::FermiDirac>(
                   OpType::h, {IndexSpace::complete}, {IndexSpace::complete})()
             : OpMaker<Statistics::FermiDirac>(
                   OpType::f, {IndexSpace::complete}, {IndexSpace::complete})();
}

ExprPtr H2() {
  return OpMaker<Statistics::FermiDirac>(
      OpType::g, {IndexSpace::complete, IndexSpace::complete},
      {IndexSpace::complete, IndexSpace::complete})();
}

ExprPtr F() {
  return OpMaker<Statistics::FermiDirac>(OpType::f, {IndexSpace::complete},
                                         {IndexSpace::complete})();
}

ExprPtr H() { return H1() + H2(); }

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  FWickTheorem wick{expr};
  wick.spinfree(false).use_topology(use_top).set_nop_connections(
      nop_connections);
  auto result = wick.compute();
  simplify(result);
  if (Logger::get_instance().wick_stats) {
    std::wcout << "WickTheorem stats: # of contractions attempted = "
               << wick.stats().num_attempted_contractions
               << " # of useful contractions = "
               << wick.stats().num_useful_contractions << std::endl;
  }
  return result;
}

namespace op {

ExprPtr H1() {
  return ex<op_t>(
      [vacuum = get_default_context().vacuum()]() -> std::wstring_view {
        return vacuum == Vacuum::Physical ? L"h" : L"f";
      },
      [=]() -> ExprPtr { return mbpt::mr::H1(); },
      [=](qnc_t& qns) {
        qns =
            combine(qnc_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}}, qns);
      });
}

ExprPtr H2() {
  return ex<op_t>(
      []() -> std::wstring_view { return L"g"; },
      [=]() -> ExprPtr { return mbpt::mr::H2(); },
      [=](qnc_t& qns) {
        qns =
            combine(qnc_t{{0, 2}, {0, 2}, {0, 2}, {0, 2}, {0, 2}, {0, 2}}, qns);
      });
}

ExprPtr H() { return H1() + H2(); }

bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t target_qns,
                    const qns_t source_qns = {}) {
  qns_t qns = source_qns;
  if (op_or_op_product.is<Product>()) {
    const auto& op_product = op_or_op_product.as<Product>();
    for (auto& op_ptr : ranges::views::reverse(op_product.factors())) {
      assert(op_ptr->template is<op_t>());
      const auto& op = op_ptr->template as<op_t>();
      qns = op(qns);
    }
    return qns.overlaps_with(target_qns);
  } else if (op_or_op_product.is<op_t>()) {
    const auto& op = op_or_op_product.as<op_t>();
    qns = op();
    return qns.overlaps_with(target_qns);
  } else
    throw std::invalid_argument(
        "sequant::mbpt::sr::contains_rank(op_or_op_product): op_or_op_product "
        "must be mbpt::sr::op_t or Product thereof");
}

bool raises_vacuum_up_to_rank(const ExprPtr& op_or_op_product,
                              const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product,
                        qns_t{{0ul, 0ul}, {0ul, k}, {0ul, k}, {0ul, 0ul}});
}

bool lowers_rank_or_lower_to_vacuum(const ExprPtr& op_or_op_product,
                                    const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{},
                        qns_t{{0ul, 0ul}, {0ul, k}, {0ul, k}, {0ul, 0ul}});
}

bool raises_vacuum_to_rank(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{0ul, k, k, 0ul});
}

bool lowers_rank_to_vacuum(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{}, qns_t{0ul, k, k, 0ul});
}

ExprPtr vac_av(
    ExprPtr expr,
    std::vector<std::pair<std::wstring, std::wstring>> op_connections) {
  // use cloned expr to avoid side effects
  expr = expr->clone();

  auto vac_av_product = [&op_connections](ExprPtr expr) {
    assert(expr.is<Product>());
    // compute connections
    std::vector<std::pair<int, int>> connections;
    {
      std::map<std::wstring, std::vector<int>>
          oplbl2pos;  // maps operator labels to the operator positions in the
                      // product
      int pos = 0;
      for (const auto& factor : expr.as<Product>()) {
        if (factor.is<op_t>()) {
          const auto& op = factor.as<op_t>();
          const std::wstring op_lbl = std::wstring(op.label());
          const auto it = oplbl2pos.find(op_lbl);
          if (it == oplbl2pos.end()) {  // new label
            oplbl2pos.emplace(op_lbl, std::vector<int>{pos});
          } else {
            it->second.emplace_back(pos);
          }
          ++pos;
        } else if (factor.is<FNOperator>() || factor.is<BNOperator>()) {
          ++pos;  // skip FNOperator and BNOperator
        }
      }

      for (const auto& [op1_lbl, op2_lbl] : op_connections) {
        auto it1 = oplbl2pos.find(op1_lbl);
        auto it2 = oplbl2pos.find(op2_lbl);
        if (it1 == oplbl2pos.end() || it2 == oplbl2pos.end())
          continue;  // one of the op labels is not present in the product
        const auto& [dummy1, op1_indices] = *it1;
        const auto& [dummy2, op2_indices] = *it2;
        for (const auto& op1_idx : op1_indices) {
          for (const auto& op2_idx : op2_indices) {
            using std::min;
            using std::max;
            connections.emplace_back(min(op1_idx, op2_idx),
                                     max(op1_idx, op2_idx));
          }
        }
      }
    }

    // lower to tensor form
    auto lower_to_tensor_form = [](ExprPtr& expr) {
      auto op_lowerer = [](ExprPtr& leaf) {
        if (leaf.is<op_t>()) leaf = leaf.as<op_t>().tensor_form();
      };
      expr->visit(op_lowerer, /* atoms only = */ true);
    };
    lower_to_tensor_form(expr);
    expr = simplify(expr);

    // compute VEVs
    return mbpt::mr::vac_av(expr, connections, /* use_topology = */ true);
  };

  ExprPtr result;
  if (expr.is<Product>()) {
    return vac_av_product(expr);
  } else if (expr.is<Sum>()) {
    result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [&op_connections](const auto& op_product) {
          return vac_av(op_product, op_connections);
        });
    return result;
  } else if (expr.is<op_t>()) {
    return ex<Constant>(
        0);  // expectation value of a normal-ordered operator is 0
  } else if (expr.is<Constant>()) {
    return expr;  // vacuum is normalized
  }
  throw std::invalid_argument(
      "mpbt::mr::op::vac_av(expr): unknown expression type");
}

}  // namespace op

}  // namespace mr

// must be defined including op.ipp since it's used there
template <>
bool is_vacuum<mr::qns_t>(mr::qns_t qns) {
  return qns == mr::qns_t{};
}

}  // namespace mbpt

template <Statistics S>
std::wstring to_latex(const mbpt::Operator<mbpt::mr::qns_t, S>& op) {
  using namespace sequant::mbpt;
  using namespace sequant::mbpt::mr;

  auto lbl = std::wstring(op.label());
  std::wstring result = L"{\\hat{" + lbl + L"}";
  auto it = label2optype.find(lbl);
  OpType optype = OpType::invalid;
  if (it != label2optype.end()) {  // handle special cases
    optype = it->second;
    if (optype == OpType::lambda) {  // λ -> \lambda
      result = L"{\\hat{\\lambda}";
    }
    if (to_class(optype) == OpClass::gen) {
      result += L"}";
      return result;
    }
  }

  // generic operator ... can only handle definite case
  const auto dN = op();
  if (!is_definite(ncre_occ(dN)) || !is_definite(nann_occ(dN)) ||
      !is_definite(ncre_uocc(dN)) || !is_definite(nann_uocc(dN))) {
    throw std::invalid_argument(
        "to_latex(const Operator<qns_t, S>& op): "
        "can only handle  generic operators with definite cre/ann numbers");
  }
  // pure quasiparticle creator/annihilator?
  const auto qprank_cre = nann_occ(dN).lower() + ncre_uocc(dN).lower();
  const auto qprank_ann = ncre_occ(dN).lower() + nann_uocc(dN).lower();
  const auto qppure = qprank_cre == 0 || qprank_ann == 0;
  auto qpaction = to_class(optype);
  if (qppure) {
    if (qprank_cre) {
      // if operator's action implied by the label and actual action agrees, use
      // subscript always
      std::wstring baseline_char = (qpaction != OpClass::deex ? L"_" : L"^");
      if (nann_occ(dN).lower() == ncre_uocc(dN).lower())
        result +=
            baseline_char + L"{" + std::to_wstring(nann_occ(dN).lower()) + L"}";
      else
        result += baseline_char + L"{" + std::to_wstring(nann_occ(dN).lower()) +
                  L"," + std::to_wstring(ncre_uocc(dN).lower()) + L"}";
    } else {
      // if operator's action implied by the label and actual action agrees, use
      // subscript always
      std::wstring baseline_char = (qpaction != OpClass::deex ? L"^" : L"_");
      if (nann_uocc(dN).lower() == ncre_occ(dN).lower()) {
        result +=
            baseline_char + L"{" + std::to_wstring(ncre_occ(dN).lower()) + L"}";
      } else
        result += baseline_char + L"{" + std::to_wstring(ncre_occ(dN).lower()) +
                  L"," + std::to_wstring(nann_uocc(dN).lower()) + L"}";
    }
  } else {  // not pure qp creator/annihilator
    result += L"_{" + std::to_wstring(nann_occ(dN).lower()) + L"," +
              std::to_wstring(ncre_uocc(dN).lower()) + L"}^{" +
              std::to_wstring(ncre_occ(dN).lower()) + L"," +
              std::to_wstring(nann_uocc(dN).lower()) + L"}";
  }
  result += L"}";
  return result;
}

}  // namespace sequant

#include "SeQuant/domain/mbpt/op.ipp"

namespace sequant {
namespace mbpt {
template class Operator<mr::qns_t, Statistics::FermiDirac>;
template class Operator<mr::qns_t, Statistics::BoseEinstein>;
}  // namespace mbpt
}  // namespace sequant
