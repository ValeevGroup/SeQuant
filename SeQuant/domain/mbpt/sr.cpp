//
// Created by Eduard Valeyev on 2019-02-19.
//

#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>

#include <range/v3/iterator/basic_iterator.hpp>
#include <range/v3/iterator/reverse_iterator.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/view.hpp>

namespace sequant {
namespace mbpt {
namespace sr {

qninterval_t ncre(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied ||
         s == IndexSpace::active_unoccupied);
  return s == IndexSpace::active_occupied ? qns[0] : qns[2];
}

qninterval_t ncre(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  return s.type() == IndexSpace::active_occupied ? qns[0] : qns[2];
}

qninterval_t ncre_occ(qns_t qns) {
  return ncre(qns, IndexSpace::active_occupied);
}

qninterval_t ncre_uocc(qns_t qns) {
  return ncre(qns, IndexSpace::active_unoccupied);
}

qninterval_t ncre(qns_t qns) { return qns[0] + qns[2]; }

qninterval_t nann(qns_t qns, const IndexSpace::Type& s) {
  assert(s == IndexSpace::active_occupied ||
         s == IndexSpace::active_unoccupied);
  return s == IndexSpace::active_occupied ? qns[1] : qns[3];
}

qninterval_t nann(qns_t qns, const IndexSpace& s) {
  assert((s.type() == IndexSpace::active_occupied ||
          s.type() == IndexSpace::active_unoccupied) &&
         s.qns() == IndexSpace::nullqns);
  return s.type() == IndexSpace::active_occupied ? qns[1] : qns[3];
}

qninterval_t nann_occ(qns_t qns) {
  return nann(qns, IndexSpace::active_occupied);
}

qninterval_t nann_uocc(qns_t qns) {
  return nann(qns, IndexSpace::active_unoccupied);
}

qninterval_t nann(qns_t qns) { return qns[1] + qns[3]; }

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

}  // namespace sr
}  // namespace mbpt

mbpt::sr::qns_t adjoint(mbpt::sr::qns_t qns) {
  return mbpt::sr::qns_t{nann(qns, IndexSpace::active_occupied),
                         ncre(qns, IndexSpace::active_occupied),
                         nann(qns, IndexSpace::active_unoccupied),
                         ncre(qns, IndexSpace::active_unoccupied)};
}

namespace mbpt {
namespace sr {

OpMaker::OpMaker(OpType op, std::size_t nbra, std::size_t nket)
    : base_type(op) {
  nket = nket == std::numeric_limits<std::size_t>::max() ? nbra : nket;
  assert(nbra > 0 || nket > 0);

  const auto unocc = IndexSpace::active_unoccupied;
  const auto occ = IndexSpace::active_occupied;
  switch (to_class(op)) {
    case OpClass::ex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, unocc);
      ket_spaces_ = decltype(ket_spaces_)(nket, occ);
      break;
    case OpClass::deex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, occ);
      ket_spaces_ = decltype(ket_spaces_)(nket, unocc);
      break;
    case OpClass::gen:
      bra_spaces_ = decltype(bra_spaces_)(nbra, IndexSpace::complete);
      ket_spaces_ = decltype(ket_spaces_)(nket, IndexSpace::complete);
      break;
  }
}

#include <SeQuant/domain/mbpt/sr/op.impl.cpp>

ExprPtr H0mp() { return F(); }

ExprPtr H1mp() {
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  return H_(2);
}

ExprPtr W() {
  assert(get_default_context().vacuum() == Vacuum::SingleProduct);
  return H1mp();
}

ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);
  switch (k) {
    case 1:
      switch (get_default_context().vacuum()) {
        case Vacuum::Physical:
          return OpMaker(OpType::h, 1)();
        case Vacuum::SingleProduct:
          return OpMaker(OpType::f, 1)();
        case Vacuum::MultiProduct:
          abort();
        default:
          abort();
      }

    case 2:
      return OpMaker(OpType::g, 2)();

    default:
      abort();
  }
}

ExprPtr H(std::size_t k) {
  assert(k > 0 && k <= 2);
  return k == 1 ? H_(1) : H_(1) + H_(2);
}

ExprPtr F(bool use_f_tensor) {
  if (use_f_tensor) return OpMaker(OpType::f, 1)();

  // add \bar{g}^{\kappa x}_{\lambda y} \gamma^y_x with x,y in occ_space_type
  auto make_g_contribution = [](const auto occ_space_type) {
    return mbpt::OpMaker<Statistics::FermiDirac>::make(
        {IndexSpace::complete}, {IndexSpace::complete},
        [=](auto braidxs, auto ketidxs, Symmetry opsymm) {
          auto m1 = Index::make_tmp_index(
              IndexSpace{occ_space_type, IndexSpace::nullqns});
          auto m2 = Index::make_tmp_index(
              IndexSpace{occ_space_type, IndexSpace::nullqns});
          assert(opsymm == Symmetry::antisymm || opsymm == Symmetry::nonsymm);
          if (opsymm == Symmetry::antisymm) {
            braidxs.push_back(m1);
            ketidxs.push_back(m2);
            return ex<Tensor>(to_wstring(mbpt::OpType::g), braidxs, ketidxs,
                              Symmetry::antisymm) *
                   ex<Tensor>(to_wstring(mbpt::OpType::δ), IndexList{m2},
                              IndexList{m1}, Symmetry::nonsymm);
          } else {  // opsymm == Symmetry::nonsymm
            auto braidx_J = braidxs;
            braidx_J.push_back(m1);
            auto ketidxs_J = ketidxs;
            ketidxs_J.push_back(m2);
            auto braidx_K = braidxs;
            braidx_K.push_back(m1);
            auto ketidxs_K = ketidxs;
            ketidxs_K.emplace(begin(ketidxs_K), m2);
            return (ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_J, ketidxs_J,
                               Symmetry::nonsymm) -
                    ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_K, ketidxs_K,
                               Symmetry::nonsymm)) *
                   ex<Tensor>(to_wstring(mbpt::OpType::δ), IndexList{m2},
                              IndexList{m1}, Symmetry::nonsymm);
          }
        });
  };

  switch (get_default_context().vacuum()) {
    case Vacuum::Physical:
      return OpMaker(OpType::h, 1)() +
             make_g_contribution(IndexSpace::occupied);  // all occupieds
    case Vacuum::SingleProduct:
      return OpMaker(OpType::f, 1)();
    case Vacuum::MultiProduct:
      abort();
    default:
      abort();
  }
}

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  FWickTheorem wick{expr};
  wick.use_topology(use_top).set_nop_connections(nop_connections);
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

ExprPtr H2_oo_vv() {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::g); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return OpMaker(
            OpType::g,
            {IndexSpace::active_occupied, IndexSpace::active_occupied},
            {IndexSpace::active_unoccupied, IndexSpace::active_unoccupied})();
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{2, 0, 0, 2}, qns);
      });
}

ExprPtr H2_vv_vv() {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::g); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return OpMaker(
            OpType::g,
            {IndexSpace::active_unoccupied, IndexSpace::active_unoccupied},
            {IndexSpace::active_unoccupied, IndexSpace::active_unoccupied})();
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{0, 0, 2, 2}, qns);
      });
}

ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);
  switch (k) {
    case 1:
      return ex<op_t>(
          [vacuum = get_default_context().vacuum()]() -> std::wstring_view {
            switch (vacuum) {
              case Vacuum::Physical:
                return optype2label.at(OpType::h);
              case Vacuum::SingleProduct:
                return optype2label.at(OpType::f);
              case Vacuum::MultiProduct:
                abort();
              default:
                abort();
            }
          },
          [=]() -> ExprPtr {
            using namespace sequant::mbpt::sr;
            return sr::H_(1);
          },
          [=](qnc_t& qns) {
            qns = combine(qnc_t{{0, 1}, {0, 1}, {0, 1}, {0, 1}}, qns);
          });

    case 2:
      return ex<op_t>(
          []() -> std::wstring_view { return optype2label.at(OpType::g); },
          [=]() -> ExprPtr {
            using namespace sequant::mbpt::sr;
            return sr::H_(2);
          },
          [=](qnc_t& qns) {
            qns = combine(qnc_t{{0, 2}, {0, 2}, {0, 2}, {0, 2}}, qns);
          });

    default:
      abort();
  }
}

ExprPtr H(std::size_t k) {
  assert(k > 0 && k <= 2);
  return k == 1 ? H_(1) : H_(1) + H_(2);
}

ExprPtr T_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::t); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::T_(K);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{0ul, K, K, 0ul}, qns);
      });
}

ExprPtr T(std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));

  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result += T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::λ); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::Λ_(K);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{K, 0ul, 0ul, K}, qns);
      });
}

ExprPtr Λ(std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));

  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result += Λ_(k);
  }
  return result;
}

ExprPtr R_(std::size_t K_occ, std::size_t K_uocc) {
  assert(K_occ > 0 || K_uocc > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::R); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::R_(K_uocc, K_occ);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{0ul, K_occ, K_uocc, 0ul}, qns);
      });
}

ExprPtr R(std::size_t K_occ, std::size_t K_uocc) {
  using boost::numeric_cast;
  assert(K_occ > 0 || K_uocc > 0);
  ExprPtr result;

  for (auto o = numeric_cast<std::int64_t>(K_occ),
            u = numeric_cast<std::int64_t>(K_uocc);
       o > 0 || u > 0; --o, --u) {
    result += R_(o, u);
  }
  return result;
}

ExprPtr L_(std::size_t K_occ, std::size_t K_uocc) {
  assert(K_occ > 0 || K_uocc > 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::L); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::L_(K_occ, K_uocc);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{K_occ, 0ul, 0ul, K_uocc}, qns);
      });
}

ExprPtr L(std::size_t K_occ, std::size_t K_uocc) {
  using boost::numeric_cast;
  assert(K_occ > 0 || K_uocc > 0);
  ExprPtr result;

  for (auto o = numeric_cast<std::int64_t>(K_occ),
            u = numeric_cast<std::int64_t>(K_uocc);
       o > 0 || u > 0; --o, --u) {
    result += L_(o, u);
  }
  return result;
}

ExprPtr A(std::int64_t Kh, std::int64_t Kp) {
  if (Kp == std::numeric_limits<std::int64_t>::max()) Kp = Kh;

  assert(!(Kh == 0 && Kp == 0));
  // if they are not zero, Kh and Kp should have the same sign
  if (Kp != 0 && Kh != 0) {
    assert((Kh > 0 && Kp > 0) || (Kh < 0 && Kp < 0));
  }
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::A); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::A(Kh, Kp);
      },
      [=](qnc_t& qns) {
        const std::size_t abs_Kh = std::abs(Kh);
        const std::size_t abs_Kp = std::abs(Kp);
        if (Kp < 0 || Kh < 0)
          qns = combine(qnc_t{abs_Kp, 0ul, 0ul, abs_Kh}, qns);
        else
          qns = combine(qnc_t{0ul, abs_Kh, abs_Kp, 0ul}, qns);
      });
}

ExprPtr S(std::int64_t K) {
  assert(K != 0);
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::S); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::S(K, K);
      },
      [=](qnc_t& qns) {
        const std::size_t abs_K = std::abs(K);
        if (K < 0)
          qns = combine(qnc_t{abs_K, 0ul, 0ul, abs_K}, qns);
        else
          qns = combine(qnc_t{0ul, abs_K, abs_K, 0ul}, qns);
      });
}

ExprPtr P(std::int64_t Kh, std::int64_t Kp) {
  if (Kp == std::numeric_limits<std::int64_t>::max()) Kp = Kh;
  if (get_default_context().spbasis() == SPBasis::spinfree) {
    assert(Kp == Kh &&
           "non-particle conserving cases does not work with spin-free basis");
    const auto K = Kh;  // K = Kp = Kh
    return S(-K);
  } else {
    assert(get_default_context().spbasis() == SPBasis::spinorbital);
    return A(-Kh, -Kp);
  }
}

ExprPtr H_pt(std::size_t order, std::size_t R) {
  assert(R > 0);
  assert(order == 1 && "only first order perturbation is supported now");
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::h_1); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::H_pt(1, R);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{{0ul, R}, {0ul, R}, {0ul, R}, {0ul, R}}, qns);
      });
}

ExprPtr T_pt_(std::size_t order, std::size_t K) {
  assert(K > 0);
  assert(order == 1 && "only first order perturbation is supported now");
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::t_1); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::T_pt_(order, K);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{0ul, K, K, 0ul}, qns);
      });
}

ExprPtr T_pt(std::size_t order, std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + T_pt_(order, k) : T_pt_(order, k);
  }
  return result;
}

ExprPtr Λ_pt_(std::size_t order, std::size_t K) {
  assert(K > 0);
  assert(order == 1 && "only first order perturbation is supported now");
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::λ_1); },
      [=]() -> ExprPtr {
        using namespace sequant::mbpt::sr;
        return sr::Λ_pt_(order, K);
      },
      [=](qnc_t& qns) {
        qns = combine(qnc_t{K, 0ul, 0ul, K}, qns);
      });
}

ExprPtr Λ_pt(std::size_t order, std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + Λ_pt_(order, k) : Λ_pt_(order, k);
  }
  return result;
}

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

using mbpt::sr::vac_av;

#include <SeQuant/domain/mbpt/vac_av.ipp>

}  // namespace op

}  // namespace sr

// must be defined including op.ipp since it's used there
template <>
bool is_vacuum<sr::qns_t>(sr::qns_t qns) {
  return qns == sr::qns_t{};
}

}  // namespace mbpt

template <Statistics S>
std::wstring to_latex(const mbpt::Operator<mbpt::sr::qns_t, S>& op) {
  using namespace sequant::mbpt;
  using namespace sequant::mbpt::sr;

  auto result = L"{\\hat{" + utf_to_latex(op.label()) + L"}";

  // check if operator has adjoint label, remove if present for base label
  auto base_lbl = sequant::to_wstring(op.label());
  if (base_lbl.back() == adjoint_label) base_lbl.pop_back();

  auto it = label2optype.find(std::wstring(base_lbl));
  OpType optype = OpType::invalid;
  if (it != label2optype.end()) {  // handle special cases
    optype = it->second;
    if (to_class(optype) == OpClass::gen) {
      result += L"}";
      return result;
      // TODO: A(2, 1) produces \hat{A}, which has no info about ranks. Fix this
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

#include <SeQuant/domain/mbpt/op.ipp>

namespace sequant {
namespace mbpt {
template class Operator<sr::qns_t, Statistics::FermiDirac>;
template class Operator<sr::qns_t, Statistics::BoseEinstein>;
}  // namespace mbpt
}  // namespace sequant
