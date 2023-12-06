#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/domain/mbpt/context.hpp"

#include "SeQuant/core/math.hpp"
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"

#include <stdexcept>

namespace sequant::mbpt {

std::vector<std::wstring> cardinal_tensor_labels() {
  return {L"κ",  L"γ",
          L"Γ",  L"A",
          L"S",  L"P",
          L"L",  L"λ",
          L"λ¹", L"h",
          L"f",  L"f̃",
          L"g",  L"t",
          L"t¹", L"R",
          L"F",  L"X",
          L"μ",  L"V",
          L"Ṽ",  L"B",
          L"U",  L"GR",
          L"C",  overlap_label(),
          L"a",  L"ã",
          L"b",  L"b̃",
          L"E"};
}

std::wstring to_wstring(OpType op) {
  auto found_it = optype2label.find(op);
  if (found_it != optype2label.end())
    return found_it->second;
  else
    throw std::invalid_argument("to_wstring(OpType op): invalid op");
}

OpClass to_class(OpType op) {
  switch (op) {
    case OpType::h:
    case OpType::f:
    case OpType::f̃:
    case OpType::g:
    case OpType::RDM:
    case OpType::RDMCumulant:
    case OpType::δ:
    case OpType::A:
    case OpType::S:
    case OpType::h_1:
      return OpClass::gen;
    case OpType::t:
    case OpType::R:
    case OpType::R12:
    case OpType::t_1:
      return OpClass::ex;
    case OpType::λ:
    case OpType::L:
    case OpType::λ_1:
      return OpClass::deex;
    default:
      throw std::invalid_argument("to_class(OpType op): invalid op");
  }
}

qninterval_t ncre(qns_t qns, IndexSpace s) {
  assert(s.type() == IndexSpace::nonnulltype && s.qns() == IndexSpace::nullqns);
  return qns[0];
}

qninterval_t ncre(qns_t qns, IndexSpace::Type s) {
  assert(s == IndexSpace::nonnulltype);
  return qns[0];
}

qninterval_t ncre(qns_t qns) { return qns[0]; }

qninterval_t nann(qns_t qns, IndexSpace s) {
  assert(s.type() == IndexSpace::nonnulltype && s.qns() == IndexSpace::nullqns);
  return qns[1];
}

qninterval_t nann(qns_t qns, IndexSpace::Type s) {
  assert(s == IndexSpace::nonnulltype);
  return qns[1];
}

qninterval_t nann(qns_t qns) { return qns[1]; }

qns_t combine(qns_t a, qns_t b) {
  const auto ncontr =
      qninterval_t{0, std::min(ncre(b).upper(), nann(a).upper())};
  const auto nc = nonnegative(ncre(a) + ncre(b) - ncontr);
  const auto na = nonnegative(nann(a) + nann(b) - ncontr);
  return qns_t{nc, na};
}

// must be defined including op.ipp since it's used there
template <>
bool is_vacuum<qns_t>(qns_t qns) {
  return qns == qns_t{};
}

}  // namespace sequant::mbpt

namespace sequant {

mbpt::qns_t adjoint(mbpt::qns_t qns) {
  return mbpt::qns_t{nann(qns), ncre(qns)};
}

template <Statistics S>
std::wstring to_latex(const mbpt::Operator<mbpt::qns_t, S>& op) {
  using namespace sequant::mbpt;

  auto result = L"{\\hat{" + utf_to_latex(op.label()) + L"}";

  // check if operator has adjoint label, remove if present for base label
  auto base_lbl = sequant::to_wstring(op.label());
  if (base_lbl.back() == adjoint_label) base_lbl.pop_back();

  auto it = label2optype.find(base_lbl);
  if (it != label2optype.end()) {  // handle special cases
    const auto optype = it->second;
    if (to_class(optype) == OpClass::gen) {
      result += L"}";
      return result;
    }
  }

  // generic operator ... can only handle definite case
  const auto dN = op();
  if (!is_definite(ncre(dN)) || !is_definite(nann(dN))) {
    throw std::invalid_argument(
        "to_latex(const Operator<qns_t, S>& op): "
        "can only handle  generic operators with definite cre/ann numbers");
  }
  const auto dN_total = ncre(dN).lower() - nann(dN).lower();
  if (dN_total == 0) {  // N-conserving
    result += L"_{" + std::to_wstring(ncre(dN).lower()) + L"}";
  } else {  // N-nonconserving
    result += L"_{" + std::to_wstring(nann(dN).lower()) + L"}^{" +
              std::to_wstring(ncre(dN).lower()) + L"}";
  }
  result += L"}";
  return result;
}

}  // namespace sequant

#include "SeQuant/domain/mbpt/op.ipp"

namespace sequant::mbpt {

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, std::initializer_list<IndexSpace::Type> bras,
                    std::initializer_list<IndexSpace::Type> kets)
    : op_(op),
      bra_spaces_(bras.begin(), bras.end()),
      ket_spaces_(kets.begin(), kets.end()) {
  assert(nbra() > 0 || nket() > 0);
}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op) : op_(op) {}

template <Statistics S>
ExprPtr OpMaker<S>::operator()(std::optional<UseDepIdx> dep,
                               std::optional<Symmetry> opsymm_opt) const {
  // if not given dep, use mbpt::Context::CSV to determine whether to use
  // dependent indices for pure (de)excitation ops
  if (!dep && get_default_formalism().csv() == mbpt::CSV::Yes) {
    if (to_class(op_) == OpClass::ex) {
      for (auto&& s : bra_spaces_) {
        assert(s == IndexSpace::complete_unoccupied ||
               s == IndexSpace::active_unoccupied);
      }
      dep = UseDepIdx::Bra;
    } else if (to_class(op_) == OpClass::deex) {
      for (auto&& s : ket_spaces_) {
        assert(s == IndexSpace::complete_unoccupied ||
               s == IndexSpace::active_unoccupied);
      }
      dep = UseDepIdx::Ket;
    } else {
      dep = UseDepIdx::None;
    }
  }

  return make(
      bra_spaces_, ket_spaces_,
      [this, opsymm_opt](const auto& braidxs, const auto& ketidxs,
                         Symmetry opsymm) {
        return ex<Tensor>(to_wstring(op_), braidxs, ketidxs,
                          opsymm_opt ? *opsymm_opt : opsymm);
      },
      dep ? *dep : UseDepIdx::None);
}

template class OpMaker<Statistics::FermiDirac>;
template class OpMaker<Statistics::BoseEinstein>;

template class Operator<qns_t, Statistics::FermiDirac>;
template class Operator<qns_t, Statistics::BoseEinstein>;

}  // namespace sequant::mbpt
