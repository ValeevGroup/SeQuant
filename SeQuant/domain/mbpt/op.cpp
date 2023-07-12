#include "SeQuant/domain/mbpt/op.hpp"
#include "SeQuant/domain/mbpt/context.hpp"

#include "SeQuant/core/math.hpp"
#include "SeQuant/core/op.hpp"
#include "SeQuant/core/tensor.hpp"

#include <stdexcept>

namespace sequant::mbpt {

std::vector<std::wstring> cardinal_tensor_labels() {
  return {L"\\lambda", L"\\gamma",
          L"\\Gamma",  L"A",
          L"S",        L"P",
          L"L",        L"λ",
          L"h",        L"f",
          L"g",        L"t",
          L"R",        L"F",
          L"X",        L"V",
          L"Ṽ",        L"B",
          L"U",        L"GR",
          L"C",        overlap_label(),
          L"a",        L"ã",
          L"b",        L"ᵬ",
          L"E"};
}

std::wstring to_wstring(OpType op) {
  switch (op) {
    case OpType::h:
      return L"h";
    case OpType::f:
      return L"f";
    case OpType::g:
      return L"g";
    case OpType::t:
      return L"t";
    case OpType::lambda:
      return L"λ";
    case OpType::A:
      return L"A";
    case OpType::L:
      return L"L";
    case OpType::R:
      return L"R";
    case OpType::R12:
      return L"F";
    case OpType::GR:
      return L"GR";
    case OpType::C:
      return L"C";
    default:
      throw std::invalid_argument("to_wstring(OpType op): invalid op");
  }
}

OpClass to_class(OpType op) {
  switch (op) {
    case OpType::h:
    case OpType::f:
    case OpType::g:
      return OpClass::gen;
    case OpType::t:
    case OpType::R:
    case OpType::R12:
      return OpClass::ex;
    case OpType::lambda:
    case OpType::A:
    case OpType::L:
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

  auto lbl = std::wstring(op.label());
  std::wstring result = L"{\\hat{" + lbl + L"}";
  auto it = label2optype.find(lbl);
  if (it != label2optype.end()) {  // handle special cases
    const auto optype = it->second;
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
ExprPtr OpMaker<S>::operator()() const {
  const bool symm = get_default_formalism().nbody_interaction_tensor_symm() ==
                    Context::NBodyInteractionTensorSymm::Yes;
  const bool csv = get_default_formalism().csv() == Context::CSV::Yes;
  bool csv_bra = false;
  bool csv_ket = false;
  if (csv) {  // validate spaces
    if (to_class(op_) == OpClass::ex) {
      for (auto&& s : bra_spaces_) {
        assert(s == IndexSpace::complete_unoccupied ||
               s == IndexSpace::active_unoccupied);
      }
      csv_bra = true;
    } else if (to_class(op_) == OpClass::deex) {
      for (auto&& s : ket_spaces_) {
        assert(s == IndexSpace::complete_unoccupied ||
               s == IndexSpace::active_unoccupied);
      }
      csv_ket = true;
    }
  }

  // not sure what it means to use nonsymmetric operator if nbra != nket
  if (!symm) assert(nbra() == nket());

  auto make_idx_vector = [](const auto& spacetypes) {
    std::vector<Index> result;
    const auto n = spacetypes.size();
    result.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      auto space = IndexSpace::instance(spacetypes[i]);
      result.push_back(Index::make_tmp_index(space));
    }
    return result;
  };

  auto make_depidx_vector = [](const auto& spacetypes, auto&& protoidxs) {
    const auto n = spacetypes.size();
    std::vector<Index> result;
    result.reserve(n);
    for (size_t i = 0; i != n; ++i) {
      auto space = IndexSpace::instance(spacetypes[i]);
      result.push_back(Index::make_tmp_index(space, protoidxs, true));
    }
    return result;
  };

  std::vector<Index> braidxs, ketidxs;
  if (csv_bra) {
    ketidxs = make_idx_vector(ket_spaces_);
    braidxs = make_depidx_vector(bra_spaces_, ketidxs);
  } else if (csv_ket) {
    braidxs = make_idx_vector(bra_spaces_);
    ketidxs = make_depidx_vector(ket_spaces_, braidxs);
  } else {
    braidxs = make_idx_vector(bra_spaces_);
    ketidxs = make_idx_vector(ket_spaces_);
  }

  const auto mult =
      symm ? factorial(nbra()) * factorial(nket()) : factorial(nbra());
  const auto opsymm =
      symm ? (S == Statistics::FermiDirac ? Symmetry::antisymm : Symmetry::symm)
           : Symmetry::nonsymm;
  return ex<Constant>(rational{1, mult}) *
         ex<Tensor>(to_wstring(op_), braidxs, ketidxs, opsymm) *
         ex<NormalOperator<S>>(/* creators */ braidxs,
                               /* annihilators */ ketidxs,
                               get_default_context().vacuum());
}

template class OpMaker<Statistics::FermiDirac>;
template class OpMaker<Statistics::BoseEinstein>;

template class Operator<qns_t, Statistics::FermiDirac>;
template class Operator<qns_t, Statistics::BoseEinstein>;

}  // namespace sequant::mbpt
