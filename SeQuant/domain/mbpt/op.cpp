#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

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

// Excitation type QNs will have quasiparticle annihilators in every space which
// intersects with the active hole space. The active particle space will have
// quasiparticle creators. The presence of additional blocks depends on whether
// the corresponding active hole or active particle space is a base space.
qns_t excitation_type_qns(std::size_t K, const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, K};
    result[1] = {0ul, K};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    // are the active qp spaces base spaces?
    bool aps_base = isr->is_base(isr->particle_space(SQN));
    bool ahs_base = isr->is_base(isr->hole_space(SQN));

    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // ex -> creators in particle space
      if (includes(isr->particle_space(SQN).type(), base_space.type())) {
        result[i * 2] = {aps_base ? K : 0ul, K};
      }
      // ex -> annihilators in hole space
      if (includes(isr->hole_space(SQN).type(), base_space.type())) {
        result[i * 2 + 1] = {ahs_base ? K : 0ul, K};
      }
    }
  }
  return result;
}

qns_t interval_excitation_type_qns(std::size_t K,
                                   const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, K};
    result[1] = {0ul, K};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();

    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // ex -> creators in particle space
      if (includes(isr->particle_space(SQN).type(), base_space.type())) {
        result[i * 2] = {0ul, K};
      }
      // ex -> annihilators in hole space
      if (includes(isr->hole_space(SQN).type(), base_space.type())) {
        result[i * 2 + 1] = {0ul, K};
      }
    }
  }
  return result;
}

qns_t deexcitation_type_qns(std::size_t K,
                            const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, K};
    result[1] = {0ul, K};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    bool aps_base = isr->is_base(isr->particle_space(SQN));
    bool ahs_base = isr->is_base(isr->hole_space(SQN));
    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // deex -> annihilators in particle space
      if (includes(isr->particle_space(SQN).type(), base_space.type())) {
        result[i * 2 + 1] = {aps_base ? K : 0ul, K};
      }
      // deex -> creators in hole space
      if (includes(isr->hole_space(SQN).type(), base_space.type())) {
        result[i * 2] = {ahs_base ? K : 0ul, K};
      }
    }
  }
  return result;
}

qns_t interval_deexcitation_type_qns(std::size_t K,
                                     const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, K};
    result[1] = {0ul, K};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // deex -> annihilators in particle space
      if (includes(isr->particle_space(SQN).type(), base_space.type())) {
        result[i * 2 + 1] = {0ul, K};
      }
      // deex -> creators in hole space
      if (includes(isr->hole_space(SQN).type(), base_space.type())) {
        result[i * 2] = {0ul, K};
      }
    }
  }
  return result;
}

qns_t general_type_qns(std::size_t K) {
  qnc_t result;
  for (int i = 0; i < result.size(); i++) {
    result[i] = {0, K};
  }
  return result;
}

qns_t generic_excitation_qns(std::size_t particle_rank, std::size_t hole_rank,
                             IndexSpace particle_space, IndexSpace hole_space,
                             const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, hole_rank};
    result[1] = {0ul, particle_rank};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    bool aps_base = isr->is_base(isr->particle_space(SQN));
    bool ahs_base = isr->is_base(isr->hole_space(SQN));
    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // ex -> creators in particle space
      if (includes(particle_space.type(), base_space.type())) {
        result[i * 2] = {aps_base ? particle_rank : 0ul,
                         particle_rank};  // creators
      }
      // ex -> annihilators in hole space
      if (includes(hole_space.type(), base_space.type())) {
        result[i * 2 + 1] = {ahs_base ? hole_rank : 0ul,
                             hole_rank};  // annihilators
      }
    }
  }
  return result;
}

qns_t generic_deexcitation_qns(std::size_t particle_rank, std::size_t hole_rank,
                               IndexSpace particle_space, IndexSpace hole_space,
                               const IndexSpace::QuantumNumbers SQN) {
  qnc_t result;
  if (get_default_context().vacuum() == Vacuum::Physical) {
    result[0] = {0ul, particle_rank};
    result[1] = {0ul, hole_rank};
  } else {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    bool aps_base = isr->is_base(isr->particle_space(SQN));
    bool ahs_base = isr->is_base(isr->hole_space(SQN));
    for (int i = 0; i < base_spaces.size(); i++) {
      const auto& base_space = base_spaces[i];
      result[i * 2] = {0ul, 0ul};
      result[i * 2 + 1] = {0ul, 0ul};
      if (base_space.qns() != SQN) continue;

      // deex -> creators in hole space
      if (includes(hole_space.type(), base_space.type())) {
        result[i * 2] = {ahs_base ? hole_rank : 0ul, hole_rank};  // creators
      }
      // deex -> annihilators in particle space
      if (includes(particle_space.type(), base_space.type())) {
        result[i * 2 + 1] = {aps_base ? particle_rank : 0ul,
                             particle_rank};  // annihilators
      }
    }
  }
  return result;
}

// By counting the number of contractions between indices of proper type, we can
// know the quantum numbers of a combined result.
qns_t combine(qns_t a, qns_t b) {
  assert(a.size() == b.size());
  qns_t result;

  if (get_default_context().vacuum() == Vacuum::Physical) {
    qns_t result;
    const auto ncontr = qninterval_t{0, std::min(b[0].upper(), a[1].upper())};
    const auto nc = nonnegative(a[0] + b[0] - ncontr);
    const auto na = nonnegative(a[1] + b[1] - ncontr);
    result[0] = nc;
    result[1] = na;
    return result;
  } else if (get_default_context().vacuum() == Vacuum::SingleProduct) {
    auto isr = get_default_context().index_space_registry();
    const auto& base_spaces = isr->base_spaces();
    for (auto i = 0; i < base_spaces.size(); i++) {
      auto cre = i * 2;
      auto ann = (i * 2) + 1;
      auto base_is_fermi_occupied = isr->is_pure_occupied(
          base_spaces[i]);  // need to distinguish particle and hole
                            // contractions.
      auto ncontr_space =
          base_is_fermi_occupied
              ? qninterval_t{0, std::min(b[ann].upper(), a[cre].upper())}
              : qninterval_t{0, std::min(b[cre].upper(), a[ann].upper())};
      auto nc_space = nonnegative(b[cre] + a[cre] - ncontr_space);
      auto na_space = nonnegative(b[ann] + a[ann] - ncontr_space);
      result[cre] = nc_space;
      result[ann] = na_space;
    }
    return result;
  } else {
    throw "unsupported vacuum context.";
  }
}

// must be defined including op.ipp since it's used there
template <>
bool is_vacuum<qns_t>(qns_t qns) {
  return qns == qns_t{};
}

}  // namespace sequant::mbpt

namespace sequant {

mbpt::qns_t adjoint(mbpt::qns_t qns) {
  mbpt::qns_t new_qnst;
  new_qnst.resize(qns.size());
  auto is_even = [](int i) { return i % 2 == 0; };
  for (int i = 0; i < qns.size(); i++) {
    if (is_even(i)) {
      new_qnst[i + 1] = qns[i];
    } else {
      new_qnst[i - 1] = qns[i];
    }
  }
  return new_qnst;
}

template <Statistics S>
std::wstring to_latex(const mbpt::Operator<mbpt::qns_t, S>& op) {
  using namespace sequant::mbpt;

  auto result = L"{\\hat{" + utf_to_latex(op.label()) + L"}";

  // check if operator has adjoint label, remove if present for base label
  auto base_lbl = sequant::to_wstring(op.label());
  bool is_adjoint = false;
  if (base_lbl.back() == adjoint_label) {
    is_adjoint = true;
    base_lbl.pop_back();
  }

  auto it = label2optype.find(base_lbl);
  OpType optype = OpType::invalid;
  if (it != label2optype.end()) {  // handle special cases
    optype = it->second;
    if (to_class(optype) == OpClass::gen) {
      result += L"}";
      return result;
    }
  }
  std::wstring baseline_char = is_adjoint ? L"^" : L"_";
  if (get_default_context().vacuum() == Vacuum::Physical) {
    if (op()[0] == op()[1]) {  // particle conserving
      result += L"_{" + std::to_wstring(op()[0].lower()) + L"}";
    } else {  // non-particle conserving
      result += L"_{" + std::to_wstring(op()[1].lower()) + L"}^{" +
                std::to_wstring(op()[0].lower()) + L"}";
    }
  } else {  // single product vacuum
    auto nann_p = is_adjoint ? op().ncre_particles() : op().nann_particles();
    auto ncre_h = is_adjoint ? op().nann_holes() : op().ncre_holes();
    auto ncre_p = is_adjoint ? op().nann_particles() : op().ncre_particles();
    auto nann_h = is_adjoint ? op().ncre_holes() : op().nann_holes();

    if (!is_definite(nann_p) || !is_definite(ncre_h) || !is_definite(ncre_p) ||
        !is_definite(nann_h)) {
      throw std::invalid_argument(
          "to_latex(const Operator<qns_t, S>& op): "
          "can only handle generic operators with definite cre/ann numbers");
    }

    // pure quasiparticle creator/annihilator?
    const auto qprank_cre = ncre_p.lower() + nann_h.lower();
    const auto qprank_ann = nann_p.lower() + ncre_h.lower();
    const auto qppure = qprank_cre == 0 || qprank_ann == 0;
    if (qppure) {
      if (qprank_cre) {
        if (ncre_p.lower() == nann_h.lower()) {  // q-particle conserving
          result +=
              baseline_char + L"{" + std::to_wstring(nann_h.lower()) + L"}";
        } else {  // particle non-conserving
          result += baseline_char + L"{" + std::to_wstring(nann_h.lower()) +
                    L"," + std::to_wstring(ncre_p.lower()) + L"}";
        }
      } else {
        if (ncre_h.lower() == nann_p.lower()) {  // q-particle conserving
          result +=
              baseline_char + L"{" + std::to_wstring(ncre_h.lower()) + L"}";
        } else {  // q-particle non-conserving
          result += baseline_char + L"{" + std::to_wstring(ncre_h.lower()) +
                    L"," + std::to_wstring(nann_p.lower()) + L"}";
        }
      }
    } else {  // not pure qp creator/annihilator
      result += L"_{" + std::to_wstring(nann_h.lower()) + L"," +
                std::to_wstring(ncre_p.lower()) + L"}^{" +
                std::to_wstring(ncre_h.lower()) + L"," +
                std::to_wstring(nann_p.lower()) + L"}";
    }
  }
  result += L"}";
  return result;
}

}  // namespace sequant

#include <SeQuant/domain/mbpt/op.ipp>

namespace sequant::mbpt {

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, std::initializer_list<IndexSpace> bras,
                    std::initializer_list<IndexSpace> kets)
    : op_(op),
      bra_spaces_(bras.begin(), bras.end()),
      ket_spaces_(kets.begin(), kets.end()) {
  assert(nbra() > 0 || nket() > 0);
}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op) : op_(op) {}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, std::size_t nbra, std::size_t nket,
                    IndexSpace particle_space, IndexSpace hole_space) {
  op_ = op;
  assert(nbra > 0 || nket > 0);
  auto isr = get_default_context().index_space_registry();
  switch (to_class(op)) {
    case OpClass::ex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, particle_space);
      ket_spaces_ = decltype(ket_spaces_)(nket, hole_space);
      break;
    case OpClass::deex:
      bra_spaces_ = decltype(bra_spaces_)(nbra, hole_space);
      ket_spaces_ = decltype(ket_spaces_)(nket, particle_space);
      break;
    case OpClass::gen:
      bra_spaces_ = decltype(bra_spaces_)(nbra, isr->complete_space(Spin::any));
      ket_spaces_ = decltype(ket_spaces_)(nket, isr->complete_space(Spin::any));
      break;
  }
}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, std::size_t nparticle) {
  op_ = op;
  auto isr = get_default_context().index_space_registry();
  auto current_context = get_default_context();
  const auto hole_space = isr->hole_space(Spin::any);
  const auto particle_space = isr->particle_space(Spin::any);
  switch (to_class(op)) {
    case OpClass::ex:
      bra_spaces_ = decltype(bra_spaces_)(nparticle, particle_space);
      ket_spaces_ = decltype(ket_spaces_)(nparticle, hole_space);
      break;
    case OpClass::deex:
      bra_spaces_ = decltype(bra_spaces_)(nparticle, hole_space);
      ket_spaces_ = decltype(ket_spaces_)(nparticle, particle_space);
      break;
    case OpClass::gen:
      bra_spaces_ =
          decltype(bra_spaces_)(nparticle, isr->complete_space(Spin::any));
      ket_spaces_ =
          decltype(ket_spaces_)(nparticle, isr->complete_space(Spin::any));
      break;
  }
}

template <Statistics S>
ExprPtr OpMaker<S>::operator()(std::optional<UseDepIdx> dep,
                               std::optional<Symmetry> opsymm_opt) const {
  auto isr = get_default_context(Statistics::FermiDirac).index_space_registry();
  // if not given dep, use mbpt::Context::CSV to determine whether to use
  // dependent indices for pure (de)excitation ops
  if (!dep && get_default_formalism().csv() == mbpt::CSV::Yes) {
    if (to_class(op_) == OpClass::ex) {
      for (auto&& s : bra_spaces_) {
        assert(isr->contains_unoccupied(s));
      }
      dep = UseDepIdx::Bra;
    } else if (to_class(op_) == OpClass::deex) {
      for (auto&& s : ket_spaces_) {
        assert(isr->contains_unoccupied(s));
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

inline namespace op {

namespace tensor {
ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);
  switch (k) {
    case 1:
      switch (get_default_context().vacuum()) {
        case Vacuum::Physical:
          return OpMaker<Statistics::FermiDirac>(OpType::h, 1)();
        case Vacuum::SingleProduct:
          return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
        case Vacuum::MultiProduct:
          return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
        default:
          abort();
      }

    case 2:
      return OpMaker<Statistics::FermiDirac>(OpType::g, 2)();

    default:
      abort();
  }
}

ExprPtr H(std::size_t k) {
  assert(k > 0 && k <= 2);
  return k == 1 ? tensor::H_(1) : tensor::H_(1) + tensor::H_(2);
}

ExprPtr F(bool use_tensor, IndexSpace reference_occupied) {
  if (use_tensor) {
    return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
  } else {                       // explicit density matrix construction
    assert(reference_occupied);  // cannot explicitly instantiate fock operator
                                 // without providing an occupied indexspace
    // add \bar{g}^{\kappa x}_{\lambda y} \gamma^y_x with x,y in occ_space_type
    auto make_g_contribution = [](const auto occ_space) {
      auto isr = get_default_context().index_space_registry();
      return mbpt::OpMaker<Statistics::FermiDirac>::make(
          {isr->complete_space(Spin::any)}, {isr->complete_space(Spin::any)},
          [=](auto braidxs, auto ketidxs, Symmetry opsymm) {
            auto m1 = Index::make_tmp_index(occ_space);
            auto m2 = Index::make_tmp_index(occ_space);
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
              using std::begin;
              ketidxs_K.emplace(begin(ketidxs_K), m2);
              return (ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_J,
                                 ketidxs_J, Symmetry::nonsymm) -
                      ex<Tensor>(to_wstring(mbpt::OpType::g), braidx_K,
                                 ketidxs_K, Symmetry::nonsymm)) *
                     ex<Tensor>(to_wstring(mbpt::OpType::δ), IndexList{m2},
                                IndexList{m1}, Symmetry::nonsymm);
            }
          });
    };
    auto isr = get_default_context().index_space_registry();
    return OpMaker<Statistics::FermiDirac>(OpType::h, 1)() +
           make_g_contribution(reference_occupied);
  }
}

ExprPtr T_(std::size_t K) {
  return OpMaker<Statistics::FermiDirac>(OpType::t, K)();
}

ExprPtr T(std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = skip1 ? 2ul : 1ul; k <= K; ++k) {
    result += tensor::T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  return OpMaker<Statistics::FermiDirac>(OpType::λ, K)();
}

ExprPtr Λ(std::size_t K) {
  assert(K > 0);

  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + tensor::Λ_(k) : tensor::Λ_(k);
  }
  return result;
}

ExprPtr R_(std::size_t nbra, std::size_t nket, IndexSpace hole_space,
           IndexSpace particle_space) {
  return OpMaker<Statistics::FermiDirac>(OpType::R, nbra, nket, particle_space,
                                         hole_space)();
}

ExprPtr L_(std::size_t nbra, std::size_t nket, IndexSpace hole_space,
           IndexSpace particle_space) {
  return OpMaker<Statistics::FermiDirac>(OpType::L, nbra, nket, particle_space,
                                         hole_space)();
}

ExprPtr P(std::int64_t K) {
  return get_default_context().spbasis() == SPBasis::spinfree ? tensor::S(-K)
                                                              : tensor::A(-K);
}

ExprPtr A(std::int64_t K) {
  auto isr = get_default_context().index_space_registry();
  assert(K != 0);
  container::svector<IndexSpace> creators;
  container::svector<IndexSpace> annihilators;
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      annihilators.emplace_back(isr->hole_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      creators.emplace_back(isr->hole_space(Spin::any));
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      creators.emplace_back(isr->particle_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      annihilators.emplace_back(isr->particle_space(Spin::any));

  std::optional<OpMaker<Statistics::FermiDirac>::UseDepIdx> dep;
  if (get_default_formalism().csv() == mbpt::CSV::Yes)
    dep = K > 0 ? OpMaker<Statistics::FermiDirac>::UseDepIdx::Bra
                : OpMaker<Statistics::FermiDirac>::UseDepIdx::Ket;
  return OpMaker<Statistics::FermiDirac>(OpType::A, creators, annihilators)(
      dep, {Symmetry::antisymm});
}

ExprPtr S(std::int64_t K) {
  auto isr = get_default_context().index_space_registry();
  assert(K != 0);
  container::svector<IndexSpace> creators;
  container::svector<IndexSpace> annihilators;
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      annihilators.emplace_back(isr->hole_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      creators.emplace_back(isr->hole_space(Spin::any));
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      creators.emplace_back(isr->particle_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      annihilators.emplace_back(isr->particle_space(Spin::any));

  std::optional<OpMaker<Statistics::FermiDirac>::UseDepIdx> dep;
  if (get_default_formalism().csv() == mbpt::CSV::Yes)
    dep = K > 0 ? OpMaker<Statistics::FermiDirac>::UseDepIdx::Bra
                : OpMaker<Statistics::FermiDirac>::UseDepIdx::Ket;
  return OpMaker<Statistics::FermiDirac>(OpType::S, creators, annihilators)(
      dep, {Symmetry::nonsymm});
}

ExprPtr H_pt(std::size_t order, std::size_t R) {
  assert(order == 1 &&
         "sequant::sr::H_pt(): only supports first order perturbation");
  assert(R > 0);
  return OpMaker<Statistics::FermiDirac>(OpType::h_1, R)();
}

ExprPtr T_pt_(std::size_t order, std::size_t K) {
  assert(order == 1 &&
         "sequant::sr::T_pt_(): only supports first order perturbation");
  return OpMaker<Statistics::FermiDirac>(OpType::t_1, order, K)();
}

ExprPtr T_pt(std::size_t order, std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + tensor::T_pt_(order, k) : tensor::T_pt_(order, k);
  }
  return result;
}

ExprPtr Λ_pt_(std::size_t order, std::size_t K) {
  assert(order == 1 &&
         "sequant::sr::Λ_pt_(): only supports first order perturbation");
  return OpMaker<Statistics::FermiDirac>(OpType::λ_1, order, K)();
}

ExprPtr Λ_pt(std::size_t order, std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + tensor::Λ_pt_(order, k) : tensor::Λ_pt_(order, k);
  }
  return result;
}

}  // namespace tensor

ExprPtr H_(std::size_t k) {
  assert(k > 0 && k <= 2);
  switch (k) {
    case 1:
      return ex<op_t>(
          [vacuum = get_default_context().vacuum()]() -> std::wstring_view {
            switch (vacuum) {
              case Vacuum::Physical:
                return L"h";
              case Vacuum::SingleProduct:
                return L"f";
              case Vacuum::MultiProduct:
                return L"f";
              default:
                abort();
            }
          },
          [=]() -> ExprPtr { return tensor::H_(1); },
          [=](qnc_t& qns) {
            qnc_t op_qnc_t = general_type_qns(1);
            qns = combine(op_qnc_t, qns);
          });

    case 2:
      return ex<op_t>([]() -> std::wstring_view { return L"g"; },
                      [=]() -> ExprPtr { return tensor::H_(2); },
                      [=](qnc_t& qns) {
                        qnc_t op_qnc_t = general_type_qns(2);
                        qns = combine(op_qnc_t, qns);
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
  return ex<op_t>([]() -> std::wstring_view { return L"t"; },
                  [=]() -> ExprPtr { return tensor::T_(K); },
                  [=](qnc_t& qns) {
                    qnc_t op_qnc_t = excitation_type_qns(K);
                    qns = combine(op_qnc_t, qns);
                  });
}

ExprPtr T(std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = skip1 ? 2ul : 1ul; k <= K; ++k) {
    result += T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  assert(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                  [=]() -> ExprPtr { return tensor::Λ_(K); },
                  [=](qnc_t& qns) {
                    qnc_t op_qnc_t = deexcitation_type_qns(K);
                    qns = combine(op_qnc_t, qns);
                  });
}

ExprPtr Λ(std::size_t K) {
  assert(K > 0);
  ExprPtr result;
  for (auto k = 1ul; k <= K; ++k) {
    result = k > 1 ? result + Λ_(k) : Λ_(k);
  }
  return result;
}

ExprPtr F(bool use_f_tensor, IndexSpace occupied_density) {
  if (use_f_tensor) {
    return ex<op_t>(
        []() -> std::wstring_view { return L"f"; },
        [=]() -> ExprPtr { return tensor::F(true, occupied_density); },
        [=](qnc_t& qns) {
          qnc_t op_qnc_t = general_type_qns(1);
          qns = combine(op_qnc_t, qns);
        });
  } else {
    throw "non-tensor use at operator level not yet supported";
  }
}

ExprPtr A(std::int64_t K) {
  auto isr = get_default_context().index_space_registry();
  assert(K != 0);
  return ex<op_t>([]() -> std::wstring_view { return L"A"; },
                  [=]() -> ExprPtr { return tensor::A(K); },
                  [=](qnc_t& qns) {
                    const std::size_t abs_K = std::abs(K);
                    if (K < 0) {
                      qnc_t op_qnc_t = deexcitation_type_qns(abs_K);
                      qns = combine(op_qnc_t, qns);
                    } else {
                      qnc_t op_qnc_t = excitation_type_qns(abs_K);
                      qns = combine(op_qnc_t, qns);
                    }
                  });
}

ExprPtr S(std::int64_t K) {
  assert(K != 0);
  return ex<op_t>([]() -> std::wstring_view { return L"S"; },
                  [=]() -> ExprPtr { return tensor::S(K); },
                  [=](qnc_t& qns) {
                    const std::size_t abs_K = std::abs(K);
                    if (K < 0) {
                      qnc_t op_qnc_t = deexcitation_type_qns(abs_K);
                      qns = combine(op_qnc_t, qns);
                    } else {
                      qnc_t op_qnc_t = excitation_type_qns(abs_K);
                      qns = combine(op_qnc_t, qns);
                    }
                  });
}

ExprPtr P(std::int64_t K) {
  return get_default_context().spbasis() == SPBasis::spinfree ? S(-K) : A(-K);
}

ExprPtr H_pt(std::size_t order, std::size_t R) {
  assert(R > 0);
  assert(order == 1 && "only first order perturbation is supported now");
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::h_1); },
      [=]() -> ExprPtr { return tensor::H_pt(order, R); },
      [=](qnc_t& qns) { qns = combine(general_type_qns(R), qns); });
}

ExprPtr T_pt_(std::size_t order, std::size_t K) {
  assert(K > 0);
  assert(order == 1 && "only first order perturbation is supported now");
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::t_1); },
      [=]() -> ExprPtr { return tensor::T_pt_(order, K); },
      [=](qnc_t& qns) { qns = combine(excitation_type_qns(K), qns); });
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
      [=]() -> ExprPtr { return tensor::Λ_pt_(order, K); },
      [=](qnc_t& qns) { qns = combine(deexcitation_type_qns(K), qns); });
}

ExprPtr Λ_pt(std::size_t order, std::size_t K, bool skip1) {
  assert(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + Λ_pt_(order, k) : Λ_pt_(order, k);
  }
  return result;
}

ExprPtr R_(std::size_t Nbra, std::size_t Nket, IndexSpace hole_space,
           IndexSpace particle_space) {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::R); },
      [=]() -> ExprPtr {
        return tensor::R_(Nbra, Nket, hole_space, particle_space);
      },
      [=](qnc_t& qns) {
        qns = combine(
            generic_excitation_qns(Nket, Nbra, particle_space, hole_space),
            qns);
      });
}

ExprPtr L_(std::size_t Nbra, std::size_t Nket, IndexSpace hole_space,
           IndexSpace particle_space) {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::L); },
      [=]() -> ExprPtr {
        return tensor::L_(Nbra, Nket, hole_space, particle_space);
      },
      [=](qnc_t& qns) {
        qns = combine(
            generic_deexcitation_qns(Nbra, Nket, particle_space, hole_space),
            qns);
      });
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

  return can_change_qns(op_or_op_product, interval_excitation_type_qns(k));
}

bool lowers_rank_or_lower_to_vacuum(const ExprPtr& op_or_op_product,
                                    const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{},
                        interval_excitation_type_qns(k));
}

bool raises_vacuum_to_rank(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, excitation_type_qns(k));
}

bool lowers_rank_to_vacuum(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  assert(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{}, excitation_type_qns(k));
}

#include <SeQuant/domain/mbpt/vac_av.ipp>

namespace tensor {

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  simplify(expr);
  auto isr = get_default_context().index_space_registry();
  const auto spinorbital =
      get_default_context().spbasis() == SPBasis::spinorbital;
  // convention is to use different label for spin-orbital and spin-free RDM
  const auto rdm_label = spinorbital ? optype2label.at(OpType::RDM) : L"Γ";

  // only need full contractions if don't have any density outside of
  // the orbitals occupied in the vacuum
  bool full_contractions =
      (isr->reference_occupied_space() == isr->vacuum_occupied_space()) ? true
                                                                        : false;
  // N.B. reference < vacuum is not yet supported
  if (isr->reference_occupied_space().intersection(
          isr->vacuum_occupied_space()) != isr->vacuum_occupied_space()) {
    throw std::invalid_argument(
        "mbpt::tensor::vac_av: vacuum occupied orbitals must be same as or "
        "subset of the reference orbital set.");
  }

  FWickTheorem wick{expr};
  wick.use_topology(use_top).set_nop_connections(nop_connections);
  wick.full_contractions(full_contractions);
  auto result = wick.compute();
  simplify(result);
  // std::wcout << "post wick: " << to_latex_align(result,20,1) << std::endl;
  if (Logger::get_instance().wick_stats) {
    std::wcout << "WickTheorem stats: # of contractions attempted = "
               << wick.stats().num_attempted_contractions
               << " # of useful contractions = "
               << wick.stats().num_useful_contractions << std::endl;
  }
  // only need to handle the special case where the dense(at least partially
  // occupied) states, contain additional functions to the vacuum_occupied.
  // including a density occupied partition using a "single-reference" method
  // will replace FNOPs with RDMs. i.e. "multi-reference" RDM replacement rules
  // work in the limit of one reference.
  if (isr->reference_occupied_space() == IndexSpace::Type{} ||
      isr->reference_occupied_space(Spin::any) ==
          isr->vacuum_occupied_space(Spin::any)) {
    return result;
  } else {
    const auto target_rdm_space_type =
        get_default_context().vacuum() == Vacuum::SingleProduct
            ? isr->intersection(isr->particle_space(Spin::any),
                                isr->hole_space(Spin::any))
            : isr->reference_occupied_space(Spin::any);

    // STEP1. replace NOPs by RDM
    auto replace_nop_with_rdm = [&rdm_label, spinorbital](ExprPtr& exptr) {
      auto replace = [&rdm_label, spinorbital](const auto& nop) -> ExprPtr {
        using index_container = container::svector<Index>;
        auto braidxs = nop.annihilators() |
                       ranges::views::transform(
                           [](const auto& op) { return op.index(); }) |
                       ranges::to<index_container>();
        auto ketidxs = nop.creators() |
                       ranges::views::transform(
                           [](const auto& op) { return op.index(); }) |
                       ranges::to<index_container>();
        assert(braidxs.size() ==
               ketidxs.size());  // need to handle particle # violating case?
        const auto rank = braidxs.size();
        return ex<Tensor>(
            rdm_label, braidxs, ketidxs,
            rank > 1 && spinorbital ? Symmetry::antisymm : Symmetry::nonsymm);
      };

      if (exptr.template is<FNOperator>()) {
        exptr = replace(exptr.template as<FNOperator>());
      } else if (exptr.template is<BNOperator>()) {
        exptr = replace(exptr.template as<BNOperator>());
      }
    };
    result->visit(replace_nop_with_rdm, true);

    // STEP 2: project RDM indices onto the target RDM subspace
    // since RDM indices only make sense within a single TN expand + flatten
    // first, then do the projection individually for each TN
    expand(result);
    // flatten(result);  // TODO where is flatten?
    auto project_rdm_indices_to_target = [&](ExprPtr& exptr) {
      auto impl_for_single_tn = [&](ProductPtr& product_ptr) {
        // enlist all indices and count their instances
        auto for_each_index_in_tn = [](const auto& product_ptr,
                                       const auto& op) {
          ranges::for_each(product_ptr->factors(), [&](auto& factor) {
            auto tensor_ptr = std::dynamic_pointer_cast<AbstractTensor>(factor);
            if (tensor_ptr) {
              ranges::for_each(tensor_ptr->_braket(),
                               [&](auto& idx) { op(idx, *tensor_ptr); });
            }
          });
        };

        // compute external indices
        container::map<Index, std::size_t> indices_w_counts;
        auto retrieve_indices_with_counts = [&indices_w_counts](const auto& idx,
                                                                auto&) {
          auto found_it = indices_w_counts.find(idx);
          if (found_it != indices_w_counts.end()) {
            found_it->second++;
          } else {
            indices_w_counts.emplace(idx, 1);
          }
        };
        for_each_index_in_tn(product_ptr, retrieve_indices_with_counts);

        container::set<Index> external_indices =
            indices_w_counts | ranges::views::filter([](auto& idx_cnt) {
              auto& [idx, cnt] = idx_cnt;
              return cnt == 1;
            }) |
            ranges::views::keys | ranges::to<container::set<Index>>;

        // extract RDM-only and all indices
        container::set<Index> rdm_indices;
        std::set<Index, Index::LabelCompare> all_indices;
        auto retrieve_rdm_and_all_indices = [&rdm_indices, &all_indices,
                                             &rdm_label](const auto& idx,
                                                         const auto& tensor) {
          all_indices.insert(idx);
          if (tensor._label() == rdm_label) {
            rdm_indices.insert(idx);
          }
        };
        for_each_index_in_tn(product_ptr, retrieve_rdm_and_all_indices);

        // compute RDM->target replacement rules
        container::map<Index, Index> replacement_rules;
        ranges::for_each(rdm_indices, [&](const Index& idx) {
          const auto target_type =
              isr->intersection(idx.space(), target_rdm_space_type);
          if (target_type) {
            Index target = Index::make_tmp_index(target_type);
            replacement_rules.emplace(idx, target);
          }
        });

        if (false) {
          std::wcout << "expr = " << product_ptr->to_latex()
                     << "\n  external_indices = ";
          ranges::for_each(external_indices, [](auto& index) {
            std::wcout << index.label() << " ";
          });
          std::wcout << "\n  replrules = ";
          ranges::for_each(replacement_rules, [](auto& index) {
            std::wcout << to_latex(index.first) << "\\to"
                       << to_latex(index.second) << "\\,";
          });
          std::wcout.flush();
        }

        if (!replacement_rules.empty()) {
          sequant::detail::apply_index_replacement_rules(
              product_ptr, replacement_rules, external_indices, all_indices,
              isr);
        }
      };

      if (exptr.template is<Product>()) {
        auto product_ptr = exptr.template as_shared_ptr<Product>();
        impl_for_single_tn(product_ptr);
        exptr = product_ptr;
      } else {
        assert(exptr.template is<Sum>());
        auto result = std::make_shared<Sum>();
        for (auto& summand : exptr.template as<Sum>().summands()) {
          assert(summand.template is<Product>());
          auto result_summand = summand.template as<Product>().clone();
          auto product_ptr = result_summand.template as_shared_ptr<Product>();
          impl_for_single_tn(product_ptr);
          result->append(product_ptr);
        }
        exptr = result;
      }
    };
    project_rdm_indices_to_target(result);

    // rename dummy indices that might have been generated by
    // project_rdm_indices_to_target
    // + may combine terms

    // TensorCanonicalizer is given a custom comparer that moves active
    // indices to the front external-vs-internal trait still takes precedence
    auto current_index_comparer =
        TensorCanonicalizer::instance()->index_comparer();
    TensorCanonicalizer::instance()->index_comparer(
        [&](const Index& idx1, const Index& idx2) -> bool {
          auto active_space = isr->intersection(isr->particle_space(Spin::any),
                                                isr->hole_space(Spin::any));
          const auto idx1_active = idx1.space().type() == active_space.type();
          const auto idx2_active = idx2.space().type() == active_space.type();
          if (idx1_active) {
            if (idx2_active)
              return current_index_comparer(idx1, idx2);
            else
              return true;
          } else {
            if (idx2_active)
              return false;
            else
              return current_index_comparer(idx1, idx2);
          }
        });
    simplify(result);
    TensorCanonicalizer::instance()->index_comparer(
        std::move(current_index_comparer));

    if (Logger::get_instance().wick_stats) {
      std::wcout << "WickTheorem stats: # of contractions attempted = "
                 << wick.stats().num_attempted_contractions
                 << " # of useful contractions = "
                 << wick.stats().num_useful_contractions << std::endl;
    }
    return result;
  }
}

}  // namespace tensor
}  // namespace op

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

}  // namespace sequant::mbpt
