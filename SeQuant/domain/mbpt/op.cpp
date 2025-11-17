#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/math.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>

#include <stdexcept>

namespace sequant::mbpt {

std::vector<std::wstring> cardinal_tensor_labels() {
  return {L"κ",
          L"γ",
          L"Γ",
          L"A",
          L"S",
          L"P",
          L"L",
          L"λ",
          L"λ¹",
          L"h",
          L"f",
          L"f̃",
          L"g",
          L"θ",
          L"t",
          L"t¹",
          L"R",
          L"F",
          L"X",
          L"μ",
          L"V",
          L"Ṽ",
          L"B",
          L"U",
          L"GR",
          L"C",
          overlap_label(),
          kronecker_label(),
          L"a",
          L"ã",
          L"b",
          L"b̃",
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
    case OpType::θ:
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
  SEQUANT_ASSERT(a.size() == b.size());
  qns_t result;

  if (get_default_context().vacuum() == Vacuum::Physical) {
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
    throw std::runtime_error("Unsupported vacuum context.");
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

  auto op_qns = op();  // operator action i.e. quantum number change

  auto it = label2optype.find(base_lbl);  // look for OpType
  const bool known_optype = it != label2optype.end();

  // special handling for general operators
  // - Ops like f and g does not need ranks, it is implied
  // - Ops like A, S, θ are general, but need rank information
  // - θ needs to be treated differently because it can have variable number of
  // quantum numbers

  auto skip_rank_info = [](const OpType& optype) {
    return to_class(optype) == OpClass::gen &&
           !(optype == OpType::θ || optype == OpType::A || optype == OpType::S);
  };

  // batch index handling
  const auto has_batching = op.batching_ordinals();
  auto add_batch_suffix = [&op](const std::wstring& inp) {
    SEQUANT_ASSERT(op.batching_ordinals() && "Batching ordinals is not set");
    std::wstring str = inp;
    using namespace ranges::views;

    const auto ordinals = op.batching_ordinals().value();
    str += L"{[";
    str += ordinals | transform([](const auto& ord) {
             return L"{z}_{" + std::to_wstring(ord) + L"}";
           }) |
           join(L',') | ranges::to<std::wstring>();
    str += L"]}";
    return str;
  };

  if (known_optype && skip_rank_info(it->second)) {
    result += L"}";  // close the brace
    return has_batching ? add_batch_suffix(result) : result;
  }
  // specially handle θ operator
  if (known_optype && it->second == OpType::θ) {
    result += L"_{" + std::to_wstring(op_qns[0].upper()) + L"}";
    result += L"}";  // close the brace
    return has_batching ? add_batch_suffix(result) : result;
  }

  if (get_default_context().vacuum() == Vacuum::Physical) {
    if (op_qns[0] == op_qns[1]) {
      // particle conserving
      result += L"_{" + std::to_wstring(op_qns[0].lower()) + L"}";
    } else {
      // non-particle conserving
      result += L"_{" + std::to_wstring(op_qns[1].lower()) + L"}^{" +
                std::to_wstring(op_qns[0].lower()) + L"}";
    }
  } else {
    // single product vacuum
    auto nann_p =
        is_adjoint ? op_qns.ncre_particles() : op_qns.nann_particles();
    auto ncre_h = is_adjoint ? op_qns.nann_holes() : op_qns.ncre_holes();
    auto ncre_p =
        is_adjoint ? op_qns.nann_particles() : op_qns.ncre_particles();
    auto nann_h = is_adjoint ? op_qns.ncre_holes() : op_qns.nann_holes();

    if (!is_definite(nann_p) || !is_definite(ncre_h) || !is_definite(ncre_p) ||
        !is_definite(nann_h)) {
      throw std::invalid_argument(
          "to_latex(const Operator<qns_t, S>& op): "
          "can only handle generic operators with definite cre/ann numbers");
    }

    // check if the Op is a projector (A or S)
    // projectors can have negative ranks, need special handling
    [[maybe_unused]] const bool is_projector =
        known_optype && (it->second == OpType::A || it->second == OpType::S);

    // pure quasiparticle creator/annihilator?
    const auto qprank_cre = ncre_p.lower() + nann_h.lower();
    const auto qprank_ann = nann_p.lower() + ncre_h.lower();
    const auto qppure = qprank_cre == 0 || qprank_ann == 0;
    if (qppure) {
      const std::wstring baseline_char = is_adjoint ? L"^" : L"_";
      if (qprank_cre) {
        // here there is no need for sign, positive ranks of projectors
        if (ncre_p.lower() == nann_h.lower()) {  // q-particle conserving
          result +=
              baseline_char + L"{" + std::to_wstring(nann_h.lower()) + L"}";
        } else {  // particle non-conserving
          result += baseline_char + L"{" + std::to_wstring(nann_h.lower()) +
                    L"," + std::to_wstring(ncre_p.lower()) + L"}";
        }
      } else {
        // if projector, add negative sign to ranks
        const std::wstring sign = is_projector ? L"-" : L"";
        if (ncre_h.lower() == nann_p.lower()) {  // q-particle conserving
          result += baseline_char + L"{" + sign +
                    std::to_wstring(ncre_h.lower()) + L"}";
        } else {  // q-particle non-conserving
          result += baseline_char + L"{" + sign +
                    std::to_wstring(ncre_h.lower()) + L"," + sign +
                    std::to_wstring(nann_p.lower()) + L"}";
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
  return has_batching ? add_batch_suffix(result) : result;
}

}  // namespace sequant

#include <SeQuant/domain/mbpt/op.ipp>

namespace {
/// @brief Make batching indices from a vector of IndexSpaces. IndexSpaces must
/// be batching spaces. Indexing starts from 1 up to the size of \p spaces.
/// @param spaces The vector of Auxiliary IndexSpaces
/// @return A vector of Index objects
sequant::container::svector<sequant::Index> make_batch_indices(
    const sequant::container::svector<sequant::IndexSpace>& spaces) {
  using namespace sequant;
  auto validator = [](const Index& idx) {
    return idx.space().base_key() ==
           L"z";  // for now only z is allowed, i.e. batching index space
  };

  IndexFactory aux_factory{validator, 1};
  return spaces |
         ranges::views::transform([&aux_factory](const IndexSpace& space) {
           return aux_factory.make(space);
         }) |
         ranges::to<container::svector<Index>>();
}
}  // namespace

namespace sequant::mbpt {

template <Statistics S>
OpMaker<S>::OpMaker(OpType op) : op_(op) {}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, ncre nc, nann na) {
  op_ = op;
  SEQUANT_ASSERT(nc > 0 || na > 0);
  switch (to_class(op)) {
    case OpClass::ex:
      cre_spaces_ = IndexSpaceContainer(nc, get_particle_space(Spin::any));
      ann_spaces_ = IndexSpaceContainer(na, get_hole_space(Spin::any));
      break;
    case OpClass::deex:
      cre_spaces_ = IndexSpaceContainer(nc, get_hole_space(Spin::any));
      ann_spaces_ = IndexSpaceContainer(na, get_particle_space(Spin::any));
      break;
    case OpClass::gen:
      cre_spaces_ = IndexSpaceContainer(nc, get_complete_space(Spin::any));
      ann_spaces_ = IndexSpaceContainer(na, get_complete_space(Spin::any));
      break;
  }
}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, std::size_t rank)
    : OpMaker(op, ncre(rank), nann(rank)) {}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, ncre nc, nann na,
                    const cre<IndexSpace>& cre_space,
                    const ann<IndexSpace>& ann_space) {
  op_ = op;
  SEQUANT_ASSERT(nc > 0 || na > 0);
  cre_spaces_ = IndexSpaceContainer(nc, cre_space);
  ann_spaces_ = IndexSpaceContainer(na, ann_space);
}

template <Statistics S>
OpMaker<S>::OpMaker(OpType op, const OpParams& params) {
  params.validate();
  if (params.nh && params.np) {
    *this = OpMaker(op, ncre(params.np.value()), nann(params.nh.value()));
  } else {
    *this = OpMaker(op, ncre(params.rank), nann(params.rank));
  }

  // Handle batching indices if specified
  if (!params.batch_ordinals.empty()) {
    SEQUANT_ASSERT(!params.batch_ordinals.empty() &&
                   "Batch ordinals cannot be empty");
    auto isr = get_default_context().index_space_registry();
    SEQUANT_ASSERT(isr->contains(L"z") &&
                   "ISR does not contain any batching space");
    const auto batch_space = isr->retrieve(L"z");
    container::svector<Index> batch_indices;
    for (const auto& ord : params.batch_ordinals) {
      auto idx = Index(batch_space, ord);
      batch_indices.push_back(idx);
    }
    batch_indices_ = std::move(batch_indices);
  } else if (params.nbatch) {
    SEQUANT_ASSERT(params.nbatch.value() > 0 &&
                   "Number of batching indices must be > 0");
    auto isr = get_default_context().index_space_registry();
    SEQUANT_ASSERT(isr->contains(L"z") &&
                   "ISR does not contain any batching space");
    const auto batch_space = isr->retrieve(L"z");
    batch_indices_ = make_batch_indices(
        IndexSpaceContainer(params.nbatch.value(), batch_space));
  }
}

template <Statistics S>
ExprPtr OpMaker<S>::operator()(std::optional<UseDepIdx> dep,
                               std::optional<Symmetry> opsymm_opt) const {
  auto isr = get_default_context(Statistics::FermiDirac).index_space_registry();
  // if not given dep, use mbpt::Context::CSV to determine whether to use
  // dependent indices for pure (de)excitation ops
  if (!dep && get_default_mbpt_context().csv() == mbpt::CSV::Yes) {
    if (to_class(op_) == OpClass::ex) {
#ifdef SEQUANT_ASSERT_ENABLED
      for (auto&& s : cre_spaces_) {
        SEQUANT_ASSERT(isr->contains_unoccupied(s));
      }
#endif
      dep = UseDepIdx::Bra;
    } else if (to_class(op_) == OpClass::deex) {
#ifdef SEQUANT_ASSERT_ENABLED
      for (auto&& s : ann_spaces_) {
        SEQUANT_ASSERT(isr->contains_unoccupied(s));
      }
#endif
      dep = UseDepIdx::Ket;
    } else {
      dep = UseDepIdx::None;
    }
  }

  // if batching indices are given, use them
  if (batch_indices_) {
    return make(
        cre_spaces_, ann_spaces_, batch_indices_.value(),
        [this, opsymm_opt](const auto& creidxs, const auto& annidxs,
                           const auto& batchidxs, Symmetry opsymm) {
          return ex<Tensor>(to_wstring(op_), bra(creidxs), ket(annidxs),
                            aux(batchidxs), opsymm_opt ? *opsymm_opt : opsymm);
        },
        dep ? *dep : UseDepIdx::None);
  }
  // else no batching
  return make(
      cre_spaces_, ann_spaces_,
      [this, opsymm_opt](const auto& creidxs, const auto& annidxs,
                         Symmetry opsymm) {
        return ex<Tensor>(to_wstring(op_), bra(creidxs), ket(annidxs),
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
  SEQUANT_ASSERT(k > 0 && k <= 2);
  switch (k) {
    case 1:
      switch (get_default_context().vacuum()) {
        case Vacuum::Physical:
          return OpMaker<Statistics::FermiDirac>(OpType::h, 1)();
        case Vacuum::SingleProduct:
          return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
        case Vacuum::MultiProduct:
          return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
      }
      SEQUANT_UNREACHABLE;

    case 2:
      return OpMaker<Statistics::FermiDirac>(OpType::g, 2)();
  }

  SEQUANT_ABORT("Unhandled k value");
}

ExprPtr H(std::size_t k) {
  SEQUANT_ASSERT(k > 0 && k <= 2);
  return k == 1 ? tensor::H_(1) : tensor::H_(1) + tensor::H_(2);
}

ExprPtr F(bool use_tensor, IndexSpace reference_occupied) {
  if (use_tensor) {
    return OpMaker<Statistics::FermiDirac>(OpType::f, 1)();
  } else {  // explicit density matrix construction
    SEQUANT_ASSERT(
        reference_occupied);  // cannot explicitly instantiate fock operator
                              // without providing an occupied indexspace
    // add \bar{g}^{\kappa x}_{\lambda y} \gamma^y_x with x,y in occ_space_type
    auto make_g_contribution = [](const auto occ_space) {
      auto isr = get_default_context().index_space_registry();
      return mbpt::OpMaker<Statistics::FermiDirac>::make(
          {isr->complete_space(Spin::any)}, {isr->complete_space(Spin::any)},
          [=](auto braidxs, auto ketidxs, Symmetry opsymm) {
            auto m1 = Index::make_tmp_index(occ_space);
            auto m2 = Index::make_tmp_index(occ_space);
            SEQUANT_ASSERT(opsymm == Symmetry::Antisymm ||
                           opsymm == Symmetry::Nonsymm);
            if (opsymm == Symmetry::Antisymm) {
              braidxs.push_back(m1);
              ketidxs.push_back(m2);
              return ex<Tensor>(to_wstring(mbpt::OpType::g),
                                bra(std::move(braidxs)),
                                ket(std::move(ketidxs)), Symmetry::Antisymm) *
                     ex<Tensor>(to_wstring(mbpt::OpType::δ), bra{m2}, ket{m1},
                                Symmetry::Nonsymm);
            } else {  // opsymm == Symmetry::Nonsymm
              auto braidx_J = braidxs;
              braidx_J.push_back(m1);
              auto ketidxs_J = ketidxs;
              ketidxs_J.push_back(m2);
              auto braidx_K = braidxs;
              braidx_K.push_back(m1);
              auto ketidxs_K = ketidxs;
              using std::begin;
              ketidxs_K.emplace(begin(ketidxs_K), m2);
              return (ex<Tensor>(to_wstring(mbpt::OpType::g),
                                 bra(std::move(braidx_J)),
                                 ket(std::move(ketidxs_J)), Symmetry::Nonsymm) -
                      ex<Tensor>(
                          to_wstring(mbpt::OpType::g), bra(std::move(braidx_K)),
                          ket(std::move(ketidxs_K)), Symmetry::Nonsymm)) *
                     ex<Tensor>(to_wstring(mbpt::OpType::δ), bra{m2}, ket{m1},
                                Symmetry::Nonsymm);
            }
          });
    };
    auto isr = get_default_context().index_space_registry();
    return OpMaker<Statistics::FermiDirac>(OpType::h, 1)() +
           make_g_contribution(reference_occupied);
  }
}

ExprPtr θ(std::size_t K) {
  return OpMaker<Statistics::FermiDirac>(OpType::θ, K)();
}

ExprPtr T_(std::size_t K) {
  return OpMaker<Statistics::FermiDirac>(OpType::t, K)();
}

ExprPtr T(std::size_t K, bool skip1) {
  SEQUANT_ASSERT(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = skip1 ? 2ul : 1ul; k <= K; ++k) {
    result += tensor::T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  return OpMaker<Statistics::FermiDirac>(OpType::λ, K)();
}

ExprPtr Λ(std::size_t K, bool skip1) {
  SEQUANT_ASSERT(K > (skip1 ? 1 : 0));

  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
    result = k > 1 ? result + tensor::Λ_(k) : tensor::Λ_(k);
  }
  return result;
}

ExprPtr R_(nann na, ncre nc, const cre<IndexSpace>& cre_space,
           const ann<IndexSpace>& ann_space) {
  return OpMaker<Statistics::FermiDirac>(OpType::R, nc, na, cre_space,
                                         ann_space)();
}
ExprPtr R_(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(np >= 0 && nh >= 0);
  return OpMaker<Statistics::FermiDirac>(OpType::R, ncre(np.value()),
                                         nann(nh.value()))();
}

ExprPtr L_(nann na, ncre nc, const cre<IndexSpace>& cre_space,
           const ann<IndexSpace>& ann_space) {
  return OpMaker<Statistics::FermiDirac>(OpType::L, nc, na, cre_space,
                                         ann_space)();
}

ExprPtr L_(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(np >= 0 && nh >= 0);
  return OpMaker<Statistics::FermiDirac>(OpType::L, ncre(nh.value()),
                                         nann(np.value()))();
}

ExprPtr P(nₚ np, nₕ nh) {
  if (np != nh)
    SEQUANT_ASSERT(
        get_default_context().spbasis() != SPBasis::Spinfree &&
        "Spinfree basis does not support non-particle conserving projectors");
  return get_default_context().spbasis() == SPBasis::Spinfree
             ? tensor::S(-nh /* nh == np */)
             : tensor::A(-np, -nh);
}

ExprPtr A(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(!(np == 0 && nh == 0));
  // if one of them is not zero, nh and np should have the same sign
  if (np != 0 && nh != 0) {
    SEQUANT_ASSERT((np > 0 && nh > 0) || (np < 0 && nh < 0));
  }

  container::svector<IndexSpace> creators;
  container::svector<IndexSpace> annihilators;
  if (nh > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, nh))
      annihilators.emplace_back(get_hole_space(Spin::any));
  else if (nh < 0)  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -nh))
      creators.emplace_back(get_hole_space(Spin::any));
  if (np > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, np))
      creators.emplace_back(get_particle_space(Spin::any));
  else if (np < 0)  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -np))
      annihilators.emplace_back(get_particle_space(Spin::any));
  // don't populate if rank is zero

  std::optional<OpMaker<Statistics::FermiDirac>::UseDepIdx> dep;
  if (get_default_mbpt_context().csv() == mbpt::CSV::Yes)
    dep = (np > 0 || nh > 0) ? OpMaker<Statistics::FermiDirac>::UseDepIdx::Bra
                             : OpMaker<Statistics::FermiDirac>::UseDepIdx::Ket;
  return OpMaker<Statistics::FermiDirac>(
      OpType::A, cre(creators), ann(annihilators))(dep, {Symmetry::Antisymm});
}

ExprPtr S(std::int64_t K) {
  SEQUANT_ASSERT(K != 0);
  container::svector<IndexSpace> creators;
  container::svector<IndexSpace> annihilators;
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      annihilators.emplace_back(get_hole_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      creators.emplace_back(get_hole_space(Spin::any));
  if (K > 0)  // ex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, K))
      creators.emplace_back(get_particle_space(Spin::any));
  else  // deex
    for ([[maybe_unused]] auto i : ranges::views::iota(0, -K))
      annihilators.emplace_back(get_particle_space(Spin::any));

  std::optional<OpMaker<Statistics::FermiDirac>::UseDepIdx> dep;
  if (get_default_mbpt_context().csv() == mbpt::CSV::Yes)
    dep = K > 0 ? OpMaker<Statistics::FermiDirac>::UseDepIdx::Bra
                : OpMaker<Statistics::FermiDirac>::UseDepIdx::Ket;
  return OpMaker<Statistics::FermiDirac>(
      OpType::S, cre(creators), ann(annihilators))(dep, {Symmetry::Nonsymm});
}

ExprPtr H_pt(std::size_t R, std::size_t order) {
  return tensor::H_pt({.rank = R, .order = order});
}

ExprPtr H_pt(const OpParams& params) {
  params.validate();
  SEQUANT_ASSERT(params.order == 1 &&
                 "sequant::mbpt::tensor::H_pt(OpParams): only supports first "
                 "order perturbation");
  SEQUANT_ASSERT(params.rank > 0 && "H_pt requires rank > 0");
  return OpMaker<Statistics::FermiDirac>(OpType::h_1, params)();
}

ExprPtr T_pt_(std::size_t K, std::size_t order) {
  return tensor::T_pt_({.rank = K, .order = order});
}

ExprPtr T_pt_(const OpParams& params) {
  params.validate();
  SEQUANT_ASSERT(params.order == 1 &&
                 "sequant::mbpt::tensor::T_pt_(OpParams): only supports first "
                 "order perturbation");
  SEQUANT_ASSERT(params.rank > 0 && "T_pt_ requires rank > 0");
  return OpMaker<Statistics::FermiDirac>(OpType::t_1, params)();
}

ExprPtr T_pt(std::size_t K, std::size_t order, bool skip1) {
  return tensor::T_pt({.rank = K, .order = order, .skip1 = skip1});
}

ExprPtr T_pt(const OpParams& params) {
  params.validate();
  if (params.skip1) SEQUANT_ASSERT(params.rank > 1);
  ExprPtr result;
  for (auto k = (params.skip1 ? 2ul : 1ul); k <= params.rank; ++k) {
    result += tensor::T_pt_({.rank = k,
                             .nh = params.nh,
                             .np = params.np,
                             .order = params.order,
                             .nbatch = params.nbatch,
                             .batch_ordinals = params.batch_ordinals,
                             .skip1 = false});
  }
  return result;
}

ExprPtr Λ_pt_(std::size_t K, std::size_t order) {
  return tensor::Λ_pt_({.rank = K, .order = order});
}

ExprPtr Λ_pt_(const OpParams& params) {
  SEQUANT_ASSERT(params.order == 1 &&
                 "sequant::mbpt::tensor::Λ_pt_(OpParams): only supports first "
                 "order perturbation");
  SEQUANT_ASSERT(params.rank > 0 && "Λ_pt_ requires rank > 0");
  return OpMaker<Statistics::FermiDirac>(OpType::λ_1, params)();
}

ExprPtr Λ_pt(std::size_t K, std::size_t order, bool skip1) {
  return tensor::Λ_pt({.rank = K, .order = order, .skip1 = skip1});
}

ExprPtr Λ_pt(const OpParams& params) {
  params.validate();
  if (params.skip1) SEQUANT_ASSERT(params.rank > 1);
  ExprPtr result;
  for (auto k = (params.skip1 ? 2ul : 1ul); k <= params.rank; ++k) {
    result += tensor::Λ_pt_({.rank = k,
                             .nh = params.nh,
                             .np = params.np,
                             .order = params.order,
                             .nbatch = params.nbatch,
                             .batch_ordinals = params.batch_ordinals,
                             .skip1 = false});
  }
  return result;
}

}  // namespace tensor

ExprPtr H_(std::size_t k) {
  SEQUANT_ASSERT(k > 0 && k <= 2);
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
            }

            SEQUANT_UNREACHABLE;
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
  }

  SEQUANT_ABORT("Unhandled k value");
}

ExprPtr H(std::size_t k) {
  SEQUANT_ASSERT(k > 0 && k <= 2);
  return k == 1 ? H_(1) : H_(1) + H_(2);
}

ExprPtr θ(std::size_t K) {
  SEQUANT_ASSERT(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"θ"; },
                  [=]() -> ExprPtr { return tensor::θ(K); },
                  [=](qnc_t& qns) {
                    qnc_t op_qnc_t = general_type_qns(K);
                    qns = combine(op_qnc_t, qns);
                  });
}

ExprPtr T_(std::size_t K) {
  SEQUANT_ASSERT(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"t"; },
                  [=]() -> ExprPtr { return tensor::T_(K); },
                  [=](qnc_t& qns) {
                    qnc_t op_qnc_t = excitation_type_qns(K);
                    qns = combine(op_qnc_t, qns);
                  });
}

ExprPtr T(std::size_t K, bool skip1) {
  SEQUANT_ASSERT(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = skip1 ? 2ul : 1ul; k <= K; ++k) {
    result += T_(k);
  }
  return result;
}

ExprPtr Λ_(std::size_t K) {
  SEQUANT_ASSERT(K > 0);
  return ex<op_t>([]() -> std::wstring_view { return L"λ"; },
                  [=]() -> ExprPtr { return tensor::Λ_(K); },
                  [=](qnc_t& qns) {
                    qnc_t op_qnc_t = deexcitation_type_qns(K);
                    qns = combine(op_qnc_t, qns);
                  });
}

ExprPtr Λ(std::size_t K, bool skip1) {
  SEQUANT_ASSERT(K > (skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (skip1 ? 2ul : 1ul); k <= K; ++k) {
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

ExprPtr A(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(!(nh == 0 && np == 0));
  // if one of them is not zero, nh and np should have the same sign
  if (nh != 0 && np != 0) {
    SEQUANT_ASSERT((nh > 0 && np > 0) || (nh < 0 && np < 0));
  }
  // if np or nh is negative, it's a deexcitation operator
  const auto deexcitation = (np < 0 || nh < 0);

  auto particle_space = get_particle_space(Spin::any);
  auto hole_space = get_hole_space(Spin::any);
  return ex<op_t>([]() -> std::wstring_view { return L"A"; },
                  [=]() -> ExprPtr { return tensor::A(np, nh); },
                  [=](qnc_t& qns) {
                    const std::size_t abs_nh = std::abs(nh);
                    const std::size_t abs_np = std::abs(np);
                    if (deexcitation) {
                      qnc_t op_qnc_t = generic_deexcitation_qns(
                          abs_np, abs_nh, particle_space, hole_space);
                      qns = combine(op_qnc_t, qns);
                    } else {
                      qnc_t op_qnc_t = generic_excitation_qns(
                          abs_np, abs_nh, particle_space, hole_space);
                      qns = combine(op_qnc_t, qns);
                    }
                  });
}

ExprPtr S(std::int64_t K) {
  SEQUANT_ASSERT(K != 0);
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

ExprPtr P(nₚ np, nₕ nh) {
  if (get_default_context().spbasis() == SPBasis::Spinfree) {
    SEQUANT_ASSERT(
        nh == np &&
        "Only particle number conserving cases are supported with spinfree "
        "basis for now");
    const auto K = np;  // K = np = nh
    return S(-K);
  } else {
    SEQUANT_ASSERT(get_default_context().spbasis() == SPBasis::Spinor);
    return A(-np, -nh);
  }
}

namespace {
/// @brief Helper to create Operator from OpParams
/// @param label_gen Function to generate operator label
/// @param tensor_gen Function to generate tensor form
/// @param qn_action Function to apply quantum number changes
/// @param params OpParams containing operator parameters (batching, etc.)
/// @return ExprPtr to the created Operator
ExprPtr make_op_from_params(std::function<std::wstring_view()> label_gen,
                            std::function<ExprPtr()> tensor_gen,
                            std::function<void(qnc_t&)> qn_action,
                            const OpParams& params) {
  params.validate();
  if (!params.batch_ordinals.empty()) {
    return ex<op_t>(label_gen, tensor_gen, qn_action, params.batch_ordinals);
  } else if (params.nbatch) {
    return ex<op_t>(label_gen, tensor_gen, qn_action, params.nbatch.value());
  } else {
    return ex<op_t>(label_gen, tensor_gen, qn_action);
  }
}
}  // anonymous namespace

ExprPtr H_pt(std::size_t R, std::size_t order) {
  return H_pt(OpParams{.rank = R, .order = order});
}

ExprPtr H_pt(const OpParams& params) {
  SEQUANT_ASSERT(params.rank > 0);
  SEQUANT_ASSERT(params.order == 1 &&
                 "only first order perturbation is supported now");

  return make_op_from_params(
      []() -> std::wstring_view { return optype2label.at(OpType::h_1); },
      [params]() -> ExprPtr { return tensor::H_pt(params); },
      [params](qnc_t& qns) {
        qns = combine(general_type_qns(params.rank), qns);
      },
      params);
}

ExprPtr T_pt_(std::size_t K, std::size_t order) {
  return op::T_pt_({.rank = K, .order = order});
}

ExprPtr T_pt_(const OpParams& params) {
  SEQUANT_ASSERT(params.rank > 0);
  SEQUANT_ASSERT(params.order == 1 &&
                 "only first order perturbation is supported now");

  return make_op_from_params(
      []() -> std::wstring_view { return optype2label.at(OpType::t_1); },
      [params]() -> ExprPtr { return tensor::T_pt_(params); },
      [params](qnc_t& qns) {
        qns = combine(excitation_type_qns(params.rank), qns);
      },
      params);
}

ExprPtr T_pt(std::size_t K, std::size_t order, bool skip1) {
  return op::T_pt({.rank = K, .order = order, .skip1 = skip1});
}

ExprPtr T_pt(const OpParams& params) {
  params.validate();
  SEQUANT_ASSERT(params.rank > (params.skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (params.skip1 ? 2ul : 1ul); k <= params.rank; ++k) {
    result += T_pt_({.rank = k,
                     .nh = params.nh,
                     .np = params.np,
                     .order = params.order,
                     .nbatch = params.nbatch,
                     .batch_ordinals = params.batch_ordinals,
                     .skip1 = false});
  }
  return result;
}

ExprPtr Λ_pt_(std::size_t K, std::size_t order) {
  return op::Λ_pt_({.rank = K, .order = order});
}

ExprPtr Λ_pt_(const OpParams& params) {
  SEQUANT_ASSERT(params.rank > 0);
  SEQUANT_ASSERT(params.order == 1 &&
                 "only first order perturbation is supported now");

  return make_op_from_params(
      []() -> std::wstring_view { return optype2label.at(OpType::λ_1); },
      [params]() -> ExprPtr { return tensor::Λ_pt_(params); },
      [params](qnc_t& qns) {
        qns = combine(deexcitation_type_qns(params.rank), qns);
      },
      params);
}

ExprPtr Λ_pt(std::size_t K, std::size_t order, bool skip1) {
  return op::Λ_pt({.rank = K, .order = order, .skip1 = skip1});
}

ExprPtr Λ_pt(const OpParams& params) {
  params.validate();
  SEQUANT_ASSERT(params.rank > (params.skip1 ? 1 : 0));
  ExprPtr result;
  for (auto k = (params.skip1 ? 2ul : 1ul); k <= params.rank; ++k) {
    result += Λ_pt_({.rank = k,
                     .nh = params.nh,
                     .np = params.np,
                     .order = params.order,
                     .nbatch = params.nbatch,
                     .batch_ordinals = params.batch_ordinals,
                     .skip1 = false});
  }
  return result;
}

ExprPtr R_(nann na, ncre nc, const cre<IndexSpace>& cre_space,
           const ann<IndexSpace>& ann_space) {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::R); },
      [=]() -> ExprPtr { return tensor::R_(na, nc, cre_space, ann_space); },
      [=](qnc_t& qns) {
        // ex -> creators in particle_space, annihilators in hole_space
        qns = combine(
            generic_excitation_qns(/*particle_rank*/ nc, /*hole_rank*/ na,
                                   cre_space, ann_space),
            qns);
      });
}

ExprPtr R_(nₚ np, nₕ nh) { return R_(nann(nh), ncre(np)); }

ExprPtr L_(nann na, ncre nc, const cre<IndexSpace>& cre_space,
           const ann<IndexSpace>& ann_space) {
  return ex<op_t>(
      []() -> std::wstring_view { return optype2label.at(OpType::L); },
      [=]() -> ExprPtr { return tensor::L_(na, nc, cre_space, ann_space); },
      [=](qnc_t& qns) {
        // deex -> creators in hole_space, annihilators in particle_space
        qns = combine(
            generic_deexcitation_qns(
                /*particle_rank*/ na, /*hole_rank*/ nc, ann_space, cre_space),
            qns);
      });
}

ExprPtr L_(nₚ np, nₕ nh) { return L_(nann(np), ncre(nh)); }

ExprPtr R(nann na, ncre nc, const cre<IndexSpace>& cre_space,
          const ann<IndexSpace>& ann_space) {
  SEQUANT_ASSERT(na > 0 || nc > 0);
  ExprPtr result;

  std::int64_t ra = na, rc = nc;
  while (ra >= 0 && rc >= 0) {
    if (ra == 0 && rc == 0) break;
    result += R_(nann(ra), ncre(rc), cre_space, ann_space);
    if (ra == 0 || rc == 0) break;
    --ra;
    --rc;
  }
  return result;
}

ExprPtr R(nₚ np, nₕ nh) { return R(nann(nh), ncre(np)); }

ExprPtr L(nann na, ncre nc, const cre<IndexSpace>& cre_space,
          const ann<IndexSpace>& ann_space) {
  SEQUANT_ASSERT(na > 0 || nc > 0);
  ExprPtr result;

  std::int64_t ra = na, rc = nc;
  while (ra >= 0 && rc >= 0) {
    if (ra == 0 && rc == 0) break;
    result += L_(nann(ra), ncre(rc), cre_space, ann_space);
    if (ra == 0 || rc == 0) break;
    --ra;
    --rc;
  }
  return result;
}

ExprPtr L(nₚ np, nₕ nh) { return L(nann(np), ncre(nh)); }

qns_t apply_to_vac(const ExprPtr& expr) {
  SEQUANT_ASSERT(expr.is<op_t>() || expr.is<Product>());
  qns_t qns;
  if (expr.is<op_t>()) {
    qns = expr.as<op_t>()();
  } else if (expr.is<Product>()) {
    const auto& op_product = expr.as<Product>();
    for (auto& op_ptr : ranges::views::reverse(op_product.factors())) {
      SEQUANT_ASSERT(op_ptr->template is<op_t>());
      const auto& op = op_ptr->template as<op_t>();
      qns = op(qns);
    }
  }
  return qns;
}

bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t& target_qns,
                    const qns_t& source_qns) {
  qns_t qns = source_qns;
  if (op_or_op_product.is<Product>() || op_or_op_product.is<op_t>()) {
    auto qnc = apply_to_vac(op_or_op_product);
    qns = combine(qnc, qns);  // apply the operator qnc on the source qns
    return qns.overlaps_with(target_qns);
  } else
    throw std::invalid_argument(
        "sequant::mbpt::sr::contains_rank(op_or_op_product): op_or_op_product "
        "must be mbpt::sr::op_t or Product thereof");
}

bool raises_vacuum_up_to_rank(const ExprPtr& op_or_op_product,
                              const unsigned long k) {
  SEQUANT_ASSERT(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());

  return can_change_qns(op_or_op_product, interval_excitation_type_qns(k));
}

bool lowers_rank_or_lower_to_vacuum(const ExprPtr& op_or_op_product,
                                    const unsigned long k) {
  SEQUANT_ASSERT(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{},
                        interval_excitation_type_qns(k));
}

bool raises_vacuum_to_rank(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  SEQUANT_ASSERT(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, excitation_type_qns(k));
}

bool lowers_rank_to_vacuum(const ExprPtr& op_or_op_product,
                           const unsigned long k) {
  SEQUANT_ASSERT(op_or_op_product.is<op_t>() || op_or_op_product.is<Product>());
  return can_change_qns(op_or_op_product, qns_t{}, excitation_type_qns(k));
}

#include <SeQuant/domain/mbpt/vac_av.ipp>

namespace tensor {

ExprPtr detail::expectation_value_impl(
    ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
    bool use_top, bool full_contractions) {
  simplify(expr);
  auto isr = get_default_context().index_space_registry();
  const auto spinor = get_default_context().spbasis() == SPBasis::Spinor;
  // convention is to use different label for spin-orbital and spin-free RDM
  const auto rdm_label = spinor ? optype2label.at(OpType::RDM) : L"Γ";

  // N.B. reference < vacuum is not yet supported
  if (isr->reference_occupied_space().intersection(
          isr->vacuum_occupied_space()) != isr->vacuum_occupied_space()) {
    throw std::invalid_argument(
        "mbpt::tensor::expectation_value_impl: vacuum occupied orbitals must "
        "be same as or "
        "subset of the reference orbital set.");
  }

  FWickTheorem wick{expr};
  wick.use_topology(use_top).set_nop_connections(nop_connections);
  wick.full_contractions(full_contractions);
  auto result = wick.compute(/* count_only = */ false,
                             /* skip_input_canonicalization? true since already
                                did simplification above */
                             true);
  simplify(result);

  if (Logger::instance().wick_stats) {
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
    auto replace_nop_with_rdm = [&rdm_label, spinor](ExprPtr& exptr) {
      auto replace = [&rdm_label, spinor](const auto& nop) -> ExprPtr {
        using index_container = container::svector<Index>;
        auto braidxs = nop.annihilators() |
                       ranges::views::transform(
                           [](const auto& op) { return op.index(); }) |
                       ranges::to<index_container>();
        auto ketidxs = nop.creators() |
                       ranges::views::transform(
                           [](const auto& op) { return op.index(); }) |
                       ranges::to<index_container>();
        SEQUANT_ASSERT(
            braidxs.size() ==
            ketidxs.size());  // need to handle particle # violating case?
        const auto rank = braidxs.size();
        return ex<Tensor>(
            rdm_label, bra(std::move(braidxs)), ket(std::move(ketidxs)),
            rank > 1 && spinor ? Symmetry::Antisymm : Symmetry::Nonsymm);
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
            std::wcout << index.full_label() << " ";
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
              product_ptr, replacement_rules, all_indices);
        }
      };

      if (exptr.template is<Product>()) {
        auto product_ptr = exptr.template as_shared_ptr<Product>();
        impl_for_single_tn(product_ptr);
        exptr = product_ptr;
      } else {
        SEQUANT_ASSERT(exptr.template is<Sum>());
        auto result = std::make_shared<Sum>();
        for (auto& summand : exptr.template as<Sum>().summands()) {
          SEQUANT_ASSERT(summand.template is<Product>());
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

    if (Logger::instance().wick_stats) {
      std::wcout << "WickTheorem stats: # of contractions attempted = "
                 << wick.stats().num_attempted_contractions
                 << " # of useful contractions = "
                 << wick.stats().num_useful_contractions << std::endl;
    }
    return result;
  }
}

ExprPtr ref_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  auto isr = get_default_context().index_space_registry();
  const bool full_contractions =
      (isr->reference_occupied_space() == isr->vacuum_occupied_space()) ? true
                                                                        : false;
  return detail::expectation_value_impl(expr, nop_connections, use_top,
                                        full_contractions);
}

ExprPtr vac_av(ExprPtr expr, std::vector<std::pair<int, int>> nop_connections,
               bool use_top) {
  return detail::expectation_value_impl(expr, nop_connections, use_top,
                                        /* full_contractions*/ true);
}

}  // namespace tensor
}  // namespace op

bool can_change_qns(const ExprPtr& op_or_op_product, const qns_t target_qns,
                    const qns_t source_qns = {}) {
  qns_t qns = source_qns;
  if (op_or_op_product.is<Product>()) {
    const auto& op_product = op_or_op_product.as<Product>();
    for (auto& op_ptr : ranges::views::reverse(op_product.factors())) {
      SEQUANT_ASSERT(op_ptr->template is<op_t>());
      const auto& op = op_ptr->template as<op_t>();
      qns = op(qns);
    }
    return qns.overlaps_with(target_qns);
  } else if (op_or_op_product.is<op_t>()) {
    const auto& op = op_or_op_product.as<op_t>();
    qns = op(qns);
    return qns.overlaps_with(target_qns);
  } else
    throw std::invalid_argument(
        "sequant::mbpt::sr::contains_rank(op_or_op_product): op_or_op_product "
        "must be mbpt::sr::op_t or Product thereof");
}

}  // namespace sequant::mbpt
