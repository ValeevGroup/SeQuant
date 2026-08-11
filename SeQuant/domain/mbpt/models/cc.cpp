#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/bernoulli.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
// alias reserved labels for readability
const auto& asymm = sequant::reserved::antisymm_label();
const auto& symm = sequant::reserved::symm_label();

/// simple commutator of A and B: [A,B] = AB - BA
inline auto commutator(const sequant::ExprPtr& A, const sequant::ExprPtr& B) {
  auto result = A * B - B * A;
  sequant::non_canon_simplify(result);
  return result;
}
}  // namespace

namespace sequant::mbpt {

CC::CC(size_t n) : CC(n, Options{}) {}

CC::CC(size_t n, const Options& opts) : N(n), opts_(opts) {
  if (opts_.hbar_expansion == HbarExpansion::Bernoulli &&
      opts_.ansatz != Ansatz::U)
    throw Exception("CC: Bernoulli expansion requires the U ansatz");
  if (unitary() && !opts_.hbar_comm_rank)
    throw Exception("CC: hbar_comm_rank is required for unitary ansatz");
  if (opts_.ansatz == Ansatz::oT || opts_.ansatz == Ansatz::oU)
    SEQUANT_ASSERT(skip_singles(),
                   "CC: skip_singles must be true for orbital-optimized "
                   "ansatz");
}

CC::Ansatz CC::ansatz() const { return opts_.ansatz; }

bool CC::unitary() const {
  return opts_.ansatz == Ansatz::U || opts_.ansatz == Ansatz::oU;
}

std::optional<size_t> CC::hbar_comm_rank() const {
  return opts_.hbar_comm_rank;
}

CC::HbarExpansion CC::hbar_expansion() const { return opts_.hbar_expansion; }

/// resolves the `skip_singles` default: on for orbital-optimized ansätze, off
/// otherwise. Kept here rather than in the ctor so `Options` round-trips.
bool CC::skip_singles() const {
  return opts_.skip_singles.value_or(opts_.ansatz == Ansatz::oT ||
                                     opts_.ansatz == Ansatz::oU);
}

bool CC::screen() const { return opts_.screen; }

bool CC::use_topology() const { return opts_.use_topology; }

ExprPtr CC::hbar(std::optional<size_t> truncation_rank) const {
  const auto truncation =
      truncation_rank.value_or(opts_.hbar_comm_rank.value_or(4));

  if (opts_.hbar_expansion == HbarExpansion::Bernoulli)
    return bernoulli::hbar(N, truncation, skip_singles());

  // for a non-unitary ansatz this is the cheaper connected-product form, which
  // is only equivalent to the commutator once the caller supplies operator
  // connectivity to ref_av (see lst_options() and the @warning on hbar())
  auto result = mbpt::lst(H(), T(N, skip_singles()), truncation, lst_options());

  // extra singles-only commutators: wrap H̄ with the t1 similarity transform
  // to order K. Same commutator form as above, since the same connectivity is
  // supplied downstream.
  if (opts_.hbar_singles_comm_rank.value_or(0) > 0) {
    auto opts = lst_options();
    opts.skip_clone = true;
    result = mbpt::lst(result, op::t(1), *opts_.hbar_singles_comm_rank, opts);
  }
  return result;
}

ExprPtr CC::energy(std::optional<size_t> comm_rank) const {
  // Bernoulli: the hbar expansion is at tensor level, call the tensor level
  // ref_av directly. No connectivity or screening.
  if (opts_.hbar_expansion == HbarExpansion::Bernoulli) {
    const auto erank = comm_rank.value_or(opts_.hbar_comm_rank.value());
    return op::tensor::ref_av(this->hbar(erank));
  }
  // <0|H̄|0>: reference expectation value of H̄ at the requested commutator
  // truncation. No projector ⇒ this is the energy. ref_av applies the
  // connectivity (empty for unitary, default otherwise).
  return this->unitary() ? this->ref_av(this->hbar(comm_rank),
                                        mbpt::OpConnections<std::wstring>{})
                         : this->ref_av(this->hbar(comm_rank));
}

std::vector<ExprPtr> CC::t(size_t pmax, size_t pmin) const {
  pmax = (pmax == std::numeric_limits<size_t>::max() ? N : pmax);
  SEQUANT_ASSERT(pmax >= pmin, "pmax should be >= pmin");

  // Bernoulli: the hbar expansion is at tensor level, project and call the
  // tensor level ref_av directly.
  if (opts_.hbar_expansion == HbarExpansion::Bernoulli) {
    const auto hbar = this->hbar();
    std::vector<ExprPtr> result(pmax + 1);
    for (std::int64_t p = pmax; p >= static_cast<std::int64_t>(pmin); --p) {
      const auto projected = (p != 0) ? op::tensor::P(nₚ(p)) * hbar : hbar;
      result.at(p) = op::tensor::ref_av(projected);
    }
    return result;
  }

  // 1. construct hbar(op) in canonical form
  auto hbar = this->hbar();

  // connectivity: empty for unitary ansatz, default otherwise
  const auto connectivity = this->unitary()
                                ? mbpt::OpConnections<std::wstring>{}
                                : default_op_connections();

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(pmax + 1);
  for (std::int64_t p = pmax; p >= static_cast<std::int64_t>(pmin); --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <p|
    std::shared_ptr<Sum>
        hbar_for_vev;  // keeps products that can produce non-zero VEV
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p

    if (opts_.screen) {  // if operator level screening is on
      for (auto& term : *hbar) {
        SEQUANT_ASSERT(term->is<Product>() || term->is<op_t>());
        if (raises_vacuum_up_to_rank(term, p)) {
          if (!hbar_le_p)
            hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_le_p->append(term);
          if (raises_vacuum_to_rank(term, p)) {
            if (!hbar_for_vev)
              hbar_for_vev = std::make_shared<Sum>(ExprPtrList{term});
            else
              hbar_for_vev->append(term);
          }
        }
      }
      hbar = hbar_le_p;
    } else {  // no screening, use full hbar
      hbar_for_vev = hbar.is<Sum>() ? hbar.as_shared_ptr<Sum>()
                                    : std::make_shared<Sum>(hbar);
    }

    // 2.b project onto <p| (i.e., multiply by P(p) if p>0) and compute VEV
    result.at(p) = this->ref_av(p != 0 ? P(nₚ(p)) * hbar_for_vev : hbar_for_vev,
                                connectivity);
  }

  return result;
}

std::vector<ExprPtr> CC::λ() const {
  SEQUANT_ASSERT(!unitary(), "there is no need for CC::λ for unitary ansatz");

  // construct hbar
  const auto commutator_rank = opts_.hbar_comm_rank.value_or(4);
  SEQUANT_ASSERT(commutator_rank >= 1, "CC::λ: hbar_comm_rank must be >= 1");
  auto hbar = this->hbar(commutator_rank -
                         1);  // -1 because of the connection with the projector

  auto lhbar = simplify((1 + Λ(N, skip_singles())) * hbar);

  const auto op_connect = concat(default_op_connections(), {{L"h", asymm},
                                                            {L"f", asymm},
                                                            {L"g", asymm},
                                                            {L"h", symm},
                                                            {L"f", symm},
                                                            {L"g", symm}});

  std::vector<ExprPtr> result(N + 1);

  // element 0: λ pseudoenergy, computed as the CC energy with T → Λ⁺.
  {
    // connected form; the ref_av below supplies the connectivity that makes it
    // equivalent to the commutator, here {h,f,f̃,g} with λ⁺ rather than with t
    const auto hbar_λ = mbpt::lst(H(), adjoint(Λ(N, skip_singles())),
                                  commutator_rank, lst_options());
    result.at(0) = this->ref_av(
        hbar_λ, {{L"h", L"λ⁺"}, {L"f", L"λ⁺"}, {L"f̃", L"λ⁺"}, {L"g", L"λ⁺"}});
  }

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  for (auto p = N; p >= 1; --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <P|
    std::shared_ptr<Sum>
        lhbar_for_vev;  // keeps products that can produce non-zero VEV
    std::shared_ptr<Sum>
        lhbar_le_p;      // keeps products that can produce excitations rank <=p
    if (opts_.screen) {  // if operator level screening is enabled
      for (auto& term : *lhbar) {  // pick terms from lhbar
        SEQUANT_ASSERT(term->is<Product>() || term->is<op_t>());

        if (lowers_rank_or_lower_to_vacuum(term, p)) {
          if (!lhbar_le_p)
            lhbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            lhbar_le_p->append(term);
          if (lowers_rank_to_vacuum(term, p)) {
            if (!lhbar_for_vev)
              lhbar_for_vev = std::make_shared<Sum>(ExprPtrList{term});
            else
              lhbar_for_vev->append(term);
          }
        }
      }
      lhbar = lhbar_le_p;
    } else {  // no screening
      lhbar_for_vev = lhbar.is<Sum>() ? lhbar.as_shared_ptr<Sum>()
                                      : std::make_shared<Sum>(lhbar);
    }

    // 2.b multiply by adjoint of P(p) (i.e., P(-p)) on the right side and
    // compute VEV
    result.at(p) = this->ref_av(lhbar_for_vev * P(nₚ(-p)), op_connect);
  }
  return result;
}

ExprPtr CC::rdm(size_t rank, std::optional<size_t> comm_rank) const {
  SEQUANT_ASSERT(opts_.hbar_expansion != HbarExpansion::Bernoulli,
                 "CC::rdm: the Bernoulli expansion is not supported yet");

  // 1. replacement operator {ã^{p_1..p_r}_{p_{r+1}..p_{2r}}} (see op::ã); its
  // indices are free, so they become the free indices of γ.
  auto replacer = op::ã(rank);

  // 2. similarity transform e^{-σ} ã e^{σ} (σ = T for traditional and σ = T−T⁺
  // for unitary). Traditional ansatz: the expansion terminates exactly, so the
  // default is the largest number of T's that can survive the VEV below: each T
  // must consume at least one of ã's 2r legs (k <= 2r), and T's 2k
  // quasi-creators must be absorbed by ã's 2r plus Λ's 2N legs (k <= r + N).
  // Unitary ansatz: T⁺ contracts with T, the expansion never terminates, so
  // there is no safe default; use the engine's hbar_comm_rank (the ctor
  // guarantees it is set). The traditional branch takes lst_options()'s
  // connected-product form; the {ã,t} connectivity handed to ref_av below is
  // what makes it equivalent to the explicit commutator.
  const auto commutator_rank = comm_rank.value_or(
      unitary() ? opts_.hbar_comm_rank.value() : std::min(2 * rank, rank + N));
  auto bar =
      mbpt::lst(replacer, T(N, skip_singles()), commutator_rank, lst_options());

  // 3. reference expectation value with the ansatz's left wavefunction:
  // γ = <0|(1+Λ) e^{-σ} ã e^{σ}|0>; the unitary ansatz drops Λ.
  auto expr = unitary() ? bar : simplify((1 + Λ(N, skip_singles())) * bar);
  const auto connect =
      unitary() ? mbpt::OpConnections<std::wstring>{}
                : mbpt::OpConnections<std::wstring>{{L"ã", L"t"}, {L"t", L"ã"}};
  auto gamma = this->ref_av(expr, connect);

  return gamma;
}

std::vector<ExprPtr> CC::tʼ(size_t rank, size_t order,
                            std::optional<size_t> nbatch) const {
  SEQUANT_ASSERT(order == 1,
                 "sequant::mbpt::CC::tʼ(): only first-order perturbation is "
                 "supported now");
  SEQUANT_ASSERT(rank == 1,
                 "sequant::mbpt::CC::tʼ(): only one-body perturbation "
                 "operator is supported now");
  if (unitary())
    SEQUANT_ASSERT(opts_.pertbar_comm_rank,
                   "pertbar_comm_rank must be specified for unitary ansatz");
  SEQUANT_ASSERT(opts_.hbar_expansion != HbarExpansion::Bernoulli,
                 "CC::tʼ: the Bernoulli expansion is not supported yet");

  // construct h1_bar
  // truncate h1_bar at rank 2 for one-body perturbation operator and at rank 4
  // for two-body perturbation operator; unless specified otherwise
  const auto h1_truncate_default = rank == 1 ? 2 : 4;
  const auto h1_truncate_at =
      opts_.pertbar_comm_rank.value_or(h1_truncate_default);
  const auto h1_bar =
      mbpt::lst(Hʼ(rank, {.order = order, .nbatch = nbatch}),
                T(N, skip_singles()), h1_truncate_at, lst_options());

  // construct [hbar, Tʼ(1)]
  const auto hbar_truncate_at = opts_.hbar_comm_rank.value_or(
      3);  // notice 3 instead of 4 here, this is because of the commutator with
           // T'(1). In case 4 is used, it will generate more terms but they
           // will not contribute.
  const auto hbar = this->hbar(hbar_truncate_at);

  ExprPtr hbar_pert;
  if (unitary()) {
    // for unitary ansatz, we need to compute the commutator [hbar, Tʼ],
    // otherwise just hbar * Tʼ is sufficient because ref_av uses connectivity
    hbar_pert = commutator(hbar, Tʼ(N, {.order = order, .nbatch = nbatch}));
  } else {
    hbar_pert = hbar * Tʼ(N, {.order = order, .nbatch = nbatch});
  }

  // [Eq. 34, WIREs Comput Mol Sci. 2019; 9:e1406]
  const auto expr = simplify(h1_bar + hbar_pert);

  // connectivity: empty for unitary ansatz, build otherwise
  OpConnections<std::wstring> op_connect;
  if (!this->unitary()) {
    // connect t and t1 with {h,f,g}
    // connect h1 with t
    op_connect =
        concat(default_op_connections(),
               {{L"h", L"t¹"}, {L"f", L"t¹"}, {L"g", L"t¹"}, {L"h¹", L"t"}});
  }

  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    const auto freq_term =
        L"ω" * P(nₚ(p)) * op::tʼ(p, {.order = order, .nbatch = nbatch});
    result.at(p) =
        this->ref_av(P(nₚ(p)) * expr, op_connect) - this->ref_av(freq_term);
  }
  return result;
}

std::vector<ExprPtr> CC::λʼ(size_t rank, size_t order,
                            std::optional<size_t> nbatch) const {
  SEQUANT_ASSERT(order == 1,
                 "sequant::mbpt::CC::λʼ(): only first-order perturbation is "
                 "supported now");
  SEQUANT_ASSERT(rank == 1,
                 "sequant::mbpt::CC::λʼ(): only one-body perturbation "
                 "operator is supported now");
  SEQUANT_ASSERT(!unitary(), "there is no need for CC::λʼ for unitary ansatz");
  SEQUANT_ASSERT(opts_.ansatz == Ansatz::T,
                 "CC::λʼ: only traditional ansatz is supported");

  // construct hbar
  const auto hbar = this->hbar();

  // construct h1_bar
  // truncate h1_bar at rank 2 for one-body perturbation operator and at rank 4
  // for two-body perturbation operator; unless specified otherwise
  const auto h1_truncate_at = (rank == 1) ? opts_.pertbar_comm_rank.value_or(2)
                                          : opts_.pertbar_comm_rank.value_or(4);
  // connected form (this path is non-unitary, see the assert above); the
  // op_connect built below is a superset of default_op_connections() and so
  // supplies the connectivity that makes it equivalent to the commutator
  const auto h1_bar =
      mbpt::lst(Hʼ(rank, {.order = order, .nbatch = nbatch}),
                T(N, skip_singles()), h1_truncate_at, lst_options());

  // construct [hbar, T(1)]
  const auto hbar_pert =
      this->hbar(3) * Tʼ(N, {.order = order, .nbatch = nbatch});

  // [Eq. 35, WIREs Comput Mol Sci. 2019; 9:e1406]
  const auto expr = simplify((1 + Λ(N, skip_singles())) * (h1_bar + hbar_pert) +
                             Λʼ(N, {.order = order, .nbatch = nbatch}) * hbar);

  // connectivity:
  // t and t1 with {h,f,g}
  // projectors with {h,f,g}
  // h1 with t
  // h1 with projectors
  const auto op_connect = concat(default_op_connections(), {{L"h", L"t¹"},
                                                            {L"f", L"t¹"},
                                                            {L"g", L"t¹"},
                                                            {L"h¹", L"t"},
                                                            {L"h", asymm},
                                                            {L"f", asymm},
                                                            {L"g", asymm},
                                                            {L"h", symm},
                                                            {L"f", symm},
                                                            {L"g", symm},
                                                            {L"h¹", asymm},
                                                            {L"h¹", symm}});

  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    const auto freq_term =
        L"ω" * op::λʼ(p, {.order = order, .nbatch = nbatch}) * P(nₚ(-p));
    result.at(p) =
        this->ref_av(expr * P(nₚ(-p)), op_connect) + this->ref_av(freq_term);
  }
  return result;
}

namespace {
// EOM eigenvector operators R and L use SquareRoot normalization
constexpr Normalization eom_norm = Normalization::SquareRoot;

/// Projection manifolds for an EOM operator with @p np particle and @p nh hole
/// counts, lowest rank first. Descends (np, nh) together, stopping before the
/// reference and, for IP/EA, once either count reaches zero.
container::svector<std::pair<std::int64_t, std::int64_t>> eom_manifolds(nₚ np,
                                                                        nₕ nh) {
  container::svector<std::pair<std::int64_t, std::int64_t>> manifolds;
  for (std::int64_t rp = np, rh = nh; rp >= 0 && rh >= 0; --rp, --rh) {
    if (rp == 0 && rh == 0) break;
    manifolds.emplace_back(rp, rh);
    if (rp == 0 || rh == 0) break;
  }
  std::ranges::reverse(manifolds);
  return manifolds;
}
}  // namespace

// UCC EOM sigma equations. For the qUCCSD ranks see
// 10.1063/5.0062090 Sec. II C, Eqs. (29)-(48); for IP/EA,
// 10.1021/acs.jctc.5c01991 Fig. 1.
//
// Eq. (10) splits the single physical E_gr from the normal-ordered components
// that build the blocks of Eq. (7). In projected-H̄ UCC assembly, each
// cumulative H̄^(k) still carries its rank-dependent scalar part, so remove
// that same scalar on the diagonal.
std::vector<ExprPtr> CC::assemble_ucc_eom(
    nₚ np, nₕ nh, const std::vector<size_t>& block_ranks,
    UCCEOMAssembly assembly) const {
  if (assembly == UCCEOMAssembly::Commutator &&
      opts_.hbar_expansion == HbarExpansion::Bernoulli)
    throw Exception("CC::eom_r: Bernoulli requires projected Hbar assembly");

  using std::min;
  if (assembly == UCCEOMAssembly::Commutator && block_ranks.empty()) {
    const auto hbar_R = commutator(hbar(), R(np, nh, eom_norm));
    std::vector<ExprPtr> result(min(np, nh) + 1);
    for (const auto& [rp, rh] : eom_manifolds(np, nh))
      result.at(min(rp, rh)) = ref_av(δl(nₚ(rp), nₕ(rh)) * hbar_R, {});
    return result;
  }

  const auto manifolds = eom_manifolds(np, nh);
  const auto K = manifolds.size();
  // `block_ranks` is read at i * K + j, so the ascending order above is what
  // makes `{2,1,1,0}` mean SS, SD, DS, DD; empty means the configured H̄ rank,
  // or the fourth commutator when no rank is configured.
  const std::vector<size_t> ranks =
      block_ranks.empty()
          ? std::vector<size_t>(K * K, hbar_comm_rank().value_or(4))
          : block_ranks;
  if (ranks.size() != K * K)
    throw Exception(
        "CC::eom_r: block_ranks must be a K x K row-major matrix, "
        "K = number of projection manifolds");

  // Bernoulli H̄ is tensor-level, BCH H̄ operator-level; the bra/ket/vev trio
  // below must match it. Connectivity is empty either way, as everywhere on the
  // unitary path; only the operator-level vev forwards screen/use_topology.
  const bool tensor_level = opts_.hbar_expansion == HbarExpansion::Bernoulli;

  // One H̄ per distinct truncation order, reduced to its R part (Bernoulli only;
  // BCH H̄ is operator-level and has no N/R split). The N part is the amplitude
  // residual <Φl|H̄|Φ0>, which Eq. (6) zeroes only at the amplitude rank and
  // Eqs. (41)-(47) carry nowhere. It shifts the manifold rank, so removing it
  // leaves diagonal blocks untouched.
  container::map<size_t, ExprPtr> hbars;
  for (const auto k : ranks) {
    auto [it, fresh] = hbars.try_emplace(k);
    if (!fresh) continue;  // deriving H̄ twice for one rank is not cheap
    it->second = hbar(k);
    if (tensor_level)
      it->second =
          bernoulli::detail::R_part(it->second, N, skip_singles() ? 2 : 1);
  }
  auto bra_of = [tensor_level](std::int64_t p, std::int64_t h) {
    return tensor_level ? op::tensor::δl(nₚ(p), nₕ(h)) : op::δl(nₚ(p), nₕ(h));
  };
  auto ket_of = [tensor_level](std::int64_t p, std::int64_t h) {
    return tensor_level ? op::tensor::r(nₚ(p), nₕ(h), eom_norm)
                        : op::r(nₚ(p), nₕ(h), eom_norm);
  };
  auto vev = [tensor_level, this](const ExprPtr& e) {
    return tensor_level ? op::tensor::ref_av(e)
                        : op::ref_av(e, {.connect = {},
                                         .screen = opts_.screen,
                                         .use_topology = opts_.use_topology});
  };

  std::vector<ExprPtr> result(min(np, nh) + 1);
  for (size_t i = 0; i < K; ++i) {
    const auto [bp, bh] = manifolds[i];
    const auto bra = bra_of(bp, bh);
    auto acc = std::make_shared<Sum>();
    for (size_t j = 0; j < K; ++j) {
      const auto [kp, kh] = manifolds[j];
      const auto& hbar_ij = hbars.at(ranks.at(i * K + j));
      const auto ket = ket_of(kp, kh);
      if (assembly == UCCEOMAssembly::Commutator) {
        acc->append(vev(bra * commutator(hbar_ij, ket)));
      } else {
        acc->append(vev(bra * hbar_ij * ket));
        // Remove the scalar part of this block's temporary H̄^(k_ii). This is
        // not a block-dependent physical E_gr: it leaves the normal-ordered
        // coefficients selected for this block in Eq. (10). Written as
        // <i|r_i H̄|0> so Wick keeps its summed indices disjoint from the
        // block's external ones.
        if (i == j) acc->append(ex<Constant>(-1) * vev(bra * ket * hbar_ij));
      }
    }
    result.at(static_cast<size_t>(min(bp, bh))) = simplify(ExprPtr{acc});
  }
  return result;
}

std::vector<ExprPtr> CC::eom_r(nₚ np, nₕ nh,
                               const std::vector<size_t>& block_ranks,
                               std::optional<UCCEOMAssembly> assembly) const {
  SEQUANT_ASSERT(np > 0 || nh > 0, "Unsupported excitation order");
  if (np != nh)
    SEQUANT_ASSERT(
        get_default_context().spbasis() != SPBasis::Spinfree,
        "spin-free basis does not yet support non particle-conserving cases");

  const auto selected_assembly = assembly.value_or(
      (!block_ranks.empty() || opts_.hbar_expansion == HbarExpansion::Bernoulli)
          ? UCCEOMAssembly::ProjectedHbar
          : UCCEOMAssembly::Commutator);
  if (unitary())
    return assemble_ucc_eom(np, nh, block_ranks, selected_assembly);

  if (!block_ranks.empty())
    throw Exception("CC::eom_r: block_ranks require a unitary ansatz");
  if (selected_assembly == UCCEOMAssembly::ProjectedHbar)
    throw Exception(
        "CC::eom_r: projected Hbar assembly requires a unitary ansatz");

  const auto hbar = this->hbar();
  const auto hbar_R = hbar * R(np, nh, eom_norm);
  const auto op_connect = concat(default_op_connections(),
                                 {{L"h", L"R"}, {L"f", L"R"}, {L"g", L"R"}});

  using std::min;
  std::vector<ExprPtr> result(min(np, nh) + 1);
  for (const auto& [rp, rh] : eom_manifolds(np, nh))
    result.at(min(rp, rh)) = ref_av(δl(nₚ(rp), nₕ(rh)) * hbar_R, op_connect);

  return result;
}

std::vector<ExprPtr> CC::eom_l(nₚ np, nₕ nh) const {
  SEQUANT_ASSERT(!unitary(),
                 "there is no need for CC::eom_l for unitary ansatz");
  SEQUANT_ASSERT(np > 0 || nh > 0, "Unsupported excitation order");

  if (np != nh)
    SEQUANT_ASSERT(
        get_default_context().spbasis() != SPBasis::Spinfree &&
        "spin-free basis does not support non particle-conserving cases");

  // construct hbar
  const auto hbar = this->hbar();

  // L * hbar
  const auto L_hbar = L(np, nh, eom_norm) * hbar;

  // connectivity:
  // default connections + connect H with projectors
  const auto op_connect = concat(default_op_connections(), {{L"h", asymm},
                                                            {L"f", asymm},
                                                            {L"g", asymm},
                                                            {L"h", symm},
                                                            {L"f", symm},
                                                            {L"g", symm}});

  using std::min;
  std::vector<ExprPtr> result(min(np, nh) + 1);  // for EE element 0 stays null
  // right project with |rp,rh> (i.e., multiply δr(rp, rh)) and compute VEV
  for (const auto& [rp, rh] : eom_manifolds(np, nh))
    result.at(min(rp, rh)) =
        this->ref_av(L_hbar * δr(nₚ(rp), nₕ(rh)), op_connect);

  return result;
}
}  // namespace sequant::mbpt
