#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>

#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>

namespace sequant::mbpt {

CC::CC(size_t n, Ansatz a, bool screen, bool use_topology)
    : N(n), ansatz_(a), screen_(screen), use_topology_(use_topology) {}

CC::Ansatz CC::ansatz() const { return ansatz_; }

bool CC::unitary() const {
  return ansatz_ == Ansatz::U || ansatz_ == Ansatz::oU;
}

bool CC::screen() const { return screen_; }

bool CC::use_topology() const { return use_topology_; }

std::vector<ExprPtr> CC::t(size_t commutator_rank, size_t pmax, size_t pmin) {
  pmax = (pmax == std::numeric_limits<size_t>::max() ? N : pmax);
  const bool skip_singles = ansatz_ == Ansatz::oT || ansatz_ == Ansatz::oU;

  SEQUANT_ASSERT(pmax >= pmin && "pmax should be >= pmin");

  // 1. construct hbar(op) in canonical form
  auto hbar = mbpt::lst(H(), T(N, skip_singles), commutator_rank, unitary());

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(pmax + 1);
  for (std::int64_t p = pmax; p >= static_cast<std::int64_t>(pmin); --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <p|
    std::shared_ptr<Sum>
        hbar_for_vev;  // keeps products that can produce non-zero VEV
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p

    if (screen_) {  // if operator level screening is on
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
    result.at(p) =
        this->ref_av(p != 0 ? P(nₚ(p)) * hbar_for_vev : hbar_for_vev);
  }

  return result;
}

std::vector<ExprPtr> CC::λ(size_t commutator_rank) {
  SEQUANT_ASSERT(commutator_rank >= 1 && "commutator rank should be >= 1");
  SEQUANT_ASSERT(!unitary() && "there is no need for CC::λ for unitary ansatz");
  const bool skip_singles = ansatz_ == Ansatz::oT || ansatz_ == Ansatz::oU;

  // construct hbar
  auto hbar =
      mbpt::lst(H(), T(N, skip_singles), commutator_rank - 1, unitary());

  const auto One = ex<Constant>(1);
  auto lhbar = simplify((One + Λ(N)) * hbar);

  const auto op_connect = concat(default_op_connections(),
                                 OpConnections<OpType>{{OpType::h, OpType::A},
                                                       {OpType::f, OpType::A},
                                                       {OpType::g, OpType::A},
                                                       {OpType::h, OpType::S},
                                                       {OpType::f, OpType::S},
                                                       {OpType::g, OpType::S}});

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <P|
    std::shared_ptr<Sum>
        lhbar_for_vev;  // keeps products that can produce non-zero VEV
    std::shared_ptr<Sum>
        lhbar_le_p;  // keeps products that can produce excitations rank <=p
    if (screen_) {   // if operator level screening is enabled
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

std::vector<ExprPtr> CC::t_pt(size_t rank, [[maybe_unused]] size_t order,
                              std::optional<size_t> nbatch) {
  SEQUANT_ASSERT(order == 1 &&
                 "sequant::mbpt::CC::t_pt(): only first-order perturbation is "
                 "supported now");
  SEQUANT_ASSERT(rank == 1 &&
                 "sequant::mbpt::CC::t_pt(): only one-body perturbation "
                 "operator is supported now");
  SEQUANT_ASSERT(ansatz_ == Ansatz::T && "unitary ansatz is not yet supported");

  // OpParams for operators
  OpParams hpt_params{.rank = rank, .order = order};
  OpParams tpt_params{.rank = N, .order = order};
  if (nbatch) {
    hpt_params.nbatch = nbatch.value();
    tpt_params.nbatch = nbatch.value();
  }

  // construct h1_bar
  // truncate h1_bar at rank 2 for one-body perturbation
  // operator and at rank 4 for two-body perturbation operator
  const auto h1_truncate_at = rank == 1 ? 2 : 4;
  const auto h1_bar = mbpt::lst(H_pt(hpt_params), T(N), h1_truncate_at);

  // construct [hbar, T(1)]
  const auto hbar_pert = mbpt::lst(H(), T(N), 3) * T_pt(tpt_params);

  // [Eq. 34, WIREs Comput Mol Sci. 2019; 9:e1406]
  const auto expr = simplify(h1_bar + hbar_pert);

  // connectivity:
  // connect t and t1 with {h,f,g}
  // connect h1 with t
  const auto op_connect =
      concat(default_op_connections(),
             OpConnections<OpType>{{OpType::h, OpType::t_1},
                                   {OpType::f, OpType::t_1},
                                   {OpType::g, OpType::t_1},
                                   {OpType::h_1, OpType::t}});

  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    OpParams tp_params{.rank = p, .order = order};
    if (nbatch) {
      tp_params.nbatch = nbatch.value();
    }
    const auto freq_term = ex<Variable>(L"ω") * P(nₚ(p)) * T_pt_(tp_params);
    result.at(p) =
        this->ref_av(P(nₚ(p)) * expr, op_connect) - this->ref_av(freq_term);
  }
  return result;
}

std::vector<ExprPtr> CC::λ_pt(size_t rank, [[maybe_unused]] size_t order,
                              std::optional<size_t> nbatch) {
  SEQUANT_ASSERT(order == 1 &&
                 "sequant::mbpt::CC::λ_pt(): only first-order perturbation is "
                 "supported now");
  SEQUANT_ASSERT(rank == 1 &&
                 "sequant::mbpt::CC::λ_pt(): only one-body perturbation "
                 "operator is supported now");
  SEQUANT_ASSERT(ansatz_ == Ansatz::T && "unitary ansatz is not yet supported");

  // OpParams for operators
  OpParams hpt_params{.rank = rank, .order = order};
  OpParams tpt_params{.rank = N, .order = order};
  OpParams lpt_params{.rank = N, .order = order};
  if (nbatch) {
    hpt_params.nbatch = nbatch.value();
    tpt_params.nbatch = nbatch.value();
    lpt_params.nbatch = nbatch.value();
  }

  // construct hbar
  const auto hbar = mbpt::lst(H(), T(N), 4);

  // construct h1_bar
  // truncate h1_bar at rank 2 for one-body perturbation
  // operator and at rank 4 for two-body perturbation operator
  const auto h1_truncate_at = rank == 1 ? 2 : 4;
  const auto h1_bar = mbpt::lst(H_pt(hpt_params), T(N), h1_truncate_at);

  // construct [hbar, T(1)]
  const auto hbar_pert = mbpt::lst(H(), T(N), 3) * T_pt(tpt_params);

  // [Eq. 35, WIREs Comput Mol Sci. 2019; 9:e1406]
  const auto One = ex<Constant>(1);
  const auto expr =
      simplify((One + Λ(N)) * (h1_bar + hbar_pert) + Λ_pt(lpt_params) * hbar);

  // connectivity:
  // t and t1 with {h,f,g}
  // projectors with {h,f,g}
  // h1 with t
  // h1 with projectors
  const auto op_connect =
      concat(default_op_connections(),
             OpConnections<OpType>{{OpType::h, OpType::t_1},
                                   {OpType::f, OpType::t_1},
                                   {OpType::g, OpType::t_1},
                                   {OpType::h_1, OpType::t},
                                   {OpType::h, OpType::A},
                                   {OpType::f, OpType::A},
                                   {OpType::g, OpType::A},
                                   {OpType::h, OpType::S},
                                   {OpType::f, OpType::S},
                                   {OpType::g, OpType::S},
                                   {OpType::h_1, OpType::A},
                                   {OpType::h_1, OpType::S}});

  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    OpParams lp_params{.rank = p, .order = order};
    if (nbatch) {
      lp_params.nbatch = nbatch.value();
    }
    const auto freq_term = ex<Variable>(L"ω") * Λ_pt_(lp_params) * P(nₚ(-p));
    result.at(p) =
        this->ref_av(expr * P(nₚ(-p)), op_connect) + this->ref_av(freq_term);
  }
  return result;
}

std::vector<ExprPtr> CC::eom_r(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(!unitary() && "Unitary ansatz is not yet supported");
  SEQUANT_ASSERT((np > 0 || nh > 0) && "Unsupported excitation order");

  if (np != nh)
    SEQUANT_ASSERT(
        get_default_context().spbasis() != SPBasis::Spinfree &&
        "spin-free basis does not yet support non particle-conserving cases");
  const bool skip_singles = ansatz_ == Ansatz::oT;

  // construct hbar
  const auto hbar = mbpt::lst(H(), T(N, skip_singles), 4);

  // hbar * R
  const auto hbar_R = hbar * R(np, nh);

  // connectivity:
  // default connections + connect R with {h,f,g}
  const auto op_connect = concat(default_op_connections(),
                                 OpConnections<OpType>{{OpType::h, OpType::R},
                                                       {OpType::f, OpType::R},
                                                       {OpType::g, OpType::R}});

  // initialize result vector
  std::vector<ExprPtr> result;
  using std::min;
  result.resize(min(np, nh) + 1);  // for EE first element will be empty

  std::int64_t rp = np, rh = nh;
  while (rp >= 0 && rh >= 0) {
    if (rp == 0 && rh == 0) break;
    // project with <rp, rh| (i.e., multiply P(rp, rh)) and compute VEV
    result.at(min(rp, rh)) =
        this->ref_av(P(nₚ(rp), nₕ(rh)) * hbar_R, op_connect);
    if (rp == 0 || rh == 0) break;
    rp--;
    rh--;
  }

  return result;
}

std::vector<ExprPtr> CC::eom_l(nₚ np, nₕ nh) {
  SEQUANT_ASSERT(!unitary() && "Unitary ansatz is not yet supported");
  SEQUANT_ASSERT((np > 0 || nh > 0) && "Unsupported excitation order");

  if (np != nh)
    SEQUANT_ASSERT(
        get_default_context().spbasis() != SPBasis::Spinfree &&
        "spin-free basis does not support non particle-conserving cases");
  const bool skip_singles = ansatz_ == Ansatz::oT;

  // construct hbar
  const auto hbar = mbpt::lst(H(), T(N, skip_singles), 4);

  // L * hbar
  const auto L_hbar = L(np, nh) * hbar;

  // connectivity:
  // default connections + connect H with projectors
  const auto op_connect = concat(default_op_connections(),
                                 OpConnections<OpType>{{OpType::h, OpType::A},
                                                       {OpType::f, OpType::A},
                                                       {OpType::g, OpType::A},
                                                       {OpType::h, OpType::S},
                                                       {OpType::f, OpType::S},
                                                       {OpType::g, OpType::S}});

  // initialize result vector
  std::vector<ExprPtr> result;
  using std::min;
  result.resize(min(nh, np) + 1);  // for EE first element will be empty

  std::int64_t rp = np, rh = nh;
  while (rp >= 0 && rh >= 0) {
    if (rp == 0 && rh == 0) break;
    // right project with |rp,rh> (i.e., multiply P(-rp, -rh)) and compute VEV
    result.at(min(rp, rh)) =
        this->ref_av(L_hbar * P(nₚ(-rp), nₕ(-rh)), op_connect);
    if (rp == 0 || rh == 0) break;
    rp--;
    rh--;
  }

  return result;
}
}  // namespace sequant::mbpt
