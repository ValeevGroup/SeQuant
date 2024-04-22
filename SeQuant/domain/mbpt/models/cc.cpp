#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/models/cc.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <cassert>
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>
#include <utility>

namespace sequant::mbpt {

CC::CC(std::size_t n, Ansatz a) : N(n), ansatz_(a) {}

CC::Ansatz CC::ansatz() const { return ansatz_; }

bool CC::unitary() const {
  return ansatz_ == Ansatz::U || ansatz_ == Ansatz::oU;
}

ExprPtr CC::sim_tr(ExprPtr expr, size_t commutator_rank) {
  const bool skip_singles = ansatz_ == Ansatz::oT || ansatz_ == Ansatz::oU;
  std::wcout << "expression: " << to_latex_align(expr, 20, 2) << std::endl;
  auto transform_op_op_pdt = [this, &commutator_rank,
                              skip_singles](const ExprPtr& expr) {
    // TODO: find the order at which the commutator expression should truncate
    // from op/op product
    assert(expr.is<op_t>() || expr.is<Product>());
    auto result = expr;
    auto op_Sk = result;
    std::wcout << "N" << N << std::endl;
    for (size_t k = 1; k <= commutator_rank; ++k) {
      ExprPtr op_Sk_comm_w_S;
      op_Sk_comm_w_S =
          op_Sk *
          op::T(N, skip_singles);  // traditional SR ansatz: [O,T] = (O T)_c
      if (unitary())  // unitary SR ansatz: [O,T-T^+] = (O T)_c + (T^+ O)_c
        op_Sk_comm_w_S += adjoint(op::T(N, skip_singles)) * op_Sk;
      std::wcout << "opskcommws presimplify "
                 << to_latex_align(op_Sk_comm_w_S, 20, 2) << std::endl;
      op_Sk = ex<Constant>(rational{1, k}) * op_Sk_comm_w_S;
      simplify(op_Sk);
      std::wcout << "op Skcomm w S " << to_latex_align(op_Sk, 20, 2)
                 << std::endl;
      result += op_Sk;
    }
    return result;
  };

  if (expr.is<op_t>()) {
    return transform_op_op_pdt(expr);
  } else if (expr.is<Product>()) {
    auto& product = expr.as<Product>();
    // Expand product as sum
    if (ranges::any_of(product.factors(), [](const auto& factor) {
          return factor.template is<Sum>();
        })) {
      expr = sequant::expand(expr);
      simplify(expr);
      return sim_tr(expr, commutator_rank);
    } else {
      return transform_op_op_pdt(expr);
    }
  } else if (expr.is<Sum>()) {
    auto result = sequant::transform_reduce(
        *expr, ex<Sum>(),
        [](const ExprPtr& running_total, const ExprPtr& summand) {
          return running_total + summand;
        },
        [=](const auto& op_product) {
          return transform_op_op_pdt(op_product);
        });
    return result;
  } else if (expr.is<Constant>() || expr.is<Variable>())
    return expr;
  else
    throw std::invalid_argument(
        "CC::sim_tr(expr): Unsupported expression type");
}

std::vector<ExprPtr> CC::t(size_t commutator_rank, size_t pmax, size_t pmin) {
  pmax = (pmax == std::numeric_limits<size_t>::max() ? N : pmax);

  assert(commutator_rank >= 1 && "commutator rank should be >= 1");
  assert(pmax >= pmin && "pmax should be >= pmin");

  // 1. construct hbar(op) in canonical form
  auto hbar = sim_tr(op::H(), commutator_rank);
  std::wcout << "runtime hbar: " << to_latex_align(hbar, 20, 2) << std::endl;

  // 2. project onto each manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> result(pmax + 1);
  for (std::int64_t p = pmax; p >= static_cast<std::int64_t>(pmin); --p) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <p|
    std::shared_ptr<Sum>
        hbar_p;  // products that can produce excitations of rank p
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p
    for (auto& term : *hbar) {
      assert(term->is<Product>() || term->is<op_t>());

      if (op::raises_vacuum_up_to_rank(term, p)) {
        if (!hbar_le_p)
          hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
        else
          hbar_le_p->append(term);
        if (op::raises_vacuum_to_rank(term, p)) {
          if (!hbar_p)
            hbar_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_p->append(term);
        }
      }
    }
    hbar = hbar_le_p;
    std::wcout << "runtime p: " << p
               << " hbar le p: " << to_latex_align(hbar_le_p, 20, 3)
               << std::endl;
    // 2.b project onto <p| (i.e., multiply by P(p) if p>0) and compute VEV
    result.at(p) = op::vac_av(p != 0 ? op::P(p) * hbar_p : hbar_p);
    std::wcout << "runtime vac av: " << to_latex_align(result.at(p), 20, 2)
               << std::endl;
  }

  return result;
}

std::vector<ExprPtr> CC::λ(std::size_t commutator_rank) {
  assert(commutator_rank >= 1 && "commutator rank should be >= 1");
  assert(!unitary() && "there is no need for CC::λ for unitary ansatz");

  // construct hbar
  auto hbar = sim_tr(op::H(), commutator_rank - 1);

  const auto One = ex<Constant>(1);
  auto lhbar = simplify((One + op::Λ(N)) * hbar);

  auto op_connect =
      op::concat(op::default_op_connections(),
                 std::vector<std::pair<mbpt::OpType, mbpt::OpType>>{
                     {OpType::h, OpType::A},
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
        hbar_p;  // products that can produce excitations of rank p
    std::shared_ptr<Sum>
        hbar_le_p;  // keeps products that can produce excitations rank <=p
    for (auto& term : *lhbar) {  // pick terms from lhbar
      assert(term->is<Product>() || term->is<op_t>());

      if (op::lowers_rank_or_lower_to_vacuum(term, p)) {
        if (!hbar_le_p)
          hbar_le_p = std::make_shared<Sum>(ExprPtrList{term});
        else
          hbar_le_p->append(term);
        if (op::lowers_rank_to_vacuum(term, p)) {
          if (!hbar_p)
            hbar_p = std::make_shared<Sum>(ExprPtrList{term});
          else
            hbar_p->append(term);
        }
      }
    }
    lhbar = hbar_le_p;

    // 2.b multiply by adjoint of P(p) (i.e., P(-p)) on the right side and
    // compute VEV
    result.at(p) = op::vac_av(hbar_p * op::P(-p), op_connect);
  }
  return result;
}

std::vector<sequant::ExprPtr> CC::t_pt(std::size_t order, std::size_t rank) {
  assert(order == 1 &&
         "sequant::mbpt::sr::CC::t_pt(): only first-order perturbation is "
         "supported now");
  assert(rank == 1 &&
         "sequant::mbpt::sr::CC::t_pt(): only one-body perturbation "
         "operator is supported now");
  assert(ansatz_ == Ansatz::T && "unitary ansatz is not yet supported");

  // construct h1_bar

  // truncate h1_bar at rank 2 for one-body perturbation
  // operator and at rank 4 for two-body perturbation operator
  const auto h1_truncate_at = rank == 1 ? 2 : 4;

  auto h1_bar = sim_tr(op::H_pt(1, rank), h1_truncate_at);

  // construct [hbar, T(1)]
  auto hbar_pert = sim_tr(op::H(), 3) * op::T_pt(order, N);

  // [Eq. 34, WIREs Comput Mol Sci. 2019; 9:e1406]
  auto expr = simplify(h1_bar + hbar_pert);

  // connectivity:
  // connect t and t1 with {h,f,g}
  // connect h1 with t
  auto op_connect =
      op::concat(op::default_op_connections(),
                 std::vector<std::pair<mbpt::OpType, mbpt::OpType>>{
                     {OpType::h, OpType::t_1},
                     {OpType::f, OpType::t_1},
                     {OpType::g, OpType::t_1},
                     {OpType::h_1, OpType::t}});

  std::vector<ExprPtr> result(N + 1);
  for (auto p = N; p >= 1; --p) {
    auto freq_term = ex<Variable>(L"ω") * op::P(p) * op::T_pt_(order, p);
    result.at(p) =
        op::vac_av(op::P(p) * expr, op_connect) - op::vac_av(freq_term);
  }
  return result;
}

std::vector<ExprPtr> CC::λ_pt(size_t order, size_t rank) {
  assert(order == 1 &&
         "sequant::mbpt::sr::CC::λ_pt(): only first-order perturbation is "
         "supported now");
  assert(rank == 1 &&
         "sequant::mbpt::sr::CC::λ_pt(): only one-body perturbation "
         "operator is supported now");
  assert(ansatz_ == Ansatz::T && "unitary ansatz is not yet supported");

  // construct hbar
  auto hbar = sim_tr(op::H(), 4);

  // construct h1_bar

  // truncate h1_bar at rank 2 for one-body perturbation
  // operator and at rank 4 for two-body perturbation operator
  const auto h1_truncate_at = rank == 1 ? 2 : 4;

  auto h1_bar = sim_tr(op::H_pt(1, rank), h1_truncate_at);
  // construct [hbar, T(1)]
  auto hbar_pert = sim_tr(op::H(), 3) * op::T_pt(order, N);

  // [Eq. 35, WIREs Comput Mol Sci. 2019; 9:e1406]
  const auto One = ex<Constant>(1);
  auto expr = simplify((One + op::Λ(N)) * (h1_bar + hbar_pert) +
                       op::Λ_pt(order, N) * hbar);

  // connectivity:
  // t and t1 with {h,f,g}
  // projectors with {h,f,g}
  // h1 with t
  // h1 with projectors
  auto op_connect =
      op::concat(op::default_op_connections(),
                 std::vector<std::pair<mbpt::OpType, mbpt::OpType>>{
                     {OpType::h, OpType::t_1},
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
    auto freq_term = ex<Variable>(L"ω") * op::Λ_pt_(order, p) * op::P(-p);
    result.at(p) =
        op::vac_av(expr * op::P(-p), op_connect) + op::vac_av(freq_term);
  }
  return result;
}

}  // namespace sequant::mbpt
