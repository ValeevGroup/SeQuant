#include <SeQuant/domain/mbpt/models/hf.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index_basis_registry.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/context.hpp>
#include <SeQuant/domain/mbpt/op.hpp>
#include <SeQuant/domain/mbpt/op_registry.hpp>
#include <SeQuant/domain/mbpt/series.hpp>
#include <SeQuant/domain/mbpt/utils.hpp>
#include <SeQuant/domain/mbpt/vac_av.hpp>

#include <cstddef>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>

namespace sequant::mbpt::hf {

namespace {

constexpr std::wstring_view rotation_label = L"κ̂", h_label = L"h";

[[noreturn]] void reject(std::string_view function, const std::string& text) {
  throw Exception("mbpt::hf::" + std::string(function) + ": " + text);
}

/// the reference-occupied space and its complement, after the checks every
/// function here makes
std::pair<IndexSpace, IndexSpace> spaces(std::string_view function) {
  if (!get_default_mbpt_context().op_registry()->contains(
          std::wstring(rotation_label)))
    reject(function, "the operator label κ̂ is not registered");
  if (get_default_context().vacuum() != Vacuum::SingleProduct ||
      !detail::reference_is_vacuum())
    reject(function,
           "the vacuum must be SingleProduct and the reference the vacuum");
  if (get_default_mbpt_context().csv() == CSV::Yes)
    reject(function,
           "orbital rotations with pair-specific orbitals are not supported");
  // the reference being the vacuum, the vacuum-unoccupied space is the
  // complement of the reference-occupied one
  const auto ibr = get_default_context().index_basis_registry();
  const auto occupied = ibr->reference_occupied_space(Spin::any);
  try {
    return {occupied, ibr->vacuum_unoccupied_space(Spin::any)};
  } catch (const Exception&) {
    reject(function, "the unoccupied space is not registered");
  }
}

ExprPtr rotated(std::size_t order, const std::pair<IndexSpace, IndexSpace>& s) {
  const auto K = ex<op_t>(
      [] { return rotation_label; },
      [s] {
        return OpMaker<Statistics::FermiDirac>(std::wstring(rotation_label),
                                               ncre(1), nann(1), cre(s.second),
                                               ann(s.first))();
      },
      [s](qnc_t& qns) {
        qns = combine(generic_excitation_qns(1, 1, s.second, s.first), qns);
      });
  return simplify(op::ref_av(lst(op::H(), K - adjoint(K), order)));
}

/// the integral tensor of the @p k -body part of op::H()
Tensor integral(std::size_t k) {
  ExprPtr result;
  [[maybe_unused]] std::size_t found = 0;
  tensor::h(k)->visit(
      [&](const ExprPtr& e) {
        if (!e->is<Tensor>()) return;
        result = e;
        ++found;
      },
      /*atoms_only=*/true);
  SEQUANT_ASSERT(found == 1, "the k-body part of op::H() must hold one tensor");
  return result->as<Tensor>();
}

/// a tensor labelled @p label over @p b and @p k , with the auxiliary
/// indices and the symmetries of @p model
ExprPtr like(const Tensor& model, std::wstring_view label,
             container::svector<Index> b, container::svector<Index> k) {
  return detail::tensor_like(model, label, std::move(b), std::move(k),
                             {model.aux().begin(), model.aux().end()});
}

}  // namespace

ExprPtr rotated_energy(std::size_t order) {
  return rotated(order, spaces("rotated_energy"));
}

ExprPtr energy(std::size_t order) {
  const auto s = spaces("energy");
  if (get_default_context().spbasis() != SPBasis::Spinor)
    reject("energy", "the basis must be the spin-orbital one");
  const auto f = integral(1), g = integral(2);
  // h over the slots of t, with its symmetries
  const auto h = [](const Tensor& t) {
    return like(t, h_label, {t.bra().begin(), t.bra().end()},
                {t.ket().begin(), t.ket().end()});
  };
  // ⟨0|H|0⟩, H written in plain creators and annihilators
  const auto& p = g.bra();
  const auto& q = g.ket();
  const auto E0 = op::tensor::ref_av(
      h(f) * fcrex(f.bra()[0]) * fannx(f.ket()[0]) +
      ex<Constant>(rational{1, 4}) * g.clone() * fcrex(p[0]) * fcrex(p[1]) *
          fannx(q[1]) * fannx(q[0]));
  // f{p;q} = h{p;q} + Σ_i g{p,i;q,i}
  const auto E = rotated(order, s);
  const auto fresh = detail::fresh_indices(E);
  return simplify(
      E0 + detail::substituted(E, [&](const Tensor& t) -> ExprPtr {
        if (t.label() != f.label()) return nullptr;
        const auto i = fresh(s.first);
        return h(t) + like(g, g.label(), {t.bra()[0], i}, {t.ket()[0], i});
      }));
}

ExprPtr orthonormalized(const ExprPtr& expr, const ExprPtr& strength,
                        std::size_t order) {
  constexpr std::string_view function = "orthonormalized";
  const auto f = integral(1), g = integral(2);
  const auto complete =
      get_default_context().index_basis_registry()->complete_space(Spin::any);
  // a null expr is check_series()'s to reject
  if (expr)
    detail::for_each_atom(expr, [&](const ExprPtr& leaf) {
      if (detail::is_operator(*leaf))
        reject(function,
               "the expression holds operators; expand its tensor form");
      const auto& b = detail::base_of(leaf);
      if (!b->is<Tensor>()) return;
      const auto& t = b->as<Tensor>();
      const auto l = strip_adjoint_label(t.label());
      // under some vacua op::H()'s one-body tensor is h itself
      if (l == f.label() && l != h_label)
        reject(function, "the expression holds f; write it in h and g");
      if (l != h_label && l != g.label()) return;
      if (l != t.label()) reject(function, "the expression holds h⁺ or g⁺");
      for (const auto& leg : t.braket_indices())
        if (leg.has_proto_indices() || leg.space().qns() != complete.qns())
          reject(function, "a leg of h or g carries protoindices or a spin");
    });
  const auto& overlap = reserved::overlap_label();
  for (const auto l : {h_label, g.label()})
    detail::check_series(expr, l, strength, order);
  detail::check_series(expr, overlap, strength, order, /*expanded=*/false);
  const auto fresh = detail::fresh_indices(expr);

  // the term of order k ≥ 1 of X from a to b: C(−1/2,m) Δ_{j_1} ⋯ Δ_{j_m}
  // over the compositions j_1 + … + j_m = k, Δ_j the term of order j of S
  const auto X = [&](std::size_t k, const Index& a, const Index& b) {
    container::svector<ExprPtr> terms;
    for (std::size_t cuts = 0; cuts != std::size_t{1} << (k - 1); ++cuts) {
      container::svector<ExprPtr> factors;
      auto from = a;
      for (std::size_t n = 1, j = 1; n <= k; ++n, ++j)
        if (n == k || (cuts >> (n - 1) & 1)) {
          const auto to = n == k ? b : fresh(complete);
          const Tensor s(overlap, bra{from}, ket{to}, Symmetry::Nonsymm,
                         f.braket_symmetry(), f.column_symmetry());
          factors.push_back(
              detail::series_term(s, overlap, strength, j, fresh));
          from = to;
          j = 0;
        }
      rational binomial = 1;  // C(−1/2, m)
      for (std::size_t i = 0; i != factors.size(); ++i)
        binomial *= (rational{-1, 2} - i) / (i + 1);
      terms.push_back(ex<Product>(binomial, std::move(factors)));
    }
    return ex<Sum>(std::move(terms));
  };

  return detail::substituted(expr, [&](const Tensor& t) -> ExprPtr {
    if (t.label() != h_label && t.label() != g.label()) return nullptr;
    // the orders of the legs, bra then ket, and of the leaf
    const auto legs = t.bra_rank() + t.ket_rank();
    container::svector<std::size_t> k(legs + 1, 0);
    container::svector<ExprPtr> terms;
    while (true) {
      if (std::accumulate(k.begin(), k.end(), std::size_t{0}) <= order) {
        container::svector<Index> b(t.bra().begin(), t.bra().end()),
            c(t.ket().begin(), t.ket().end());
        container::svector<ExprPtr> factors;
        for (std::size_t j = 0; j != legs; ++j) {
          if (k[j] == 0) continue;
          const auto n = fresh(complete);
          auto& leg = j < b.size() ? b[j] : c[j - b.size()];
          factors.push_back(j < b.size() ? X(k[j], leg, n) : X(k[j], n, leg));
          leg = n;
        }
        factors.push_back(
            detail::series_term(like(t, t.label(), b, c)->as<Tensor>(),
                                t.label(), strength, k[legs], fresh));
        terms.push_back(ex<Product>(std::move(factors)));
      }
      std::size_t j = 0;
      while (j != k.size() && k[j] == order) k[j++] = 0;
      if (j == k.size()) break;
      ++k[j];
    }
    return ex<Sum>(std::move(terms));
  });
}

}  // namespace sequant::mbpt::hf
