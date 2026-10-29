#include <SeQuant/domain/mbpt/rules/pair_basis.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/domain/mbpt/rules/pno.hpp>

#include <range/v3/algorithm/equal.hpp>

#include <optional>

namespace sequant::mbpt {

namespace {

/// a one-bra, one-ket tensor with the given label and aux rank, no
/// K-conjugation mark
bool is_rank11(Tensor const& t, std::wstring const& label,
               std::size_t aux_rank) {
  return t.label() == label && t.bra_rank() == 1 && t.ket_rank() == 1 &&
         t.aux_rank() == aux_rank && !t.kconjugated();
}

/// @p factors[pos] as a Tensor, or null if it is not one
Tensor const* tensor_at(container::svector<ExprPtr> const& factors,
                        std::size_t pos) {
  return factors[pos]->is<Tensor>() ? &factors[pos]->as<Tensor>() : nullptr;
}

ExprPtr fold_expr(ExprPtr const& x, PairBasisFoldOptions const& opts,
                  PairBasisFoldReport& report);

/// @p prod with its pair-basis triples folded; @p x itself if none
ExprPtr fold_product(ExprPtr const& x, PairBasisFoldOptions const& opts,
                     PairBasisFoldReport& report) {
  auto const& prod = x->as<Product>();
  // the integral and its coefficients must be siblings: flatten nested
  // products first, and fold inside sum factors on their own
  Product flat(prod.scalar(), prod.factors(), Product::Flatten::Recursively);
  container::svector<ExprPtr> factors(flat.begin(), flat.end());
  bool changed = false;
  for (auto& f : factors) {
    if (!f->is<Sum>()) continue;
    auto folded = fold_expr(f, opts, report);
    changed = changed || folded != f;
    f = std::move(folded);
  }
  container::svector<bool> consumed(factors.size(), false);
  for (std::size_t gi = 0; gi < factors.size(); ++gi) {
    Tensor const* g = tensor_at(factors, gi);
    if (!g || !is_rank11(*g, opts.integral_label, 1) || g->adjointed())
      continue;
    Index const& a = g->bra().at(0);
    Index const& b = g->ket().at(0);
    if (a.has_proto_indices() || b.has_proto_indices()) continue;
    if (is_occupied(a) || is_occupied(b)) continue;
    // the conjugate coefficient on the bra leg a: the proto index in the
    // bra, spelled by its slots alone (the CSV transform's `C{x<P>;a}`) or
    // with the adjoint mark (`C⁺{x<P>;a}`), both conj(C{a;x<P>}); the stored
    // one (proto in the ket, unmarked) on the ket leg b
    std::optional<std::size_t> conj_pos, stored_pos, mirror_conj, mirror_stored;
    for (std::size_t ci = 0; ci < factors.size(); ++ci) {
      if (ci == gi || consumed[ci]) continue;
      Tensor const* c = tensor_at(factors, ci);
      if (!c || !is_rank11(*c, opts.coeff_label, 0)) continue;
      Index const& cb = c->bra().at(0);
      Index const& ck = c->ket().at(0);
      if (cb.has_proto_indices() && !ck.has_proto_indices()) {
        if (ck == a) conj_pos = ci;
        if (ck == b) mirror_conj = ci;
      } else if (ck.has_proto_indices() && !cb.has_proto_indices() &&
                 !c->adjointed()) {
        if (cb == b) stored_pos = ci;
        if (cb == a) mirror_stored = ci;
      }
    }
    if (!conj_pos || !stored_pos) {
      if (mirror_conj && mirror_stored) ++report.mirrored;
      continue;
    }
    Index const& xidx = tensor_at(factors, *conj_pos)->bra().at(0);
    Index const& yidx = tensor_at(factors, *stored_pos)->ket().at(0);
    if (!ranges::equal(xidx.proto_indices(), yidx.proto_indices())) continue;
    Tensor leaf(opts.integral_label, bra({xidx}), ket({yidx}),
                aux({g->aux().at(0)}),
                TensorSymmetries{.perm = Symmetry::Nonsymm,
                                 .hermiticity = Hermiticity::NonHermitian,
                                 .conjugation_parity = g->conjugation_parity(),
                                 .column = ColumnSymmetry::Nonsymm,
                                 .kramers = g->kramers_symmetry()});
    factors[gi] = ex<Tensor>(std::move(leaf));
    consumed[*conj_pos] = consumed[*stored_pos] = true;
    ++report.folded;
    changed = true;
  }
  if (!changed) return x;
  auto result = std::make_shared<Product>();
  result->scale(flat.scalar());
  for (std::size_t i = 0; i < factors.size(); ++i)
    if (!consumed[i]) result->append(1, factors[i], Product::Flatten::No);
  return result;
}

ExprPtr fold_expr(ExprPtr const& x, PairBasisFoldOptions const& opts,
                  PairBasisFoldReport& report) {
  if (x->is<Product>()) return fold_product(x, opts, report);
  if (!x->is<Sum>()) return x;
  container::svector<ExprPtr> summands;
  summands.reserve(x->size());
  bool changed = false;
  for (auto const& s : *x) {
    summands.push_back(fold_expr(s, opts, report));
    changed = changed || summands.back() != s;
  }
  return changed ? ex<Sum>(std::move(summands)) : x;
}

}  // namespace

ExprPtr fold_pair_basis_integrals(ExprPtr const& expr,
                                  PairBasisFoldOptions const& opts,
                                  PairBasisFoldReport* report) {
  if (!expr) return expr;
  PairBasisFoldReport local;
  auto result = fold_expr(expr, opts, report ? *report : local);
  return result;
}

}  // namespace sequant::mbpt
