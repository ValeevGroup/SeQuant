//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/csv.hpp>
#include <SeQuant/domain/mbpt/spinor.hpp>

#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/view/transform.hpp>

#include <string>

namespace sequant::mbpt {

namespace {

/// @return whether @p idx carries a definite spin (Kramers) label, ↑ or ↓;
///         a spin-free index (Spin::any) has none
bool spin_labeled(const Index& idx) {
  const auto qns = idx.space().qns();
  if ((bitset_t(qns) & mask_v<Spin>) == 0) return false;
  const auto s = to_spin(qns);
  return s == Spin::alpha || s == Spin::beta;
}

/// a fresh expansion dummy in the @p s-flavored counterpart of @p csv_basis:
/// the basis of a probe index in @p csv_basis respelled with the spin case
/// (make_index_with_spincase keeps the basis instance and resolves the
/// registered flavored basis). N.B. make_spinalpha(make_tmp_index(csv_basis))
/// would feed the tmp-range ordinal back through the non-tmp Index ctor and
/// throw, hence the probe.
Index make_flavored_csv_index(const IndexBasis& csv_basis, Spin s) {
  const Index probe(csv_basis, 1);
  return Index::make_tmp_index(
      (s == Spin::alpha ? make_spinalpha(probe) : make_spinbeta(probe))
          .basis());
}

}  // namespace

/// expands CSVs in a tensor in terms of a basis (standard unoccupieds, PAOs,
/// AOs, etc.)
/// @param tnsr a Tensor object
/// @param csv_basis the basis in terms of which the CSVs are expanded
/// @param orthonormal whether @p csv_basis is orthonormal
/// @param kramers if true, spin(Kramers)-labeled CSV indices expand over BOTH
///        spin flavors of the expansion dummy (see csv_transform)
/// @param kramers_union (Kramers mode only) if true, a spin-labeled CSV index
///        expands over ONE spin-free dummy spanning both flavors
ExprPtr csv_transform_impl(Tensor const& tnsr, const IndexBasis& csv_basis,
                           bool orthonormal,
                           std::wstring_view coeff_tensor_label, bool kramers,
                           bool kramers_union) {
  using ranges::views::transform;
  using sequant::reserved::overlap_label;

  if (ranges::none_of(tnsr.const_braket_indices(), &Index::has_proto_indices))
    return nullptr;

  SEQUANT_ASSERT(ranges::none_of(tnsr.aux(), &Index::has_proto_indices));
  SEQUANT_ASSERT(get_default_context().index_basis_registry());
  SEQUANT_ASSERT(
      get_default_context().index_basis_registry()->contains(csv_basis));

  // a basis instance is registered in one space, so an index minted in it
  // cannot keep a leg's space (e.g. its spin)
  auto require_in_basis_space = [&csv_basis](const Index& leg) {
    if (csv_basis.has_basis_instance() && leg.space() != csv_basis.space())
      throw Exception("csv_transform: a CSV leg is in space " +
                      toUtf8(leg.space().base_key()) +
                      ", not in the target basis instance's space " +
                      toUtf8(csv_basis.space().base_key()));
  };

  // the CSV coefficient enters the product unconjugated with the default
  // (non-Hermitian) traits; a conjugated C is spelled through the core
  // states (see set_states below) rather than a pinned braket symmetry,
  // which the core model refuses over a real basis
  auto csv_coeff = [&](Index a, Index b) {
    return ex<Tensor>(coeff_tensor_label, bra({std::move(a)}),
                      ket({std::move(b)}), aux({}));
  };

  // The expansion dummy of one CSV-projected leg in the default
  // (ms-conserving) mode. A spin-labeled leg keeps its Kramers flavor on the
  // fresh index: the relativistic (Kramers-traced) residual carries ↑/↓ on
  // its virtuals, and the downstream Kramers reconstruction / leaf yielder
  // must tell which Kramers block the CSV coefficient belongs to. A
  // spin-free leg (the non-relativistic path) is minted in csv_basis itself.
  auto make_csv_index = [&](const Index& leg) -> Index {
    if (spin_labeled(leg))
      return make_flavored_csv_index(csv_basis, to_spin(leg.space().qns()));
    require_in_basis_space(leg);
    return Index::make_tmp_index(csv_basis);
  };
  // The expansion-dummy candidates of one leg: the single dummy above, or in
  // Kramers mode for a spin-labeled leg BOTH flavors of csv_basis (the CSV
  // vectors have support on both Kramers row-blocks, so the ms-conserving
  // single-flavor expansion would silently drop the Kramers-off-diagonal
  // coefficient blocks), or (kramers_union) ONE spin-free dummy spanning the
  // union of both flavors, which the consumer serves as the union leaf.
  auto candidates = [&](const Index& leg) -> container::svector<Index> {
    if (kramers && spin_labeled(leg)) {
      if (kramers_union) return {Index::make_tmp_index(csv_basis)};
      return {make_flavored_csv_index(csv_basis, Spin::alpha),
              make_flavored_csv_index(csv_basis, Spin::beta)};
    }
    return {make_csv_index(leg)};
  };

  // shortcut for the CSV overlap if csv_basis is orthonormal
  if (orthonormal && tnsr.label() == overlap_label()) {
    SEQUANT_ASSERT(tnsr.bra_rank() == 1     //
                   && tnsr.ket_rank() == 1  //
                   && tnsr.aux_rank() == 0);

    auto&& bra_idx = tnsr.bra().at(0);
    auto&& ket_idx = tnsr.ket().at(0);
    [[maybe_unused]] const auto bra_has_proto_indices =
        bra_idx.has_proto_indices();
    [[maybe_unused]] const auto ket_has_proto_indices =
        ket_idx.has_proto_indices();
    SEQUANT_ASSERT(bra_has_proto_indices || ket_has_proto_indices);

    if (bra_has_proto_indices && ket_has_proto_indices) {
      // the contracted index of the C†C expansion is a FRESH dummy (another
      // overlap may share either leg) with the Kramers flavor of the leg that
      // compares first by ordinal: in Kramers mode both flavors (or their
      // union) are summed over, see candidates()
      const Index& end = ordinal_compare(bra_idx, ket_idx) ? bra_idx : ket_idx;
      auto cc = [&](const Index& dummy) {
        return ex<Product>(1, ExprPtrList{csv_coeff(bra_idx, dummy),
                                          csv_coeff(dummy, ket_idx)});
      };
      if (kramers && spin_labeled(end)) {
        Sum out;
        for (const auto& dummy : candidates(end)) out.append(cc(dummy));
        if (out.size() == 1) return out.summand(0);
        return ex<Sum>(std::move(out));
      }
      if (spin_labeled(end)) return cc(make_csv_index(end));
      require_in_basis_space(end);
      return cc(Index::make_tmp_index(csv_basis.has_basis_instance()
                                          ? csv_basis
                                          : IndexBasis(end.space())));
    } else if ((bra_has_proto_indices ? ket_idx : bra_idx).basis() ==
               csv_basis) {
      // the other leg is in the target basis
      return ex<Product>(
          1,
          ExprPtrList{ex<Tensor>(coeff_tensor_label,  //
                                 bra({bra_idx}), ket({ket_idx}))});
    }
  }

  // Gather the CSV-projected (proto-indexed) slots and their expansion-dummy
  // candidates, then emit the sum over all candidate combinations (one term
  // outside Kramers mode).
  const container::svector<Index> bras(tnsr.bra().begin(), tnsr.bra().end());
  const container::svector<Index> kets(tnsr.ket().begin(), tnsr.ket().end());
  struct Slot {
    bool in_bra;
    std::size_t pos;
  };
  container::svector<Slot> slots;
  container::svector<container::svector<Index>> cands;
  auto gather = [&](container::svector<Index> const& idxs, bool in_bra) {
    for (std::size_t p = 0; p < idxs.size(); ++p) {
      const auto& idx = idxs[p];
      if (!idx.has_proto_indices()) continue;
      slots.push_back(Slot{in_bra, p});
      cands.push_back(candidates(idx));
    }
  };
  gather(bras, true);
  gather(kets, false);

  std::size_t ncombo = 1;
  for (auto const& c : cands) ncombo *= c.size();

  Sum out;
  for (std::size_t combo = 0; combo < ncombo; ++combo) {
    Product result;
    container::svector<Index> rbra(bras), rket(kets);
    std::size_t rem = combo;
    for (std::size_t sl = 0; sl < slots.size(); ++sl) {
      const Index& x = cands[sl][rem % cands[sl].size()];
      rem /= cands[sl].size();
      const auto& orig =
          slots[sl].in_bra ? bras[slots[sl].pos] : kets[slots[sl].pos];
      if (slots[sl].in_bra) {
        result.append(1, csv_coeff(orig, x));
        rbra[slots[sl].pos] = x;
      } else {
        result.append(1, csv_coeff(x, orig));
        rket[slots[sl].pos] = x;
      }
    }
    auto xtnsr = ex<Tensor>(tnsr.label(), bra(rbra), ket(rket), tnsr.aux(),
                            tnsr.symmetries());
    // the source's traits and its two core states ride onto the csv slots.
    // Where those slots are over a different field the normalization can
    // trade the adjointed state for the K-conjugated one -- exchanging the
    // bundles back -- and clear that one against an odd conjugation parity,
    // consuming a sign. A Tensor holds no sign, so it goes to this Product's
    // scalar. N.B. carrying a state across a change of field presumes the
    // coefficient tensors are real: each C straddles a complex source index
    // and a real csv index and enters the product unconjugated, so it
    // commutes with the conjugation the state denotes.
    const auto states_sign =
        xtnsr->as<Tensor>().set_states(tnsr.adjointed(), tnsr.kconjugated());
    result.prepend(states_sign, std::move(xtnsr));
    out.append(ex<Product>(std::move(result)));
  }
  if (ncombo == 1) return out.summand(0);
  return ex<Sum>(std::move(out));
}

namespace {

ExprPtr csv_transform_rec(ExprPtr const& expr, const IndexBasis& csv_basis,
                          bool orthonormal,
                          std::wstring const& coeff_tensor_label,
                          container::svector<std::wstring> const& tensor_labels,
                          bool kramers, bool kramers_union) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(
        *expr                                                       //
        | transform([&csv_basis, orthonormal, &coeff_tensor_label,  //
                     &tensor_labels, kramers, kramers_union](auto&& x) {
            return csv_transform_rec(x, csv_basis, orthonormal,
                                     coeff_tensor_label, tensor_labels, kramers,
                                     kramers_union);
          }));
  else if (expr->is<Tensor>()) {
    auto const& tnsr = expr->as<Tensor>();
    if (!ranges::contains(tensor_labels, tnsr.label())) return expr;
    if (ranges::none_of(tnsr.indices(), &Index::has_proto_indices)) return expr;
    auto out = csv_transform_impl(tnsr, csv_basis, orthonormal,
                                  coeff_tensor_label, kramers, kramers_union);
    // the CSV coefficients (and the rebuilt tensor) inherit the source
    // tensor's Kramers symmetry: C_μ^{a↓<..>} and C_μ^{a↑<..>} are
    // time-reversal partners exactly when the tensor they expand is
    if (tnsr.kramers_symmetry() != KramersSymmetry::Nonsymm)
      out = mark_kramers_symmetric(out, tnsr.kramers_symmetry());
    return out;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());

    for (auto&& f : prod.factors()) {
      auto trans =
          csv_transform_rec(f, csv_basis, orthonormal, coeff_tensor_label,
                            tensor_labels, kramers, kramers_union);
      // N.B. do not flatten the product to ensure that CSV transform of
      // each factor is performed before assembling the final product
      // this way for DF-factorized integrals each DF factor is transformed
      // to CSV basis before DF reconstruction
      result.append(1, trans ? trans : f, Product::Flatten::No);
    }

    return ex<Product>(std::move(result));

  } else
    return expr;
}

}  // namespace

ExprPtr csv_transform(ExprPtr const& expr, const IndexBasis& csv_basis,
                      std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& tensor_labels,
                      bool kramers, bool kramers_union) {
  if (!csv_basis.has_basis_instance())
    return csv_transform_rec(expr, csv_basis, /*orthonormal=*/true,
                             coeff_tensor_label, tensor_labels, kramers,
                             kramers_union);
  const auto& registry = get_default_context().index_basis_registry();
  if (!registry || !registry->basis_label(csv_basis))
    throw Exception(
        "csv_transform: the target basis instance " +
        toUtf8(csv_basis.space().base_key()) + ";" +
        std::to_string(*csv_basis.basis_instance()) +
        " is not registered under a name; indices minted in it would print "
        "and be keyed as the space's own basis");
  const auto resolved = registry->resolve(csv_basis);
  return csv_transform_rec(
      expr, resolved, resolved.metric() == IndexBasisMetric::Unit,
      coeff_tensor_label, tensor_labels, kramers, kramers_union);
}

}  // namespace sequant::mbpt
