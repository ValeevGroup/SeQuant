//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/csv.hpp>

#include <SeQuant/domain/mbpt/space_qns.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/view/transform.hpp>

namespace sequant::mbpt {

/// expands CSVs in a tensor in terms of a basis (standard unoccupieds, PAOs,
/// AOs, etc.)
/// @param tnsr a Tensor object
/// @param csv_basis the basis in terms of which the CSVs are expanded
ExprPtr csv_transform_impl(Tensor const& tnsr, const IndexSpace& csv_basis,
                           std::wstring_view coeff_tensor_label) {
  using ranges::views::transform;
  using sequant::reserved::overlap_label;

  if (ranges::none_of(tnsr.const_braket_indices(), &Index::has_proto_indices))
    return nullptr;

  SEQUANT_ASSERT(ranges::none_of(tnsr.aux(), &Index::has_proto_indices));
  SEQUANT_ASSERT(get_default_context().index_space_registry());
  SEQUANT_ASSERT(
      get_default_context().index_space_registry()->contains(csv_basis));

  // shortcut for the CSV overlap if csv_basis is orthonormal (assume LCAO bases
  // are orthonormal!)
  const bool csv_basis_is_ao =
      bitset_t(csv_basis.qns()) & bitset_t(LCAOQNS::ao);
  const bool csv_basis_is_pao =
      bitset_t(csv_basis.qns()) & bitset_t(LCAOQNS::pao);
  const bool csv_basis_is_orthonormal = !(csv_basis_is_ao || csv_basis_is_pao);
  if (csv_basis_is_orthonormal && tnsr.label() == overlap_label()) {
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
      // The contracted index of the C†C overlap expansion must be a FRESH index
      // in the (Kramers-labeled) CSV basis. drop_proto_indices() keeps the base
      // ordinal (e.g. a↑_1), so two overlaps that share a base virtual collide
      // on the same dummy -- that index then appears >2 times across the term
      // and the eval-time tensor network is rejected as invalid. Take the
      // spin-labeled base space from drop_proto_indices() (which correctly
      // preserves the ↑/↓ Kramers label) but mint a unique tmp index in it.
      const auto dummy_space = (ordinal_compare(bra_idx, ket_idx)  //
                                    ? bra_idx.drop_proto_indices()
                                    : ket_idx.drop_proto_indices())
                                   .space();
      auto dummy_idx = Index::make_tmp_index(dummy_space);

      return ex<Product>(
          1,
          ExprPtrList{ex<Tensor>(coeff_tensor_label,                 //
                                 bra({bra_idx}), ket({dummy_idx})),  //
                      ex<Tensor>(coeff_tensor_label,                 //
                                 bra({dummy_idx}), ket({ket_idx}))});
    } else {
      return ex<Product>(
          1,
          ExprPtrList{ex<Tensor>(coeff_tensor_label,  //
                                 bra({bra_idx}), ket({ket_idx}))});
    }
  }

  Product result;
  container::svector<Index> rbra, rket;

  // Preserve the Kramers spincase (↑/↓) of the CSV-projected index on the fresh
  // CSV-basis index. The relativistic (Kramers-traced) residual carries α/β
  // labels on its virtuals; make_tmp_index(csv_basis) alone would drop them, so
  // the downstream Kramers reconstruction / leaf yielder could not tell which
  // Kramers block (↑ vs ↓) the CSV coefficient belongs to. Spin-free indices
  // are left untouched (non-relativistic path).
  auto make_csv_index = [&csv_basis](const Index& src) -> Index {
    const auto qns = src.space().qns();
    if ((bitset_t(qns) & mask_v<Spin>) != 0) {
      const auto s = to_spin(qns);
      // Make the fresh CSV index directly in the spin-labeled csv_basis space.
      // NB: make_spinalpha(make_tmp_index(csv_basis)) would rebuild the index
      // via make_index_with_spincase, which feeds the tmp-range ordinal back
      // through the non-tmp Index ctor and throws ("ordinal must be less than
      // min_tmp_index"). Derive the spin-labeled space from a throwaway non-tmp
      // probe, then make a fresh tmp index in that space.
      const Index probe(csv_basis, 1);
      if (s == Spin::alpha)
        return Index::make_tmp_index(make_spinalpha(probe).space());
      if (s == Spin::beta)
        return Index::make_tmp_index(make_spinbeta(probe).space());
    }
    return Index::make_tmp_index(csv_basis);
  };

  // The CSV coefficient's bra<->ket exchange is a complex conjugation (C_μ^{a}
  // vs C^μ_{a}), so mark it BraKetSymmetry::Conjugate rather than letting the
  // ctor derive it from the (over a complex/Kramers field, unreliable)
  // CSV-basis index field. This lets eval-time exploit_conjugate fold the two
  // orientations onto one cached value and serve conj(C) via an Adjoint
  // (elementwise conj) on retrieval; without it a conj(C) leaf is silently
  // served as bare C, biasing a complex (Kramers) result with a spurious
  // imaginary part. No-op over a real field (conjugation is the identity and
  // exploit_conjugate is off); Conjugate is a no-swap symmetry to the
  // canonicalizer, like the previous default.
  auto csv_coeff = [&](Index a, Index b) {
    return ex<Tensor>(coeff_tensor_label, bra({std::move(a)}),
                      ket({std::move(b)}), aux({}), Symmetry::Nonsymm,
                      BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm);
  };

  rbra.reserve(tnsr.bra_rank());
  for (auto&& idx : tnsr.bra()) {
    if (idx.has_proto_indices()) {
      Index xidx = make_csv_index(idx);
      result.append(1, csv_coeff(idx, xidx));
      rbra.emplace_back(std::move(xidx));
    } else
      rbra.emplace_back(idx);
  }

  rket.reserve(tnsr.ket_rank());
  for (auto&& idx : tnsr.ket()) {
    if (idx.has_proto_indices()) {
      Index xidx = make_csv_index(idx);
      result.append(1, csv_coeff(xidx, idx));
      rket.emplace_back(std::move(xidx));
    } else
      rket.emplace_back(idx);
  }

  auto xtnsr = ex<Tensor>(tnsr.label(), bra(rbra), ket(rket), tnsr.aux(),
                          tnsr.symmetry(), tnsr.braket_symmetry(),
                          tnsr.column_symmetry());
  result.prepend(1, std::move(xtnsr));

  return ex<Product>(std::move(result));
}

ExprPtr csv_transform(ExprPtr const& expr, const IndexSpace& csv_basis,
                      std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& tensor_labels) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(*expr                                          //
                   | transform([&csv_basis, &coeff_tensor_label,  //
                                &tensor_labels](auto&& x) {
                       return csv_transform(x, csv_basis, coeff_tensor_label,
                                            tensor_labels);
                     }));
  else if (expr->is<Tensor>()) {
    auto const& tnsr = expr->as<Tensor>();
    if (!ranges::contains(tensor_labels, tnsr.label())) return expr;
    if (ranges::none_of(tnsr.indices(), &Index::has_proto_indices)) return expr;
    return csv_transform_impl(tnsr, csv_basis, coeff_tensor_label);
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());

    for (auto&& f : prod.factors()) {
      auto trans =
          csv_transform(f, csv_basis, coeff_tensor_label, tensor_labels);
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

}  // namespace sequant::mbpt
