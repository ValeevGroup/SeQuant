//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/csv.hpp>

#include <SeQuant/domain/mbpt/space_qns.hpp>

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

/// expands CSVs in a tensor in terms of a basis (standard unoccupieds, PAOs,
/// AOs, etc.)
/// @param tnsr a Tensor object
/// @param csv_basis the basis in terms of which the CSVs are expanded
/// @param orthonormal whether @p csv_basis is orthonormal
ExprPtr csv_transform_impl(Tensor const& tnsr, const IndexBasis& csv_basis,
                           bool orthonormal,
                           std::wstring_view coeff_tensor_label) {
  using ranges::views::transform;
  using sequant::reserved::overlap_label;

  if (ranges::none_of(tnsr.const_braket_indices(), &Index::has_proto_indices))
    return nullptr;

  SEQUANT_ASSERT(ranges::none_of(tnsr.aux(), &Index::has_proto_indices));
  SEQUANT_ASSERT(get_default_context().index_space_registry());
  SEQUANT_ASSERT(
      get_default_context().index_space_registry()->contains(csv_basis));

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
      const Index& end = ordinal_compare(bra_idx, ket_idx) ? bra_idx : ket_idx;
      auto dummy_idx = end.drop_proto_indices().replace_basis_instance({});

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

  rbra.reserve(tnsr.bra_rank());
  for (auto&& idx : tnsr.bra()) {
    if (idx.has_proto_indices()) {
      Index xidx = Index::make_tmp_index(csv_basis);
      result.append(
          1, ex<Tensor>(coeff_tensor_label, bra({idx}), ket({xidx}), aux({})));
      rbra.emplace_back(std::move(xidx));
    } else
      rbra.emplace_back(idx);
  }

  rket.reserve(tnsr.ket_rank());
  for (auto&& idx : tnsr.ket()) {
    if (idx.has_proto_indices()) {
      Index xidx = Index::make_tmp_index(csv_basis);
      result.append(
          1, ex<Tensor>(coeff_tensor_label, bra({xidx}), ket({idx}), aux({})));
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

namespace {

ExprPtr csv_transform_rec(
    ExprPtr const& expr, const IndexBasis& csv_basis, bool orthonormal,
    std::wstring const& coeff_tensor_label,
    container::svector<std::wstring> const& tensor_labels) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(
        *expr                                                       //
        | transform([&csv_basis, orthonormal, &coeff_tensor_label,  //
                     &tensor_labels](auto&& x) {
            return csv_transform_rec(x, csv_basis, orthonormal,
                                     coeff_tensor_label, tensor_labels);
          }));
  else if (expr->is<Tensor>()) {
    auto const& tnsr = expr->as<Tensor>();
    if (!ranges::contains(tensor_labels, tnsr.label())) return expr;
    if (ranges::none_of(tnsr.indices(), &Index::has_proto_indices)) return expr;
    return csv_transform_impl(tnsr, csv_basis, orthonormal, coeff_tensor_label);
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());

    for (auto&& f : prod.factors()) {
      auto trans = csv_transform_rec(f, csv_basis, orthonormal,
                                     coeff_tensor_label, tensor_labels);
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
                      bool orthonormal, std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& tensor_labels) {
  if (csv_basis.has_basis_instance()) {
    const auto& registry = get_default_context().index_space_registry();
    if (!registry || !registry->basis_label(csv_basis))
      throw Exception(
          "csv_transform: the target basis instance " +
          toUtf8(csv_basis.space().base_key()) + ";" +
          std::to_string(*csv_basis.basis_instance()) +
          " is not registered under a name; indices minted in it would print "
          "and be keyed as the space's own basis");
    return csv_transform_rec(expr, registry->resolve(csv_basis), orthonormal,
                             coeff_tensor_label, tensor_labels);
  }
  return csv_transform_rec(expr, csv_basis, orthonormal, coeff_tensor_label,
                           tensor_labels);
}

ExprPtr csv_transform(ExprPtr const& expr, const IndexSpace& csv_basis,
                      std::wstring const& coeff_tensor_label,
                      container::svector<std::wstring> const& tensor_labels) {
  const bool is_ao = bitset_t(csv_basis.qns()) & bitset_t(LCAOQNS::ao);
  const bool is_pao = bitset_t(csv_basis.qns()) & bitset_t(LCAOQNS::pao);
  return csv_transform(expr, IndexBasis{csv_basis},
                       /*orthonormal=*/!(is_ao || is_pao), coeff_tensor_label,
                       tensor_labels);
}

}  // namespace sequant::mbpt
