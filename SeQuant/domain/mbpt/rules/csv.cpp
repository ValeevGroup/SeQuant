//
// Created by Eduard Valeyev on 3/8/25.
//

#include "SeQuant/domain/mbpt/rules/csv.hpp"

#include "SeQuant/domain/mbpt/lcao.hpp"

#include "SeQuant/core/tensor.hpp"
#include "SeQuant/core/utility/indices.hpp"

namespace sequant::mbpt {

/// expands CSVs in a tensor in terms of a basis (standard unoccupieds, PAOs,
/// AOs, etc.)
/// @param tnsr a Tensor object
/// @param csv_basis the basis in terms of which the CSVs are expanded
ExprPtr csv_transform_impl(Tensor const& tnsr, const IndexSpace& csv_basis,
                           std::wstring_view coeff_tensor_label) {
  using ranges::views::transform;

  if (ranges::none_of(tnsr.const_braket(), &Index::has_proto_indices))
    return nullptr;

  assert(ranges::none_of(tnsr.aux(), &Index::has_proto_indices));
  assert(get_default_context().index_space_registry());
  assert(get_default_context().index_space_registry()->contains(csv_basis));

  // shortcut for the CSV overlap if csv_basis is orthonormal (assume LCAO bases
  // are orthonormal!)
  const bool csv_basis_is_orthonormal =
      bitset_t(csv_basis.qns()) & bitset_t(LCAOQNS::lcao);
  if (csv_basis_is_orthonormal && tnsr.label() == overlap_label()) {
    assert(tnsr.bra_rank() == 1     //
           && tnsr.ket_rank() == 1  //
           && tnsr.aux_rank() == 0);

    auto&& bra_idx = tnsr.bra().at(0);
    auto&& ket_idx = tnsr.ket().at(0);

    auto dummy_idx = suffix_compare(bra_idx, ket_idx)    //
                         ? bra_idx.drop_proto_indices()  //
                         : ket_idx.drop_proto_indices();

    return ex<Product>(
        1,
        ExprPtrList{ex<Tensor>(coeff_tensor_label,                 //
                               bra({bra_idx}), ket({dummy_idx})),  //
                    ex<Tensor>(coeff_tensor_label,                 //
                               bra({dummy_idx}), ket({ket_idx}))});
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

  auto xtnsr =
      ex<Tensor>(tnsr.label(), bra(rbra), ket(rket), aux(), tnsr.symmetry(),
                 tnsr.braket_symmetry(), tnsr.particle_symmetry());
  result.append(1, std::move(xtnsr));

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
      result.append(1, trans ? trans : f,
                    (f->is<Product>() || f->is<Sum>()) ? Product::Flatten::No
                                                       : Product::Flatten::Yes);
    }

    return ex<Product>(std::move(result));

  } else
    return expr;
}

}  // namespace sequant::mbpt
