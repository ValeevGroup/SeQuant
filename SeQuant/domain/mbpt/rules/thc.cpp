//
// Created by Oliver Backhouse on 22/11/2025.
//

#include <SeQuant/domain/mbpt/rules/thc.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/view.hpp>

#include <string_view>

namespace sequant::mbpt {

ExprPtr tensor_hypercontract_impl(Tensor const& tnsr, Index const& aux_idx_1,
                                  Index const& aux_idx_2,
                                  std::wstring_view factor_label,
                                  std::wstring_view aux_label) {
  SEQUANT_ASSERT(tnsr.bra_rank() == 2     //
                 && tnsr.ket_rank() == 2  //
                 && tnsr.aux_rank() == 0);

  auto t1 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}), ket(),
                       aux({aux_idx_1}));
  auto t2 = ex<Tensor>(factor_label, bra(), ket({ranges::front(tnsr.ket())}),
                       aux({aux_idx_1}));
  auto t3 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}), ket(),
                       aux({aux_idx_2}));
  auto t4 = ex<Tensor>(factor_label, bra(), ket({ranges::back(tnsr.ket())}),
                       aux({aux_idx_2}));
  auto z = ex<Tensor>(aux_label, bra(), ket(), aux({aux_idx_1, aux_idx_2}));

  if (tnsr.symmetry() == Symmetry::Antisymm) {
    auto t1a = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}), ket(),
                          aux({aux_idx_1}));
    auto t3a = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}), ket(),
                          aux({aux_idx_2}));

    return (t1 * t2 * z * t3 * t4) - (t1a * t2 * z * t3a * t4);
  }

  return t1 * t2 * z * t3 * t4;
}

ExprPtr tensor_hypercontract(ExprPtr const& expr, IndexSpace aux_space,
                             std::wstring_view tensor_label,
                             std::wstring_view factor_label,
                             std::wstring_view aux_label) {
  using ranges::views::transform;

  if (expr->is<Sum>())
    return ex<Sum>(*expr | transform([&](auto&& x) {
      return tensor_hypercontract(x, aux_space, tensor_label, factor_label,
                                  aux_label);
    }));

  else if (expr->is<Tensor>()) {
    auto const& tensor = expr->as<Tensor>();
    if (tensor.label() == tensor_label  //
        && tensor.bra_rank() == 2       //
        && tensor.ket_rank() == 2)
      return tensor_hypercontract_impl(tensor, Index(aux_space, 1),
                                       Index(aux_space, 2), factor_label,
                                       aux_label);
    else
      return expr;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());
    size_t aux_ix = 0;
    for (auto&& f : prod.factors()) {
      if (f->is<Tensor>() && f->as<Tensor>().label() == tensor_label) {
        auto const& g = f->as<Tensor>();
        auto index1 = Index(aux_space, ++aux_ix);
        auto index2 = Index(aux_space, ++aux_ix);
        auto g_thc = tensor_hypercontract_impl(g, index1, index2, factor_label,
                                               aux_label);
        result.append(1, std::move(g_thc), Product::Flatten::Yes);
      } else {
        result.append(1, f, Product::Flatten::No);
      }
    }
    return ex<Product>(std::move(result));
  } else {
    return expr;
  }
}

}  // namespace sequant::mbpt
