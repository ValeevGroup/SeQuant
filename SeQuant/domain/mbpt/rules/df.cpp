//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/df.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#include <string_view>

namespace sequant::mbpt {

ExprPtr density_fit_impl(Tensor const& tnsr, Index const& aux_idx,
                         std::wstring_view factor_label) {
  assert(tnsr.bra_rank() == 2     //
         && tnsr.ket_rank() == 2  //
         && tnsr.aux_rank() == 0);

  auto t1 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}),
                       ket({ranges::front(tnsr.ket())}), aux({aux_idx}));

  auto t2 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}),
                       ket({ranges::back(tnsr.ket())}), aux({aux_idx}));

  return ex<Product>(1, ExprPtrList{t1, t2});
}

ExprPtr density_fit(ExprPtr const& expr, IndexSpace aux_space,
                    std::wstring_view tensor_label,
                    std::wstring_view factor_label) {
  using ranges::views::transform;
  if (expr->is<Sum>())
    return ex<Sum>(*expr | transform([&](auto&& x) {
      return density_fit(x, aux_space, tensor_label, factor_label);
    }));

  else if (expr->is<Tensor>()) {
    auto const& tensor = expr->as<Tensor>();
    if (tensor.label() == tensor_label  //
        && tensor.bra_rank() == 2       //
        && tensor.ket_rank() == 2       //
        && ranges::none_of(tensor.indices(), &Index::has_proto_indices))
      return density_fit_impl(expr->as<Tensor>(), Index(L"1", aux_space),
                              factor_label);
    else
      return expr;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());
    size_t aux_ix = 0;
    for (auto&& f : prod.factors())
      if (f.is<Tensor>() && f.as<Tensor>().label() == L"g") {
        auto const& g = f->as<Tensor>();
        auto g_df = density_fit_impl(
            g, Index(std::to_wstring(++aux_ix), aux_space), factor_label);
        result.append(1, std::move(g_df), Product::Flatten::Yes);
      } else {
        result.append(1, f, Product::Flatten::No);
      }
    return ex<Product>(std::move(result));
  } else
    return expr;
}

}  // namespace sequant::mbpt
