//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/df.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/operations.hpp>
#include <range/v3/view/transform.hpp>

#include <string_view>

namespace sequant::mbpt {

ExprPtr density_fit_impl(Tensor const& tnsr, Index const& aux_idx,
                         std::wstring_view factor_label) {
  SEQUANT_ASSERT(tnsr.bra_rank() == 2     //
                 && tnsr.ket_rank() == 2  //
                 && tnsr.aux_rank() == 0);

  // The 3-center DF factor (pq|X) is a matrix element of a (real, symmetric)
  // Coulomb metric and is therefore Hermitian in its p<->q (bra<->ket) pair --
  // (pq|X) = (qp|X) -- regardless of the spaces of p, q. Declaring it Hermitian
  // (rather than the default non-Hermitian) lets a real computation treat it as
  // bra<->ket symmetric, so e.g. (pq|X) C^p and (pq|X) C^q collapse to one
  // intermediate. The concrete BraKetSymmetry (Symm vs Conjugate) is derived
  // from the bra/ket indices' IndexSpace::field() (see sequant::base_field)
  // when the Tensor is built.
  auto t1 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}),
                       ket({ranges::front(tnsr.ket())}), aux({aux_idx}),
                       Symmetry::Nonsymm, Hermiticity::Hermitian);

  auto t2 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}),
                       ket({ranges::back(tnsr.ket())}), aux({aux_idx}),
                       Symmetry::Nonsymm, Hermiticity::Hermitian);

  if (tnsr.symmetry() == Symmetry::Antisymm) {
    auto t3 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}),
                         ket({ranges::front(tnsr.ket())}), aux({aux_idx}),
                         Symmetry::Nonsymm, Hermiticity::Hermitian);

    auto t4 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}),
                         ket({ranges::back(tnsr.ket())}), aux({aux_idx}),
                         Symmetry::Nonsymm, Hermiticity::Hermitian);
    return t1 * t2 - t3 * t4;
  }

  return t1 * t2;
}

ExprPtr density_fit(ExprPtr const& expr, IndexSpace aux_space,
                    std::wstring_view tensor_label,
                    std::wstring_view factor_label) {
  using ranges::views::transform;

  auto process_tensor = [&](const Tensor& tensor,
                            std::size_t idx_ordinal) -> ExprPtr {
    if (tensor.label() == tensor_label && tensor.bra_net_rank() == 2 &&
        tensor.ket_net_rank() == 2 && tensor.aux_rank() == 0) {
      return density_fit_impl(tensor, Index(aux_space, idx_ordinal),
                              factor_label);
    }

    return nullptr;
  };

  if (expr->is<Sum>()) {
    return ex<Sum>(*expr | transform([&](auto&& x) {
      return density_fit(x, aux_space, tensor_label, factor_label);
    }));
  } else if (expr->is<Tensor>()) {
    if (auto factorized = process_tensor(expr->as<Tensor>(), 1); factorized) {
      return factorized;
    }
    return expr;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());
    size_t aux_ix = 0;
    for (auto&& f : prod.factors())
      if (f.is<Tensor>()) {
        if (auto factorized = process_tensor(f->as<Tensor>(), ++aux_ix);
            factorized) {
          result.append(1, std::move(factorized), Product::Flatten::Yes);
        } else {
          result.append(1, f, Product::Flatten::No);
        }
      } else {
        result.append(1, f, Product::Flatten::No);
      }
    return ex<Product>(std::move(result));
  } else
    return expr;
}

}  // namespace sequant::mbpt
