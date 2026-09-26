//
// Created by Eduard Valeyev on 3/8/25.
//

#include <SeQuant/domain/mbpt/rules/df.hpp>

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/operations.hpp>
#include <range/v3/view/transform.hpp>

#include <algorithm>
#include <string_view>

namespace sequant::mbpt {

/// \param aux_idx the auxiliary index of the first factor (and of the second
/// one in the folded form)
/// \param metric_label if non-empty, the factors are raw and the inverse
/// metric `metric_label{aux_idx;aux_idx2}` sits between them
/// \param aux_idx2 the auxiliary index of the second factor in the metric
/// form (ignored otherwise)
ExprPtr density_fit_impl(Tensor const& tnsr_in, Index const& aux_idx,
                         std::wstring_view factor_label,
                         std::wstring_view metric_label = {},
                         Index const& aux_idx2 = {}) {
  // Normalize to the VALUE orientation first: a marker-conjugated (folded)
  // tensor spells conj(bra<->ket-swapped); rebuilding from its raw slot
  // layout would silently drop the conjugation (see sequant::value_oriented).
  const Tensor tnsr = value_oriented(tnsr_in);
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
  // metric form: the second factor carries its own auxiliary index and the
  // (real symmetric, hence Hermitian) inverse metric sits between the two
  const bool with_metric = !metric_label.empty();
  const Index& aux_idx_2nd = with_metric ? aux_idx2 : aux_idx;
  auto metric = [&]() -> ExprPtr {
    return ex<Tensor>(metric_label, bra({aux_idx}), ket({aux_idx2}),
                      Symmetry::Nonsymm, Hermiticity::Hermitian,
                      ColumnSymmetry::Symm);
  };
  auto join = [&](ExprPtr const& a, ExprPtr const& b) -> ExprPtr {
    return with_metric ? a * metric() * b : a * b;
  };

  auto t1 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}),
                       ket({ranges::front(tnsr.ket())}), aux({aux_idx}),
                       Symmetry::Nonsymm, Hermiticity::Hermitian,
                       ColumnSymmetry::Symm, tnsr.kramers_symmetry());

  auto t2 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}),
                       ket({ranges::back(tnsr.ket())}), aux({aux_idx_2nd}),
                       Symmetry::Nonsymm, Hermiticity::Hermitian,
                       ColumnSymmetry::Symm, tnsr.kramers_symmetry());

  if (tnsr.symmetry() == Symmetry::Antisymm) {
    auto t3 = ex<Tensor>(factor_label, bra({ranges::back(tnsr.bra())}),
                         ket({ranges::front(tnsr.ket())}), aux({aux_idx}),
                         Symmetry::Nonsymm, Hermiticity::Hermitian,
                         ColumnSymmetry::Symm, tnsr.kramers_symmetry());

    auto t4 = ex<Tensor>(factor_label, bra({ranges::front(tnsr.bra())}),
                         ket({ranges::back(tnsr.ket())}), aux({aux_idx_2nd}),
                         Symmetry::Nonsymm, Hermiticity::Hermitian,
                         ColumnSymmetry::Symm, tnsr.kramers_symmetry());
    return join(t1, t2) - join(t3, t4);
  }

  return join(t1, t2);
}

namespace {

/// recursive worker of density_fit. \p aux_ix is a running counter of the aux
/// indices of ONE term: a decomposed tensor takes the next value, and a nested
/// factor (a Sum left behind by a flavor or CSV expansion, or a sub-product)
/// continues from where its siblings stopped, so no two decomposed tensors of
/// the same term share an aux label. Sibling SUMMANDS restart from the value
/// their parent passed in -- they are alternatives, not co-existing factors,
/// and equal spellings across summands are what downstream common-subexpression
/// elimination keys on.
ExprPtr density_fit_rec(ExprPtr const& expr, IndexSpace const& aux_space,
                        std::wstring_view tensor_label,
                        std::wstring_view factor_label,
                        std::function<bool(Tensor const&)> const& should_split,
                        std::wstring_view metric_label, std::size_t& aux_ix) {
  // auxiliary indices drawn per decomposed tensor: one (folded metric) or two
  // (explicit metric leaf between raw factors)
  const std::size_t n_aux = metric_label.empty() ? 1 : 2;
  auto process_tensor = [&](const Tensor& tensor,
                            std::size_t idx_ordinal) -> ExprPtr {
    if (tensor.label() == tensor_label && tensor.bra_net_rank() == 2 &&
        tensor.ket_net_rank() == 2 && tensor.aux_rank() == 0 &&
        (!should_split || should_split(tensor))) {
      return density_fit_impl(tensor, Index(aux_space, idx_ordinal),
                              factor_label, metric_label,
                              Index(aux_space, idx_ordinal + 1));
    }

    return nullptr;
  };

  if (expr->is<Sum>()) {
    const std::size_t aux_ix_in = aux_ix;
    std::size_t aux_ix_max = aux_ix;
    auto out = ex<Sum>(*expr | ranges::views::transform([&](auto&& x) {
      std::size_t aux_ix_summand = aux_ix_in;
      auto res = density_fit_rec(x, aux_space, tensor_label, factor_label,
                                 should_split, metric_label, aux_ix_summand);
      aux_ix_max = std::max(aux_ix_max, aux_ix_summand);
      return res;
    }));
    aux_ix = aux_ix_max;
    return out;
  } else if (expr->is<Tensor>()) {
    if (auto factorized = process_tensor(expr->as<Tensor>(), aux_ix + 1);
        factorized) {
      aux_ix += n_aux;
      return factorized;
    }
    return expr;
  } else if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();

    Product result;
    result.scale(prod.scalar());
    for (auto&& f : prod.factors())
      if (f.is<Tensor>()) {
        if (auto factorized = process_tensor(f->as<Tensor>(), aux_ix + 1);
            factorized) {
          aux_ix += n_aux;
          result.append(1, std::move(factorized), Product::Flatten::Yes);
        } else {
          result.append(1, f, Product::Flatten::No);
        }
      } else {
        // a nested factor -- a Sum left behind by a flavor or CSV expansion, or
        // a sub-product: decompose inside it, continuing this term's aux
        // numbering
        result.append(1,
                      density_fit_rec(f, aux_space, tensor_label, factor_label,
                                      should_split, metric_label, aux_ix),
                      Product::Flatten::No);
      }
    return ex<Product>(std::move(result));
  } else
    return expr;
}

}  // namespace

ExprPtr density_fit(ExprPtr const& expr, IndexSpace aux_space,
                    std::wstring_view tensor_label,
                    std::wstring_view factor_label,
                    std::function<bool(Tensor const&)> const& should_split,
                    std::wstring_view metric_label) {
  std::size_t aux_ix = 0;
  return density_fit_rec(expr, aux_space, tensor_label, factor_label,
                         should_split, metric_label, aux_ix);
}

}  // namespace sequant::mbpt
