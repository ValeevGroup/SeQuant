#ifndef SEQUANT_CORE_DENSITY_HPP
#define SEQUANT_CORE_DENSITY_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include <string>
#include <string_view>

namespace sequant::density {

/// labels of the reference-state density tensors produced by the extended
/// Wick theorem and consumed by mbpt::decompositions
inline std::wstring rdm_label() { return L"γ"; }
inline std::wstring hole_rdm_label() { return L"η"; }
inline std::wstring cumulant_label() { return L"κ"; }

/// densities and cumulants describe indistinguishable particles, hence are
/// column symmetric, and are Hermitian by definition; both take part in the
/// tensor hash, so every producer must spell them identically
inline constexpr TensorSymmetries rdm_symmetries{
    .hermiticity = Hermiticity::Hermitian, .column = ColumnSymmetry::Symm};
inline constexpr TensorSymmetries cumulant_symmetries{
    .perm = Symmetry::Antisymm,
    .hermiticity = Hermiticity::Hermitian,
    .column = ColumnSymmetry::Symm};

/// @return γ with bra = @p ann and ket = @p cre, i.e. ⟨a†_cre a_ann⟩
inline ExprPtr make_rdm(const Index &ann, const Index &cre) {
  return ex<Tensor>(rdm_label(), bra{ann}, ket{cre}, rdm_symmetries);
}

/// @return η with bra = @p ann and ket = @p cre, i.e. ⟨a_ann a†_cre⟩ = δ - γ
inline ExprPtr make_hole_rdm(const Index &ann, const Index &cre) {
  return ex<Tensor>(hole_rdm_label(), bra{ann}, ket{cre}, rdm_symmetries);
}

/// @return the reference expectation value of @p nop as a tensor with the
/// given label: bra = annihilator indices, ket = creator indices, both in
/// particle order
template <Statistics S>
ExprPtr rdm_from_nop(const NormalOperator<S> &nop, std::wstring_view label,
                     TensorSymmetries syms) {
  using index_container = container::svector<Index>;
  auto braidxs =
      nop.annihilators() |
      ranges::views::transform([](const auto &op) { return op.index(); }) |
      ranges::to<index_container>();
  auto ketidxs =
      nop.creators() |
      ranges::views::transform([](const auto &op) { return op.index(); }) |
      ranges::to<index_container>();
  SEQUANT_ASSERT(braidxs.size() == ketidxs.size());
  return ex<Tensor>(std::wstring(label), bra(std::move(braidxs)),
                    ket(std::move(ketidxs)), syms);
}

/// @return κ_k for the k-body @p nop (k ≥ 2)
template <Statistics S>
ExprPtr make_cumulant(const NormalOperator<S> &nop) {
  SEQUANT_ASSERT(nop.rank() >= 2);
  return rdm_from_nop(nop, cumulant_label(), cumulant_symmetries);
}

}  // namespace sequant::density

#endif  // SEQUANT_CORE_DENSITY_HPP
