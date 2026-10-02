#ifndef SEQUANT_CORE_DENSITY_HPP
#define SEQUANT_CORE_DENSITY_HPP

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/algorithm/contains.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include <cstddef>
#include <string>
#include <string_view>

namespace sequant::density {

/// labels of the reference-state density tensors produced by the extended
/// Wick theorem and consumed by mbpt::decompositions. They are reserved (see
/// reserved::density_labels()): a density's symmetries take part in the tensor
/// hash, so a γ built by hand with other symmetries would not merge with the
/// γ the engine builds. Construct densities with the factories below.
inline const std::wstring &rdm_label() { return reserved::rdm_label(); }
inline const std::wstring &hole_rdm_label() {
  return reserved::hole_rdm_label();
}
inline const std::wstring &cumulant_label() {
  return reserved::cumulant_label();
}
inline const std::wstring &spinfree_rdm_label() {
  return reserved::spinfree_rdm_label();
}

/// densities and cumulants describe indistinguishable particles, hence are
/// column symmetric, and are Hermitian by definition; both take part in the
/// tensor hash, so every producer must spell them identically
inline constexpr TensorSymmetries rdm_symmetries{
    .hermiticity = Hermiticity::Hermitian, .column = ColumnSymmetry::Symm};
inline constexpr TensorSymmetries cumulant_symmetries{
    .perm = Symmetry::Antisymm,
    .hermiticity = Hermiticity::Hermitian,
    .column = ColumnSymmetry::Symm};

/// @return the defining symmetries of the density tensor @p label of rank
/// @p rank: a multi-body spin-orbital density (γ, η, κ) is antisymmetric,
/// a spin-free one (Γ) and a 1-body one are not
/// @note a multi-body spin-orbital density may also be perm-nonsymmetric: that
///       is a spin component of it, as spin tracing produces
inline TensorSymmetries symmetries(std::wstring_view label, std::size_t rank) {
  return rank > 1 && label != spinfree_rdm_label() ? cumulant_symmetries
                                                   : rdm_symmetries;
}

/// @return the density tensor @p label (one of reserved::density_labels())
/// with the given bra (annihilator) and ket (creator) indices and its
/// defining symmetries()
template <range_of_castables_to_index IndexRange1,
          range_of_castables_to_index IndexRange2>
ExprPtr make_density(std::wstring_view label,
                     const bra<IndexRange1> &bra_indices,
                     const ket<IndexRange2> &ket_indices) {
  SEQUANT_ASSERT(ranges::contains(reserved::density_labels(), label));
  SEQUANT_ASSERT(bra_indices.size() == ket_indices.size());
  return ex<Tensor>(std::wstring(label), bra_indices, ket_indices,
                    symmetries(label, bra_indices.size()));
}

/// @return γ with bra = @p ann and ket = @p cre, i.e. ⟨a†_cre a_ann⟩
inline ExprPtr make_rdm(const Index &ann, const Index &cre) {
  return make_density(rdm_label(), bra{ann}, ket{cre});
}

/// @return the multi-body γ with bra = @p ann and ket = @p cre
template <range_of_castables_to_index IndexRange1,
          range_of_castables_to_index IndexRange2>
ExprPtr make_rdm(const bra<IndexRange1> &ann, const ket<IndexRange2> &cre) {
  return make_density(rdm_label(), ann, cre);
}

/// @return η with bra = @p ann and ket = @p cre, i.e. ⟨a_ann a†_cre⟩ = δ - γ
inline ExprPtr make_hole_rdm(const Index &ann, const Index &cre) {
  return make_density(hole_rdm_label(), bra{ann}, ket{cre});
}

/// @return the reference expectation value of @p nop as the density tensor
/// @p label: bra = annihilator indices, ket = creator indices, both in
/// particle order
template <Statistics S>
ExprPtr rdm_from_nop(const NormalOperator<S> &nop, std::wstring_view label) {
  using index_container = container::svector<Index>;
  auto braidxs =
      nop.annihilators() |
      ranges::views::transform([](const auto &op) { return op.index(); }) |
      ranges::to<index_container>();
  auto ketidxs =
      nop.creators() |
      ranges::views::transform([](const auto &op) { return op.index(); }) |
      ranges::to<index_container>();
  return make_density(label, bra(std::move(braidxs)), ket(std::move(ketidxs)));
}

/// @return κ with bra = @p ann and ket = @p cre
template <range_of_castables_to_index IndexRange1,
          range_of_castables_to_index IndexRange2>
ExprPtr make_cumulant(const bra<IndexRange1> &ann,
                      const ket<IndexRange2> &cre) {
  return make_density(cumulant_label(), ann, cre);
}

/// @return κ_k for the k-body @p nop (k ≥ 2)
template <Statistics S>
ExprPtr make_cumulant(const NormalOperator<S> &nop) {
  SEQUANT_ASSERT(nop.rank() >= 2);
  return rdm_from_nop(nop, cumulant_label());
}

}  // namespace sequant::density

#endif  // SEQUANT_CORE_DENSITY_HPP
