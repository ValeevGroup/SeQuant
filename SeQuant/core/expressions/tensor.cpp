//
// Created by Eduard Valeyev on 2019-01-30.
//

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <stdexcept>

#include <range/v3/algorithm/contains.hpp>

namespace sequant {

Tensor::~Tensor() = default;

void Tensor::assert_nonreserved_label(
    [[maybe_unused]] std::wstring_view label) const {
  SEQUANT_ASSERT(!ranges::contains(FNOperator::labels(), label) &&
                 !ranges::contains(BNOperator::labels(), label));
}

void Tensor::adjoint() {
  // _swap_bra_ket() swaps bra<->ket *and* the derived net ranks, then
  // re-canonicalizes slots (needed when empty slots are present) and resets the
  // hash; a bare std::swap of the index containers would leave the net ranks
  // and slot order inconsistent
  _swap_bra_ket();

  // adjointness is tracked solely by the label marker, for Nonsymm braket
  if (braket_symmetry() == BraKetSymmetry::Nonsymm) {
    if (!label_.empty() && label_.back() == sequant::adjoint_label)
      label_.pop_back();
    else
      label_.push_back(sequant::adjoint_label);
  }

  reset_hash_value();
}

ExprPtr Tensor::canonicalize(CanonicalizeOptions) {
  return TensorCanonicalizer::instance()->apply(*this);
}

Tensor value_oriented(Tensor const &t) {
  if (!t.conjugated()) return t;
  if (t.braket_symmetry() == BraKetSymmetry::Nonsymm)
    throw std::logic_error(
        "sequant::value_oriented: an elementwise-conjugated "
        "BraKetSymmetry::Nonsymm tensor has no value-oriented slot spelling "
        "(the conjugation cannot be consumed into slots)");
  Tensor bare{t};
  bare.conjugate();
  if (t.braket_symmetry() == BraKetSymmetry::Conjugate)
    bare.adjoint();  // pure bra<->ket swap: undoes the fold
  // Symm: conj is the identity in value -- clearing the marker suffices
  return bare;
}

}  // namespace sequant
