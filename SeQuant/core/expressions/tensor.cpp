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

#include <range/v3/algorithm/contains.hpp>

namespace sequant {

Tensor::~Tensor() = default;

void Tensor::assert_nonreserved_label(
    [[maybe_unused]] std::wstring_view label) const {
  SEQUANT_ASSERT(!ranges::contains(FNOperator::labels(), label) &&
                 !ranges::contains(BNOperator::labels(), label));
}

void Tensor::adjoint() {
  std::swap(bra_.value(), ket_.value());

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
  auto canonicalizer_ptr = TensorCanonicalizer::instance_ptr(label_);
  return canonicalizer_ptr ? canonicalizer_ptr->apply(*this) : ExprPtr{};
}

}  // namespace sequant
