//
// Created by Eduard Valeyev on 2019-01-30.
//

#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

namespace sequant {

Tensor::~Tensor() = default;

void Tensor::assert_nonreserved_label(std::wstring_view label) const {
  assert(!ranges::contains(FNOperator::labels(), label) &&
         !ranges::contains(BNOperator::labels(), label));
}

void Tensor::adjoint() {
  std::swap(bra_.value(), ket_.value());

  // only track adjointness for BraKetSymmetry::nonsymm cases
  if (braket_symmetry() == BraKetSymmetry::nonsymm) {
    if (label_.back() == sequant::adjoint_label) {
      assert(is_adjoint_);
      label_.pop_back();
    } else {
      assert(!is_adjoint_);
      label_.push_back(sequant::adjoint_label);
    }
    is_adjoint_ = !is_adjoint_;
  }

  reset_hash_value();
}

ExprPtr Tensor::canonicalize(CanonicalizeOptions) {
  auto canonicalizer_ptr = TensorCanonicalizer::instance_ptr(label_);
  return canonicalizer_ptr ? canonicalizer_ptr->apply(*this) : ExprPtr{};
}

}  // namespace sequant
