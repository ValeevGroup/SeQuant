//
// Created by Eduard Valeyev on 2019-01-30.
//

#include "tensor.hpp"

namespace sequant {

Tensor::~Tensor() = default;

void Tensor::assert_nonreserved_label(std::wstring_view label) const {
  // assert(label != overlap_label());
}

void Tensor::adjoint() {
  std::swap(bra_, ket_);

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

ExprPtr Tensor::canonicalize() {
  auto canonicalizer_ptr = TensorCanonicalizer::instance_ptr(label_);
  return canonicalizer_ptr ? canonicalizer_ptr->apply(*this) : ExprPtr{};
}

}  // namespace sequant
