//
// Created by Eduard Valeyev on 2019-01-30.
//

#include "tensor.hpp"

namespace sequant {

Tensor::~Tensor() = default;

void Tensor::assert_nonreserved_label(std::wstring_view label) const {
  // assert(label != overlap_label());
}

void Tensor::toggle_adjoint_label() {
  if (is_adjoint_) {
    if (label_.back() != sequant::adjoint_label)
      label_.push_back(sequant::adjoint_label);
  } else {
    if (label_.back() == sequant::adjoint_label) label_.pop_back();
  }
}

void Tensor::adjoint() {
  // toggle adjoint flag and update label
  is_adjoint_ = !is_adjoint_;
  toggle_adjoint_label();
  std::swap(bra_, ket_);
  reset_hash_value();
}

ExprPtr Tensor::canonicalize() {
  auto canonicalizer_ptr = TensorCanonicalizer::instance_ptr(label_);
  return canonicalizer_ptr ? canonicalizer_ptr->apply(*this) : ExprPtr{};
}

}  // namespace sequant
