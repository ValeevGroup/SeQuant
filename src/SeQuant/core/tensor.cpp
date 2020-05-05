//
// Created by Eduard Valeyev on 2019-01-30.
//

#include "tensor.hpp"

namespace sequant {

Tensor::~Tensor() = default;

void
Tensor::adjoint() {
  std::swap(bra_, ket_);
  reset_hash_value();
}

ExprPtr Tensor::canonicalize() {
  const auto &canonicalizer = TensorCanonicalizer::instance(label_);
  return canonicalizer->apply(*this);
}

}  // namespace sequant