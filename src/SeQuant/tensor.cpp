//
// Created by Eduard Valeyev on 2019-01-30.
//

#include "tensor.hpp"

namespace sequant {

Tensor::~Tensor() = default;

std::shared_ptr<Expr> Tensor::canonicalize() {
  const auto &canonicalizer = TensorCanonicalizer::instance(label_);
  return canonicalizer->apply(*this);
}

}  // namespace sequant