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

ExprPtr Tensor::canonicalize(CanonicalizeOptions) {
  return TensorCanonicalizer::instance()->apply(*this);
}

Tensor value_oriented(Tensor const &t) {
  switch (t.value_modifier()) {
    case ValueModifier::None:
    case ValueModifier::Adjoint:
      // t⁺ names a distinct array whose slots are as written; every consumer
      // has always treated the '⁺' spelling that way
      return t;
    case ValueModifier::Transpose: {
      // T^T{q;p} = T{p;q}: a pure respelling (Nonsymm only; the other
      // symmetries normalize the transposition away)
      Tensor r{t};
      r.transpose();
      return r;
    }
    case ValueModifier::Conjugate:
      if (t.braket_symmetry() == BraKetSymmetry::Conjugate) {
        // T^*{q;p} = T{p;q}: transpose() folds into, and so clears, the
        // conjugation bit
        Tensor r{t};
        r.transpose();
        return r;
      }
      throw Exception(
          "sequant::value_oriented: an elementwise-conjugated "
          "BraKetSymmetry::Nonsymm tensor has no slot spelling of its value "
          "(the conjugation cannot be consumed into slots)");
  }
  SEQUANT_UNREACHABLE;
}

}  // namespace sequant
