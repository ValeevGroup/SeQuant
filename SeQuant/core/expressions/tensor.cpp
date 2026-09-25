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

ValueOriented value_oriented(Tensor const &t) {
  switch (t.value_modifier()) {
    case ValueModifier::None:
    case ValueModifier::Adjoint:
      // t⁺ names a distinct array whose slots are as written; every consumer
      // treats the '⁺' spelling that way
      return {t, 1};
    case ValueModifier::Transpose: {
      // T^T{q;p} = T{p;q}: a pure respelling (only a tensor without an
      // exchange relation carries this state; the others normalize the
      // transposition away)
      Tensor r{t};
      const auto sign = r.transpose();
      return {std::move(r), sign};
    }
    case ValueModifier::Conjugate:
      if (braket_conjugate_swap_sign(t.braket_symmetry())) {
        // T^*{q;p} = s T{p;q}: transpose() folds into, and so clears, the
        // conjugation bit, contributing the relation's sign
        Tensor r{t};
        const auto sign = r.transpose();
        return {std::move(r), sign};
      }
      throw Exception(
          "sequant::value_oriented: an elementwise-conjugated tensor without "
          "an adjoint relation has no slot spelling of its value "
          "(the conjugation cannot be consumed into slots)");
  }
  SEQUANT_UNREACHABLE;
}

}  // namespace sequant
