//
// Created by Eduard Valeyev on 2019-01-30.
//

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expressions/abstract_tensor.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/utility/exception.hpp>
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
  // _swap_bra_ket() swaps bra<->ket *and* the derived net ranks, then
  // re-canonicalizes slots (needed when empty slots are present) and resets the
  // hash; a bare std::swap of the index containers would leave the net ranks
  // and slot order inconsistent
  _swap_bra_ket();

  // adjointness is tracked solely by the label marker, for Nonsymm braket
  if (braket_symmetry() == BraKetSymmetry::Nonsymm) {
    toggle_adjoint_label(label_);
  }

  reset_hash_value();
}

namespace {

/// @return whether an index occurs more than once in the slots of @p t,
/// protoindices included
bool has_repeated_index(const Tensor &t) {
  container::set<Index> seen;
  for (const auto &idx : t.const_slots()) {
    if (!seen.insert(idx).second) return true;
    for (const auto &proto : idx.proto_indices())
      if (!seen.insert(proto).second) return true;
  }
  return false;
}

}  // namespace

ExprPtr Tensor::canonicalize(CanonicalizeOptions opts) {
  if (is_canonical(opts)) return {};
  const auto contexts_version = current_contexts_version();
  const auto ctx = get_default_context_snapshot();
  ExprPtr byproduct;
  if (has_repeated_index(*this)) {
    // a repeated index is contracted, so this is a tensor network and is
    // canonicalized as one, like a Product of tensors
    TensorNetwork tn(static_cast<const Expr &>(*this));
    byproduct = tn.canonicalize(ctx.cardinal_tensor_labels(), opts);
    SEQUANT_ASSERT(tn.tensors().size() == 1);
    const auto canonical = std::dynamic_pointer_cast<Tensor>(tn.tensors()[0]);
    SEQUANT_ASSERT(canonical);
    *this = *canonical;
  } else {
    const auto canonicalizer = ctx.tensor_canonicalizer_ptr(L"");
    if (!canonicalizer)
      throw Exception(
          "Tensor::canonicalize: the current context has no default tensor "
          "canonicalizer");
    byproduct = canonicalizer->apply(*this);
  }
  mark_canonical(opts, contexts_version);
  return byproduct;
}

}  // namespace sequant
