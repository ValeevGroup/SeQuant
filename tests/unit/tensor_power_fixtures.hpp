//
// A Hadamard product with an element-wise power below a batched contraction,
// E{i_2;i_3} = L · (M ⊙ pow(D,-1)), shared by the dry-run, ordered-executor and
// legality unit tests of batched evaluation through a Power node.
//
// Lives in a named namespace (not an anonymous one) so it can be included by
// several test files under a unity build without ODR collisions.
//

#ifndef SEQUANT_TESTS_UNIT_TENSOR_POWER_FIXTURES_HPP
#define SEQUANT_TESTS_UNIT_TENSOR_POWER_FIXTURES_HPP

#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/result_expr.hpp>
#include <SeQuant/core/index.hpp>

#include <cstddef>
#include <string_view>
#include <utility>

namespace sequant::tests::tensor_power {

/// the default context with AssertStrictBraKetSymmetry::No: i_1 sits in a bra
/// slot of both M and D, which binarize rejects under the strict knob
inline Context lax_braket_context() {
  auto ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  return ctx;
}

/// the symmetries of the power's base: totally symmetric
inline TensorSymmetries symmetric() {
  return {.perm = Symmetry::Symm,
          .braket = BraKetSymmetry::Symm,
          .column = ColumnSymmetry::Symm};
}

/// the symmetries of the other factors: antisymmetric
inline TensorSymmetries antisymmetric() {
  return {.perm = Symmetry::Antisymm,
          .braket = BraKetSymmetry::Symm,
          .column = ColumnSymmetry::Symm};
}

/// D{i_1,i_2;a_1,a_2}, totally symmetric
inline ExprPtr power_base() {
  return ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                    symmetric());
}

/// an antisymmetric tensor \p label{bra;ket}
template <typename Bra, typename Ket>
ExprPtr antisymm(std::wstring_view label, Bra b, Ket k) {
  return ex<Tensor>(label, std::move(b), std::move(k), antisymmetric());
}

/// L{a_1,a_2;i_1,i_3}: contracts i_1, a_1 and a_2 with D's slots and leaves
/// i_3 open
inline ExprPtr contractor() {
  return antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_3"});
}

/// \p label{i_2;i_3}, the head of either term
inline Tensor head(std::wstring_view label) {
  return Tensor(label, bra{L"i_2"}, ket{L"i_3"},
                TensorSymmetries{.perm = Symmetry::Nonsymm});
}

/// E{i_2;i_3} = L{a_1,a_2;i_1,i_3} · (M{i_1,i_2;a_1,a_2} ⊙ pow(D,-1)), the
/// Hadamard product kept as the right operand of the contraction; the head
/// keeps i_2, which only M and D carry, open through the Hadamard product
inline ResultExpr hadamard_power_term() {
  auto e = ex<Product>(ExprPtrList{});
  e->as<Product>().append(1, contractor(), Product::Flatten::No);
  e->as<Product>().append(
      1,
      ex<Product>(
          ExprPtrList{antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
                      pow(power_base(), rational{-1})}),
      Product::Flatten::No);
  return ResultExpr(head(L"E"), e);
}

/// F{i_2;i_3} = L{a_1,a_2;i_1,i_3} · D{i_1,i_2;a_1,a_2}: reads the base of
/// the power directly
inline ResultExpr base_term() {
  auto e = ex<Product>(ExprPtrList{});
  e->as<Product>().append(1, contractor(), Product::Flatten::No);
  e->as<Product>().append(1, power_base(), Product::Flatten::No);
  return ResultExpr(head(L"F"), e);
}

/// the occupied index i_1 the root of either term contracts
inline Index const& batched_occ() {
  static Index const ix{L"i_1"};
  return ix;
}

/// stamps \p root as the node that contracts batched_occ() in batches, the
/// annotation the optimizer gives the contraction that opens the loop
template <typename Node>
void stamp_occ_loop(Node& root) {
  root->set_node_slice_mask({{batched_occ(), BatchModeType::Contracted}});
  root->set_batch_loops_opened_here(
      {{batched_occ(), BatchModeType::Contracted}});
}

/// occupied indices batchable as contracted modes, one element per batch
inline BatchPolicy occ_contracted_policy() {
  BatchPolicy policy;
  policy.is_batchable_contracted_index = [](Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_target_size = [](Index const&) -> std::size_t { return 1; };
  policy.peak_threshold = 1;
  return policy;
}

}  // namespace sequant::tests::tensor_power

#endif  // SEQUANT_TESTS_UNIT_TENSOR_POWER_FIXTURES_HPP
