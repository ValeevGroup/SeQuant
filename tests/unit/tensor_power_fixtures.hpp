//
// A Hadamard product with an element-wise power below a batched contraction,
// E = L · (M ⊙ pow(D,-1)), shared by the dry-run, ordered-executor and
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

/// D{i_1,i_2;a_1,a_2}, totally symmetric
inline ExprPtr power_base() {
  return ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                    TensorSymmetries{.perm = Symmetry::Symm,
                                     .braket = BraKetSymmetry::Symm,
                                     .column = ColumnSymmetry::Symm});
}

/// an antisymmetric tensor \p label{bra;ket}
template <typename Bra, typename Ket>
ExprPtr antisymm(std::wstring_view label, Bra b, Ket k) {
  return ex<Tensor>(label, std::move(b), std::move(k),
                    TensorSymmetries{.perm = Symmetry::Antisymm,
                                     .braket = BraKetSymmetry::Symm,
                                     .column = ColumnSymmetry::Symm});
}

/// L{a_1,a_2;i_1,i_2} · (M{i_1,i_2;a_1,a_2} ⊙ pow(D,-1)), the Hadamard
/// product kept as the right operand of the full contraction
inline ExprPtr hadamard_power_term() {
  auto e = ex<Product>(ExprPtrList{});
  e->as<Product>().append(
      1, antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}),
      Product::Flatten::No);
  e->as<Product>().append(
      1,
      ex<Product>(
          ExprPtrList{antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
                      pow(power_base(), rational{-1})}),
      Product::Flatten::No);
  return e;
}

/// L{a_1,a_2;i_1,i_2} · D{i_1,i_2;a_1,a_2}: reads the base of the power
/// directly
inline ExprPtr base_term() {
  auto e = ex<Product>(ExprPtrList{});
  e->as<Product>().append(
      1, antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}),
      Product::Flatten::No);
  e->as<Product>().append(1, power_base(), Product::Flatten::No);
  return e;
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
