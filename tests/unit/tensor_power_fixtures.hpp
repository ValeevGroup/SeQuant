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

#include <cmath>
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

/// \p t with its bra and ket replaced, slot by slot, by \p b and \p k; label,
/// aux and symmetries kept
template <typename Bra, typename Ket>
Tensor relabelled(Tensor const& t, Bra b, Ket k) {
  return Tensor(t.label(), std::move(b), std::move(k), aux(t.aux()),
                TensorSymmetries{.perm = t.symmetry(),
                                 .braket = t.braket_symmetry(),
                                 .column = t.column_symmetry()});
}

/// ℓ{i_1,i_2;a_1,a_2} = M{i_1,i_2;a_1,a_2} ⊙ pow(D,-1), the Hadamard factor
/// of hadamard_power_term() as its own equation; antisymmetric like M
inline ResultExpr ell_definition() {
  return ResultExpr(
      Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, antisymmetric()),
      ex<Product>(
          ExprPtrList{antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
                      pow(power_base(), rational{-1})}));
}

/// E{i_2;i_3} = L{a_1,a_2;i_1,i_3} · ℓ{i_1,i_2;a_1,a_2}: hadamard_power_term()
/// with the Hadamard factor read through the head of ell_definition()
inline ResultExpr consumer_E() {
  auto e = ex<Product>(ExprPtrList{});
  e->as<Product>().append(1, contractor(), Product::Flatten::No);
  e->as<Product>().append(
      1,
      ex<Tensor>(relabelled(ell_definition().result_as_tensor(),
                            bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"})),
      Product::Flatten::No);
  return ResultExpr(head(L"E"), e);
}

/// D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}, the base of the
/// power as a broadcast Sum of orbital energies
inline ResultExpr D_broadcast() {
  auto eps = [](auto b, auto k) {
    return ex<Tensor>(L"ε", std::move(b), std::move(k));
  };
  return ResultExpr(power_base()->as<Tensor>(),
                    eps(bra{}, ket{L"a_1"}) + eps(bra{}, ket{L"a_2"}) -
                        eps(bra{L"i_1"}, ket{}) - eps(bra{L"i_2"}, ket{}));
}

/// leaf data for numerical tests of the inlined system, chosen so that the
/// leaves honor their declared symmetries (random data would not, and the
/// evaluator relies on them when it reorders slots) and every denominator is
/// bounded away from zero. Arguments are slot positions (0-based ordinals).
namespace data {

/// ε{i;}: occupied orbital energies, all negative
inline double eps_occ(std::size_t i) { return -0.5 - 0.37 * i; }

/// ε{;a}: virtual orbital energies, all positive
inline double eps_virt(std::size_t a) { return 0.8 + 0.53 * a; }

/// D{i_1,i_2;a_1,a_2} as D_broadcast() defines it; totally symmetric
inline double D(std::size_t i1, std::size_t i2, std::size_t a1,
                std::size_t a2) {
  return eps_virt(a1) + eps_virt(a2) - eps_occ(i1) - eps_occ(i2);
}

/// an element of an antisymmetric 4-index tensor X{p,q;r,s}, \p seed telling
/// tensors apart
inline double antisymm(double seed, std::size_t p, std::size_t q, std::size_t r,
                       std::size_t s) {
  auto const g = [seed](std::size_t p, std::size_t q, std::size_t r,
                        std::size_t s) {
    return std::sin(seed + 1.3 * p + 2.9 * q * q + 0.7 * r + 1.1 * s * s +
                    0.37 * p * s);
  };
  return g(p, q, r, s) - g(q, p, r, s) - g(p, q, s, r) + g(q, p, s, r);
}

inline double M(std::size_t i1, std::size_t i2, std::size_t a1,
                std::size_t a2) {
  return antisymm(0.3, i1, i2, a1, a2);
}

inline double L(std::size_t a1, std::size_t a2, std::size_t i1,
                std::size_t i3) {
  return antisymm(1.7, a1, a2, i1, i3);
}

inline double K(std::size_t a1, std::size_t a2, std::size_t i1,
                std::size_t i2) {
  return antisymm(2.9, a1, a2, i1, i2);
}

/// t{i;a}, nonsymmetric
inline double t(std::size_t i, std::size_t a) {
  return std::cos(0.4 + 1.3 * i + 0.8 * a * a);
}

/// ℓ = M ⊙ pow(D,-1), per ell_definition()
inline double ell(std::size_t i1, std::size_t i2, std::size_t a1,
                  std::size_t a2) {
  return M(i1, i2, a1, a2) / D(i1, i2, a1, a2);
}

/// Σ_{i_1,a_1,a_2} X{a_1,a_2;i_1,x} ℓ{i_1,y;a_1,a_2} over \p nocc occupied and
/// \p nvirt virtual orbitals; E{i_2;i_3} of consumer_E() is
/// contract_with_ell(L, i_3, i_2)
template <typename X>
double contract_with_ell(X const& x, std::size_t xi, std::size_t yi,
                         std::size_t nocc, std::size_t nvirt) {
  double result = 0;
  for (std::size_t i1 = 0; i1 < nocc; ++i1)
    for (std::size_t a1 = 0; a1 < nvirt; ++a1)
      for (std::size_t a2 = 0; a2 < nvirt; ++a2)
        result += x(a1, a2, i1, xi) * ell(i1, yi, a1, a2);
  return result;
}

}  // namespace data

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
