//
// Created by Kshitij Surjuse on 2026-08-31.
//

// The conjugation case catalogue: every symbolic identity involving complex
// conjugation that SeQuant is expected to honor, one TEST_CASE per identity
// family. Canonicalization- and network-level cases live at the bottom and
// grow with the conjugation-symbolic work.

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

using namespace sequant;

namespace {
using C = Complex<rational>;
const auto i_unit = C{0, 1};  // the imaginary unit as a Constant value
}  // namespace

TEST_CASE("conj_constant_involution", "[conjugation]") {
  // conj on a complex Constant conjugates the value; twice restores it
  auto c = ex<Constant>(C{1, 2});
  auto cc = conjugate(c);
  REQUIRE(cc->as<Constant>().value() == (C{1, -2}));
  REQUIRE(*conjugate(cc) == *c);
}

TEST_CASE("conj_variable_marker_hash_reset", "[conjugation]") {
  // Variable::conjugate() must reset the memoized hash: compute the hash
  // FIRST, then conjugate, then verify the hash actually changed and that
  // toggling back restores it
  auto v = Variable(L"x");
  const auto h0 = v.hash_value();  // memoize
  v.conjugate();
  REQUIRE(v.conjugated());
  REQUIRE(v.hash_value() != h0);
  v.conjugate();
  REQUIRE(v.hash_value() == h0);
}

TEST_CASE("conj_tensor_marker_roundtrip", "[conjugation]") {
  // (A*)* = A bit-for-bit, hash included; the marker is pure elementwise
  // conjugation: slots are untouched
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  auto t = ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm);
  const auto h0 = t->hash_value();
  auto tc = conjugate(t);
  REQUIRE(tc->as<Tensor>().conjugated());
  REQUIRE(tc->hash_value() != h0);
  // slots untouched
  REQUIRE(tc->as<Tensor>().bra()[0] == t->as<Tensor>().bra()[0]);
  auto tcc = conjugate(tc);
  REQUIRE(*tcc == *t);
  REQUIRE(tcc->hash_value() == h0);
}

TEST_CASE("conjugate_free_function_total", "[conjugation]") {
  // sequant::conjugate dispatches over every scalar node kind and is an
  // involution on each; operator-valued content is rejected loudly
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  // Sum + Product distribution: (c A B)* = conj(c) A* B*, NO reversal
  auto A = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm);
  auto B = ex<Tensor>(L"B", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm);
  auto prod = ex<Constant>(i_unit) * A->clone() * B->clone();
  auto pc = conjugate(prod);
  const auto& p = pc->as<Product>();
  REQUIRE(p.scalar() == (C{0, -1}));
  REQUIRE(p.factors().size() == 2);
  // factor ORDER preserved (contrast adjoint, which reverses)
  REQUIRE(p.factors()[0]->as<Tensor>().label() == L"A");
  REQUIRE(p.factors()[0]->as<Tensor>().conjugated());
  REQUIRE(p.factors()[1]->as<Tensor>().label() == L"B");
  REQUIRE(p.factors()[1]->as<Tensor>().conjugated());
  REQUIRE(*conjugate(pc) == *prod);

  auto sum = A->clone() + B->clone();
  auto sc = conjugate(sum);
  for (const auto& s : *sc) REQUIRE(s->as<Tensor>().conjugated());
  REQUIRE(*conjugate(sc) == *sum);

  // Re/Im are real-valued: conj is the identity on them
  auto re = real_part(A->clone() * B->clone());
  REQUIRE(*conjugate(re) == *re);
}

TEST_CASE("re_im_composition_table", "[conjugation]") {
  // Re/Im are real-valued, so the four compositions collapse:
  //   Re(Re x) = Re x,  Re(Im x) = Im x,  Im(Re x) = 0,  Im(Im x) = 0
  auto x = ex<Variable>(L"x");
  auto re = real_part(x->clone());
  auto im = imaginary_part(x->clone());
  REQUIRE(*real_part(re->clone()) == *re);
  REQUIRE(*real_part(im->clone()) == *im);
  REQUIRE(imaginary_part(re->clone())->as<Constant>().is_zero());
  REQUIRE(imaginary_part(im->clone())->as<Constant>().is_zero());
}

TEST_CASE("re_im_constant_evaluation", "[conjugation]") {
  // Re/Im of a complex Constant evaluate exactly through the Complex ring
  auto c = ex<Constant>(C{3, -4});
  REQUIRE(real_part(c->clone())->as<Constant>().value() == (C{3, 0}));
  REQUIRE(imaginary_part(c->clone())->as<Constant>().value() == (C{-4, 0}));
}

TEST_CASE("braket_foldable_predicates", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // Conjugate c-number: value fold applies
  Tensor h(L"h", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  REQUIRE(braket_conjugate_foldable(h));
  REQUIRE(braket_foldable(h));

  // Nonsymm: no fold of any kind
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  REQUIRE_FALSE(braket_conjugate_foldable(t));
  REQUIRE_FALSE(braket_foldable(t));

  // Symm: free swap, not the Conjugate value fold
  Tensor s(L"s", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  REQUIRE_FALSE(braket_conjugate_foldable(s));
  REQUIRE(braket_foldable(s));
}

TEST_CASE("conjugate_braket_fold_per_tensor", "[conjugation]") {
  // The Conjugate bra<->ket fold engages in per-tensor canonicalization:
  // both orientations of a c-number Conjugate tensor land on ONE canonical
  // spelling, the originally-swapped input acquiring the
  // elementwise-conjugation marker so the represented value is unchanged.
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // occ (i) vs virt (a) bundles: space-decidable orientation
  Tensor A(L"h", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor B(L"h", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  DefaultTensorCanonicalizer::canonicalize_braket(A);
  DefaultTensorCanonicalizer::canonicalize_braket(B);
  // both spell the same slot order; exactly one carries the marker
  REQUIRE(A.bra()[0].label() == B.bra()[0].label());
  REQUIRE(A.conjugated() != B.conjugated());

  // full space tie (same-space named indices): label tie-break folds too
  Tensor C1(L"T", bra{L"p_1"}, ket{L"p_2"}, Symmetry::Nonsymm,
            BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor C2(L"T", bra{L"p_2"}, ket{L"p_1"}, Symmetry::Nonsymm,
            BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  DefaultTensorCanonicalizer::canonicalize_braket(C1);
  DefaultTensorCanonicalizer::canonicalize_braket(C2);
  REQUIRE(C1.bra()[0].label() == C2.bra()[0].label());
  REQUIRE(C1.conjugated() != C2.conjugated());

  // identical bundles (diagonal): never swapped, never marked
  Tensor D(L"T", bra{L"p_1"}, ket{L"p_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  DefaultTensorCanonicalizer::canonicalize_braket(D);
  REQUIRE_FALSE(D.conjugated());
}

TEST_CASE("with_slots_carries_attributes", "[conjugation]") {
  // Tensor::with_slots rebuilds the slots and carries label, symmetries,
  // hermiticity, and the conjugation marker -- the sanctioned rebuild API
  // for transforms (rebuilding through a plain ctor drops the marker)
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  Tensor t(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, Symmetry::Antisymm,
           Hermiticity::Hermitian);
  t.conjugate();
  using ixvec = container::svector<Index>;
  auto r = t.with_slots(bra<ixvec>{ixvec{Index{L"i_3"}, Index{L"i_4"}}},
                        ket<ixvec>{ixvec{Index{L"a_3"}, Index{L"a_4"}}},
                        aux<ixvec>{});
  REQUIRE(r.label() == t.label());
  REQUIRE(r.symmetry() == t.symmetry());
  REQUIRE(r.braket_symmetry() == t.braket_symmetry());
  REQUIRE(r.hermiticity() == t.hermiticity());
  REQUIRE(r.column_symmetry() == t.column_symmetry());
  REQUIRE(r.conjugated());
  REQUIRE(r.bra()[0].label() == L"i_3");
}

TEST_CASE("fold_conjugate_pairs", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // a fully-contracted (energy-like) summand of BraKetSymmetry::Conjugate
  // tensors, and its conjugate written independently: real scalar kept,
  // factor order reversed, bra<->ket swapped, dummies renamed
  auto term = deserialize(L"1/2 h{i_1;a_1}:N-C-S t{a_1;i_1}:N-C-S");
  auto term_adj = deserialize(L"1/2 t{i_2;a_2}:N-C-S h{a_2;i_2}:N-C-S");
  // a manifestly real (self-conjugate) summand: BraKetSymmetry::Symm tensors
  auto self_adj = deserialize(L"1/4 f{i_1;a_1}:N-S-S u{a_1;i_1}:N-S-S");

  SECTION("sum pair emits 2 Re(A)") {
    auto folded = fold_conjugate_pairs(term->clone() + term_adj->clone());
    auto expected = ex<Constant>(2) * real_part(term->clone());
    REQUIRE(*folded == *expected);
  }

  SECTION("difference pair emits 2i Im(A)") {
    auto folded = fold_conjugate_pairs(term->clone() +
                                       ex<Constant>(-1) * term_adj->clone());
    auto expected =
        ex<Constant>(Complex<rational>(0, 2)) * imaginary_part(term->clone());
    REQUIRE(*folded == *expected);
  }

  SECTION("self-conjugate and unpaired summands stay untouched") {
    auto folded = fold_conjugate_pairs(self_adj->clone() + term->clone());
    auto expected = self_adj->clone() + term->clone();
    simplify(folded, SimplifyOptions::default_options().copy_and_set(
                         SimplifyOptions::FoldConjugatePairs::No));
    simplify(expected, SimplifyOptions::default_options().copy_and_set(
                           SimplifyOptions::FoldConjugatePairs::No));
    REQUIRE(*folded == *expected);
  }

  SECTION("mixed sum: pair folds, bystander survives") {
    auto folded = fold_conjugate_pairs(term->clone() + self_adj->clone() +
                                       term_adj->clone());
    REQUIRE(folded->is<Sum>());
    REQUIRE(folded->as<Sum>().summands().size() == 2);
    bool have_re = false;
    for (auto&& sm : *folded)
      if (sm->is<Product>())
        for (auto&& f : sm->as<Product>().factors())
          if (f->is<RealPart>()) have_re = true;
    REQUIRE(have_re);
  }

  SECTION("back-compat real-sum fold emits 2 A") {
    auto folded =
        fold_conjugate_pairs_of_real_sum(term->clone() + term_adj->clone());
    auto expected = ex<Constant>(2) * term->clone();
    simplify(folded);
    simplify(expected);
    REQUIRE(*folded == *expected);
  }

  SECTION("custom conjugate_op recognizes relabeling-based pairs") {
    // spin-annotated spaces for the ↑/↓ labels
    auto spin_ctx = get_default_context();
    spin_ctx.set(mbpt::make_min_sr_spaces());
    auto spin_resetter = set_scoped_default_context(spin_ctx);

    auto term_up = deserialize(L"1/2 h{i↑_1;a↑_1}:N-C-S t{a↑_1;i↑_1}:N-C-S");
    auto term_dn = deserialize(L"1/2 h{i↓_1;a↓_1}:N-C-S t{a↓_1;i↓_1}:N-C-S");
    {  // default (adjoint) pairing finds nothing: label flip, not swap
      auto folded = fold_conjugate_pairs(term_up->clone() + term_dn->clone());
      REQUIRE(folded->is<Sum>());
      REQUIRE(folded->as<Sum>().summands().size() == 2);
    }
    {  // with the label-flip map the pair folds
      auto folded = fold_conjugate_pairs(
          term_up->clone() + term_dn->clone(),
          CanonicalizeOptions::default_options(),
          [](ExprPtr const& sm) { return mbpt::swap_spin(sm); });
      auto expected = ex<Constant>(2) * real_part(term_up->clone());
      REQUIRE(*folded == *expected);
    }
  }

  SECTION("opt-in fold in simplify, complex field") {
    auto sum = term->clone() + term_adj->clone();
    auto folded = sum->clone();
    simplify(folded, SimplifyOptions::default_options().copy_and_set(
                         SimplifyOptions::FoldConjugatePairs::Yes));
    // min_sr spaces are complex-field: the pair folds to 2 Re(...)
    bool have_re = false;
    folded->visit(
        [&](ExprPtr const& node) {
          if (node->is<RealPart>()) have_re = true;
        },
        /*atoms_only=*/false);
    if (folded->is<RealPart>()) have_re = true;
    REQUIRE(have_re);

    // default is No until evaluation understands RealPart/ImagPart nodes
    auto unfolded = sum->clone();
    simplify(unfolded);
    REQUIRE(unfolded->is<Sum>());
    REQUIRE(unfolded->as<Sum>().summands().size() == 2);
  }
}

TEST_CASE("hermitian_network_recognition", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // closed |C|^2 network: C conj(C), fully contracted -> real
  REQUIRE(is_hermitian_network(
      deserialize(L"C{a_1;i_1}:N-C-S C^*{a_1;i_1}:N-C-S")));
  // closed C*C without the conjugation: a generically complex scalar
  REQUIRE_FALSE(
      is_hermitian_network(deserialize(L"C{a_1;i_1}:N-N-S C{a_1;i_1}:N-N-S")));
  // a declared-Hermitian energy-like scalar: h t + t* h* is self-adjoint
  auto term = deserialize(L"h{i_1;a_1}:N-C-S t{a_1;i_1}:N-C-S");
  auto sum = term->clone() + conjugate(term);
  REQUIRE(is_hermitian_network(sum));
}

TEST_CASE("swap_bra_ket_carries_marker_and_aux", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto t = deserialize(L"C{a_1;i_1;p_5}:N-C-S");
  t->as<Tensor>().conjugate();
  auto sw = mbpt::swap_bra_ket(t);
  auto const& st = sw->as<Tensor>();
  REQUIRE(st.conjugated());
  REQUIRE(st.aux().size() == 1);
  REQUIRE(st.bra()[0].label() == L"i_1");
  REQUIRE(st.ket()[0].label() == L"a_1");
}
