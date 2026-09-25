//
// Created by Kshitij Surjuse on 2026-08-31.
//

// The conjugation case catalogue: every symbolic identity involving complex
// conjugation that SeQuant is expected to honor, one TEST_CASE per identity
// family. Canonicalization- and network-level cases live at the bottom and
// grow with the conjugation-symbolic work.

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_algorithms.hpp>
#include <SeQuant/core/expressions/power.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/sum.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/expressions/variable.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/v3.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

using namespace sequant;

namespace {
using C = Complex<rational>;
const auto i_unit = C{0, 1};  // the imaginary unit as a Constant value

/// @return an Index named @p label whose space carries the given @p field
///         (IndexSpace's default field is Complex)
Index idx(std::wstring_view label, Field field) {
  Index i(label);
  IndexSpace sp = i.space();
  sp.field(field);
  return Index(label, sp);
}
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

  // Antisymm (anti-Hermitian over a real basis): a swap carrying -1, not the
  // conjugate value fold
  Tensor n(L"n", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
           TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                            .column = ColumnSymmetry::Symm});
  REQUIRE(n.braket_symmetry() == BraKetSymmetry::Antisymm);
  REQUIRE_FALSE(braket_conjugate_foldable(n));
  REQUIRE(braket_foldable(n));

  // AntiConjugate (anti-Hermitian over the complex basis): the conjugate
  // value fold applies, at -1
  Tensor d(L"d", bra{L"i_1"}, ket{L"a_1"},
           TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                            .column = ColumnSymmetry::Symm});
  REQUIRE(d.braket_symmetry() == BraKetSymmetry::AntiConjugate);
  REQUIRE(braket_conjugate_foldable(d));
  REQUIRE(braket_foldable(d));

  // operator-valued: reorienting would exchange creators and annihilators, so
  // no fold applies however symmetric the bra/ket exchange looks (this one's
  // bra and ket agree, hence Hermitian, hence Conjugate over the complex
  // basis)
  FNOperator op(cre({L"i_1"}), ann({L"i_1"}));
  REQUIRE(braket_symmetry(op) == BraKetSymmetry::Conjugate);
  REQUIRE_FALSE(braket_conjugate_foldable(op));
  REQUIRE_FALSE(braket_foldable(op));
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
  REQUIRE(t.conjugate() == 1);
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
  REQUIRE(t->as<Tensor>().conjugate() == 1);
  auto sw = mbpt::swap_bra_ket(t);
  auto const& st = sw->as<Tensor>();
  REQUIRE(st.conjugated());
  REQUIRE(st.aux().size() == 1);
  REQUIRE(st.bra()[0].label() == L"i_1");
  REQUIRE(st.ket()[0].label() == L"a_1");
}

TEST_CASE("re_im_scalar_rules", "[conjugation]") {
  auto x = ex<Variable>(L"x");
  auto E = ex<Variable>(L"E");

  // real scalar hoists: Re(2E) = 2 Re(E), Im(2E) = 2 Im(E)
  {
    auto re = real_part(ex<Constant>(2) * E->clone());
    REQUIRE(re->is<Product>());
    REQUIRE(re->as<Product>().scalar() == (C{2, 0}));
    REQUIRE(re->as<Product>().factors()[0]->is<RealPart>());
  }
  // i-rotation: Re(i A) = -Im(A), Im(i A) = Re(A)
  {
    auto re = real_part(ex<Constant>(i_unit) * x->clone());
    REQUIRE(re->is<Product>());
    REQUIRE(re->as<Product>().scalar() == (C{-1, 0}));
    REQUIRE(re->as<Product>().factors()[0]->is<ImagPart>());
    auto im = imaginary_part(ex<Constant>(i_unit) * x->clone());
    REQUIRE(im->as<Product>().scalar() == (C{1, 0}));
    REQUIRE(im->as<Product>().factors()[0]->is<RealPart>());
  }
  // general complex scalar stays wrapped (recognized, not auto-expanded)
  {
    auto re = real_part(ex<Constant>(C{1, 1}) * x->clone());
    REQUIRE(re->is<RealPart>());
  }
}

TEST_CASE("adjoint_conjugate_transpose_relations", "[conjugation]") {
  // the Klein four-group {id, conj, transpose, adjoint}: each op is an
  // involution, adjoint = transpose o conj = conj o transpose, and the
  // modifier bits record exactly which group element was applied
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto t0 = deserialize(L"t{a_1;i_1}:N-N-S");
  auto& T0 = t0->as<Tensor>();

  SECTION("involutions") {
    for (auto op : {&Tensor::conjugate, &Tensor::transpose, &Tensor::adjoint}) {
      Tensor t = T0;
      REQUIRE((t.*op)() == 1);
      REQUIRE(t != T0);
      REQUIRE((t.*op)() == 1);
      REQUIRE(t == T0);
      REQUIRE(t.hash_value() == T0.hash_value());
    }
  }

  SECTION("transpose swaps slots and sets the bit") {
    Tensor t = T0;
    REQUIRE(t.transpose() == 1);
    REQUIRE(t.value_modifier() == ValueModifier::Transpose);
    REQUIRE(t.bra()[0].label() == L"i_1");
    REQUIRE(t.ket()[0].label() == L"a_1");
    REQUIRE(serialize(ex<Tensor>(t), {.annot_symm = true}) ==
            L"t^T{i_1;a_1}:N-N-S");
  }

  SECTION("adjoint = transpose o conjugate = conjugate o transpose") {
    Tensor a = T0;
    REQUIRE(a.adjoint() == 1);
    Tensor tc = T0;
    REQUIRE(tc.transpose() == 1);
    REQUIRE(tc.conjugate() == 1);
    Tensor ct = T0;
    REQUIRE(ct.conjugate() == 1);
    REQUIRE(ct.transpose() == 1);
    REQUIRE(a == tc);
    REQUIRE(a == ct);
    REQUIRE(a.value_modifier() == ValueModifier::Adjoint);
    Tensor ctt = T0;
    REQUIRE(ctt.conjugate_transpose() == 1);
    REQUIRE(a == ctt);
  }

  SECTION("conj(adjoint(t)) is the transpose") {
    Tensor t = T0;
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.conjugate() == 1);
    REQUIRE(t.value_modifier() == ValueModifier::Transpose);
  }
}

TEST_CASE("value_modifier_normalization", "[conjugation]") {
  // the bits are normalized against the braket symmetry: Symm clears both,
  // Conjugate folds the transposition into the conjugation, Nonsymm keeps both
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  Tensor g(L"g", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor s(L"s", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Symm, ColumnSymmetry::Symm);

  SECTION("Conjugate: transpose() is the starred swapped spelling") {
    Tensor gt = g;
    REQUIRE(gt.transpose() == 1);
    REQUIRE(gt.value_modifier() == ValueModifier::Conjugate);
    REQUIRE(gt.bra()[0].label() == L"a_1");
    REQUIRE(serialize(ex<Tensor>(gt), {.annot_symm = true}) ==
            L"g^*{a_1;i_1}:N-C-S");
    // and adjoint() is a pure swap
    Tensor ga = g;
    REQUIRE(ga.adjoint() == 1);
    REQUIRE(ga.value_modifier() == ValueModifier::None);
    REQUIRE(ga.bra()[0].label() == L"a_1");
    // transpose() again unfolds
    REQUIRE(gt.transpose() == 1);
    REQUIRE(gt == g);
  }

  SECTION("Symm: every modifier is the identity") {
    Tensor sc = s;
    REQUIRE(sc.conjugate() == 1);
    REQUIRE(sc == s);
    Tensor st = s;
    REQUIRE(st.transpose() == 1);
    REQUIRE(st.value_modifier() == ValueModifier::None);
    REQUIRE(st.bra()[0].label() == L"a_1");  // slots did swap
  }

  SECTION("a ⁺ label on a Hermitian tensor normalizes away") {
    Tensor g_adj(L"g⁺", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                 BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    Tensor g_swapped(L"g", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    REQUIRE(g_adj == g_swapped);
    REQUIRE(g_adj.value_modifier() == ValueModifier::None);
  }

  SECTION("set_value_modifier normalizes too") {
    Tensor gt = g;
    REQUIRE(gt.set_value_modifier(ValueModifier::Transpose) == 1);
    REQUIRE(gt.value_modifier() == ValueModifier::Conjugate);
    Tensor sa = s;
    REQUIRE(sa.set_value_modifier(ValueModifier::Adjoint) == 1);
    REQUIRE(sa.value_modifier() == ValueModifier::None);
  }
}

TEST_CASE("conj_serialization_roundtrip", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto t = deserialize(L"g{i_1,i_2;a_1,a_2}:A-C-S");
  REQUIRE(t->as<Tensor>().conjugate() == 1);
  auto rt = deserialize(serialize(t));
  REQUIRE(rt->as<Tensor>().conjugated());
  REQUIRE(*rt == *t);
}

TEST_CASE("conjugation_parity_serialization", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);
  // fourth letter: parity; absent means Even
  // a trait letter back-fills the default parity, Even, so an Even tensor
  // spelled by one needs no fourth letter
  auto d = deserialize(L"d{i_1;i_2}:N-A-N");  // anti-Hermitian trait letter
  REQUIRE(d->as<Tensor>().hermiticity() == Hermiticity::AntiHermitian);
  REQUIRE(d->as<Tensor>().braket_symmetry() == BraKetSymmetry::AntiConjugate);
  REQUIRE(serialize(d, {.annot_symm = true}) == L"d{i_1;i_2}:N-A-N");
  auto p = deserialize(L"p{i_1;i_2}:N-H-N-O");
  REQUIRE(p->as<Tensor>().conjugation_parity() == ConjugationParity::Odd);
  // over the default complex field a Hermitian tensor's observable braket
  // symmetry is Conjugate, spelled 'C'; the parity letter survives
  REQUIRE(serialize(p, {.annot_symm = true}) == L"p{i_1;i_2}:N-C-N-O");
  // an unsigned braket letter (Nonsymm) together with an explicit Even
  // parity back-fills NonHermitian; re-serializing drops the parity letter
  // again since Even is never spelled out
  auto t = deserialize(L"t{i_1;i_2}:N-N-N-E");
  REQUIRE(t->as<Tensor>().conjugation_parity() == ConjugationParity::Even);
  REQUIRE(serialize(t, {.annot_symm = true}) == L"t{i_1;i_2}:N-N-N");
}

TEST_CASE("conjugation_parity_serialization_real_field", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  // a real-field Odd Hermitian tensor: over Field::Real this observable
  // braket symmetry is a plain (anti)symmetry, not a conjugation, so the
  // trait letter 'H' (not 'C') is what survives serialization
  Index i1 = idx(L"i_1", Field::Real);
  Index i2 = idx(L"i_2", Field::Real);
  auto t = ex<Tensor>(
      L"t", bra{i1}, ket{i2},
      TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                       .conjugation_parity = ConjugationParity::Odd});
  REQUIRE(t->as<Tensor>().braket_symmetry() == BraKetSymmetry::Antisymm);
  auto str = serialize(t, {.annot_symm = true});
  REQUIRE(str == L"t{i_1;i_2}:N-H-N-O");

  // point the deserializer at spaces that are Real too, and check the
  // parity and braket symmetry round-trip
  auto sr_reg = std::make_shared<IndexSpaceRegistry>(
      get_default_context().index_space_registry()->clone());
  std::vector<std::wstring> keys;
  for (const auto& s : *sr_reg) keys.push_back(s.base_key());
  for (const auto& k : keys)
    if (auto* sp = sr_reg->retrieve_ptr(k)) sp->field(Field::Real);
  auto real_resetter =
      set_scoped_default_context(Context(get_default_context()).set(sr_reg));

  auto rt = deserialize(str);
  REQUIRE(rt->as<Tensor>().conjugation_parity() == ConjugationParity::Odd);
  REQUIRE(rt->as<Tensor>().braket_symmetry() == BraKetSymmetry::Antisymm);

  // the fourth letter is emitted only when the parser could not back-fill the
  // parity from the bra/ket letter. Over a real field a pinned Conjugate *is*
  // what parity None spells, so `:A-C-S` needs no `-N` and round-trips
  // verbatim; the trait letters back-fill Even, so `:N-A-N` needs nothing
  // either while an Odd parity under `H` still spells itself out
  REQUIRE(serialize(deserialize(L"t{i_1;i_2}:A-C-S"), {.annot_symm = true}) ==
          L"t{i_1;i_2}:A-C-S");
  REQUIRE(serialize(deserialize(L"d{i_1;i_2}:N-A-N"), {.annot_symm = true}) ==
          L"d{i_1;i_2}:N-A-N");
  REQUIRE(serialize(deserialize(L"p{i_1;i_2}:N-H-N-O"), {.annot_symm = true}) ==
          L"p{i_1;i_2}:N-H-N-O");

  // and the expression round-trip the other way round: what the letters leave
  // out is exactly what the parser puts back
  auto pinned = ex<Tensor>(L"c", bra{i1}, ket{i2},
                           TensorSymmetries{.braket = BraKetSymmetry::Conjugate,
                                            .column = ColumnSymmetry::Symm});
  REQUIRE(pinned->as<Tensor>().conjugation_parity() == ConjugationParity::None);
  REQUIRE(serialize(pinned, {.annot_symm = true}) == L"c{i_1;i_2}:N-C-S");
  REQUIRE(*deserialize(serialize(pinned, {.annot_symm = true})) == *pinned);
  REQUIRE(*deserialize(str) == *t);
}

TEST_CASE("conj_power_roundtrip", "[conjugation]") {
  auto p = ex<Power>(ex<Variable>(L"x"), 2);
  auto pc = conjugate(p);
  REQUIRE(pc->hash_value() != p->hash_value());
  REQUIRE(*conjugate(pc) == *p);
}

TEST_CASE("tn_slots_determinism", "[conjugation]") {
  // T17: canonicalize_slots is presentation-independent for a
  // conjugate-marked network -- both factor orders and both conj
  // placements of C(x;m) C*(y;m) land on one graph/hash family with a
  // consistent per-tensor conj report
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto C = [](const wchar_t* ext) {
    return ex<Tensor>(L"C", bra{Index{ext}}, ket{Index{L"p_1"}},
                      Symmetry::Nonsymm, BraKetSymmetry::Conjugate,
                      ColumnSymmetry::Symm);
  };
  auto Cstar = [&](const wchar_t* ext) {
    auto t = C(ext);
    REQUIRE(t->as<Tensor>().conjugate() == 1);
    return t;
  };
  auto md = [](ExprPtr e) {
    TensorNetworkV3 tn(e);
    return tn.canonicalize_slots(TensorNetworkV3::CanonicalizeSlotsOptions{});
  };
  auto m1 = md(C(L"a_1") * Cstar(L"a_2"));
  auto m2 = md(Cstar(L"a_2") * C(L"a_1"));  // factor order flipped
  REQUIRE(m1.hash_value() == m2.hash_value());
  REQUIRE(m1.graph->cmp(*m2.graph) == 0);
  // the conj-swapped spelling is a DIFFERENT value and keeps its own slot
  auto m3 = md(Cstar(L"a_1") * C(L"a_2"));
  REQUIRE(m3.hash_value() == m1.hash_value());  // one shared graph family
}

TEST_CASE("conjugate_fold_skips_reserved", "[conjugation]") {
  // reserved bookkeeping operators ((anti)symmetrizers) never reorient and
  // never acquire the marker
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  auto e = deserialize(L"Â{i_1,i_2;a_1,a_2}:A g{a_1,a_2;i_1,i_2}:A-C-S");
  canonicalize(e);
  bool found_A = false;
  e->visit(
      [&](ExprPtr const& node) {
        if (node->is<Tensor>() &&
            node->as<Tensor>().label() == reserved::antisymm_label()) {
          found_A = true;
          REQUIRE_FALSE(node->as<Tensor>().conjugated());
          REQUIRE(node->as<Tensor>().bra()[0].space() == Index(L"i_1").space());
        }
      },
      /*atoms_only=*/true);
  REQUIRE(found_A);
}

TEST_CASE("sum_merge_conjugate_marked_terms", "[conjugation]") {
  // identically-marked summands merge; a marked and an unmarked spelling of
  // DIFFERENT values do not
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto t = deserialize(L"t{a_1;i_1}:N-N-S");
  auto tstar = conjugate(t);
  auto sum = tstar->clone() + tstar->clone();
  simplify(sum);
  REQUIRE(sum->is<Product>());
  REQUIRE(sum->as<Product>().scalar() == (C{2, 0}));

  auto mixed = t->clone() + tstar->clone();
  simplify(mixed);
  REQUIRE(mixed->is<Sum>());
  REQUIRE(mixed->as<Sum>().summands().size() == 2);
}

TEST_CASE("eval_tot_leaf_named_index_comparator", "[conjugation]") {
  // a proto-indexed (ToT) leaf's canon_indices puts occupieds first: the
  // DECLARED default comparator (default_idxptr_slottype_lesscompare) orders
  // by proto-index count before space -- the layout downstream
  // coefficient-shape detectors rely on
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto C = deserialize(L"C{a_1<i_1>;i_2}:N-C-S")->as<Tensor>();
  EvalExpr leaf{C};
  auto const& ci = leaf.canon_indices();
  // named indices: the proto i_1, the ket i_2, and the ToT virtual a_1<i_1>
  REQUIRE(ci.size() == 3);
  // proto-free indices precede proto-indexed ones (the comparator orders by
  // proto-index count before space)
  REQUIRE_FALSE(ci[0].has_proto_indices());
  REQUIRE_FALSE(ci[1].has_proto_indices());
  REQUIRE(ci[2].has_proto_indices());
}

TEST_CASE("value_modifier_encoding", "[conjugation]") {
  // The adjoint mark is a value modifier, not a label character: a '⁺'
  // arriving in a label is adopted into the bits, label() is bare, and
  // decorated_label() reproduces that spelling for printing.
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  REQUIRE(t.value_modifier() == ValueModifier::None);
  REQUIRE(t.decorated_label() == L"t");

  SECTION("adjoint() sets both bits and keeps the ⁺ spelling") {
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    REQUIRE(ta.label() == L"t");
    REQUIRE(ta.value_modifier() == ValueModifier::Adjoint);
    REQUIRE(ta.conjugated());
    REQUIRE(ta.transposed());
    REQUIRE(ta.decorated_label() == L"t⁺");
    REQUIRE(to_latex(ta) == L"{t⁺^{{a_1}}_{{i_1}}}");
    REQUIRE(serialize(ex<Tensor>(ta), {.annot_symm = true}) ==
            L"t⁺{i_1;a_1}:N-N-N");
  }

  SECTION("a ⁺ in the label is adopted into the bits") {
    Tensor from_label(L"t⁺", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    REQUIRE(from_label.label() == L"t");
    REQUIRE(from_label.value_modifier() == ValueModifier::Adjoint);
    REQUIRE(from_label == ta);
    REQUIRE(from_label.hash_value() == ta.hash_value());
    // same through set_label
    Tensor relabeled = t;
    relabeled.set_label(L"t⁺");
    REQUIRE(relabeled.label() == L"t");
    REQUIRE(relabeled.value_modifier() == ValueModifier::Adjoint);

    // a Hermitian tensor's adjoint is itself: the mark is dropped, no bits
    Tensor g_adj(L"g⁺", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                 BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    REQUIRE(g_adj.label() == L"g");
    REQUIRE(g_adj.value_modifier() == ValueModifier::None);
    REQUIRE(g_adj.decorated_label() == L"g");
  }

  SECTION(
      "every value modifier enters the hash the same way, on the bare "
      "label") {
    // t, t^*, t^T and t⁺ share one label; each modifier must still give the
    // tensor a distinct hash, since a colliding pair would alias distinct
    // values onto one cache slot.
    Tensor t_star = t;
    REQUIRE(t_star.conjugate() == 1);
    Tensor t_transposed = t;
    REQUIRE(t_transposed.transpose() == 1);
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);

    const auto h = t.hash_value();
    const auto h_star = t_star.hash_value();
    const auto h_transposed = t_transposed.hash_value();
    const auto h_adj = t_adj.hash_value();
    REQUIRE(h != h_star);
    REQUIRE(h != h_transposed);
    REQUIRE(h != h_adj);
    REQUIRE(h_star != h_transposed);
    REQUIRE(h_star != h_adj);
    REQUIRE(h_transposed != h_adj);

    // the same holds one level up, for the EvalExpr leaves binarize serves a
    // marked tensor as (hash_terminal_tensor's invariant)
    REQUIRE(EvalExpr{t_adj}.hash_value() != EvalExpr{t_star}.hash_value());
  }

  SECTION("deserializer: ⁺, ^* and ^T round-trip") {
    for (auto spelling :
         {L"t⁺{i_1;a_1}:N-N-N", L"t^*{a_1;i_1}:N-N-N", L"t^T{i_1;a_1}:N-N-N"}) {
      auto e = deserialize(spelling);
      REQUIRE(e->is<Tensor>());
      REQUIRE(e->as<Tensor>().label() == L"t");
      REQUIRE(serialize(e, {.annot_symm = true}) == spelling);
    }
    REQUIRE(deserialize(L"t^T{i_1;a_1}:N-N-N")->as<Tensor>().value_modifier() ==
            ValueModifier::Transpose);
    REQUIRE(to_latex(deserialize(L"t^T{i_1;a_1}:N-N-N")) ==
            L"{{t^T}^{{a_1}}_{{i_1}}}");

    // a hand-written ⁺ followed by ^* composes to the transpose
    REQUIRE(
        deserialize(L"t⁺^*{i_1;a_1}:N-N-N")->as<Tensor>().value_modifier() ==
        ValueModifier::Transpose);

    // a ⁺ on an anti-Hermitian tensor names minus the bare tensor: the sign
    // goes to a Product, as it does for a ^* or ^T that costs one
    auto z_adj = deserialize(L"z⁺{i_1;a_1}:N-A-N");
    REQUIRE(z_adj->is<Product>());
    REQUIRE(z_adj->as<Product>().scalar() == -1);
    REQUIRE(*z_adj == *adjoint(deserialize(L"z{a_1;i_1}:N-A-N")));
  }

  SECTION("set_value_modifier copies bits without touching slots") {
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    Tensor rebuilt(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                   BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    REQUIRE(rebuilt.set_value_modifier(ta.value_modifier()) == 1);
    REQUIRE(rebuilt == ta);
    REQUIRE(rebuilt.hash_value() == ta.hash_value());
  }

  SECTION("with_slots carries both bits") {
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    using ixvec = container::svector<Index>;
    auto w = ta.with_slots(bra<ixvec>{ixvec{Index{L"i_2"}}},
                           ket<ixvec>{ixvec{Index{L"a_2"}}}, aux<ixvec>{});
    REQUIRE(w.value_modifier() == ValueModifier::Adjoint);

    Tensor tt = t;
    REQUIRE(tt.adjoint() == 1);
    REQUIRE(tt.conjugate() == 1);
    REQUIRE(tt.value_modifier() == ValueModifier::Transpose);
    REQUIRE(tt.with_slots(bra<ixvec>{ixvec{Index{L"i_2"}}},
                          ket<ixvec>{ixvec{Index{L"a_2"}}}, aux<ixvec>{})
                .value_modifier() == ValueModifier::Transpose);
  }

  SECTION("ordering: t < t^* < t^T < t⁺, then by slots") {
    Tensor tc = t;
    REQUIRE(tc.conjugate() == 1);
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    Tensor tt = t;
    REQUIRE(tt.adjoint() == 1);
    REQUIRE(tt.conjugate() == 1);
    REQUIRE(tt.value_modifier() == ValueModifier::Transpose);
    REQUIRE(t < tc);
    REQUIRE(tc < tt);
    REQUIRE(tt < ta);

    // slot tie-break: same label and modifier, different bra index
    Tensor t2(L"t", bra{L"a_2"}, ket{L"i_1"}, Symmetry::Nonsymm,
              BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    REQUIRE(t < t2);
  }
}

TEST_CASE("value_oriented_totality", "[conjugation]") {
  // value_oriented() returns the spelling whose slot layout denotes the
  // value directly, for every state that has one
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  Tensor g(L"g", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);

  // Conjugate: the folded (starred + swapped) spelling unfolds back
  Tensor folded = g;
  REQUIRE(folded.transpose() == 1);
  REQUIRE(folded.value_modifier() == ValueModifier::Conjugate);
  {
    auto [vo, sign] = value_oriented(folded);
    REQUIRE(sign == 1);  // Conjugate: the fold and its undo are both free
    REQUIRE(vo == g);
  }
  {
    auto [vo, sign] = value_oriented(g);
    REQUIRE(sign == 1);
    REQUIRE(vo == g);
  }

  // Nonsymm transpose: a pure respelling
  Tensor tt = t;
  REQUIRE(tt.transpose() == 1);
  {
    auto [vo, sign] = value_oriented(tt);
    REQUIRE(sign == 1);
    REQUIRE(vo == t);
  }

  // Nonsymm adjoint: a distinct array, slots as written -- unchanged
  Tensor ta = t;
  REQUIRE(ta.adjoint() == 1);
  {
    auto [vo, sign] = value_oriented(ta);
    REQUIRE(sign == 1);
    REQUIRE(vo == ta);
  }

  // Nonsymm conjugate: no slot spelling -> refuse loudly
  Tensor tc = t;
  REQUIRE(tc.conjugate() == 1);
  REQUIRE_THROWS_AS(value_oriented(tc), sequant::Exception);
}

TEST_CASE("canonicalize_marked_nonsymm_network", "[conjugation]") {
  // A Conjugate/Transpose modifier on a Nonsymm tensor is part of the
  // graph colouring in every canonicalization, so networks that differ
  // only in which factor carries the mark stay distinguishable, and
  // canonicalization is idempotent on them.
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto e1 = deserialize(L"t^*{a_1;i_1}:N-N-N u{i_1;a_1}:N-N-N");
  auto e2 = deserialize(L"t{a_1;i_1}:N-N-N u^*{i_1;a_1}:N-N-N");
  auto c1 = canonicalize(e1->clone());
  auto c2 = canonicalize(e2->clone());
  REQUIRE(*c1 != *c2);
  REQUIRE(*canonicalize(c1->clone()) == *c1);
  REQUIRE(*canonicalize(c2->clone()) == *c2);
}

TEST_CASE("value_modifier_group", "[conjugation]") {
  // ValueModifier is the direct product of its two Z2 factors, with `*` the
  // group operation on all three types
  using CM = ConjugateModifier;
  using TM = TransposeModifier;
  using VM = ValueModifier;
  STATIC_REQUIRE(CM::Yes * CM::Yes == CM::No);
  STATIC_REQUIRE(TM::Yes * TM::Yes == TM::No);
  STATIC_REQUIRE(CM::No * TM::No == VM::None);
  STATIC_REQUIRE(CM::Yes * TM::No == VM::Conjugate);
  STATIC_REQUIRE(CM::No * TM::Yes == VM::Transpose);
  STATIC_REQUIRE(CM::Yes * TM::Yes == VM::Adjoint);
  STATIC_REQUIRE(TM::Yes * CM::Yes == VM::Adjoint);
  STATIC_REQUIRE(VM::Adjoint * VM::Conjugate == VM::Transpose);
  STATIC_REQUIRE(VM::Transpose * VM::Transpose == VM::None);
  STATIC_REQUIRE(conjugate_modifier(VM::Adjoint) == CM::Yes);
  STATIC_REQUIRE(transpose_modifier(VM::Conjugate) == TM::No);

  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  REQUIRE(t.adjoint() == 1);
  REQUIRE(t.conjugate_modifier() == CM::Yes);
  REQUIRE(t.transpose_modifier() == TM::Yes);
  REQUIRE(t.value_modifier() ==
          t.conjugate_modifier() * t.transpose_modifier());
}

TEST_CASE("conjugation_parity_trait", "[conjugation]") {
  // the parity is a stored trait; the observable symmetries are derived from
  // it, the hermiticity and the base field at construction
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);
  SECTION("defaults") {
    Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    REQUIRE(t.conjugation_parity() == ConjugationParity::Even);
    REQUIRE(t.base_field() == Field::Complex);
    REQUIRE(t.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);
  }

  SECTION("odd parity over a real field: imaginary antisymmetric array") {
    Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(p.conjugation_parity() == ConjugationParity::Odd);
    REQUIRE(p.conjugation_symmetry() == ConjugationSymmetry::Antisymm);
    REQUIRE(p.braket_symmetry() == BraKetSymmetry::Antisymm);
    REQUIRE(p.hermiticity() == Hermiticity::Hermitian);
  }

  SECTION("anti-Hermitian over the complex field") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"i_2"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(d.braket_symmetry() == BraKetSymmetry::AntiConjugate);
  }

  SECTION("explicit braket contradicting the traits throws") {
    REQUIRE_THROWS_AS(
        Tensor(L"d", bra{idx(L"i_1", Field::Real)},
               ket{idx(L"i_2", Field::Real)},
               TensorSymmetries{.braket = BraKetSymmetry::Symm,
                                .hermiticity = Hermiticity::AntiHermitian}),
        sequant::Exception);
  }

  SECTION("with_slots carries the parity") {
    Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    using ixvec = container::svector<Index>;
    auto q =
        p.with_slots(bra<ixvec>{ixvec{idx(L"i_3", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"i_4", Field::Real)}}, aux<ixvec>{});
    REQUIRE(q.conjugation_parity() == ConjugationParity::Odd);
    REQUIRE(q.braket_symmetry() == BraKetSymmetry::Antisymm);
  }

  SECTION("aux slots do not enter the conjugation field") {
    // the elementwise conjugation relation is stated over the bra/ket basis;
    // aux slots are array-like and pair nothing, so a tensor carrying aux
    // slots alone has no basis to state the relation over and asserts none
    Tensor w(L"w", bra{}, ket{}, aux{L"p_1"}, Symmetry::Nonsymm);
    REQUIRE(w.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);
  }

  SECTION(
      "all slots (bra, ket, aux) real-field assert Symm conjugation "
      "symmetry") {
    Tensor t(L"t", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             aux{idx(L"p_1", Field::Real)}, Symmetry::Nonsymm);
    REQUIRE(t.conjugation_symmetry() == ConjugationSymmetry::Symm);
  }

  SECTION("the conjugation field is the bra/ket field, not the aux field") {
    // an aux slot's own field must not feed base_field(): a real-field
    // bra/ket pair alongside a complex-field aux index still reads Symm
    // under Even parity, as if the aux slot were not there
    Tensor t(L"t", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             aux{idx(L"p_1", Field::Complex)}, Symmetry::Nonsymm);
    REQUIRE(t.conjugation_symmetry() == ConjugationSymmetry::Symm);

    // twin: swap the ket to complex-field. base_field() is the OR of the
    // bra/ket fields, so it is now Complex and Even parity reads Nonsymm --
    // the aux slot's (still complex) field plays no part in the change.
    Tensor u(L"u", bra{idx(L"i_1", Field::Real)},
             ket{idx(L"i_2", Field::Complex)}, aux{idx(L"p_1", Field::Complex)},
             Symmetry::Nonsymm);
    REQUIRE(u.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);
  }

  SECTION("with_slots re-derives the field-dependent symmetries") {
    // real-field, Odd parity, Hermitian: braket Antisymm, conjugation Antisymm
    Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(p.braket_symmetry() == BraKetSymmetry::Antisymm);
    REQUIRE(p.conjugation_symmetry() == ConjugationSymmetry::Antisymm);

    // rebuild onto default (complex) indices: the traits (hermiticity,
    // parity) are carried, but the field-dependent symmetries must be
    // re-derived against the new (complex) field, not copied verbatim
    using ixvec = container::svector<Index>;
    auto q = p.with_slots(bra<ixvec>{ixvec{Index{L"i_3"}}},
                          ket<ixvec>{ixvec{Index{L"i_4"}}}, aux<ixvec>{});
    REQUIRE(q.conjugation_parity() == ConjugationParity::Odd);
    REQUIRE(q.hermiticity() == Hermiticity::Hermitian);
    REQUIRE(q.braket_symmetry() == BraKetSymmetry::Conjugate);
    REQUIRE(q.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);
  }

  SECTION("explicit braket consistent with the traits constructs") {
    Tensor d(L"d", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.braket = BraKetSymmetry::Antisymm,
                              .hermiticity = Hermiticity::Hermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(d.hermiticity() == Hermiticity::Hermitian);
  }

  SECTION("tensor with no slots at all asserts no conjugation symmetry") {
    Tensor c(L"c", bra{}, ket{}, aux{});
    REQUIRE(c.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);
  }

  SECTION(
      "an explicit Conjugate pin over a real field back-fills parity None") {
    // the pin asserts the adjoint relation and no reality, so the parity it
    // back-fills is None, and re-deriving the exchange symmetry from the
    // traits reproduces the pin
    Tensor t(L"t", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
             TensorSymmetries{.braket = BraKetSymmetry::Conjugate});
    REQUIRE(t.conjugation_parity() == ConjugationParity::None);
    REQUIRE(t.hermiticity() == Hermiticity::Hermitian);
    REQUIRE(t.braket_symmetry() == BraKetSymmetry::Conjugate);
    REQUIRE(t.conjugation_symmetry() == ConjugationSymmetry::Nonsymm);

    using ixvec = container::svector<Index>;
    auto u =
        t.with_slots(bra<ixvec>{ixvec{idx(L"i_2", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"a_2", Field::Real)}}, aux<ixvec>{});
    REQUIRE(u.braket_symmetry() == BraKetSymmetry::Conjugate);
  }

  SECTION("with_slots normalizes the modifier bits against the new field") {
    // complex-field Hermitian (Conjugate): the starred spelling is a state of
    // its own, the folded orientation
    Tensor h(L"h", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    REQUIRE(h.conjugate() == 1);
    REQUIRE(h.value_modifier() == ValueModifier::Conjugate);
    // rebuilt onto real-field slots the tensor is Symm, where a set bit
    // denotes nothing: it normalizes away, as it does in a fresh construction
    using ixvec = container::svector<Index>;
    auto r =
        h.with_slots(bra<ixvec>{ixvec{idx(L"i_2", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"a_2", Field::Real)}}, aux<ixvec>{});
    REQUIRE(r.braket_symmetry() == BraKetSymmetry::Symm);
    REQUIRE(r.value_modifier() == ValueModifier::None);
    Tensor fresh(L"h", bra{idx(L"i_2", Field::Real)},
                 ket{idx(L"a_2", Field::Real)},
                 TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    REQUIRE(r == fresh);
    REQUIRE(r.hash_value() == fresh.hash_value());
  }

  SECTION("with_slots refuses a rebuild whose normalization costs a sign") {
    // complex-field anti-Hermitian of odd parity (AntiConjugate): the starred
    // spelling is a state. Over a real field the odd parity makes
    // conj(T) = -T, so the bit would be consumed at -1, which no Tensor can
    // hold
    Tensor z(L"z", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(z.conjugate() == 1);
    REQUIRE(z.value_modifier() == ValueModifier::Conjugate);
    using ixvec = container::svector<Index>;
    REQUIRE_THROWS_AS(
        z.with_slots(bra<ixvec>{ixvec{idx(L"i_2", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"a_2", Field::Real)}}, aux<ixvec>{}),
        Exception);
  }
}

TEST_CASE("signed_normalization", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  SECTION(
      "anti-Hermitian over the complex field: adjoint is minus the tensor") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"i_2"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(d.braket_symmetry() == BraKetSymmetry::AntiConjugate);
    Tensor dt = d;
    // T^T{q;p} = T{p;q} is a respelling, but the only spelling of the
    // transpose available to an AntiConjugate tensor is the starred one,
    // and `T{q;p} = -conj(T{p;q})` prices that fold at -1
    REQUIRE(dt.transpose() == -1);
    REQUIRE(dt.value_modifier() == ValueModifier::Conjugate);  // folded
    Tensor dc = d;
    REQUIRE(dc.conjugate_transpose() == -1);
    REQUIRE(dc.value_modifier() == ValueModifier::None);
    REQUIRE(dc.bra()[0].label() == L"i_2");  // swapped
    // Expr::adjoint returns the sign as its byproduct
    Tensor da = d;
    REQUIRE(da.adjoint() == -1);
    REQUIRE(da == dc);
    // the free function absorbs it into a scalar
    auto adj = adjoint(ex<Tensor>(d));
    REQUIRE(adj->is<Product>());
    REQUIRE(adj->as<Product>().scalar() == -1);
    REQUIRE(adj->as<Product>().factors().size() == 1);
    REQUIRE(adj->as<Product>().factor(0)->as<Tensor>().value_modifier() ==
            ValueModifier::None);
    // through a Product the sign lands in the product's scalar
    auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_2"}, Symmetry::Nonsymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    auto padj = adjoint(ex<Tensor>(d) * u);
    REQUIRE(padj->is<Product>());
    REQUIRE(padj->as<Product>().scalar() == -1);
    REQUIRE(padj->as<Product>().factors().size() == 2);
    // and through a Sum it wraps the affected summand
    auto sadj = adjoint(ex<Tensor>(d) + u);
    REQUIRE(sadj->is<Sum>());
    REQUIRE(sadj->as<Sum>().summand(0)->is<Product>());
    REQUIRE(sadj->as<Sum>().summand(0)->as<Product>().scalar() == -1);
  }

  SECTION("odd parity over a real field: conjugation is a sign") {
    Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::NonHermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(p.conjugation_symmetry() == ConjugationSymmetry::Antisymm);
    Tensor pc = p;
    REQUIRE(pc.conjugate() == -1);
    REQUIRE(pc.value_modifier() == ValueModifier::None);
    // transpose of an array with known conjugation is represented as the
    // adjoint (the ⁺ spelling), with the conjugation's sign
    Tensor pt = p;
    REQUIRE(pt.transpose() == -1);
    REQUIRE(pt.value_modifier() == ValueModifier::Adjoint);
    auto c = conjugate(ex<Tensor>(p));
    REQUIRE(c->is<Product>());
    REQUIRE(c->as<Product>().scalar() == -1);
  }

  SECTION("even parity over a real field: conjugation is the identity") {
    Tensor t(L"t", bra{idx(L"a_1", Field::Real)}, ket{idx(L"i_1", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::NonHermitian});
    Tensor tc = t;
    REQUIRE(tc.conjugate() == 1);
    REQUIRE(tc == t);
    Tensor ta = t;
    REQUIRE(ta.conjugate_transpose() == 1);
    REQUIRE(ta.value_modifier() == ValueModifier::Adjoint);
    REQUIRE(ta.decorated_label() == L"t⁺");
  }

  SECTION("imaginary Hermitian over a real field: antisymmetric array") {
    Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                              .conjugation_parity = ConjugationParity::Odd});
    REQUIRE(p.braket_symmetry() == BraKetSymmetry::Antisymm);
    Tensor pt = p;
    REQUIRE(pt.transpose() == -1);
    REQUIRE(pt.value_modifier() == ValueModifier::None);
    Tensor pa = p;
    REQUIRE(pa.conjugate_transpose() == 1);  // Hermitian: adjoint is itself
    REQUIRE(pa.value_modifier() == ValueModifier::None);
  }

  SECTION("a ⁺ label on an anti-Hermitian tensor cannot be adopted") {
    REQUIRE_THROWS_AS(
        Tensor(L"d⁺", bra{L"i_1"}, ket{L"i_2"},
               TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian}),
        sequant::Exception);
  }

  SECTION("value_oriented reports the sign") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"i_2"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    // the folded spelling d^*{i_2;i_1} is -d{i_1;i_2}: unfolding it back to
    // the value orientation costs the anti-Hermitian sign, exactly as the
    // fold did
    Tensor folded = d;
    REQUIRE(folded.transpose() == -1);
    auto [vo, sign] = value_oriented(folded);
    REQUIRE(sign == -1);
    REQUIRE(vo == d);
    // deserialized d^* over the complex field is conj(d) = -d^T: unfolding
    // the conjugation costs the anti-Hermitian sign
    Tensor dstar = d;
    REQUIRE(dstar.conjugate() == 1);
    auto [vo2, sign2] = value_oriented(dstar);
    REQUIRE(sign2 == -1);
    REQUIRE(vo2.value_modifier() == ValueModifier::None);
    REQUIRE(vo2.bra()[0].label() == L"i_2");
  }
}

TEST_CASE("canonicalize_signed_braket", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // N.B. the tensors below carry the default ColumnSymmetry::Nonsymm:
  // TensorNetworkV3's graph-dictated bra<->ket reorientation folds the two
  // bundles of every braket_foldable() tensor, column-symmetric or not (only
  // permuting slots *within* a bundle needs the column symmetry)

  SECTION("anti-Hermitian: the two orientations differ by a sign") {
    // d{i;a} u{a;i} and d{a;i} u{a;i}: d{a;i} = -conj(d{i;a}), so the
    // canonical forms differ by -1 and a conjugation marker on d
    auto d = [](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(
          L"d", bra{b}, ket{k},
          TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    };
    auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"});
    auto e1 = d(L"i_1", L"a_1") * u;
    auto e2 = d(L"a_1", L"i_1") * u;
    auto c1 = canonicalize(e1->clone());
    auto c2 = canonicalize(e2->clone());
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    // canonical forms are idempotent
    REQUIRE(*canonicalize(c1->clone()) == *c1);
    REQUIRE(*canonicalize(c2->clone()) == *c2);
    // exactly one of the two carries a conjugation marker on d, and the
    // product scalars differ by the sign of the fold
    auto d_of = [](const ExprPtr& p) {
      for (auto& f : p->as<Product>().factors())
        if (f->as<Tensor>().label() == L"d") return f->as<Tensor>();
      throw Exception("test: no factor labelled d");
    };
    REQUIRE(d_of(c1).conjugated() != d_of(c2).conjugated());
    REQUIRE(c1->as<Product>().scalar() == -c2->as<Product>().scalar());
  }

  SECTION("anti-Hermitian, Complete method: the lexicographic refold signs") {
    // the test binary pins Topological; run the same check under the
    // library's default Complete so the post-relabel refold loop in
    // TensorNetworkV3::canonicalize, whose sign reaches the byproduct
    // separately from canonicalize_graph's, is exercised too
    const CanonicalizeOptions opts{.method = CanonicalizationMethod::Complete};
    auto d = [](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(
          L"d", bra{b}, ket{k},
          TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    };
    auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"});
    auto c1 = canonicalize(d(L"i_1", L"a_1") * u, opts);
    auto c2 = canonicalize(d(L"a_1", L"i_1") * u, opts);
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    REQUIRE(*canonicalize(c1->clone(), opts) == *c1);
    REQUIRE(*canonicalize(c2->clone(), opts) == *c2);
    auto d_of = [](const ExprPtr& p) {
      for (auto& f : p->as<Product>().factors())
        if (f->as<Tensor>().label() == L"d") return f->as<Tensor>();
      throw Exception("test: no factor labelled d");
    };
    REQUIRE(d_of(c1).conjugated() != d_of(c2).conjugated());
    REQUIRE(c1->as<Product>().scalar() == -c2->as<Product>().scalar());
  }

  SECTION("Hermitian, default column symmetry: one canonical form") {
    // g{i;a} u{a;i} and g{a;i} u{a;i}: g{a;i} = conj(g{i;a}), so the two
    // spellings are one value up to a conjugation marker on g. The
    // graph-dictated reorientation folds g's bundles although g is not
    // column-symmetric; only permuting slots *within* a bundle needs that
    auto g = [](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(
          L"g", bra{b}, ket{k},
          TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    };
    REQUIRE(g(L"i_1", L"a_1")->as<Tensor>().column_symmetry() ==
            ColumnSymmetry::Nonsymm);
    auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"});
    auto c1 = canonicalize(g(L"i_1", L"a_1") * u);
    auto c2 = canonicalize(g(L"a_1", L"i_1") * u);
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    auto g_of = [](const ExprPtr& p) {
      for (auto& f : p->as<Product>().factors())
        if (f->as<Tensor>().label() == L"g") return f->as<Tensor>();
      throw Exception("test: no factor labelled g");
    };
    // one slot spelling of g, exactly one marker, the same scalar
    REQUIRE(g_of(c1).bra()[0].label() == g_of(c2).bra()[0].label());
    REQUIRE(g_of(c1).ket()[0].label() == g_of(c2).ket()[0].label());
    REQUIRE(g_of(c1).conjugated() != g_of(c2).conjugated());
    REQUIRE(c1->as<Product>().scalar() == c2->as<Product>().scalar());
  }

  SECTION("real-field odd-parity Hermitian tensor: antisymmetric, no marker") {
    auto p = [&](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(
          L"p", bra{idx(b, Field::Real)}, ket{idx(k, Field::Real)},
          TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                           .conjugation_parity = ConjugationParity::Odd});
    };
    auto u = [&](std::wstring_view b, std::wstring_view k) {
      return ex<Tensor>(L"u", bra{idx(b, Field::Real)},
                        ket{idx(k, Field::Real)});
    };
    auto c1 = canonicalize(p(L"i_1", L"a_1") * u(L"a_1", L"i_1"));
    auto c2 = canonicalize(p(L"a_1", L"i_1") * u(L"a_1", L"i_1"));
    REQUIRE(c1->as<Product>().scalar() == -c2->as<Product>().scalar());
    for (auto& c : {c1, c2})
      for (auto& f : c->as<Product>().factors())
        REQUIRE(f->as<Tensor>().value_modifier() == ValueModifier::None);
    REQUIRE(*canonicalize(c1->clone()) == *c1);
  }
}

TEST_CASE("symmetries_carry_through_slot_rebuilds", "[conjugation]") {
  // Tensor::symmetries() hands the field-agnostic traits to a rebuild, which
  // derives the field-dependent symmetries from them again; a rebuild that
  // forwards the observable braket symmetry instead loses the parity, and
  // with it the array's elementwise conjugation relation
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // real-field, Hermitian, odd parity: an imaginary Hermitian array, whose
  // exchange symmetry and conjugation symmetry are both Antisymm
  Tensor p(L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
           TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                            .conjugation_parity = ConjugationParity::Odd});
  REQUIRE(p.braket_symmetry() == BraKetSymmetry::Antisymm);
  REQUIRE(p.conjugation_symmetry() == ConjugationSymmetry::Antisymm);

  SECTION("the pack carries the traits and pins no exchange symmetry") {
    const auto syms = p.symmetries();
    REQUIRE(syms.perm == p.symmetry());
    REQUIRE(syms.hermiticity == Hermiticity::Hermitian);
    REQUIRE(syms.conjugation_parity == ConjugationParity::Odd);
    REQUIRE(syms.column == p.column_symmetry());
    REQUIRE_FALSE(syms.braket.has_value());
  }

  SECTION("expand_antisymm keeps the parity") {
    auto expanded = mbpt::expand_antisymm(p);
    REQUIRE(expanded->is<Tensor>());
    const auto& q = expanded->as<Tensor>();
    REQUIRE(q.symmetry() == Symmetry::Nonsymm);  // the attribute it changes
    REQUIRE(q.hermiticity() == Hermiticity::Hermitian);
    REQUIRE(q.conjugation_parity() == ConjugationParity::Odd);
    REQUIRE(q.conjugation_symmetry() == ConjugationSymmetry::Antisymm);
    REQUIRE(q.braket_symmetry() == BraKetSymmetry::Antisymm);
  }

  SECTION("a pinned exchange symmetry the traits do not derive is carried") {
    // the Symm pin over the complex basis: the traits it back-fills
    // (Hermitian, Even) derive Conjugate over that basis, so the pack carries
    // the pin itself, or a rebuild would change the tensor's symmetry
    Tensor g(L"g", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
             TensorSymmetries{.braket = BraKetSymmetry::Symm});
    REQUIRE(g.braket_symmetry() == BraKetSymmetry::Symm);
    const auto syms = g.symmetries();
    REQUIRE(syms.braket == BraKetSymmetry::Symm);
    REQUIRE(syms.conjugation_parity == g.conjugation_parity());

    using ixvec = container::svector<Index>;
    auto r = g.with_slots(bra<ixvec>{ixvec{Index{L"i_2"}}},
                          ket<ixvec>{ixvec{Index{L"a_2"}}}, aux<ixvec>{});
    REQUIRE(r.braket_symmetry() == BraKetSymmetry::Symm);
    REQUIRE(r.hermiticity() == g.hermiticity());
    REQUIRE(r.conjugation_parity() == g.conjugation_parity());

    // and through a domain rebuild
    auto expanded = mbpt::expand_antisymm(
        deserialize(L"g{i_1,i_2;a_1,a_2}:A-S-S")->as<Tensor>());
    std::size_t n = 0;
    expanded->visit(
        [&n](const ExprPtr& e) {
          if (!e->is<Tensor>()) return;
          ++n;
          CHECK(e->as<Tensor>().braket_symmetry() == BraKetSymmetry::Symm);
        },
        /* atoms_only = */ true);
    REQUIRE(n > 1);
  }
}

TEST_CASE("signed_eval_boundary", "[conjugation]") {
  // the sign the value respelling consumes must reach the eval tree: a leaf
  // whose spelling is minus its value orientation lowers to that orientation
  // and a Constant(-1) factor, the only channel that reaches a value (a node
  // phase is a cache-orientation round trip)
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  SECTION("an AntiConjugate leaf's starred spelling carries the fold's sign") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                              .column = ColumnSymmetry::Symm});
    REQUIRE(d.braket_symmetry() == BraKetSymmetry::AntiConjugate);
    Tensor dstar = d;
    REQUIRE(dstar.conjugate() == 1);
    REQUIRE(dstar.value_modifier() == ValueModifier::Conjugate);
    // d^* = -d^T: unfolding the marker back to the value orientation costs
    // the anti-Hermitian sign
    auto [vo, vo_sign] = value_oriented(dstar);
    REQUIRE(vo_sign == -1);
    REQUIRE(vo.value_modifier() == ValueModifier::None);

    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(dstar));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE_FALSE(tree.leaf());
    REQUIRE(tree->op_type() == EvalOp::Product);
    REQUIRE(tree.right()->is_constant());
    REQUIRE(tree.right()->as_constant().value<int>() == -1);
    // the operand is the value orientation, a plain leaf; its block
    // canonicalization contributes no phase of its own
    REQUIRE(tree.left().leaf());
    REQUIRE(tree.left()->canon_phase() == 1);
    REQUIRE_FALSE(tree.left()->as_tensor().conjugated());
    REQUIRE(tree.left()->as_tensor().bra()[0].label() == L"a_1");
  }

  SECTION("inside a product the marked leaf's sign joins the scalar") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    Tensor dstar = d;
    REQUIRE(dstar.conjugate() == 1);
    auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"});

    // -1 d^*{i;a} u{a;i}: the scalar the canonicalizer carries beside the
    // marked spelling cancels the leaf's sign, so nothing is scaled at run
    // time; the tree is the bare contraction of the value orientation
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Product>(-1, ExprPtrList{ex<Tensor>(dstar), u}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree->op_type() == EvalOp::Product);
    REQUIRE(tree.left().leaf());
    REQUIRE(tree.right().leaf());
    REQUIRE_FALSE(tree.left()->as_tensor().conjugated());
    REQUIRE(tree.left()->as_tensor().bra()[0].label() == L"a_1");

    // d^*{i;a} u{a;i} alone: one scale, by the product scalar, over the
    // contraction, not one per marked factor
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree2 = binarize(ex<Product>(ExprPtrList{ex<Tensor>(dstar), u}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree2->op_type() == EvalOp::Product);
    REQUIRE(tree2.right()->is_constant());
    REQUIRE(tree2.right()->as_constant().value<int>() == -1);
    REQUIRE(tree2.left()->op_type() == EvalOp::Product);
    REQUIRE(tree2.left().left().leaf());
    REQUIRE(tree2.left().right().leaf());
  }

  SECTION("a plain Adjoint-state Nonsymm leaf carries no sign") {
    // nothing folds for a Nonsymm tensor, so there is no relation to spend
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);
    REQUIRE(t_adj.value_modifier() == ValueModifier::Adjoint);

    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree->op_type() == EvalOp::Adjoint);
    REQUIRE(tree->canon_phase() == 1);
  }
}

TEST_CASE("adjoint_sign_is_absolute", "[conjugation]") {
  // the sign an adjoint produces belongs to the expression, not to the
  // spelling it was taken in: adjoining a canonical form and canonicalizing
  // an adjoint must land on the same signed expression
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  auto d =
      ex<Tensor>(L"d", bra{L"i_1"}, ket{L"a_1"},
                 TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                                  .column = ColumnSymmetry::Symm});
  auto u = ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"});
  auto x = d * u;

  auto adj_then_canon = canonicalize(adjoint(x->clone()));
  auto canon_then_adj = canonicalize(adjoint(canonicalize(x->clone())));
  REQUIRE(*adj_then_canon == *canon_then_adj);
}

TEST_CASE("bad_conjugation_parity_letter", "[conjugation]") {
  // the fourth symmetry letter of an annotation names a ConjugationParity
  // (E, O or N); anything else is a deserialization error
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  REQUIRE_NOTHROW(deserialize(L"t{i_1;i_2}:N-N-N-O"));
  REQUIRE_THROWS_AS(deserialize(L"t{i_1;i_2}:N-N-N-X"),
                    io::serialization::SerializationError);
}

TEST_CASE("painter_colours_by_conjugation_symmetry", "[conjugation]") {
  // the graph painter's shade must follow the observable
  // ConjugationSymmetry -- the property Tensor::static_equal compares -- and
  // not the ConjugationParity trait it is derived from, so that tensors which
  // compare equal also colour, and hence canonicalize, equally
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  SECTION("parities that resolve to the same symmetry colour alike") {
    // over the complex field neither Even nor None exposes an elementwise
    // conjugation relation, so both resolve to ConjugationSymmetry::Nonsymm
    auto c = [](ConjugationParity parity) {
      return ex<Tensor>(L"c", bra{L"i_1"}, ket{L"a_1"},
                        TensorSymmetries{.braket = BraKetSymmetry::Conjugate,
                                         .conjugation_parity = parity,
                                         .column = ColumnSymmetry::Symm});
    };
    auto even = c(ConjugationParity::Even);
    auto none = c(ConjugationParity::None);
    REQUIRE(even->as<Tensor>().conjugation_parity() == ConjugationParity::Even);
    REQUIRE(none->as<Tensor>().conjugation_parity() == ConjugationParity::None);
    REQUIRE(even->as<Tensor>().conjugation_symmetry() ==
            ConjugationSymmetry::Nonsymm);
    REQUIRE(none->as<Tensor>().conjugation_symmetry() ==
            ConjugationSymmetry::Nonsymm);
    REQUIRE(*even == *none);
    REQUIRE(even->hash_value() == none->hash_value());
    auto u = [] { return ex<Tensor>(L"u", bra{L"a_1"}, ket{L"i_1"}); };
    REQUIRE(*canonicalize(even->clone() * u()) ==
            *canonicalize(none->clone() * u()));
  }

  SECTION("a real-field Odd tensor colours apart from its Even twin") {
    // the two differ in nothing but the parity, which over a real field is an
    // observable elementwise relation (Antisymm vs Symm)
    auto p = [](ConjugationParity parity) {
      return ex<Tensor>(
          L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
          TensorSymmetries{.hermiticity = Hermiticity::NonHermitian,
                           .conjugation_parity = parity,
                           .column = ColumnSymmetry::Symm});
    };
    auto u = [] {
      return ex<Tensor>(L"u", bra{idx(L"a_1", Field::Real)},
                        ket{idx(L"i_1", Field::Real)});
    };
    REQUIRE(p(ConjugationParity::Odd)->as<Tensor>().braket_symmetry() ==
            p(ConjugationParity::Even)->as<Tensor>().braket_symmetry());
    auto odd = canonicalize(p(ConjugationParity::Odd) * u());
    auto even = canonicalize(p(ConjugationParity::Even) * u());
    REQUIRE_FALSE(*odd == *even);
  }

  SECTION("a real-field None-parity tensor colours apart from its Even twin") {
    // over a real field Even resolves to Symm and None to Nonsymm: the two
    // differ in the observable conjugation symmetry although neither exposes
    // a sign, so the shade must separate them as static_equal does
    auto p = [](ConjugationParity parity) {
      return ex<Tensor>(
          L"p", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
          TensorSymmetries{.hermiticity = Hermiticity::NonHermitian,
                           .conjugation_parity = parity,
                           .column = ColumnSymmetry::Symm});
    };
    auto u = [] {
      return ex<Tensor>(L"u", bra{idx(L"a_1", Field::Real)},
                        ket{idx(L"i_1", Field::Real)});
    };
    REQUIRE(p(ConjugationParity::None)->as<Tensor>().conjugation_symmetry() ==
            ConjugationSymmetry::Nonsymm);
    REQUIRE(p(ConjugationParity::Even)->as<Tensor>().conjugation_symmetry() ==
            ConjugationSymmetry::Symm);
    REQUIRE(p(ConjugationParity::None)->as<Tensor>().braket_symmetry() ==
            p(ConjugationParity::Even)->as<Tensor>().braket_symmetry());
    auto none = canonicalize(p(ConjugationParity::None) * u());
    auto even = canonicalize(p(ConjugationParity::Even) * u());
    REQUIRE_FALSE(*none == *even);
  }
}
