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
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/v3.hpp>

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
      (t.*op)();
      REQUIRE(t != T0);
      (t.*op)();
      REQUIRE(t == T0);
      REQUIRE(t.hash_value() == T0.hash_value());
    }
  }

  SECTION("transpose swaps slots and sets the bit") {
    Tensor t = T0;
    t.transpose();
    REQUIRE(t.value_modifier() == ValueModifier::Transpose);
    REQUIRE(t.bra()[0].label() == L"i_1");
    REQUIRE(t.ket()[0].label() == L"a_1");
    REQUIRE(serialize(ex<Tensor>(t), {.annot_symm = true}) ==
            L"t^T{i_1;a_1}:N-N-S");
  }

  SECTION("adjoint = transpose o conjugate = conjugate o transpose") {
    Tensor a = T0;
    a.adjoint();
    Tensor tc = T0;
    tc.transpose();
    tc.conjugate();
    Tensor ct = T0;
    ct.conjugate();
    ct.transpose();
    REQUIRE(a == tc);
    REQUIRE(a == ct);
    REQUIRE(a.value_modifier() == ValueModifier::Adjoint);
    Tensor ctt = T0;
    ctt.conjugate_transpose();
    REQUIRE(a == ctt);
  }

  SECTION("conj(adjoint(t)) is the transpose") {
    Tensor t = T0;
    t.adjoint();
    t.conjugate();
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
    gt.transpose();
    REQUIRE(gt.value_modifier() == ValueModifier::Conjugate);
    REQUIRE(gt.bra()[0].label() == L"a_1");
    REQUIRE(serialize(ex<Tensor>(gt), {.annot_symm = true}) ==
            L"g^*{a_1;i_1}:N-C-S");
    // and adjoint() is a pure swap
    Tensor ga = g;
    ga.adjoint();
    REQUIRE(ga.value_modifier() == ValueModifier::None);
    REQUIRE(ga.bra()[0].label() == L"a_1");
    // transpose() again unfolds
    gt.transpose();
    REQUIRE(gt == g);
  }

  SECTION("Symm: every modifier is the identity") {
    Tensor sc = s;
    sc.conjugate();
    REQUIRE(sc == s);
    Tensor st = s;
    st.transpose();
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
    gt.set_value_modifier(ValueModifier::Transpose);
    REQUIRE(gt.value_modifier() == ValueModifier::Conjugate);
    Tensor sa = s;
    sa.set_value_modifier(ValueModifier::Adjoint);
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
  t->as<Tensor>().conjugate();
  auto rt = deserialize(serialize(t));
  REQUIRE(rt->as<Tensor>().conjugated());
  REQUIRE(*rt == *t);
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
    t->as<Tensor>().conjugate();
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
  // decorated_label() reproduces that spelling for printing/hashing.
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
    ta.adjoint();
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
    ta.adjoint();
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
  }

  SECTION("set_value_modifier copies bits without touching slots") {
    Tensor ta = t;
    ta.adjoint();
    Tensor rebuilt(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                   BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    rebuilt.set_value_modifier(ta.value_modifier());
    REQUIRE(rebuilt == ta);
    REQUIRE(rebuilt.hash_value() == ta.hash_value());
  }

  SECTION("with_slots carries both bits") {
    Tensor ta = t;
    ta.adjoint();
    using ixvec = container::svector<Index>;
    auto w = ta.with_slots(bra<ixvec>{ixvec{Index{L"i_2"}}},
                           ket<ixvec>{ixvec{Index{L"a_2"}}}, aux<ixvec>{});
    REQUIRE(w.value_modifier() == ValueModifier::Adjoint);

    Tensor tt = t;
    tt.adjoint();
    tt.conjugate();
    REQUIRE(tt.value_modifier() == ValueModifier::Transpose);
    REQUIRE(tt.with_slots(bra<ixvec>{ixvec{Index{L"i_2"}}},
                          ket<ixvec>{ixvec{Index{L"a_2"}}}, aux<ixvec>{})
                .value_modifier() == ValueModifier::Transpose);
  }

  SECTION("ordering: t < t^* < t^T < t⁺, then by slots") {
    Tensor tc = t;
    tc.conjugate();
    Tensor ta = t;
    ta.adjoint();
    Tensor tt = t;
    tt.adjoint();
    tt.conjugate();
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
  folded.transpose();
  REQUIRE(folded.value_modifier() == ValueModifier::Conjugate);
  REQUIRE(value_oriented(folded) == g);
  REQUIRE(value_oriented(g) == g);

  // Nonsymm transpose: a pure respelling
  Tensor tt = t;
  tt.transpose();
  REQUIRE(value_oriented(tt) == t);

  // Nonsymm adjoint: a distinct array, slots as written -- unchanged
  Tensor ta = t;
  ta.adjoint();
  REQUIRE(value_oriented(ta) == ta);

  // Nonsymm conjugate: no slot spelling -> refuse loudly
  Tensor tc = t;
  tc.conjugate();
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
  t.adjoint();
  REQUIRE(t.conjugate_modifier() == CM::Yes);
  REQUIRE(t.transpose_modifier() == TM::Yes);
  REQUIRE(t.value_modifier() ==
          t.conjugate_modifier() * t.transpose_modifier());
}
