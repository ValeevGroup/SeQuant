//
// Created by Kshitij Surjuse on 2026-08-31.
//

// The conjugation case catalogue: every symbolic identity involving complex
// conjugation that SeQuant is expected to honor, one TEST_CASE per identity
// family. Canonicalization- and network-level cases live at the bottom and
// grow with the conjugation-symbolic work.

#include <SeQuant/core/utility/macros.hpp>
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
#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/reserved.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network/v3.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/core/utility/string.hpp>

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
  // _first_, then conjugate, then verify the hash actually changed and that
  // toggling back restores it
  auto v = Variable(L"x");
  const auto h0 = v.hash_value();  // memoize
  v.conjugate();
  REQUIRE(v.conjugated());
  REQUIRE(v.hash_value() != h0);
  v.conjugate();
  REQUIRE(v.hash_value() == h0);
}

TEST_CASE("conjugate_free_function_total", "[conjugation]") {
  // sequant::conjugate is the conjugate of the _value_: the adjoint on c-number
  // content, an involution on each node kind; operator-valued content, which
  // has no value, is rejected loudly
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  // Sum + Product distribution: (c A B)* = conj(c) B⁺ A⁺
  auto A = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm);
  auto B = ex<Tensor>(L"B", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm);
  auto prod = ex<Constant>(i_unit) * A->clone() * B->clone();
  auto pc = conjugate(prod);
  const auto& p = pc->as<Product>();
  REQUIRE(p.scalar() == (C{0, -1}));
  REQUIRE(p.factors().size() == 2);
  // the factors are reversed, as the adjoint of a product reverses them
  REQUIRE(p.factors()[0]->as<Tensor>().label() == L"B");
  REQUIRE(p.factors()[0]->as<Tensor>().adjointed());
  REQUIRE(p.factors()[1]->as<Tensor>().label() == L"A");
  REQUIRE(p.factors()[1]->as<Tensor>().adjointed());
  REQUIRE(*conjugate(pc) == *prod);

  auto sum = A->clone() + B->clone();
  auto sc = conjugate(sum);
  for (const auto& s : *sc) REQUIRE(s->as<Tensor>().adjointed());
  REQUIRE(*conjugate(sc) == *sum);

  // Re/Im are real-valued: conj is the identity on them
  auto re = real_part(A->clone() * B->clone());
  REQUIRE(*conjugate(re) == *re);

  // an operator string carries no value: sequant::kconjugate is its conjugate
  REQUIRE_THROWS_AS(conjugate(ex<FNOperator>(cre({L"i_1"}), ann({L"a_1"}))),
                    Exception);
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

  // Conjugate: the two orientations are two values (T{q;p} = conj(T{p;q})),
  // so no fold applies
  Tensor h(L"h", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  REQUIRE_FALSE(braket_foldable(h));

  // Nonsymm: no fold of any kind
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  REQUIRE_FALSE(braket_foldable(t));

  // Symm (Hermitian over a real basis): the exchange is a free respelling
  Tensor s(L"s", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
           Symmetry::Nonsymm, BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  REQUIRE(braket_foldable(s));

  // Antisymm (anti-Hermitian over a real basis): a respelling carrying -1
  Tensor n(L"n", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
           TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                            .column = ColumnSymmetry::Symm});
  REQUIRE(n.braket_symmetry() == BraKetSymmetry::Antisymm);
  REQUIRE(braket_foldable(n));

  // AntiConjugate (anti-Hermitian over the complex basis): the two
  // orientations are two values as well, so no fold applies
  Tensor d(L"d", bra{L"i_1"}, ket{L"a_1"},
           TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                            .column = ColumnSymmetry::Symm});
  REQUIRE(d.braket_symmetry() == BraKetSymmetry::AntiConjugate);
  REQUIRE_FALSE(braket_foldable(d));

  // operator-valued whose braket symmetry is Conjugate (this one's creator
  // and annihilator index multisets agree, hence Hermitian, hence Conjugate
  // over the complex basis): not foldable, like every other Conjugate
  FNOperator op(cre({L"i_1"}), ann({L"i_1"}));
  REQUIRE(braket_symmetry(op) == BraKetSymmetry::Conjugate);
  REQUIRE_FALSE(braket_foldable(op));

  // operator-valued over a real basis: Hermitian (equal creator and
  // annihilator index multisets) derives Symm, and the free bra<->ket
  // exchange is a symmetry of the operator, so it is foldable; the exchange
  // is NormalOperator::_swap_bra_ket
  FNOperator rop(cre({idx(L"p_1", Field::Real), idx(L"p_2", Field::Real)}),
                 ann({idx(L"p_1", Field::Real), idx(L"p_2", Field::Real)}));
  REQUIRE(braket_symmetry(rop) == BraKetSymmetry::Symm);
  REQUIRE(braket_foldable(rop));
  // a non-Hermitian one is not
  FNOperator nop(cre({idx(L"p_1", Field::Real), idx(L"p_2", Field::Real)}),
                 ann({idx(L"p_3", Field::Real), idx(L"p_4", Field::Real)}));
  REQUIRE(braket_symmetry(nop) == BraKetSymmetry::Nonsymm);
  REQUIRE_FALSE(braket_foldable(nop));

  // a reserved bookkeeping operator's orientation defines its external
  // indices, so it is pinned and never folds; the Tensor constructor also
  // refuses to give it a braket symmetry in the first place
  Tensor A(reserved::antisymm_label(), bra{idx(L"i_1", Field::Real)},
           ket{idx(L"a_1", Field::Real)}, Symmetry::Antisymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  REQUIRE(braket_orientation_pinned(A));
  REQUIRE_FALSE(braket_foldable(A));
  REQUIRE_THROWS_AS(
      Tensor(reserved::antisymm_label(), bra{idx(L"i_1", Field::Real)},
             ket{idx(L"a_1", Field::Real)}, Symmetry::Antisymm,
             BraKetSymmetry::Symm, ColumnSymmetry::Symm),
      Exception);
}

TEST_CASE("normal_operator_swap_bra_ket", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  SECTION("the exchange swaps the bundles, each kept in particle order") {
    FNOperator op(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"}));
    static_cast<AbstractTensor&>(op)._swap_bra_ket();
    REQUIRE(op.ncreators() == 2);
    REQUIRE(op.nannihilators() == 2);
    const auto cre_labels = op.creators() |
                            ranges::views::transform([](auto const& o) {
                              return o.index().label();
                            }) |
                            ranges::to_vector;
    const auto ann_labels = op.annihilators() |
                            ranges::views::transform([](auto const& o) {
                              return o.index().label();
                            }) |
                            ranges::to_vector;
    REQUIRE(cre_labels == decltype(cre_labels){L"p_3", L"p_4"});
    REQUIRE(ann_labels == decltype(ann_labels){L"p_1", L"p_2"});
    // twice is the identity
    static_cast<AbstractTensor&>(op)._swap_bra_ket();
    REQUIRE(op == FNOperator(cre({L"p_1", L"p_2"}), ann({L"p_3", L"p_4"})));
  }

  SECTION("a Hermitian operator over a real basis canonicalizes") {
    // the graph verdict may exchange the operator's bundles; with the
    // exchange implemented that is a respelling of the same operator
    auto ridx = [](std::wstring_view l) { return idx(l, Field::Real); };
    auto make = [&](std::wstring_view c1, std::wstring_view c2,
                    std::wstring_view a1, std::wstring_view a2) {
      return ex<Tensor>(L"h", bra{ridx(L"p_5")}, ket{ridx(L"p_6")}) *
             ex<FNOperator>(cre({ridx(c1), ridx(c2)}),
                            ann({ridx(a1), ridx(a2)}));
    };
    auto e1 = make(L"p_1", L"p_2", L"p_1", L"p_2");
    auto e2 = make(L"p_2", L"p_1", L"p_2", L"p_1");  // same operator
    REQUIRE_NOTHROW(canonicalize(e1));
    REQUIRE_NOTHROW(canonicalize(e2));
    REQUIRE(*e1 == *e2);
    auto e3 = e1->clone();
    canonicalize(e3);
    REQUIRE(*e3 == *e1);
  }
}

TEST_CASE("conjugate_braket_fold_per_tensor", "[conjugation]") {
  // Per-tensor canonicalization exchanges the two bra/ket bundles only where
  // the exchange is a respelling of one value: Symm swaps freely, Antisymm
  // swaps at -1. A Conjugate/AntiConjugate tensor's two orientations are two
  // values, so both are left exactly as written and neither picks up a state.
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // occ (i) vs virt (a) bundles: space-decidable, yet Conjugate does not fold
  Tensor A(L"h", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor B(L"h", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(A) == 1);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(B) == 1);
  REQUIRE(A.bra()[0].label() == L"a_1");
  REQUIRE(A.ket()[0].label() == L"i_1");
  REQUIRE(B.bra()[0].label() == L"i_1");
  REQUIRE(B.ket()[0].label() == L"a_1");
  for (const Tensor& t : {A, B}) {
    REQUIRE_FALSE(t.adjointed());
    REQUIRE_FALSE(t.kconjugated());
  }

  // AntiConjugate: likewise untouched, and no sign is consumed
  Tensor D1(L"d", bra{L"a_1"}, ket{L"i_1"},
            TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                             .column = ColumnSymmetry::Symm});
  Tensor D2(L"d", bra{L"i_1"}, ket{L"a_1"},
            TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                             .column = ColumnSymmetry::Symm});
  REQUIRE(D1.braket_symmetry() == BraKetSymmetry::AntiConjugate);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(D1) == 1);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(D2) == 1);
  REQUIRE(D1.bra()[0].label() == L"a_1");
  REQUIRE(D2.bra()[0].label() == L"i_1");

  // Symm over a real basis: the two orientations land on one spelling, free
  auto real = [](std::wstring_view l) { return idx(l, Field::Real); };
  Tensor S1(L"s", bra{real(L"a_1")}, ket{real(L"i_1")}, Symmetry::Nonsymm,
            BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  Tensor S2(L"s", bra{real(L"i_1")}, ket{real(L"a_1")}, Symmetry::Nonsymm,
            BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(S1) == 1);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(S2) == 1);
  REQUIRE(S1.bra()[0].label() == S2.bra()[0].label());
  REQUIRE(S1.ket()[0].label() == S2.ket()[0].label());
  REQUIRE_FALSE(S1.kconjugated());
  REQUIRE_FALSE(S2.kconjugated());

  // Antisymm over a real basis: one spelling as well, the swapped one at -1
  auto anti = [&real](std::wstring_view b, std::wstring_view k) {
    return Tensor(L"n", bra{real(b)}, ket{real(k)},
                  TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                                   .column = ColumnSymmetry::Symm});
  };
  Tensor N1 = anti(L"a_1", L"i_1");
  Tensor N2 = anti(L"i_1", L"a_1");
  REQUIRE(N1.braket_symmetry() == BraKetSymmetry::Antisymm);
  const auto n1 = DefaultTensorCanonicalizer::canonicalize_braket(N1);
  const auto n2 = DefaultTensorCanonicalizer::canonicalize_braket(N2);
  REQUIRE(N1.bra()[0].label() == N2.bra()[0].label());
  REQUIRE(n1 * n2 == -1);  // exactly one was swapped, at -1

  // with fold_signed off the signed exchange is declined, so the two
  // orientations stay apart
  Tensor N3 = anti(L"a_1", L"i_1");
  Tensor N4 = anti(L"i_1", L"a_1");
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(
              N3, /*fold_signed=*/false) == 1);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(
              N4, /*fold_signed=*/false) == 1);
  REQUIRE(N3.bra()[0].label() == L"a_1");
  REQUIRE(N4.bra()[0].label() == L"i_1");

  // identical bundles (diagonal): never swapped
  Tensor E(L"s", bra{real(L"p_1")}, ket{real(L"p_1")}, Symmetry::Nonsymm,
           BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  REQUIRE(DefaultTensorCanonicalizer::canonicalize_braket(E) == 1);
  REQUIRE(E.bra()[0].label() == L"p_1");
}

TEST_CASE("with_slots_carries_attributes", "[conjugation]") {
  // Tensor::with_slots rebuilds the slots and carries label, symmetries,
  // hermiticity, and the states -- the sanctioned rebuild API for transforms
  // (rebuilding through a plain ctor drops the states)
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  Tensor t(L"g", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
           TensorSymmetries{.perm = Symmetry::Antisymm,
                            .hermiticity = Hermiticity::Hermitian,
                            .conjugation_parity = ConjugationParity::None});
  REQUIRE(t.kconjugate() == 1);
  using ixvec = container::svector<Index>;
  auto r = t.with_slots(bra<ixvec>{ixvec{Index{L"i_3"}, Index{L"i_4"}}},
                        ket<ixvec>{ixvec{Index{L"a_3"}, Index{L"a_4"}}},
                        aux<ixvec>{});
  REQUIRE(r.label() == t.label());
  REQUIRE(r.symmetry() == t.symmetry());
  REQUIRE(r.braket_symmetry() == t.braket_symmetry());
  REQUIRE(r.hermiticity() == t.hermiticity());
  REQUIRE(r.column_symmetry() == t.column_symmetry());
  REQUIRE(r.kconjugated());
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
  // a manifestly real (self-conjugate) summand: BraKetSymmetry::Symm tensors,
  // derivable only over a real basis
  auto self_adj = [] {
    auto real_basis = tests::scoped_real_basis();
    return deserialize(L"1/4 f{i_1;a_1}:N-S-S u{a_1;i_1}:N-S-S");
  }();

  SECTION("sum pair emits 2 Re(A)") {
    auto folded = fold_conjugate_pairs(term->clone() + term_adj->clone());
    auto expected = ex<Constant>(2) * real_part(term->clone());
    // the fold emits canonical-representative inners (so wrappers from
    // different passes merge and cancel); compare canonically
    auto cf = folded->clone();
    canonicalize(cf, CanonicalizeOptions::default_options());
    auto ce = expected->clone();
    canonicalize(ce, CanonicalizeOptions::default_options());
    REQUIRE(*cf == *ce);
  }

  SECTION("difference pair emits 2i Im(A)") {
    auto folded = fold_conjugate_pairs(term->clone() +
                                       ex<Constant>(-1) * term_adj->clone());
    auto expected =
        ex<Constant>(Complex<rational>(0, 2)) * imaginary_part(term->clone());
    auto cf = folded->clone();
    canonicalize(cf, CanonicalizeOptions::default_options());
    auto ce = expected->clone();
    canonicalize(ce, CanonicalizeOptions::default_options());
    REQUIRE(*cf == *ce);
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
    SEQUANT_PRAGMA_CLANG(diagnostic push)
    SEQUANT_PRAGMA_CLANG(diagnostic ignored "-Wdeprecated-declarations")
    SEQUANT_PRAGMA_GCC(diagnostic push)
    SEQUANT_PRAGMA_GCC(diagnostic ignored "-Wdeprecated-declarations")
    auto folded =
        fold_conjugate_pairs_of_real_sum(term->clone() + term_adj->clone());
    SEQUANT_PRAGMA_GCC(diagnostic pop)
    SEQUANT_PRAGMA_CLANG(diagnostic pop)
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
      // the fold emits a canonical-representative inner; the up spelling
      // and its (value-equal, per swap_spin) down partner are both valid
      auto expected_up = ex<Constant>(2) * real_part(term_up->clone());
      auto expected_dn =
          ex<Constant>(2) * real_part(mbpt::swap_spin(term_up->clone()));
      auto cf = folded->clone();
      canonicalize(cf, CanonicalizeOptions::default_options());
      auto cu = expected_up->clone();
      canonicalize(cu, CanonicalizeOptions::default_options());
      auto cd = expected_dn->clone();
      canonicalize(cd, CanonicalizeOptions::default_options());
      REQUIRE((*cf == *cu || *cf == *cd));
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

    // the default is Yes (evaluation ingests RealPart/ImagPart nodes);
    // opting out preserves the unfolded sum
    auto unfolded = sum->clone();
    simplify(unfolded, SimplifyOptions::default_options().copy_and_set(
                           SimplifyOptions::FoldConjugatePairs::No));
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

  // closed |C|^2 network: C conj(C), fully contracted -> real. The conjugate
  // of a value is the adjoint, so the second factor is the ⁺ spelling
  auto C = deserialize(L"C{a_1;i_1}:N-N-S");
  REQUIRE(is_hermitian_network(C->clone() * conjugate(C->clone())));
  REQUIRE(*conjugate(C->clone()) == *deserialize(L"C⁺{i_1;a_1}:N-N-S"));
  // closed C*C without the conjugation: a generically complex scalar
  REQUIRE_FALSE(
      is_hermitian_network(deserialize(L"C{a_1;i_1}:N-N-S C{a_1;i_1}:N-N-S")));
  // a declared-Hermitian energy-like scalar: h t + t* h* is self-adjoint
  auto term = deserialize(L"h{i_1;a_1}:N-C-S t{a_1;i_1}:N-C-S");
  auto sum = term->clone() + conjugate(term);
  REQUIRE(is_hermitian_network(sum));
}

TEST_CASE("swap_bra_ket_carries_states_and_aux", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // _swap_bra_ket exchanges the slot bundles only: both states, and the aux
  // slots, ride along
  auto t = ex<Tensor>(
      L"C", bra{L"a_1"}, ket{L"i_1"}, aux{L"p_5"},
      TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  REQUIRE(t->as<Tensor>().set_states(true, true) == 1);
  auto sw = mbpt::swap_bra_ket(t);
  auto const& st = sw->as<Tensor>();
  REQUIRE(st.adjointed());
  REQUIRE(st.kconjugated());
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

TEST_CASE("core_states", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  SECTION("NonHermitian, parity None: both states are distinct atoms") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    REQUIRE(ta.adjointed());
    REQUIRE_FALSE(ta.kconjugated());
    REQUIRE(ta.bra()[0].label() == L"i_1");  // slots exchanged
    REQUIRE(ta.decorated_label() == L"t⁺");
    Tensor tk = t;
    REQUIRE(tk.kconjugate() == 1);
    REQUIRE(tk.kconjugated());
    REQUIRE(tk.bra()[0].label() == L"a_1");  // slots in place
    REQUIRE(tk.decorated_label() == L"t꙳");
    Tensor tak = ta;
    REQUIRE(tak.kconjugate() == 1);
    REQUIRE(tak.decorated_label() == L"t⁺꙳");
    // involutions, and the two commute
    REQUIRE(tak.adjoint() == 1);
    REQUIRE(tak.decorated_label() == L"t꙳");
    REQUIRE(tak.kconjugate() == 1);
    REQUIRE(tak == t);
    // four distinct values: label, then states, order t < t⁺ < t꙳ < t⁺꙳
    REQUIRE(t < ta);
    REQUIRE(ta < tk);
    REQUIRE(t.hash_value() != ta.hash_value());
    REQUIRE(ta.hash_value() != tk.hash_value());
  }

  SECTION("Hermitian: the adjoint normalizes away, at +1") {
    Tensor h(L"h", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    REQUIRE(h.adjoint() == 1);
    REQUIRE_FALSE(h.adjointed());
    REQUIRE(h.bra()[0].label() == L"a_1");
  }

  SECTION("AntiHermitian: the adjoint normalizes away, at -1") {
    Tensor d(L"d", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(d.adjoint() == -1);
    REQUIRE_FALSE(d.adjointed());
    // a label mark that would carry -1 is refused
    REQUIRE_THROWS_AS(
        Tensor(L"d⁺", bra{L"i_1"}, ket{L"a_1"},
               TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian}),
        Exception);
  }

  SECTION("parity Even/Odd: the K state normalizes away") {
    Tensor e(L"e", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::Even});
    REQUIRE(e.kconjugate() == 1);
    REQUIRE_FALSE(e.kconjugated());
    Tensor o(L"o", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::Odd});
    REQUIRE(o.kconjugate() == -1);
    REQUIRE_FALSE(o.kconjugated());
    REQUIRE_THROWS_AS(
        Tensor(L"o꙳", bra{L"i_1"}, ket{L"a_1"},
               TensorSymmetries{.conjugation_parity = ConjugationParity::Odd}),
        Exception);
  }

  SECTION("real basis, both traits indefinite: the adjoint is spelled as K") {
    Tensor t(L"t", bra{idx(L"a_1", Field::Real)}, ket{idx(L"i_1", Field::Real)},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor ta = t;
    REQUIRE(ta.adjoint() == 1);
    REQUIRE_FALSE(ta.adjointed());
    REQUIRE(ta.kconjugated());
    REQUIRE(ta.bra()[0].label() == L"a_1");  // t⁺{i;a} = t꙳{a;i}: in place
    REQUIRE(ta.decorated_label() == L"t꙳");
  }

  SECTION("bra/ket-less: the adjoint is spelled as K") {
    Tensor w(L"w", bra{}, ket{}, aux{L"p_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor wa = w;
    REQUIRE(wa.adjoint() == 1);
    REQUIRE_FALSE(wa.adjointed());
    REQUIRE(wa.kconjugated());
    // declared real (Hermitian): the conjugate is itself
    Tensor r(L"r", bra{}, ket{}, aux{L"p_1"},
             TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    REQUIRE(r.kconjugate() == 1);
    REQUIRE_FALSE(r.kconjugated());
  }

  SECTION("marks in the label are adopted, in either order, at most once") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor a(L"t⁺꙳", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor b(L"t꙳⁺", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(a == b);
    REQUIRE(a.label() == L"t");
    REQUIRE(a.adjointed());
    REQUIRE(a.kconjugated());
    REQUIRE_THROWS_AS(Tensor(L"t⁺⁺", bra{L"i_1"}, ket{L"a_1"}), Exception);
    REQUIRE_THROWS_AS(Tensor(L"t꙳꙳", bra{L"i_1"}, ket{L"a_1"}), Exception);
  }

  SECTION("set_states does not touch the slots and normalizes") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.set_states(true, true) == 1);
    REQUIRE(t.bra()[0].label() == L"a_1");
    REQUIRE(t.decorated_label() == L"t⁺꙳");
    Tensor h(L"h", bra{L"i_1"}, ket{L"a_1"},
             TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(h.set_states(true, false) == -1);
    REQUIRE_FALSE(h.adjointed());
  }

  SECTION("with_slots carries both states and re-normalizes") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.set_states(true, true) == 1);
    using ixvec = container::svector<Index>;
    auto r = t.with_slots(bra<ixvec>{ixvec{Index{L"a_2"}}},
                          ket<ixvec>{ixvec{Index{L"i_2"}}}, aux<ixvec>{});
    REQUIRE(r.adjointed());
    REQUIRE(r.kconjugated());
  }

  SECTION("real basis, a definite hermiticity: the K state reduces") {
    const Index p1 = idx(L"p_1", Field::Real);
    const Index p2 = idx(L"p_2", Field::Real);
    auto herm = [&](const Index& b, const Index& k) {
      return ex<Tensor>(
          L"h", bra{b}, ket{k},
          TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                           .conjugation_parity = ConjugationParity::None});
    };
    // h꙳{p_1;p_2} = conj h{p_1;p_2} = h⁺{p_2;p_1} = h{p_2;p_1}
    auto hk = kconjugate(herm(p1, p2));
    REQUIRE(hk->is<Tensor>());
    REQUIRE_FALSE(hk->as<Tensor>().kconjugated());
    REQUIRE(*hk == *herm(p2, p1));
    REQUIRE(*hk == *conjugate(herm(p1, p2)));
    auto diff = herm(p2, p1) - kconjugate(herm(p1, p2));
    simplify(diff);
    REQUIRE(diff->is<Constant>());
    REQUIRE(diff->as<Constant>().value() == 0);

    // anti-Hermitian: d꙳{p_1;p_2} = -d{p_2;p_1}, the sign onto the Product
    auto anti = [&](const Index& b, const Index& k) {
      return ex<Tensor>(
          L"d", bra{b}, ket{k},
          TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian,
                           .conjugation_parity = ConjugationParity::None});
    };
    auto dk = kconjugate(anti(p1, p2));
    REQUIRE(dk->is<Product>());
    REQUIRE(dk->as<Product>().scalar() == -1);
    REQUIRE(dk->as<Product>().factors().size() == 1);
    REQUIRE(*dk->as<Product>().factors()[0] == *anti(p2, p1));
    auto asum = anti(p2, p1) + kconjugate(anti(p1, p2));
    simplify(asum);
    REQUIRE(asum->is<Constant>());
    REQUIRE(asum->as<Constant>().value() == 0);

    // an indefinite hermiticity keeps the mark, slots as written
    Tensor t(L"t", bra{p1}, ket{p2},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.kconjugate() == 1);
    REQUIRE(t.kconjugated());
    REQUIRE(t.bra()[0].label() == L"p_1");

    // over a complex basis the conjugation is not the adjoint: mark kept
    Tensor hc(L"h", bra{L"p_1"}, ket{L"p_2"},
              TensorSymmetries{.hermiticity = Hermiticity::Hermitian,
                               .conjugation_parity = ConjugationParity::None});
    REQUIRE(hc.kconjugate() == 1);
    REQUIRE(hc.kconjugated());
    REQUIRE(hc.bra()[0].label() == L"p_1");
  }
}

TEST_CASE("mark_serialization_roundtrip", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);
  for (auto spelling : {L"t⁺{i_1;a_1}:N-N-N-N", L"t꙳{a_1;i_1}:N-N-N-N",
                        L"t⁺꙳{i_1;a_1}:N-N-N-N"}) {
    auto e = deserialize(spelling);
    REQUIRE(e->is<Tensor>());
    REQUIRE(e->as<Tensor>().label() == L"t");
    REQUIRE(serialize(e, {.annot_symm = true}) == spelling);
  }
  // either mark order parses to one tensor
  REQUIRE(*deserialize(L"t꙳⁺{i_1;a_1}:N-N-N-N") ==
          *deserialize(L"t⁺꙳{i_1;a_1}:N-N-N-N"));
  // a mark carrying a sign becomes a scalar factor
  auto z = deserialize(L"z⁺{i_1;a_1}:N-A-N");
  REQUIRE(z->is<Product>());
  REQUIRE(z->as<Product>().scalar() == -1);
  auto o = deserialize(L"o꙳{i_1;a_1}:N-N-N-O");
  REQUIRE(o->is<Product>());
  REQUIRE(o->as<Product>().scalar() == -1);
  // variables
  auto v = deserialize(L"x꙳");
  REQUIRE(v->as<Variable>().conjugated());
  REQUIRE(serialize(v) == L"x꙳");
  // the caret spellings are errors
  using io::serialization::SerializationError;
  REQUIRE_THROWS_AS(deserialize(L"t^*{i_1;a_1}"), SerializationError);
  REQUIRE_THROWS_AS(deserialize(L"t^T{i_1;a_1}"), SerializationError);
  REQUIRE_THROWS_AS(deserialize(L"x^*"), SerializationError);
  // marks on an operator name are errors, not asserts
  REQUIRE_THROWS_AS(deserialize(L"ã꙳{a_1;i_1}"), SerializationError);
  REQUIRE_THROWS_AS(deserialize(L"a⁺{i_1;a_1}"), SerializationError);
  // a repeated mark names no state
  REQUIRE_THROWS_AS(deserialize(L"t⁺⁺{i_1;a_1}"), SerializationError);
  REQUIRE_THROWS_AS(deserialize(L"x꙳꙳"), SerializationError);
  // a variable has the one conjugated state, so an adjoint mark is an error
  REQUIRE_THROWS_AS(deserialize(L"x⁺"), SerializationError);
  // a mark on a reserved label constructs: δ is Hermitian by definition and
  // Even, so the state normalizes away and the defining symmetries stand
  auto k = deserialize(L"δ꙳{i_1;a_1}");
  REQUIRE(k->is<Tensor>());
  REQUIRE(k->as<Tensor>().label() == reserved::kronecker_label());
  REQUIRE_FALSE(k->as<Tensor>().kconjugated());
  REQUIRE(k->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(*k == *deserialize(L"δ{i_1;a_1}"));
  // LaTeX renders the states as superscripts, never as a raw mark
  REQUIRE(to_latex(deserialize(L"t⁺{i_1;a_1}:N-N-N-N")) ==
          L"{{t^{\\dagger}}^{{a_1}}_{{i_1}}}");
  REQUIRE(to_latex(deserialize(L"t꙳{a_1;i_1}:N-N-N-N")) ==
          L"{{t^{*}}^{{i_1}}_{{a_1}}}");
  REQUIRE(to_latex(deserialize(L"t⁺꙳{i_1;a_1}:N-N-N-N")) ==
          L"{{t^{\\dagger *}}^{{a_1}}_{{i_1}}}");
  REQUIRE(to_latex(v) == L"{{x}^{*}}");
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
  // whatever the observable braket symmetry the traits derive (Conjugate, over
  // the default complex field), the definite hermiticity is spelled with its
  // trait letter; the parity letter survives
  REQUIRE(p->as<Tensor>().braket_symmetry() == BraKetSymmetry::Conjugate);
  REQUIRE(serialize(p, {.annot_symm = true}) == L"p{i_1;i_2}:N-H-N-O");
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
  // parity from the bra/ket letter. The trait letters back-fill Even, so
  // `:N-A-N` needs nothing while an Odd parity under `H` spells itself out;
  // a pinned `C` over a real field stands for the traits Hermitian and None,
  // which is what the re-spelling names, the parity included
  REQUIRE(serialize(deserialize(L"t{i_1;i_2}:A-C-S"), {.annot_symm = true}) ==
          L"t{i_1;i_2}:A-H-S-N");
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
  REQUIRE(serialize(pinned, {.annot_symm = true}) == L"c{i_1;i_2}:N-H-S-N");
  REQUIRE(*deserialize(serialize(pinned, {.annot_symm = true})) == *pinned);
  REQUIRE(*deserialize(str) == *t);
}

TEST_CASE("hermiticity_trait_serialization", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);

  Index i1 = idx(L"i_1", Field::Real);
  Index i2 = idx(L"i_2", Field::Real);

  // over a real basis the exchange symmetry a definite hermiticity derives is
  // a plain (anti)symmetry, yet the letter spelled is the trait's, so the
  // spelling says what the tensor is rather than how its basis reads it
  auto h = ex<Tensor>(L"h", bra{i1}, ket{i2},
                      TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
  REQUIRE(h->as<Tensor>().braket_symmetry() == BraKetSymmetry::Symm);
  REQUIRE(serialize(h, {.annot_symm = true}) == L"h{i_1;i_2}:N-H-N");
  auto d =
      ex<Tensor>(L"d", bra{i1}, ket{i2},
                 TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
  REQUIRE(d->as<Tensor>().braket_symmetry() == BraKetSymmetry::Antisymm);
  REQUIRE(serialize(d, {.annot_symm = true}) == L"d{i_1;i_2}:N-A-N");

  // an exchange symmetry pinned at construction is a spelling of the traits
  // that derive it, so it too comes back as the trait letter
  auto pinned = ex<Tensor>(L"h", bra{i1}, ket{i2},
                           TensorSymmetries{.braket = BraKetSymmetry::Symm});
  REQUIRE(pinned->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(serialize(pinned, {.annot_symm = true}) == L"h{i_1;i_2}:N-H-N");
  REQUIRE(*pinned == *h);

  // the trait letter is field-agnostic, so the spelling reads back over a
  // basis of either field: over the real one it reproduces the tensor, over
  // the complex one the same traits derive the exchange symmetry that basis
  // states (Conjugate). The exchange symmetry itself has no such spelling:
  // `S` names a relation a complex basis does not have.
  {
    auto real_basis = tests::scoped_real_basis();
    REQUIRE(*deserialize(L"h{i_1;i_2}:N-H-N") == *h);
    REQUIRE(*deserialize(L"d{i_1;i_2}:N-A-N") == *d);
  }
  auto over_complex = deserialize(L"h{i_1;i_2}:N-H-N");
  REQUIRE(over_complex->as<Tensor>().hermiticity() == Hermiticity::Hermitian);
  REQUIRE(over_complex->as<Tensor>().braket_symmetry() ==
          BraKetSymmetry::Conjugate);
  using io::serialization::SerializationError;
  REQUIRE_THROWS_AS(deserialize(L"h{i_1;i_2}:N-S-N"), SerializationError);
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
    REQUIRE(t->as<Tensor>().kconjugate() == 1);
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
  // the conj-swapped spelling is a _different_ value and keeps its own slot
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
          REQUIRE_FALSE(node->as<Tensor>().kconjugated());
          REQUIRE(node->as<Tensor>().bra()[0].space() == Index(L"i_1").space());
        }
      },
      /*atoms_only=*/true);
  REQUIRE(found_A);
}

TEST_CASE("sum_merge_conjugate_marked_terms", "[conjugation]") {
  // identically-marked summands merge; a marked and an unmarked spelling of
  // _different_ values do not
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
  // _declared_ default comparator (default_idxptr_slottype_lesscompare) orders
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
  // proto-free indices precede proto-indexed ones (Nested outer;inner order)
  REQUIRE_FALSE(ci[0].has_proto_indices());
  REQUIRE_FALSE(ci[1].has_proto_indices());
  REQUIRE(ci[2].has_proto_indices());
}

TEST_CASE("canonicalize_marked_nonsymm_network", "[conjugation]") {
  // The K-conjugated state of a Nonsymm tensor is part of the graph
  // colouring in every canonicalization, so networks that differ only in
  // which factor carries the mark stay distinguishable, and canonicalization
  // is idempotent on them.
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // parity None (the fourth annotation letter): with the default parity Even
  // the ꙳ state would normalize away and there would be nothing to colour
  auto e1 = deserialize(L"t꙳{a_1;i_1}:N-N-N-N u{i_1;a_1}:N-N-N-N");
  auto e2 = deserialize(L"t{a_1;i_1}:N-N-N-N u꙳{i_1;a_1}:N-N-N-N");
  auto c1 = canonicalize(e1->clone());
  auto c2 = canonicalize(e2->clone());
  REQUIRE(*c1 != *c2);
  REQUIRE(*canonicalize(c1->clone()) == *c1);
  REQUIRE(*canonicalize(c2->clone()) == *c2);
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

  SECTION("explicit braket with a fixed trait: the parity is searched") {
    // Symm with AntiHermitian over a real basis is an imaginary symmetric
    // array: odd parity
    Tensor d(L"d", bra{idx(L"i_1", Field::Real)}, ket{idx(L"i_2", Field::Real)},
             TensorSymmetries{.braket = BraKetSymmetry::Symm,
                              .hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(d.conjugation_parity() == ConjugationParity::Odd);
    REQUIRE(d.braket_symmetry() == BraKetSymmetry::Symm);
    // no parity reproduces Conjugate from an anti-Hermitian over a real basis
    REQUIRE_THROWS_AS(
        Tensor(L"d", bra{idx(L"i_1", Field::Real)},
               ket{idx(L"i_2", Field::Real)},
               TensorSymmetries{.braket = BraKetSymmetry::Conjugate,
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

  SECTION("with_slots normalizes the states against the new field") {
    // complex-field NonHermitian: the adjoint is a state of its own
    Tensor t(L"t", bra{Index{L"i_1"}}, ket{Index{L"a_1"}});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    // rebuilt onto real-field slots the coset rule spells the adjoint
    // t⁺{a;i} as the K state with the slots exchanged, t꙳{i;a}, which the
    // even parity clears: t{i;a}, as a fresh construction would give
    using ixvec = container::svector<Index>;
    auto r =
        t.with_slots(bra<ixvec>{ixvec{idx(L"a_2", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"i_2", Field::Real)}}, aux<ixvec>{});
    REQUIRE_FALSE(r.adjointed());
    REQUIRE_FALSE(r.kconjugated());
    REQUIRE(r.bra()[0].label() == L"i_2");
    Tensor fresh(L"t", bra{idx(L"i_2", Field::Real)},
                 ket{idx(L"a_2", Field::Real)});
    REQUIRE(r == fresh);
    REQUIRE(r.hash_value() == fresh.hash_value());
  }

  SECTION("with_slots refuses a rebuild whose normalization costs a sign") {
    // complex-field NonHermitian of odd parity: the adjoint is a state. Over
    // a real field the coset rule spells it as the K state, which the odd
    // parity consumes at -1, which no Tensor can hold
    Tensor z(L"z", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
             TensorSymmetries{.conjugation_parity = ConjugationParity::Odd});
    REQUIRE(z.adjoint() == 1);
    REQUIRE(z.adjointed());
    using ixvec = container::svector<Index>;
    REQUIRE_THROWS_AS(
        z.with_slots(bra<ixvec>{ixvec{idx(L"a_2", Field::Real)}},
                     ket<ixvec>{ixvec{idx(L"i_2", Field::Real)}}, aux<ixvec>{}),
        Exception);
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

  // a fully-contracted summand and its adjoint, over a complex basis: the
  // adjoint reverses the factors and reorients each of them, so the
  // orientation of the (Anti)Conjugate factor is what tells the two apart
  auto d = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(
        L"d", bra{b}, ket{k},
        TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
  };
  auto g = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(L"g", bra{b}, ket{k},
                      TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
  };
  auto u = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(L"u", bra{b}, ket{k});
  };
  // the u of the adjoint summand: NonHermitian, so its adjoint keeps the ⁺
  auto u_adj = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(L"u⁺", bra{b}, ket{k});
  };
  auto factor_of = [](const ExprPtr& p, std::wstring_view label) {
    for (auto& f : p->as<Product>().factors())
      if (f->as<Tensor>().label() == label) return f->as<Tensor>();
    throw Exception("test: no factor with the requested label");
  };

  SECTION("anti-Hermitian: the two orientations are two values") {
    // e1 = d{i;a} u{a;i} and its adjoint e2 = d{a;i} u⁺{i;a}: since
    // d{a;i} = -conj(d{i;a}), e2 is -conj(e1), so the two canonical forms
    // are distinct and fold_conjugate_pairs emits 2i Im(e1)
    auto e1 = d(L"i_1", L"a_1") * u(L"a_1", L"i_1");
    auto e2 = d(L"a_1", L"i_1") * u_adj(L"i_1", L"a_1");
    // the adjoint reverses the factor order, so compare canonical forms
    REQUIRE(*canonicalize(adjoint(e1->clone())) ==
            *canonicalize(ex<Constant>(-1) * e2->clone()));
    auto c1 = canonicalize(e1->clone());
    auto c2 = canonicalize(e2->clone());
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    // canonical forms are idempotent
    REQUIRE(*canonicalize(c1->clone()) == *c1);
    REQUIRE(*canonicalize(c2->clone()) == *c2);
    // d keeps the orientation it was written in, and neither form acquires a
    // state: the two spellings are two values, not one respelt
    REQUIRE(factor_of(c1, L"d").bra()[0].label() !=
            factor_of(c2, L"d").bra()[0].label());
    for (auto& c : {c1, c2})
      for (auto& f : c->as<Product>().factors())
        REQUIRE_FALSE(f->as<Tensor>().kconjugated());
    REQUIRE(*c1 != *c2);
    // the pair is still a conjugate pair of the value: s + (-s*) = 2i Im(s)
    auto folded = fold_conjugate_pairs(e1->clone() + e2->clone());
    INFO(toUtf8(to_latex(folded)));
    REQUIRE_FALSE(folded->is<Sum>());
    REQUIRE(folded->is<Product>());
    REQUIRE(folded->as<Product>().scalar() == (C{0, 2}));
    REQUIRE(folded->as<Product>().factors().size() == 1);
    auto const& wrapped = folded->as<Product>().factors()[0];
    REQUIRE(wrapped->is<ImagPart>());
    // the wrapper carries one of the pair's two canonical forms
    auto const inner = canonicalize(wrapped->as<ImagPart>().inner()->clone());
    REQUIRE((*inner == *c1 || *inner == *c2));
  }

  SECTION("anti-Hermitian, Complete method: the lexicographic pass agrees") {
    // the test binary pins Topological; run the same check under the
    // library's default Complete so the lexicographic relabeling pass of
    // TensorNetworkV3::canonicalize is exercised on these tensors too
    const CanonicalizeOptions opts{.method = CanonicalizationMethod::Complete};
    auto c1 = canonicalize(d(L"i_1", L"a_1") * u(L"a_1", L"i_1"), opts);
    auto c2 = canonicalize(d(L"a_1", L"i_1") * u_adj(L"i_1", L"a_1"), opts);
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    REQUIRE(*canonicalize(c1->clone(), opts) == *c1);
    REQUIRE(*canonicalize(c2->clone(), opts) == *c2);
    REQUIRE(factor_of(c1, L"d").bra()[0].label() !=
            factor_of(c2, L"d").bra()[0].label());
    REQUIRE(*c1 != *c2);
  }

  SECTION("Hermitian, default column symmetry: two canonical forms") {
    // e1 = g{i;a} u{a;i} and its adjoint e2 = g{a;i} u⁺{i;a}: g{a;i} is
    // conj(g{i;a}), a different array, so canonicalization keeps the two
    // apart and fold_conjugate_pairs emits 2 Re(e1)
    REQUIRE(g(L"i_1", L"a_1")->as<Tensor>().column_symmetry() ==
            ColumnSymmetry::Nonsymm);
    auto e1 = g(L"i_1", L"a_1") * u(L"a_1", L"i_1");
    auto e2 = g(L"a_1", L"i_1") * u_adj(L"i_1", L"a_1");
    REQUIRE(*canonicalize(adjoint(e1->clone())) == *canonicalize(e2->clone()));
    auto c1 = canonicalize(e1->clone());
    auto c2 = canonicalize(e2->clone());
    REQUIRE(c1->is<Product>());
    REQUIRE(c2->is<Product>());
    REQUIRE(factor_of(c1, L"g").bra()[0].label() !=
            factor_of(c2, L"g").bra()[0].label());
    for (auto& c : {c1, c2})
      for (auto& f : c->as<Product>().factors())
        REQUIRE_FALSE(f->as<Tensor>().kconjugated());
    REQUIRE(*c1 != *c2);
    auto folded = fold_conjugate_pairs(e1->clone() + e2->clone());
    INFO(toUtf8(to_latex(folded)));
    REQUIRE_FALSE(folded->is<Sum>());
    REQUIRE(folded->is<Product>());
    REQUIRE(folded->as<Product>().scalar() == (C{2, 0}));
    REQUIRE(folded->as<Product>().factors().size() == 1);
    auto const& wrapped = folded->as<Product>().factors()[0];
    REQUIRE(wrapped->is<RealPart>());
    // the wrapper carries one of the pair's two canonical forms
    auto const inner = canonicalize(wrapped->as<RealPart>().inner()->clone());
    REQUIRE((*inner == *c1 || *inner == *c2));
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
      for (auto& f : c->as<Product>().factors()) {
        REQUIRE_FALSE(f->as<Tensor>().adjointed());
        REQUIRE_FALSE(f->as<Tensor>().kconjugated());
      }
    REQUIRE(*canonicalize(c1->clone()) == *c1);
  }
}

TEST_CASE("hermitian_cycle_relabeling_invariance", "[conjugation]") {
  // A Conjugate-braket tensor keeps the orientation it was written in, so a
  // canonical form must be decided by the network's shape alone. These
  // networks are Hermitian cycles: every relabeling and factor reordering of
  // one spells the same value, and a graph automorphism that exchanged a
  // tensor's bundles would let the labels leak into the verdict.
  auto ctx = get_default_context();
  ctx.set(mbpt::make_min_sr_spaces(mbpt::SpinConvention::None));
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto _ = set_scoped_default_context(ctx);
  auto canon = [](std::wstring s) {
    auto e = deserialize(s);
    canonicalize(e);
    return to_latex(e);
  };
  for (auto [X, Y] :
       std::vector<std::pair<std::wstring, std::wstring>>{{L"h", L"γ"},
                                                          {L"f", L"d"},
                                                          {L"X", L"Y"},
                                                          {L"u", L"v"},
                                                          {L"A1", L"B7"},
                                                          {L"g", L"Γ"}}) {
    std::vector<std::wstring> sp = {
        X + L"{p_1;p_2}:N-C-S * " + Y + L"{p_2;p_1}:N-C-S",
        X + L"{p_2;p_1}:N-C-S * " + Y + L"{p_1;p_2}:N-C-S",
        Y + L"{p_2;p_1}:N-C-S * " + X + L"{p_1;p_2}:N-C-S",
        X + L"{p_3;p_7}:N-C-S * " + Y + L"{p_7;p_3}:N-C-S",
        X + L"{p_7;p_3}:N-C-S * " + Y + L"{p_3;p_7}:N-C-S"};
    auto ref = canon(sp[0]);
    for (auto& s : sp) {
      INFO(toUtf8(s));
      CHECK(canon(s) == ref);
    }
    auto diff = deserialize(sp[0]) - deserialize(sp[1]);
    simplify(diff);
    INFO(toUtf8(to_latex(diff)));
    CHECK((diff->is<Constant>() && diff->as<Constant>().is_zero()));
    // 3-cycle, all Hermitian (reversal automorphism)
    auto s1 = X + L"{p_1;p_2}:N-C-S * " + Y +
              L"{p_2;p_3}:N-C-S * "
              L"Z{p_3;p_1}:N-C-S";
    auto s2 = X + L"{p_5;p_4}:N-C-S * " + Y +
              L"{p_4;p_9}:N-C-S * "
              L"Z{p_9;p_5}:N-C-S";
    auto s3 = L"Z{p_3;p_1}:N-C-S * " + Y + L"{p_2;p_3}:N-C-S * " + X +
              L"{p_1;p_2}:N-C-S";
    CHECK(canon(s1) == canon(s2));
    CHECK(canon(s1) == canon(s3));
  }
  // column-symmetric 4-index Hermitian blocks, both spellings
  CHECK(canon(L"g{p_1,p_2;p_2,p_1}:N-C-S") ==
        canon(L"g{p_2,p_1;p_1,p_2}:N-C-S"));
  CHECK(canon(L"g{p_1,p_2;p_3,p_4}:N-C-S * Γ{p_3,p_4;p_1,p_2}:N-C-S") ==
        canon(L"g{p_3,p_4;p_1,p_2}:N-C-S * Γ{p_1,p_2;p_3,p_4}:N-C-S"));
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

  SECTION("an exchange symmetry the traits do not derive is refused") {
    // Symm and Antisymm relate an array to itself under the bra<->ket
    // exchange; over a complex basis the exchange relates an operator's
    // matrix to its conjugate, so no trait combination derives them there.
    // Such an array is real, which is declared through the basis field.
    REQUIRE_THROWS_AS(Tensor(L"g", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
                             TensorSymmetries{.braket = BraKetSymmetry::Symm}),
                      Exception);
    REQUIRE_THROWS_AS(
        Tensor(L"g", bra{Index{L"i_1"}}, ket{Index{L"a_1"}},
               TensorSymmetries{.braket = BraKetSymmetry::Antisymm}),
        Exception);
    REQUIRE_THROWS_AS(deserialize(L"g{i_1;a_1}:N-S-N"),
                      io::serialization::SerializationError);
    // over a real basis both derive
    Tensor gs(L"g", bra{idx(L"i_1", Field::Real)},
              ket{idx(L"a_1", Field::Real)},
              TensorSymmetries{.braket = BraKetSymmetry::Symm});
    REQUIRE(gs.hermiticity() == Hermiticity::Hermitian);
    REQUIRE(gs.conjugation_parity() == ConjugationParity::Even);
    Tensor ga(L"g", bra{idx(L"i_1", Field::Real)},
              ket{idx(L"a_1", Field::Real)},
              TensorSymmetries{.braket = BraKetSymmetry::Antisymm});
    REQUIRE(ga.hermiticity() == Hermiticity::AntiHermitian);
    REQUIRE(ga.conjugation_parity() == ConjugationParity::Even);
    // the given traits are fixed and the missing one is searched: Antisymm
    // with Hermitian is an odd-parity (imaginary) array
    Tensor go(L"g", bra{idx(L"i_1", Field::Real)},
              ket{idx(L"a_1", Field::Real)},
              TensorSymmetries{.braket = BraKetSymmetry::Antisymm,
                               .hermiticity = Hermiticity::Hermitian});
    REQUIRE(go.conjugation_parity() == ConjugationParity::Odd);
    // and a pin the given traits contradict is refused
    REQUIRE_THROWS_AS(
        Tensor(L"g", bra{idx(L"i_1", Field::Real)},
               ket{idx(L"a_1", Field::Real)},
               TensorSymmetries{.braket = BraKetSymmetry::Conjugate,
                                .conjugation_parity = ConjugationParity::Odd}),
        Exception);
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

TEST_CASE("generic_adjoint_keeps_no_sign", "[conjugation]") {
  // sequant::adjoint(T) returns a T, which cannot hold the sign an
  // anti-Hermitian tensor's adjoint carries: it throws, and the ExprPtr
  // overload is the way to get that adjoint with its sign
  Tensor d(L"d", bra{L"i_1"}, ket{L"i_2"},
           TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
  REQUIRE_THROWS_AS(sequant::adjoint(d), Exception);
  auto viaExpr = sequant::adjoint(ex<Tensor>(d));
  REQUIRE(viaExpr->is<Product>());
  REQUIRE(viaExpr->as<Product>().scalar() == -1);

  // a Hermitian tensor's adjoint is itself, respelled
  Tensor h(L"h", bra{L"i_1"}, ket{L"i_2"},
           TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
  Tensor h_adj = sequant::adjoint(h);
  REQUIRE(h_adj.bra()[0].label() == L"i_2");
  REQUIRE_FALSE(h_adj.adjointed());
  REQUIRE_FALSE(h_adj.kconjugated());
  // and a Nonsymm one is the ⁺ state
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
  REQUIRE(sequant::adjoint(t).adjointed());
}

TEST_CASE("conjugate_is_the_value_conjugate", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  SECTION("a c-number tensor: the adjoint spelling") {
    auto t = ex<Tensor>(
        L"t", bra{L"a_1"}, ket{L"i_1"},
        TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    auto c = conjugate(t);
    REQUIRE(*c == *sequant::adjoint(t));
    REQUIRE(c->as<Tensor>().adjointed());
    REQUIRE(c->as<Tensor>().bra()[0].label() == L"i_1");
    // Hermitian: the other orientation, no mark
    auto h =
        ex<Tensor>(L"h", bra{L"i_1"}, ket{L"a_1"},
                   TensorSymmetries{.hermiticity = Hermiticity::Hermitian});
    REQUIRE(
        *conjugate(h) ==
        *ex<Tensor>(L"h", bra{L"a_1"}, ket{L"i_1"},
                    TensorSymmetries{.hermiticity = Hermiticity::Hermitian}));
    // anti-Hermitian: minus the other orientation
    auto d =
        ex<Tensor>(L"d", bra{L"i_1"}, ket{L"a_1"},
                   TensorSymmetries{.hermiticity = Hermiticity::AntiHermitian});
    REQUIRE(conjugate(d)->is<Product>());
    REQUIRE(conjugate(d)->as<Product>().scalar() == -1);
  }
  SECTION(
      "over a real basis the spelling is the star with the slots in place") {
    auto t = ex<Tensor>(
        L"t", bra{idx(L"a_1", Field::Real)}, ket{idx(L"i_1", Field::Real)},
        TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    auto c = conjugate(t);
    REQUIRE(c->as<Tensor>().kconjugated());
    REQUIRE_FALSE(c->as<Tensor>().adjointed());
    REQUIRE(c->as<Tensor>().bra()[0].label() == L"a_1");
    REQUIRE(*c == *kconjugate(t));
  }
  SECTION("a product: factorwise, the scalar conjugated") {
    auto e = ex<Constant>(Constant::scalar_type{0, 1}) *
             ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"}) * ex<Variable>(L"z");
    auto c = conjugate(e);
    REQUIRE(c->is<Product>());
    REQUIRE(c->as<Product>().scalar() == Constant::scalar_type{0, -1});
    // conjugate is the adjoint on c-number content, and Product::adjoint
    // reverses the factors, so each factor is located by its kind
    auto factor_of = [&c](auto pred) {
      for (const auto& f : c->as<Product>().factors())
        if (pred(f)) return f;
      return ExprPtr{};
    };
    auto tf = factor_of([](const ExprPtr& f) { return f->is<Tensor>(); });
    auto vf = factor_of([](const ExprPtr& f) { return f->is<Variable>(); });
    REQUIRE(tf);
    REQUIRE(vf);
    REQUIRE(tf->as<Tensor>().adjointed());
    REQUIRE(vf->as<Variable>().conjugated());
  }
  SECTION(
      "operators have no value: conjugate throws, kconjugate is the "
      "identity on the string") {
    auto e = ex<Tensor>(L"h", bra{idx(L"p_1", Field::Real)},
                        ket{idx(L"p_2", Field::Real)}) *
             ex<FNOperator>(cre({idx(L"p_1", Field::Real)}),
                            ann({idx(L"p_2", Field::Real)}));
    REQUIRE_THROWS_AS(conjugate(e), Exception);
    auto k = kconjugate(e);
    REQUIRE(k->as<Product>().factor(1)->is<FNOperator>());
    REQUIRE(*k->as<Product>().factor(1) == *e->as<Product>().factor(1));
    // over a complex basis the string has no K-closure: refused
    auto ec = ex<Tensor>(L"h", bra{L"p_1"}, ket{L"p_2"}) *
              ex<FNOperator>(cre({L"p_1"}), ann({L"p_2"}));
    REQUIRE_THROWS_AS(kconjugate(ec), Exception);
    // an index is K-closed only with its proto-index closure: a real-space
    // index over a complex-space proto index is refused as well
    const Index a1_over_real{idx(L"a_1", Field::Real).space(), 1,
                             IndexList{idx(L"i_1", Field::Real)}};
    const Index a1_over_complex{idx(L"a_1", Field::Real).space(), 1,
                                IndexList{Index(L"i_1")}};
    REQUIRE_NOTHROW(kconjugate(
        ex<FNOperator>(cre({a1_over_real}), ann({idx(L"i_1", Field::Real)}))));
    REQUIRE_THROWS_AS(
        kconjugate(ex<FNOperator>(cre({a1_over_complex}),
                                  ann({idx(L"i_1", Field::Real)}))),
        Exception);
  }
}

TEST_CASE("fold_conjugate_pairs_is_order_independent", "[conjugation]") {
  auto sr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  Context ctx = get_default_context();
  ctx.set(sr);
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);
  auto term = deserialize(L"1/2 h{i_1;a_1}:N-C-S t{a_1;i_1}:N-C-S");
  auto term_adj = deserialize(L"1/2 t{i_2;a_2}:N-C-S h{a_2;i_2}:N-C-S");
  auto f1 = fold_conjugate_pairs(term->clone() + term_adj->clone());
  auto f2 = fold_conjugate_pairs(term_adj->clone() + term->clone());
  simplify(f1);
  simplify(f2);
  REQUIRE(*f1 == *f2);
}

TEST_CASE("slot_mutation_renormalizes_states", "[conjugation]") {
  // A slot mutation can put a tensor on a basis of another field, where a set
  // state denotes a different array: over a real basis the coset rule trades a
  // `⁺` for a `꙳` with the bundles exchanged back, and the parity then keeps
  // or consumes that star. transform_indices(), set_bra() and set_ket()
  // therefore normalize the states again, and refuse a normalization that
  // costs a sign.
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto resetter = set_scoped_default_context(ctx);

  // `t⁺{i_1;a_1}`: NonHermitian over a complex basis, so the mark stands
  auto adjointed_t = [](ConjugationParity parity) {
    Tensor t(L"t", bra{Index{L"a_1"}}, ket{Index{L"i_1"}},
             TensorSymmetries{.conjugation_parity = parity});
    REQUIRE(t.base_field() == Field::Complex);
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    REQUIRE(t.bra()[0].label() == L"i_1");
    return t;
  };
  const container::map<Index, Index> to_real{
      {Index{L"a_1"}, idx(L"a_1", Field::Real)},
      {Index{L"i_1"}, idx(L"i_1", Field::Real)}};

  SECTION("onto a real basis: the coset spelling, parity None keeps the star") {
    auto t = adjointed_t(ConjugationParity::None);
    REQUIRE(t.transform_indices(to_real));
    t.reset_tags();  // transform_indices tags the replaced indices
    REQUIRE(t.base_field() == Field::Real);
    REQUIRE_FALSE(t.adjointed());
    REQUIRE(t.kconjugated());
    // the coset rule exchanged the bundles back, so the slots stand as written
    REQUIRE(t.bra()[0].label() == L"a_1");
    REQUIRE(t.ket()[0].label() == L"i_1");
    REQUIRE(t.decorated_label() == L"t꙳");
  }

  SECTION("onto a real basis: the default parity consumes the star") {
    auto t = adjointed_t(ConjugationParity::Even);
    REQUIRE(t.transform_indices(to_real));
    t.reset_tags();
    REQUIRE(t.base_field() == Field::Real);
    REQUIRE_FALSE(t.adjointed());
    REQUIRE_FALSE(t.kconjugated());
    REQUIRE(t.bra()[0].label() == L"a_1");
    REQUIRE(t.decorated_label() == L"t");
  }

  SECTION("an odd-parity array is refused: the mutation would cost a sign") {
    auto t = adjointed_t(ConjugationParity::Odd);
    const Tensor before = t;
    REQUIRE_THROWS_AS(t.transform_indices(to_real), Exception);
    // the refused mutation leaves the tensor as the call found it, down to the
    // slots' field, which Index equality does not see
    REQUIRE(t.base_field() == Field::Complex);
    REQUIRE(t.decorated_label() == L"t⁺");
    REQUIRE(t.adjointed());
    REQUIRE_FALSE(t.kconjugated());
    REQUIRE(t.bra()[0].label() == L"i_1");
    REQUIRE(t.braket_symmetry() == before.braket_symmetry());
    REQUIRE(t.conjugation_symmetry() == before.conjugation_symmetry());
    REQUIRE(t.hash_value() == before.hash_value());
  }

  SECTION("a relabeling within one field leaves the states as they are") {
    auto t = adjointed_t(ConjugationParity::None);
    const container::map<Index, Index> rename{{Index{L"i_1"}, Index{L"i_2"}}};
    REQUIRE(t.transform_indices(rename));
    t.reset_tags();
    REQUIRE(t.base_field() == Field::Complex);
    REQUIRE(t.adjointed());
    REQUIRE_FALSE(t.kconjugated());
    REQUIRE(t.bra()[0].label() == L"i_2");
  }

  SECTION("set_bra and set_ket normalize the same way") {
    auto t = adjointed_t(ConjugationParity::None);
    t.set_bra({idx(L"i_1", Field::Real)});
    // the bra alone is real, so the base field is still Complex: nothing moved
    REQUIRE(t.base_field() == Field::Complex);
    REQUIRE(t.adjointed());
    t.set_ket({idx(L"a_1", Field::Real)});
    REQUIRE(t.base_field() == Field::Real);
    REQUIRE_FALSE(t.adjointed());
    REQUIRE(t.kconjugated());
    // the coset rule exchanged the bundles, so the real a_1 is now the bra
    REQUIRE(t.bra()[0].label() == L"a_1");
    REQUIRE(t.ket()[0].label() == L"i_1");

    auto odd = adjointed_t(ConjugationParity::Odd);
    odd.set_bra({idx(L"i_1", Field::Real)});
    REQUIRE(odd.base_field() == Field::Complex);
    REQUIRE_THROWS_AS(odd.set_ket({idx(L"a_1", Field::Real)}), Exception);
    // the refused setter puts the complex ket back
    REQUIRE(odd.base_field() == Field::Complex);
    REQUIRE(odd.decorated_label() == L"t⁺");
    REQUIRE(odd.adjointed());
    REQUIRE(odd.bra()[0].label() == L"i_1");
    REQUIRE(odd.ket()[0].space().field() == Field::Complex);
  }

  SECTION("a moved tensor is the tensor a fresh construction gives") {
    // the exchange symmetry is derived from the traits and the slots' field, so
    // a mutation that moves the tensor onto another basis derives it again: a
    // Hermitian, even-parity array is Conjugate over a complex basis and Symm
    // over a real one
    const TensorSymmetries hermitian{.hermiticity = Hermiticity::Hermitian};
    Tensor moved(L"g", bra{Index{L"i_1"}}, ket{Index{L"a_1"}}, hermitian);
    REQUIRE(moved.braket_symmetry() == BraKetSymmetry::Conjugate);
    const container::map<Index, Index> g_to_real{
        {Index{L"i_1"}, idx(L"i_1", Field::Real)},
        {Index{L"a_1"}, idx(L"a_1", Field::Real)}};
    REQUIRE(moved.transform_indices(g_to_real));
    moved.reset_tags();

    const Tensor fresh(L"g", bra{idx(L"i_1", Field::Real)},
                       ket{idx(L"a_1", Field::Real)}, hermitian);
    REQUIRE(moved.base_field() == Field::Real);
    REQUIRE(moved.braket_symmetry() == fresh.braket_symmetry());
    REQUIRE(moved.braket_symmetry() == BraKetSymmetry::Symm);
    REQUIRE(moved.conjugation_symmetry() == fresh.conjugation_symmetry());
    REQUIRE(moved.hermiticity() == fresh.hermiticity());
    REQUIRE(moved.conjugation_parity() == fresh.conjugation_parity());
    REQUIRE(moved.hash_value() == fresh.hash_value());
    REQUIRE(moved == fresh);
  }
}

TEST_CASE("has_tensor_sees_through_re_im", "[conjugation]") {
  // the conjugate-pair fold (default-on in a complex field) wraps folded
  // summands in RealPart/ImagPart; tensor queries must see through them
  using namespace sequant;
  auto sr = mbpt::make_min_sr_spaces();
  Context ctx = get_default_context();
  ctx.set(sr);
  auto resetter = set_scoped_default_context(ctx);
  auto t = deserialize(L"t{a_1;i_1}:N-C-S");
  auto f = deserialize(L"f{i_1;a_1}:N-C-S");
  REQUIRE(f->as<Tensor>().kconjugate() == 1);  // f꙳
  auto summand = ex<Constant>(2) * real_part(t->clone() * f->clone());
  INFO("summand = " << toUtf8(to_latex(summand)));
  REQUIRE(has_tensor(summand, L"f"));
  REQUIRE_FALSE(has_tensor(summand, L"g"));
  auto sum = summand->clone() + deserialize(L"g{i_1;a_1}:N-C-S");
  REQUIRE(has_tensor(sum, L"f"));
  REQUIRE(has_tensor(sum, L"g"));
  auto im = ex<Constant>(Constant::scalar_type(0, 2)) *
            imaginary_part(t->clone() * f->clone());
  REQUIRE(has_tensor(im, L"f"));
}
