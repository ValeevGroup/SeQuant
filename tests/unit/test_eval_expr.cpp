#include <SeQuant/core/expressions/complex.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"
#include "tensor_power_fixtures.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/eval_node_compare.hpp>
#include <SeQuant/core/eval/kramers_blind.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/equal.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

namespace sequant {
Tensor parse_tensor(
    std::wstring_view tnsr,
    const io::serialization::DeserializationOptions& options = {}) {
  return deserialize(tnsr, options)->as<Tensor>();
}

Constant parse_constant(std::wstring_view c) {
  return deserialize(c)->as<Constant>();
}

EvalExpr result_expr(EvalExpr const& left, EvalExpr const& right, EvalOp op) {
  SEQUANT_ASSERT(op == EvalOp::Product || op == EvalOp::Sum);
  auto xpr = op == EvalOp::Product ? left.expr() * right.expr()
                                   : left.expr() + right.expr();
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  return *binarize(xpr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

}  // namespace sequant

namespace {
/// an Index whose space carries the given field (the default field is Complex)
sequant::Index idx(std::wstring_view label, sequant::Field field) {
  sequant::Index i(label);
  sequant::IndexSpace sp = i.space();
  sp.field(field);
  return sequant::Index(label, sp);
}
}  // namespace

// a summand with an empty slot: the Sum's placeholder keeps the slot (and
// hence the summand's layout) rather than dropping it with the null index
TEST_CASE("eval_expr_sum_keeps_empty_slots", "[EvalExpr]") {
  using namespace sequant;
  auto make = [](std::wstring_view label) {
    return ex<Tensor>(label, bra{Index{L"a_1"}, Index{L"a_2"}},
                      ket{Index{L"i_1"}, Index{}}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = binarize(make(L"X") + make(L"Y"));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(node->is_sum());
  auto const& t = node->as_tensor();
  REQUIRE(t.bra().size() == 2);
  REQUIRE(t.ket().size() == 2);
  REQUIRE(t.ket_net_rank() == 1);
}

// Re/Im are evaluated for scalar results only, so a wrapper over a
// tensor-valued expression is refused where the eval node is built
TEST_CASE("eval_expr_re_im_tensor_valued", "[EvalExpr][re-im]") {
  using namespace sequant;
  auto t = ex<Tensor>(L"t", bra{L"a_1"}, ket{L"i_1"});
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  REQUIRE_THROWS_AS(binarize(real_part(t->clone())), Exception);
  REQUIRE_THROWS_AS(binarize(imaginary_part(t->clone())), Exception);
  // a scalar-valued one is fine
  auto s = ex<Tensor>(L"g", bra{L"i_1"}, ket{L"a_1"}) * t->clone();
  REQUIRE_NOTHROW(binarize(real_part(s)));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
}

TEST_CASE("eval_expr", "[EvalExpr]") {
  using namespace std::string_literals;
  using sequant::EvalExpr;
  using namespace sequant;

  SECTION("Constructors") {
    auto t1 = parse_tensor(L"t_{i1, i2}^{a1, a2}");
    REQUIRE_NOTHROW(EvalExpr{t1});
    REQUIRE(!EvalExpr{t1}.op_type());

    auto p1 = deserialize(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");
    REQUIRE(!EvalExpr{p1->at(0)->as<Tensor>()}.op_type());
    REQUIRE(!EvalExpr{p1->at(1)->as<Tensor>()}.op_type());

    REQUIRE(!EvalExpr{Variable{L"λ"}}.op_type());

    REQUIRE(!EvalExpr{Constant{1}}.op_type());
    REQUIRE(!EvalExpr{deserialize(L"1/2")->as<Constant>()}.op_type());

    REQUIRE(!EvalExpr{Power(ex<Constant>(rational{1, 2}), rational{1, 2})}
                 .op_type());
    REQUIRE(!EvalExpr{Power(ex<Variable>(L"x"), rational{3, 1})}.op_type());
  }

  SECTION("ResultType types") {
    auto T = [](std::wstring_view xpr) { return EvalExpr{parse_tensor(xpr)}; };

    auto C = [](std::wstring_view xpr) {
      return EvalExpr{parse_constant(xpr)};
    };

    auto result_type = [](EvalExpr const& left,   //
                          EvalExpr const& right,  //
                          EvalOp op) -> ResultType {
      return result_expr(left, right, op).result_type();
    };

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                T(L"Y{i1;a1}"),  //
                EvalOp::Sum      //
                ) == ResultType::Tensor);

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                T(L"Y{a1;i1}"),  //
                EvalOp::Product  //
                ) == ResultType::Scalar);

    REQUIRE(result_type(                //
                T(L"X{i1,i2; a3,a4}"),  //
                T(L"Y{a3,a4; a1,a2}"),  //
                EvalOp::Product         //
                ) == ResultType::Tensor);

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                C(L"2.5"),       //
                EvalOp::Product  //
                ) == ResultType::Tensor);

    REQUIRE(result_type(         //
                C(L"1.5"),       //
                C(L"2.5"),       //
                EvalOp::Product  //
                ) == ResultType::Scalar);

    REQUIRE(result_type(     //
                C(L"1.5"),   //
                C(L"2.5"),   //
                EvalOp::Sum  //
                ) == ResultType::Scalar);
  }

  SECTION("result expr") {
    ExprPtr expr = deserialize(L"2 var");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    ExprPtr root_expr = binarize(expr)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr->is<Variable>());
    REQUIRE(*root_expr != *expr);

    expr = deserialize(L"2 t{a1;i1}");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(expr)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr->is<Tensor>());
    REQUIRE(*root_expr != *expr);

    // The binarized tree shall respect the label of the ResultExpr
    ResultExpr res =
        deserialize<ResultExpr>(L"E = g{i1,i2;a1,a2} t{a1,a2;i1,i2}");
    root_expr = binarize(res)->expr();
    REQUIRE(root_expr.is<Variable>());
    REQUIRE(root_expr.as<Variable>().label() == L"E");

    // The binarized tree shall respect the indexing of the ResultExpr
    // (the result's `S` braket letter is derivable only over a real basis)
    auto real_basis = sequant::tests::scoped_real_basis();
    res = deserialize<ResultExpr>(
        L"Result{a2;i2}:A-S-S = g{i1,i2;a1,a2} t{a1;i1}");
    root_expr = binarize(res)->expr();
    REQUIRE(root_expr.is<Tensor>());
    REQUIRE(root_expr.as<Tensor>() ==
            Tensor(L"Result", bra(IndexList{L"a_2"}), ket(IndexList{L"i_2"}),
                   Symmetry::Antisymm, BraKetSymmetry::Symm,
                   ColumnSymmetry::Symm));

    // continued ->  check that changing indexing in result changes indexing in
    // tree
    res = deserialize<ResultExpr>(
        L"Result{i2;a2}:A-S-S = g{i1,i2;a1,a2} t{a1;i1}");
    root_expr = binarize(res)->expr();
    REQUIRE(root_expr.is<Tensor>());
    REQUIRE(root_expr.as<Tensor>() ==
            Tensor(L"Result", bra(IndexList{L"i_2"}), ket(IndexList{L"a_2"}),
                   Symmetry::Antisymm, BraKetSymmetry::Symm,
                   ColumnSymmetry::Symm));

    // The name-respecting property shall also hold for terminals
    res = deserialize<ResultExpr>(L"Other = Var");
    root_expr = binarize(res)->expr();
    REQUIRE(root_expr.is<Variable>());
    REQUIRE(root_expr.as<Variable>().label() == L"Other");

    res = deserialize<ResultExpr>(L"Amplitude{i1;a1} = t{a1;i1}");
    root_expr = binarize(res)->expr();
    REQUIRE(root_expr.is<Tensor>());
    // the deserialized ResultExpr's Amplitude picks up the Context's column
    // symmetry (Symm), so the programmatic reference must request it too --
    // programmatic ctors are Context-independent (see Tensor::Defaults)
    REQUIRE(root_expr.as<Tensor>() ==
            Tensor(L"Amplitude", bra(IndexList{L"i_1"}), ket(IndexList{L"a_1"}),
                   TensorSymmetries{.column = ColumnSymmetry::Symm}));
  }

  SECTION("scalar-tensor-canon-phase") {
    // Regression: binarize(Product) hardcoded canon_phase=1 for scalar*tensor
    // nodes instead of inheriting the tensor operand's real phase.
    // Subexpression reuse relies on that phase to reconcile sign-differing
    // duplicates, so the wrong phase silently produced a wrong result.
    auto check_phase_inheritance = [](std::wstring_view expr_str) {
      auto res = deserialize<ResultExpr>(
          std::wstring{L"Result{a3;i1,i2} = "} + std::wstring{expr_str},
          {.def_perm_symm = Symmetry::Antisymm});
      auto root = binarize(res);
      auto const& tensor_operand =
          root.left()->is_tensor() ? root.left() : root.right();
      // Sanity: this contraction must exercise phase=-1, else the CHECK
      // below would pass trivially even with the bug reintroduced.
      REQUIRE((int)tensor_operand->canon_phase() == -1);
      CHECK((int)root->canon_phase() == (int)tensor_operand->canon_phase());
    };
    check_phase_inheritance(
        L"1/2 R{a2;i1,i3} g{i3,a3;i2,a2}");  // numeric scalar
    check_phase_inheritance(
        L"R{a2;i1,i3} g{i3,a3;i2,a2} λ");  // Variable scalar
  }

  SECTION("K-conjugated non-Hermitian leaves over a complex basis") {
    // The parity trait normalizes a '꙳' basis-independently: under Even it
    // clears, under None it stays. A kept '꙳' over a complex basis is a
    // different operator's matrix, its own array: a plain leaf whose hash
    // differs from the unmarked twin's.
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.hermiticity = Hermiticity::NonHermitian,
                              .conjugation_parity = ConjugationParity::None});
    REQUIRE(t.base_field() == Field::Complex);

    // Even parity: the star leaves nothing behind, the leaf is the bare one
    Tensor e(L"t", bra{L"a_1"}, ket{L"i_1"});
    Tensor e_star = e;
    REQUIRE(e_star.kconjugate() == 1);
    REQUIRE_FALSE(e_star.kconjugated());
    REQUIRE(EvalExpr{e}.hash_value() == EvalExpr{e_star}.hash_value());

    Tensor t_star = t;
    REQUIRE(t_star.kconjugate() == 1);
    REQUIRE(t_star.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto t_star_tree = binarize(ex<Tensor>(t_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(t_star_tree.leaf());
    REQUIRE(t_star_tree->as_tensor().kconjugated());

    // a directly-constructed EvalExpr must not alias t and t꙳ onto one cache
    // slot
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_star}.hash_value());

    // t and t⁺ share the bare array's slot: '⁺' decodes to {conj,
    // braket_swap} over the unmarked spelling, one provider array for both
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);
    REQUIRE(t_adj.adjointed());
    REQUIRE(EvalExpr{t}.hash_value() == EvalExpr{t_adj}.hash_value());
    Tensor t_adj_star = t;
    REQUIRE(t_adj_star.adjoint() == 1);
    REQUIRE(t_adj_star.kconjugate() == 1);
    REQUIRE(t_adj_star.adjointed());
    REQUIRE(t_adj_star.kconjugated());
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_adj_star}.hash_value());
    REQUIRE(EvalExpr{t_adj}.hash_value() != EvalExpr{t_adj_star}.hash_value());

    // both states: the '꙳' leaf is its own array and the adjoint channel
    // rides its transform
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto both = binarize(ex<Tensor>(t_adj_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(both.leaf());
    REQUIRE_FALSE(both->as_tensor().adjointed());
    REQUIRE(both->as_tensor().kconjugated());
    REQUIRE(both->as_tensor().bra().at(0).label() == L"a_1");
    REQUIRE(both->hash_value() == t_star_tree->hash_value());
  }
  SECTION("an adjointed leaf in a binarized term") {
    // Regression: a tensor leaf can carry the '⁺' state without having been
    // produced by Tensor::adjoint() — e.g. when built from a label string
    // ending in '⁺' (the constructor adopts the mark into the bits).
    //
    // The leaf ctor keys off that state to strip it into the bare spelling
    // plus a {conj, braket_swap} transform, and that path must tolerate a
    // leaf that never went through Tensor::adjoint().
    auto expr = deserialize(L"1/2 g{i_1,i_2;a_1,a_2} t⁺{a_1;i_1} t{a_2;i_2}",
                            {.def_braket_symm = BraKetSymmetry::Nonsymm});
    REQUIRE(expr->is<Product>());

    bool has_adjointed_leaf = false;
    for (auto const& factor : expr->as<Product>().factors())
      has_adjointed_leaf |=
          factor->is<Tensor>() && factor->as<Tensor>().adjointed();
    REQUIRE(has_adjointed_leaf);

    // binarize() must not throw on this term:
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    REQUIRE_NOTHROW(binarize(expr));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("Power op") {
    // a real-basis registry: the tensors here pin BraKetSymmetry::Symm, which
    // over the (default) complex basis no trait derives
    [[maybe_unused]] auto const real_basis_scope = scoped_real_field_context();
    const TensorSymmetries dsym{.perm = Symmetry::Symm,
                                .braket = BraKetSymmetry::Symm,
                                .column = ColumnSymmetry::Symm};
    auto D = ex<Tensor>(L"D", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, dsym);
    auto p = pow(D, rational{-1});
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(p);
    auto base_tree = binarize(D);
    auto p2_tree = binarize(pow(D, rational{-2}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE_FALSE(tree.leaf());
    REQUIRE(tree->op_type() == EvalOp::Power);
    REQUIRE(tree->is_elementwise_power());
    REQUIRE(tree->is_tensor());
    REQUIRE(tree.left().leaf());
    REQUIRE(tree.left()->is_tensor());
    REQUIRE(tree.left()->as_tensor().label() == L"D");
    REQUIRE(tree.right().leaf());
    REQUIRE(tree.right()->is_constant());
    REQUIRE(power_exponent(tree) == rational{-1});
    CHECK(ranges::equal(tree->canon_indices(), tree.left()->canon_indices()));
    CHECK(tree->hash_value() != base_tree->hash_value());
    CHECK(tree->hash_value() != p2_tree->hash_value());
    // round trip, modulo the block canonicalization of the D leaf (its Symm
    // bra-ket symmetry lets it swap bra and ket)
    CHECK(*to_expr(tree) == *pow(to_expr(base_tree), rational{-1}));

    // power of an adjoint-marked base: Power over the bare leaf (no Adjoint
    // node, see EvalExpr(Tensor const&)); over this real basis the adjoint is
    // the exchange of the slots, so the leaf spells it as written and the
    // round trip below is exact
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor t_adj = t;
    (void)t_adj.adjoint();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto pa = binarize(pow(ex<Tensor>(t_adj), rational{2}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(pa->op_type() == EvalOp::Power);
    REQUIRE(pa.left().leaf());
    CHECK(*to_expr(pa) == *pow(ex<Tensor>(t_adj), rational{2}));

    // a product containing a power binarizes; the power's slots are the base's
    auto M = ex<Tensor>(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                        TensorSymmetries{.perm = Symmetry::Antisymm});
    auto ell = ResultExpr(Tensor(L"ℓ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"},
                                 TensorSymmetries{.perm = Symmetry::Antisymm}),
                          M * p);
    REQUIRE_NOTHROW(binarize(ell));
    auto ell_tree = binarize(ell);
    REQUIRE(ell_tree->op_type() == EvalOp::Product);
    CHECK((ell_tree.left()->is_elementwise_power() ||
           ell_tree.right()->is_elementwise_power()));

    // the node phase is the base's canonicalization phase raised to the
    // exponent
    auto const asym = TensorSymmetries{.perm = Symmetry::Antisymm};
    auto A = ex<Tensor>(L"D", bra{L"a_1", L"i_1"}, ket{L"i_2", L"a_2"}, asym);
    auto A_sorted =
        ex<Tensor>(L"D", bra{L"i_1", L"a_1"}, ket{L"i_2", L"a_2"}, asym);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto a2 = binarize(pow(A, rational{2}));
    auto a2_sorted = binarize(pow(A_sorted, rational{2}));
    auto am1 = binarize(pow(A, rational{-1}));
    auto half = binarize(pow(D, rational{1, 2}));
    auto third = binarize(pow(D, rational{1, 3}));
    auto zeroth = binarize(pow(D, rational{0}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(a2.left()->canon_phase() == -1);
    CHECK(a2->canon_phase() == 1);
    CHECK(am1->canon_phase() == -1);
    CHECK(a2->hash_value() == a2_sorted->hash_value());
    CHECK(a2_sorted->canon_phase() == 1);

    // fractional exponents are told apart; a zero exponent is still a node
    CHECK(half->hash_value() != third->hash_value());
    REQUIRE(zeroth->op_type() == EvalOp::Power);
    CHECK(power_exponent(zeroth) == rational{0});
  }

  SECTION("external hyperindices") {
    // t{i1,i2;a1,a3} T2{a2,a3;i1,i2}: i1,i2 appear in bra of t and ket of
    // T2 (multiply-appearing), a3 also multiply-appearing
    auto expr = deserialize(L"t{i1,i2;a1,a3} T2{a2,a3;i1,i2}");

    SECTION("without external indices") {
      // i1,i2,a3 are all contracted; result has only {a1,a2}
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto tree = binarize(expr);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
      REQUIRE(tree->is_tensor());
      auto const& ixs = tree->as_tensor().const_braket() |
                        ranges::views::transform(&Index::label) |
                        ranges::to<container::set<std::wstring_view>>;
      auto expected = std::initializer_list<std::wstring_view>{L"a_1", L"a_2"} |
                      ranges::to<container::set<std::wstring_view>>;
      REQUIRE(ixs == expected);
    }

    SECTION("with external indices") {
      // external={i1,i2}: only a3 is contracted; result has {a1,a2,i1,i2}
      // with i1,i2 in aux
      IndexSet ext;
      ext.emplace(Index{L"i_1"});
      ext.emplace(Index{L"i_2"});
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto tree = binarize(expr, ext);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
      REQUIRE(tree->is_tensor());
      auto const& t = tree->as_tensor();
      auto all_labels = t.const_indices() |
                        ranges::views::transform(&Index::label) |
                        ranges::to<container::set<std::wstring_view>>;
      auto expected = std::initializer_list<std::wstring_view>{L"a_1", L"a_2",
                                                               L"i_1", L"i_2"} |
                      ranges::to<container::set<std::wstring_view>>;
      REQUIRE(all_labels == expected);
      // hyperindices should be in aux (they appear in multiple slots)
      auto aux_labels = t.aux() | ranges::views::transform(&Index::label) |
                        ranges::to<container::set<std::wstring_view>>;
      REQUIRE(aux_labels.contains(L"i_1"));
      REQUIRE(aux_labels.contains(L"i_2"));
    }
  }

  SECTION("Sequant expression") {
    const auto& str_t1 = L"g_{a1,a2}^{a3,a4}";
    const auto& str_t2 = L"t_{a3,a4}^{i1,i2}";
    const auto& t1 = deserialize(str_t1);

    const auto& t2 = deserialize(str_t2);

    const auto& x1 = EvalExpr{t1->as<Tensor>()};
    const auto& x2 = EvalExpr{t2->as<Tensor>()};

    REQUIRE(*t1 == x1.expr()->as<Tensor>());
    REQUIRE(*t2 == x2.expr()->as<Tensor>());

    const auto& x3 = result_expr(x1, x2, EvalOp::Product);

    REQUIRE_NOTHROW(x3.expr()->as<Tensor>());

    const auto& prod_indices =
        x3.expr()->as<Tensor>().const_braket() |
        ranges::views::transform([](const auto& x) { return x.label(); }) |
        ranges::to<container::set<std::wstring_view>>;

    const auto& expected_indices =
        std::initializer_list<std::wstring_view>{L"i_1", L"i_2", L"a_1",
                                                 L"a_2"} |
        ranges::to<container::set<std::wstring_view>>;

    REQUIRE(x3.op_type() == EvalOp::Product);

    REQUIRE(prod_indices == expected_indices);

    const auto t4 = parse_tensor(L"g_{i3,i4}^{a3,a4}");
    const auto t5 = parse_tensor(L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}");

    const auto& x45 = result_expr(EvalExpr{t4}, EvalExpr{t5}, EvalOp::Product);
    const auto& x54 = result_expr(EvalExpr{t5}, EvalExpr{t4}, EvalOp::Product);

    REQUIRE(x45.to_latex() == deserialize(L"I_{a1,a2}^{i1,i2}")->to_latex());
    REQUIRE(x45.to_latex() == x54.to_latex());
  }

  SECTION("Hash value") {
    const auto t1 =
        parse_tensor(L"t_{i1}^{a1}", {.def_perm_symm = Symmetry::Antisymm});
    const auto t2 =
        parse_tensor(L"t_{i2}^{a2}", {.def_perm_symm = Symmetry::Antisymm});
    const auto t3 = parse_tensor(L"t_{i1,i2}^{a1,a2}",
                                 {.def_perm_symm = Symmetry::Antisymm});

    const auto& x1 = EvalExpr{t1};
    const auto& x2 = EvalExpr{t2};

    const auto& x12 = result_expr(x1, x2, EvalOp::Product);
    const auto& x21 = result_expr(x2, x1, EvalOp::Product);

    REQUIRE(x1.hash_value() == x2.hash_value());
    REQUIRE(x12.hash_value() == x21.hash_value());

    const auto& x3 = EvalExpr{t3};

    REQUIRE_FALSE(x1.hash_value() == x3.hash_value());
    REQUIRE_FALSE(x12.hash_value() == x3.hash_value());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree1 = binarize(deserialize(L"A C"));
    auto tree2 = binarize(deserialize(L"A t{a1;i1}"));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    REQUIRE(tree1->hash_value() != tree2->hash_value());
  }

  SECTION("Symmetry of product") {
    // whole bra <-> ket contraction between two antisymmetric tensors
    const auto t1 = parse_tensor(L"g_{i3,i4}^{i1,i2}",
                                 {.def_perm_symm = Symmetry::Antisymm});
    const auto t2 = parse_tensor(L"t_{a1,a2}^{i3,i4}",
                                 {.def_perm_symm = Symmetry::Antisymm});

    const auto x12 = result_expr(EvalExpr{t1}, EvalExpr{t2}, EvalOp::Product);

    // todo:
    // REQUIRE(x12.expr()->as<Tensor>().symmetry() == Symmetry::Antisymm);
    REQUIRE(x12.expr()->as<Tensor>().symmetry() == Symmetry::Nonsymm);

    // whole bra <-> ket contraction between two symmetric tensors
    const auto t3 =
        deserialize(L"g_{i3,i4}^{i1,i2}", {.def_perm_symm = Symmetry::Symm})
            ->as<Tensor>();
    const auto t4 =
        deserialize(L"t_{a1,a2}^{i3,i4}", {.def_perm_symm = Symmetry::Symm})
            ->as<Tensor>();

    const auto x34 = result_expr(EvalExpr{t3}, EvalExpr{t4}, EvalOp::Product);

    // todo:
    // REQUIRE(x34.expr()->as<Tensor>().symmetry() == Symmetry::Symm);
    REQUIRE(x34.expr()->as<Tensor>().symmetry() == Symmetry::Nonsymm);

    // outer product of the same tensor
    const auto t5 =
        deserialize(L"f_{i1}^{a1}", {.def_perm_symm = Symmetry::Nonsymm})
            ->as<Tensor>();
    const auto t6 =
        deserialize(L"f_{i2}^{a2}", {.def_perm_symm = Symmetry::Nonsymm})
            ->as<Tensor>();

    const auto& x56 = result_expr(EvalExpr{t5}, EvalExpr{t6}, EvalOp::Product);

    // todo:
    // REQUIRE(x56.expr()->as<Tensor>().symmetry() == Symmetry::Antisymm);
    REQUIRE(x56.expr()->as<Tensor>().symmetry() == Symmetry::Nonsymm);

    // contraction of some indices from a bra to a ket
    const auto t7 = parse_tensor(L"g_{a1,a2}^{i1,a3}",
                                 {.def_perm_symm = Symmetry::Antisymm});
    const auto t8 =
        parse_tensor(L"t_{a3}^{i2}", {.def_perm_symm = Symmetry::Antisymm});

    const auto x78 = result_expr(EvalExpr{t7}, EvalExpr{t8}, EvalOp::Product);
    REQUIRE(x78.expr()->as<Tensor>().symmetry() == Symmetry::Nonsymm);

    // whole bra <-> ket contraction between symmetric and antisymmetric tensors
    auto const t9 =
        deserialize(L"g_{a1,a2}^{a3,a4}", {.def_perm_symm = Symmetry::Antisymm})
            ->as<Tensor>();
    auto const t10 =
        deserialize(L"t_{a3,a4}^{i1,i2}", {.def_perm_symm = Symmetry::Symm})
            ->as<Tensor>();
    auto const x910 = result_expr(EvalExpr{t9}, EvalExpr{t10}, EvalOp::Product);
    // todo:
    // REQUIRE(x910.expr()->as<Tensor>().symmetry() == Symmetry::Symm);
    REQUIRE(x910.expr()->as<Tensor>().symmetry() == Symmetry::Nonsymm);
  }

  SECTION("Debug") {
    auto t1 = EvalExpr{deserialize(L"O{a_1<i_1,i_2>;a_1<i_3,i_2>}",
                                   {.def_perm_symm = Symmetry::Nonsymm})
                           ->as<Tensor>()};
    auto t2 = EvalExpr{deserialize(L"O{a_2<i_1,i_2>;a_2<i_3,i_2>}",
                                   {.def_perm_symm = Symmetry::Nonsymm})
                           ->as<Tensor>()};

    REQUIRE_NOTHROW(result_expr(t1, t2, EvalOp::Product));
  }
}

TEST_CASE("leaf lowering of the core states", "[eval_expr]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  SECTION("adjointed leaf: {conj, braket_swap} over the bare array") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});  // NonHermitian by default
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);
    REQUIRE(t_adj.adjointed());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t_adj));
    auto bare = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree.leaf());
    REQUIRE_FALSE(tree->as_tensor().adjointed());  // the bare array
    REQUIRE(tree->as_tensor().bra()[0].label() == L"a_1");
    REQUIRE(tree->canon_transform().conj);
    REQUIRE(tree->canon_transform().braket_swap);
    REQUIRE(tree->hash_value() == bare->hash_value());  // one slot
    REQUIRE(bare->canon_transform().trivial());
    // the decoder is an involution on a leaf: the denoted spelling is what
    // was written
    REQUIRE(*tree->denoted_expr() == t_adj);
  }
  SECTION("adjointed leaf moved onto a real basis: it arrives normalized") {
    // transform_indices re-normalizes the states over the new slots' field, so
    // the decoder never sees a '⁺' over a real basis: the coset rule has
    // already traded it for a '꙳' on the slots as written, which the parity
    // keeps (None) or consumes (the default)
    container::map<Index, Index> to_real{
        {Index{L"a_1"}, idx(L"a_1", Field::Real)},
        {Index{L"i_1"}, idx(L"i_1", Field::Real)}};
    auto moved_onto_real = [&to_real](ConjugationParity parity) {
      Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
               TensorSymmetries{.conjugation_parity = parity});
      REQUIRE(t.adjoint() == 1);
      REQUIRE(t.adjointed());
      REQUIRE(t.transform_indices(to_real));
      t.reset_tags();  // transform_indices tags the replaced indices
      REQUIRE(t.base_field() == Field::Real);
      REQUIRE_FALSE(t.adjointed());
      REQUIRE(t.bra()[0].label() == L"a_1");
      return t;
    };

    // the default parity consumes the star: a stateless leaf, no transform
    Tensor even = moved_onto_real(ConjugationParity::Even);
    REQUIRE_FALSE(even.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto even_tree = binarize(ex<Tensor>(even));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(even_tree.leaf());
    REQUIRE_FALSE(even_tree->op_type().has_value());
    REQUIRE(even_tree->canon_transform().trivial());
    REQUIRE(even_tree->as_tensor().bra()[0].label() == L"a_1");

    // an indefinite parity keeps it: the conjugation rides the transform,
    // over the identity layout
    Tensor none = moved_onto_real(ConjugationParity::None);
    REQUIRE(none.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto none_tree = binarize(ex<Tensor>(none));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(none_tree.leaf());
    REQUIRE(none_tree->canon_transform().conj);
    REQUIRE_FALSE(none_tree->canon_transform().braket_swap);
    REQUIRE_FALSE(none_tree->as_tensor().adjointed());
    REQUIRE_FALSE(none_tree->as_tensor().kconjugated());
    REQUIRE(none_tree->as_tensor().bra()[0].label() == L"a_1");
  }
  SECTION("K-conjugated leaf over a real basis: {conj}, identity layout") {
    // a real-basis '꙳' survives only where the hermiticity is indefinite:
    // a definite one reduces it through the coset rule
    Index a = idx(L"a_1", Field::Real), i = idx(L"i_1", Field::Real);
    Tensor t(L"t", bra{a}, ket{i},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor tk = t;
    REQUIRE(tk.kconjugate() == 1);
    REQUIRE(tk.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(tk));
    auto bare = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree.leaf());
    REQUIRE_FALSE(tree->as_tensor().kconjugated());
    REQUIRE(tree->canon_transform() == CanonTransform{.conj = true});
    REQUIRE(tree->hash_value() == bare->hash_value());
    REQUIRE(tree->canon_indices() == bare->canon_indices());  // identity layout
    REQUIRE(*tree->denoted_expr() == tk);
  }
  SECTION("K-conjugated leaf over a complex basis: its own array") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor tk = t;
    REQUIRE(tk.kconjugate() == 1);
    REQUIRE(tk.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto a = binarize(ex<Tensor>(t)), b = binarize(ex<Tensor>(tk));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(b.leaf());
    REQUIRE(b->as_tensor().kconjugated());    // the state stays on the array
    REQUIRE(b->canon_transform().trivial());  // no transform serves it
    REQUIRE(a->hash_value() != b->hash_value());
    REQUIRE(*b->denoted_expr() == tk);
  }
  SECTION("adjointed and K-conjugated over a complex basis compose") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    Tensor tk = t;
    REQUIRE(tk.kconjugate() == 1);
    Tensor both = tk;
    REQUIRE(both.adjoint() == 1);
    REQUIRE(both.adjointed());
    REQUIRE(both.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(both)), star = binarize(ex<Tensor>(tk));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree.leaf());
    // the array is t꙳, the transform the adjoint channel
    REQUIRE(tree->as_tensor().kconjugated());
    REQUIRE_FALSE(tree->as_tensor().adjointed());
    REQUIRE(tree->canon_transform().conj);
    REQUIRE(tree->canon_transform().braket_swap);
    REQUIRE(tree->hash_value() == star->hash_value());
    REQUIRE(*tree->denoted_expr() == both);
  }
  SECTION("the parity keys an imaginary leaf apart over a real basis") {
    // Over a real basis an Odd-parity array is imaginary while the Even one is
    // real, so the two are different arrays and the leaf hash separates them.
    // Even and None both leave the hash alone, so they key alike.
    Index i = idx(L"i_1", Field::Real), j = idx(L"i_2", Field::Real);
    auto leaf_hash = [&i, &j](ConjugationParity parity) {
      Tensor p(L"p", bra{i}, ket{j},
               TensorSymmetries{.conjugation_parity = parity});
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto tree = binarize(ex<Tensor>(p));
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
      REQUIRE(tree.leaf());
      return tree->hash_value();
    };
    REQUIRE(leaf_hash(ConjugationParity::Even) !=
            leaf_hash(ConjugationParity::Odd));
    REQUIRE(leaf_hash(ConjugationParity::Even) ==
            leaf_hash(ConjugationParity::None));
  }
}

TEST_CASE("no eval op performs a conjugation", "[eval_expr]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  // The IR's conjugation-bearing ops are the projections Re and Im, which are
  // not invertible and so cannot ride in a CanonTransform. Every invertible
  // conjugation channel does ride there: an adjointed factor lowers to a leaf
  // over the bare array carrying {conj, braket_swap}, not to an op of its own.
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
  Tensor t_adj = t;
  REQUIRE(t_adj.adjoint() == 1);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node =
      binarize(ex<Tensor>(t_adj) * ex<Tensor>(L"w", bra{L"a_1"}, ket{L"i_1"}));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(node->op_type() == EvalOp::Product);
  std::size_t leaves = 0;
  std::size_t conjugating_leaves = 0;
  node.visit([&leaves, &conjugating_leaves](auto const& n) {
    if (n.leaf()) {
      ++leaves;
      if (n->canon_transform().conj) {
        ++conjugating_leaves;
        // the adjoint channel: a bra-ket exchange over an array the provider
        // serves unconjugated
        REQUIRE(n->canon_transform().braket_swap);
        REQUIRE_FALSE(n->as_tensor().adjointed());
        REQUIRE_FALSE(n->as_tensor().kconjugated());
      }
      return;
    }
    REQUIRE((n->op_type() == EvalOp::Product || n->op_type() == EvalOp::Sum ||
             n->op_type() == EvalOp::RealPart ||
             n->op_type() == EvalOp::ImagPart));
  });
  REQUIRE(leaves == 2);  // no sentinel child
  REQUIRE(conjugating_leaves == 1);
}

TEST_CASE("eval_expr_node_slice_mask_typed", "[EvalExpr][batched-here]") {
  using namespace sequant;
  auto const tnsr =
      parse_tensor(L"g{i_1,a_1;i_2,a_2}", {.def_perm_symm = Symmetry::Nonsymm});
  EvalExpr node{tnsr};
  container::svector<std::pair<Index, BatchModeType>> modes{
      {Index{L"a_1"}, BatchModeType::Contracted},
      {Index{L"i_1"}, BatchModeType::External}};
  node.set_node_slice_mask(modes);
  REQUIRE(node.node_slice_mask().size() == 2);
  REQUIRE(node.node_slice_mask()[0].second == BatchModeType::Contracted);
  REQUIRE(node.node_slice_mask()[1].second == BatchModeType::External);
}

// The cases below build eval trees straight from expressions: the head layout
// is irrelevant to what they check (slot identity, transforms, phases), so
// the deprecated binarize(ExprPtr) is used on purpose.
SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
TEST_CASE("eval_expr_conjugation_state_identity",
          "[EvalExpr][conjugate-state]") {
  // The Gram overlap of a Hermitian C is spelled with the adjoint (a
  // Conjugate tensor's bundles are not interchangeable, so a ket-ket
  // contraction is not a spelling of it): C{a_1;p} C⁺{p;a_2} and
  // C⁺{p;a_1} C{a_2;p} are one tensor S up to the named relabeling
  // a_1 <-> a_2 (its slots are "index of the unconjugated factor" and "index
  // of the conjugated factor"), so they may share one eval-node hash and
  // cache slot -- but only if their canonical layouts put the unconjugated
  // factor's index in the same slot. Were the two spellings the same coloured
  // graph with an automorphism exchanging the factors, bliss would pin the
  // slot order by the index labels alone, and the same buffer would be read
  // as S by one occurrence and as S^T* by the other (identical only for real
  // C). Regression for the Kramers-restricted CSV inter-pair overlap.
  using namespace sequant;

  auto C = [](std::wstring_view ext) {
    return ex<Tensor>(L"C", bra{Index{ext}}, ket{Index{L"p_1"}},
                      Symmetry::Nonsymm, BraKetSymmetry::Conjugate,
                      ColumnSymmetry::Symm);
  };
  // the adjoint of a Hermitian C is the bra/ket exchange, no state left
  auto Cadj = [&C](std::wstring_view ext) {
    auto t = C(ext);
    REQUIRE(t->as<Tensor>().adjoint() == 1);
    REQUIRE_FALSE(t->as<Tensor>().adjointed());
    REQUIRE(t->as<Tensor>().bra()[0].label() == L"p_1");
    return t;
  };

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto A = binarize(C(L"a_1") * Cadj(L"a_2"));  // S(a_1, a_2)
  auto B = binarize(Cadj(L"a_1") * C(L"a_2"));  // S(a_2, a_1)
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(!A.leaf());
  REQUIRE(!B.leaf());
  REQUIRE(A->canon_indices().size() == 2);
  REQUIRE(B->canon_indices().size() == 2);

  // slot of the unconjugated factor's index in the canonical layout
  auto unconj_slot = [](auto const& node, Index const& unconj_idx) {
    auto const& ci = node->canon_indices();
    return std::distance(ci.begin(),
                         std::find(ci.begin(), ci.end(), unconj_idx));
  };
  auto const slot_A = unconj_slot(A, Index{L"a_1"});
  auto const slot_B = unconj_slot(B, Index{L"a_2"});
  REQUIRE(slot_A < 2);
  REQUIRE(slot_B < 2);

  // same value up to relabeling: one identity ...
  REQUIRE(A->hash_value() == B->hash_value());
  // ... and a layout that agrees on which slot is the unconjugated one
  REQUIRE(slot_A == slot_B);

  // a genuinely different tensor is kept apart: for a non-Hermitian C the
  // '⁺' stays, and C{a_1;p} C⁺{p;a_2} (an adjointed leaf inside the product)
  // is not C{a_1;p} C{p;a_2}
  auto N = [](std::wstring_view b, std::wstring_view k) {
    return ex<Tensor>(L"C", bra{Index{b}}, ket{Index{k}},
                      TensorSymmetries{.hermiticity = Hermiticity::NonHermitian,
                                       .column = ColumnSymmetry::Symm});
  };
  auto Nadj = N(L"a_2", L"p_1");
  REQUIRE(Nadj->as<Tensor>().adjoint() == 1);
  REQUIRE(Nadj->as<Tensor>().adjointed());
  REQUIRE(Nadj->as<Tensor>().bra()[0].label() == L"p_1");
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto D = binarize(N(L"a_1", L"p_1") * Nadj);
  auto E = binarize(N(L"a_1", L"p_1") * N(L"p_1", L"a_2"));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(D.right().leaf());
  REQUIRE(D->hash_value() != E->hash_value());
}

TEST_CASE("kramers_leaf_slot_identity", "[eval_expr][kramers]") {
  // a down-first Kramers leaf occupies the SAME eval slot as its up-first
  // partner; the difference rides the CanonTransform as {conj, phase}
  using namespace sequant;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto mk = [](std::wstring_view b, std::wstring_view k, BraKetSymmetry bks) {
    return Tensor(L"C", bra{Index(b)}, ket{Index(k)}, Symmetry::Nonsymm, bks,
                  ColumnSymmetry::Nonsymm, KramersSymmetry::TimeReversal);
  };
  // OFF by default: the leaf provider serves as_tensor()'s spelling and the
  // transform's labels must match the sibling's, which a flavor flip breaks
  // unless the provider itself aliases the partner block (serving-level
  // aliasing); so the leaf fold is an explicit opt-in
  {
    auto resetter = set_scoped_default_context(ctx);
    EvalExpr eu{mk(L"a↑_1", L"a↓_2", BraKetSymmetry::Nonsymm)};
    EvalExpr ed{mk(L"a↓_1", L"a↑_2", BraKetSymmetry::Nonsymm)};
    REQUIRE(eu.hash_value() != ed.hash_value());
    REQUIRE(ed.canon_transform().trivial());
  }
  ctx.set(
      CanonicalizeOptions{.fold_kramers_eval_leaves =
                              CanonicalizeOptions::FoldKramersEvalLeaves::Yes});
  auto resetter = set_scoped_default_context(ctx);
  for (auto bks : {BraKetSymmetry::Nonsymm, BraKetSymmetry::Conjugate}) {
    // C(down,up) = -conj(C(up,down)): one slot flipped from down. For a
    // Hermitian C the up-row spelling is reached by the braket move instead
    // (C{a↓;a↑} = conj C{a↑;a↓} swapped): conj, no phase
    EvalExpr eu{mk(L"a↑_1", L"a↓_2", bks)};
    EvalExpr ed{mk(L"a↓_1", L"a↑_2", bks)};
    REQUIRE(ed.canon_transform().conj != eu.canon_transform().conj);
    if (bks == BraKetSymmetry::Nonsymm) {
      REQUIRE(eu.hash_value() == ed.hash_value());
      REQUIRE(ed.canon_transform().phase == -eu.canon_transform().phase);
    } else {
      REQUIRE(ed.canon_transform().braket_swap);
      REQUIRE(ed.canon_transform().phase == eu.canon_transform().phase);
    }
    // C(down,down) = +conj(C(up,up)): two slots flipped from down
    EvalExpr euu{mk(L"a↑_1", L"a↑_2", bks)};
    EvalExpr edd{mk(L"a↓_1", L"a↓_2", bks)};
    REQUIRE(euu.hash_value() == edd.hash_value());
    REQUIRE(edd.canon_transform().conj != euu.canon_transform().conj);
    REQUIRE(edd.canon_transform().phase == euu.canon_transform().phase);
    // distinct families stay distinct
    REQUIRE(eu.hash_value() != euu.hash_value());
    // T19 layer 2 contract: expr() is the FOLDED (up-row) spelling -- what a
    // leaf provider fetches -- while canon_indices() carries the as-written
    // labels (the parent contracts by label; TA matches annotations, not
    // spellings), so the served up block + {conj, phase} denotes the
    // as-written value. (For a Conjugate tensor the braket fold may swap
    // bra and ket, so label SETS are the invariant, not slot positions.)
    auto const& ed_t = ed.expr()->as<Tensor>();
    REQUIRE(!ed_t.kconjugated());
    REQUIRE(!ed_t.adjointed());
    auto has_space = [&](auto rng, std::wstring_view sp) {
      return ranges::any_of(
          rng, [&](Index const& i) { return i.space() == isr->retrieve(sp); });
    };
    REQUIRE(has_space(ed_t.const_indices(), L"a↑"));
    // the served spelling is the Kramers-canonical one: bra slot up
    REQUIRE(ed_t.bra()[0].space() == isr->retrieve_space(L"a↑"));
    auto has_label = [&](std::wstring_view lbl) {
      return ranges::any_of(ed.canon_indices(), [&](Index const& i) {
        return i.full_label() == lbl;
      });
    };
    REQUIRE(has_label(L"a↓_1"));
    REQUIRE(has_label(L"a↑_2"));
  }
  // a product over a down-flavored dummy still contracts: the folded leaf
  // shares the up-row slot while its labels keep matching the sibling
  {
    auto f_dn = Tensor(L"f", bra{Index(L"a↓_1")}, ket{Index(L"i↑_1")},
                       Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                       ColumnSymmetry::Nonsymm, KramersSymmetry::TimeReversal);
    auto f_up = Tensor(L"f", bra{Index(L"a↑_1")}, ket{Index(L"i↓_1")},
                       Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                       ColumnSymmetry::Nonsymm, KramersSymmetry::TimeReversal);
    auto t = Tensor(L"t", bra{Index(L"i↑_1")}, ket{Index(L"a↓_1")},
                    Symmetry::Nonsymm, BraKetSymmetry::Nonsymm,
                    ColumnSymmetry::Nonsymm, KramersSymmetry::TimeReversal);
    auto root = binarize(ex<Tensor>(f_dn) * ex<Tensor>(t));
    REQUIRE(root->result_type() == ResultType::Scalar);
    REQUIRE(root->canon_indices().empty());
    auto const& f_leaf =
        root.left()->as_tensor().label() == L"f" ? root.left() : root.right();
    REQUIRE(f_leaf->hash_value() == EvalExpr{f_up}.hash_value());
    REQUIRE(f_leaf->canon_transform().conj);
    REQUIRE(f_leaf->canon_transform().phase == -1);
  }
}

TEST_CASE("sum_placeholder_is_spelled_in_its_layout", "[eval_expr][sum]") {
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  // A sum hands up its FIRST summand's canonical layout, and the tensor it
  // spells for its value (expr()) is what an enclosing network sees for the
  // opaque node -- so that placeholder must be spelled in that layout. Two
  // relabeled copies of one sum share a slot (their hashes are label-blind)
  // while their values are transposes of each other: here the column
  // symmetry of g lets a_1 and a_2 trade bra slots, and the ket slots (a_3 vs
  // p_1) pin which column is which. Spelled label-sorted instead, the two
  // placeholders coincide, an enclosing product gets ONE slot with ONE
  // layout for both, and a cache serves one value for the other untransposed
  // (HSeOH PNS-CCD, (vv|vv) ladder kept 4-center, 2026-09-04).
  using namespace sequant;
  auto const bracket = [](std::wstring const& cs) {
    return ex<Sum>(ExprPtrList{deserialize(L"g{a_5,a_6;a_7,a_8}:N-C-S" + cs),
                               deserialize(L"h{a_5,a_6;a_7,a_8}:N-C-S" + cs)});
  };
  auto const product = [](ExprPtr const& sum) {
    return ex<Product>(
        ExprPtrList{sum, deserialize(L"t{a_3,p_1;i_1,i_2}:A-N-S")});
  };
  auto const PA = binarize(product(
      bracket(L" * C{a_5;a_1} * C{a_6;a_2} * C{a_7;a_3} * C{a_8;p_1}")));
  auto const PB = binarize(product(
      bracket(L" * C{a_5;a_2} * C{a_6;a_1} * C{a_7;a_3} * C{a_8;p_1}")));
  auto const swap12 = [](Index::index_vector v) {
    for (auto& ix : v) {
      if (ix.label() == L"a_1")
        ix = Index(L"a_2");
      else if (ix.label() == L"a_2")
        ix = Index(L"a_1");
    }
    return v;
  };
  auto const slots = [](EvalExpr const& e) {
    auto const& t = e.expr()->as<Tensor>();
    Index::index_vector v;
    for (auto const& ix : t.bra()) v.push_back(ix);
    for (auto const& ix : t.ket()) v.push_back(ix);
    for (auto const& ix : t.aux()) v.push_back(ix);
    return v;
  };
  auto const labels = [](Index::index_vector const& v) {
    std::wstring out;
    for (auto const& ix : v) out += std::wstring(ix.full_label()) + L" ";
    return toUtf8(out);
  };
  auto const& SA = *PA.left();
  auto const& SB = *PB.left();
  REQUIRE(SA.op_type() == EvalOp::Sum);
  REQUIRE(SB.op_type() == EvalOp::Sum);
  INFO("SA layout " << labels(SA.canon_indices()) << " placeholder "
                    << toUtf8(to_latex(SA.expr())));
  INFO("SB layout " << labels(SB.canon_indices()) << " placeholder "
                    << toUtf8(to_latex(SB.expr())));
  // relabeled copies of one sum: one slot, transposed layouts
  REQUIRE(SA.hash_value() == SB.hash_value());
  REQUIRE(SB.canon_indices() == swap12(SA.canon_indices()));
  // the placeholders spell the layouts
  CHECK(slots(SB) == swap12(slots(SA)));
  // ... so the enclosing products share a slot with transposed layouts too
  INFO("PA layout " << labels(PA->canon_indices()));
  INFO("PB layout " << labels(PB->canon_indices()));
  REQUIRE(PA->hash_value() == PB->hash_value());
  CHECK(PB->canon_indices() == swap12(PA->canon_indices()));
  // The two products ARE one cache slot: the buffer of one, read under the
  // other's labels, is the other's value (plain externals -- the layout
  // fingerprint is relabeling-invariant on purpose)
  using node_t = std::remove_cvref_t<decltype(PA)>;
  TreeNodeEqualityComparator<node_t> same;
  REQUIRE(same(PA, PA));
  REQUIRE(same(PA, PB));
}

TEST_CASE("twins_whose_composites_carry_the_swapped_externals_are_two_slots",
          "[eval_expr][cache]") {
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  // CSV/PNS composites carry their pair as proto indices, so the inner tile
  // at outer position (p,q) is the pair-(p,q) block: two relabeled twins that
  // lay the pair out as (i_1,i_2) and (i_2,i_1) are NOT value-compatible --
  // served for each other, one gets the pair-(q,p) blocks. The layout
  // fingerprint (ids of the externals AND of every composite's protos, in
  // layout order) tells them apart: binarize folds it into the node id, so
  // the twins are two slots outright, and the comparator refuses them too.
  // Measured 2026-09-05 (HSeOH PNS-MP1, residual block 3 after the brackets
  // were optimized): C†.(g.C) laid out (i↑_1,i↑_2;..) was served to its twin
  // laid out (i↑_2,i↑_1;..): |R| 0.579 instead of 0.293, E 7 % off.
  using namespace sequant;
  auto const X = binarize(ex<Product>(
      ExprPtrList{deserialize(L"f{i_1;i_3}:N-N-N"),
                  deserialize(L"t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_3,i_2}:N-N-N")}));
  auto const Y = binarize(ex<Product>(
      ExprPtrList{deserialize(L"f{i_2;i_3}:N-N-N"),
                  deserialize(L"t{a_1<i_1,i_2>,a_2<i_1,i_2>;i_3,i_1}:N-N-N")}));
  auto const labels = [](Index::index_vector const& v) {
    std::wstring out;
    for (auto const& ix : v) out += std::wstring(ix.full_label()) + L" ";
    return toUtf8(out);
  };
  INFO("X layout " << labels(X->canon_indices()));
  INFO("Y layout " << labels(Y->canon_indices()));
  REQUIRE(X->canon_indices() != Y->canon_indices());
  REQUIRE(X->layout_fingerprint() != Y->layout_fingerprint());
  REQUIRE(X->hash_value() != Y->hash_value());  // two slots
  using node_t = std::remove_cvref_t<decltype(X)>;
  TreeNodeEqualityComparator<node_t> same;
  REQUIRE(same(X, X));
  REQUIRE_FALSE(same(X, Y));
}

SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

// binarize(Sum const&, ...)'s make_sum
// lambda used to capture its prefix-hash range (imed_hashes(hvals)) as a
// LAZY, stateful view; ranges::at(hs, ++i) re-begin()s that view on every
// access, which re-drives inits' internal mutable `++n` counter and
// silently drops the LAST summand from the running hash -- so two Sums
// differing only in their last summand collided on hash_value(). The
// Product path in this same file already materializes its prefix-hash
// range eagerly (`auto const hs = imed_hashes(hvals) | ranges::to_vector;`)
// and was unaffected.
TEST_CASE("Sum-node hash is sensitive to every summand",
          "[eval][binarize][hash]") {
  using namespace sequant;

  auto const root = [](std::wstring_view s) {
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    return binarize(deserialize(s));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  };

  auto const last_a = root(L"(a * b) + c");  // differ in LAST summand
  auto const last_b = root(L"(a * b) - c");
  auto const first_a = root(L"c + (a * b)");  // differ in FIRST summand
  auto const first_b = root(L"d + (a * b)");  // (already worked pre-fix)

  CHECK(last_a->hash_value() != last_b->hash_value());
  CHECK(first_a->hash_value() != first_b->hash_value());
  // (a*b)+c and c+(a*b) are the same multiset of summands -> same hash.
  CHECK(last_a->hash_value() == first_a->hash_value());
}

TEST_CASE("proto-index-liveness", "[eval]") {
  using namespace sequant;

  // the `S` braket default below is derivable only over a real basis
  auto real_basis = tests::scoped_real_basis();
  // an energy term with proto-carrying indices: past g, i_2 appears only in
  // the proto tuples of a_3 and a_4, so it must stay live on the g * (...) node
  auto const expr = deserialize<ExprPtr>(
      L"(((("
      L"g{a5,a6;i1,i2}"
      L" * (C{a3<i1,i2>;a5}"
      L"    * ((C{a1<i2>;a7} * t{a1<i2>;i2}) * C{a7;a3<i1,i2>})))"
      L" * C{a4<i1,i2>;a6})"
      L" * C{a8;a4<i1,i2>})"
      L" * C{a2<i1>;a8})"
      L" * t{a2<i1>;i1}",
      // fixture contracts bra-with-bra, so tensors must be bra-ket symmetric
      {.def_braket_symm = BraKetSymmetry::Symm});

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto const root = sequant::binarize(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  // root = (((((g * B) * C) * C) * C) * t); descend to the g * B node.
  auto const& g_node = root.left().left().left().left();
  REQUIRE(g_node.left()->is_tensor());
  CHECK(g_node.left()->as_tensor().label() == L"g");

  auto const has_label = [](auto const& ixs, std::wstring_view lbl) {
    return ranges::any_of(ixs,
                          [lbl](auto const& ix) { return ix.label() == lbl; });
  };
  // i_1 is a slot index of the trailing t; i_2 is live only via proto tuples.
  CHECK(has_label(g_node->canon_indices(), L"i_1"));
  CHECK(has_label(g_node->canon_indices(), L"i_2"));
}

TEST_CASE("canon_transform_algebra", "[EvalExpr][conj-transform]") {
  using sequant::CanonTransform;
  CanonTransform id{};
  REQUIRE(id.trivial());
  REQUIRE(id.phase == 1);
  REQUIRE_FALSE(id.conj);
  REQUIRE_FALSE(id.braket_swap);

  CanonTransform c{.phase = 1, .conj = true, .braket_swap = false};
  CanonTransform s{.phase = -1, .conj = false, .braket_swap = true};
  REQUIRE_FALSE(c.trivial());

  // composition: phases multiply, conj/swap are Z2 (xor)
  auto cs = compose(c, s);
  REQUIRE(cs.phase == -1);
  REQUIRE(cs.conj);
  REQUIRE(cs.braket_swap);
  REQUIRE(compose(c, c).trivial());  // involution
  // structural salt: conj/swap enter, phase does _not_ (hoistable)
  REQUIRE(CanonTransform{.phase = -1}.structural_salt() ==
          CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != s.structural_salt());
  REQUIRE(c.structural_salt() != cs.structural_salt());
}

TEST_CASE("eval_expr_carries_canon_transform", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  EvalExpr ee{t};
  REQUIRE(ee.canon_transform().trivial());
  REQUIRE(ee.canon_phase() == ee.canon_transform().phase);  // compat accessor
}

TEST_CASE("leaf_transform_channels", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  // The channel table of the leaf decoder, one block per arm.

  // A Hermitian core has no state to decode: its two orientations are stored
  // as written and no transform serves either. They share one slot, because
  // one provider array serves both -- a leaf's slot identity is label-blind,
  // and it is the annotation that wires each mode to its contraction
  // partner. The spellings they denote differ, and so do the annotations.
  Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor g_swapped = g;
  REQUIRE(g_swapped.kconjugate() == 1);  // the Even parity consumes the star
  REQUIRE_FALSE(g_swapped.kconjugated());
  REQUIRE(g_swapped.adjoint() == 1);  // a Hermitian core: a bundle exchange
  REQUIRE_FALSE(g_swapped.adjointed());
  REQUIRE(g_swapped.bra()[0].label() == L"p_3");
  EvalExpr eg{g}, egs{g_swapped};
  REQUIRE(eg.canon_transform().trivial());
  REQUIRE(egs.canon_transform().trivial());
  REQUIRE(egs.as_tensor().bra()[0].label() == L"p_3");  // stored as written
  REQUIRE(eg.hash_value() == egs.hash_value());         // one slot
  REQUIRE_FALSE(eg.denoted_expr()->as<Tensor>() ==
                egs.denoted_expr()->as<Tensor>());
  REQUIRE(eg.indices_annot() != egs.indices_annot());

  // '⁺' over a NonHermitian core: the bare array plus {conj, swap}
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor t_adj = t;
  REQUIRE(t_adj.adjoint() == 1);
  EvalExpr et{t}, eta{t_adj};
  REQUIRE(et.hash_value() == eta.hash_value());
  REQUIRE(eta.canon_transform().conj);
  REQUIRE(eta.canon_transform().braket_swap);
  REQUIRE(eta.as_tensor().decorated_label() == L"t");  // bare spelling

  // '⁺' over a kept '꙳': an indefinite parity leaves the star on the array
  // over a complex basis, and the adjoint channel rides the transform
  Tensor n(L"n", bra{L"a_1"}, ket{L"i_1"},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  REQUIRE(n.base_field() == Field::Complex);
  Tensor n_star = n;
  REQUIRE(n_star.kconjugate() == 1);
  REQUIRE(n_star.kconjugated());
  Tensor n_adj_star = n_star;
  REQUIRE(n_adj_star.adjoint() == 1);
  REQUIRE(n_adj_star.adjointed());
  EvalExpr en{n}, ens{n_star}, enas{n_adj_star};
  REQUIRE(ens.hash_value() != en.hash_value());    // an array of its own
  REQUIRE(enas.hash_value() == ens.hash_value());  // served from that array
  REQUIRE(ens.canon_transform().trivial());
  REQUIRE(enas.as_tensor().kconjugated());
  REQUIRE_FALSE(enas.as_tensor().adjointed());
  REQUIRE(enas.canon_transform().conj);
  REQUIRE(enas.canon_transform().braket_swap);

  // an Even-parity star is consumed by the trait, so the decoder never meets
  // it: the starred spelling is the plain one
  Tensor s(L"s", bra{L"i_1"}, ket{L"a_1"},
           TensorSymmetries{.conjugation_parity = ConjugationParity::Even});
  Tensor s_star = s;
  REQUIRE(s_star.kconjugate() == 1);
  REQUIRE_FALSE(s_star.kconjugated());
  REQUIRE(EvalExpr{s_star}.canon_transform() == EvalExpr{s}.canon_transform());

  // ToT: the same decoder runs before canonicalize_slots
  auto ct = deserialize(L"C{a_1<i_1>;i_1}:N-N-S")->as<Tensor>();
  REQUIRE(ranges::any_of(ct.const_indices(), &Index::has_proto_indices));
  Tensor ct_adj = ct;
  REQUIRE(ct_adj.adjoint() == 1);
  REQUIRE(ct_adj.adjointed());
  EvalExpr ec{ct}, eca{ct_adj};
  REQUIRE(ec.hash_value() == eca.hash_value());  // one slot
  REQUIRE_FALSE(eca.as_tensor().adjointed());
  REQUIRE(eca.canon_transform().conj);
  REQUIRE(eca.canon_transform().braket_swap);
}

// The hash contract stated directly: a leaf's slot hash is the hash of the
// array it stores -- label, slots, the state byte of the stored spelling, and
// the AntiSymm conjugation-symmetry term, nothing else (hash_terminal_tensor).
// leaf_transform_channels above exercises the same channels through
// binarize(); this pins the contract at EvalExpr construction directly.
TEST_CASE("leaf_slot_identity_is_the_stored_array",
          "[EvalExpr][conj-transform]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  auto hash_of = [](Tensor const& t) { return EvalExpr{t}.hash_value(); };

  // the two transform cases share the bare array's slot
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
  Tensor t_adj = t;
  REQUIRE(t_adj.adjoint() == 1);
  REQUIRE(hash_of(t_adj) == hash_of(t));

  Index ra = idx(L"a_1", Field::Real), ri = idx(L"i_1", Field::Real);
  Tensor r(L"r", bra{ra}, ket{ri},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  Tensor rk = r;
  REQUIRE(rk.kconjugate() == 1);
  REQUIRE(hash_of(rk) == hash_of(r));

  // the complex-basis K-conjugate is its own array
  Tensor c(L"c", bra{L"a_1"}, ket{L"i_1"},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  Tensor ck = c;
  REQUIRE(ck.kconjugate() == 1);
  REQUIRE(hash_of(ck) != hash_of(c));

  // the two orientations of a Conjugate (Hermitian) tensor sharing one
  // provider array's slot is pinned by leaf_transform_channels above.

  // an Odd-parity real-basis array keying apart from the Even one is pinned
  // by "the parity keys an imaginary leaf apart over a real basis" above.
}

TEST_CASE("conj_hoisting_structural_identity", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  // over a real basis with an indefinite parity the conjugate of a matrix
  // element is the K-conjugate with the slots in place, so every factor
  // decodes to {conj} and the conjugation is a whole-node transform
  auto rx = [](std::wstring_view l) { return idx(l, Field::Real); };
  auto rt = [&rx](std::wstring_view lbl, std::wstring_view b,
                  std::wstring_view k) {
    return ex<Tensor>(
        lbl, bra{rx(b)}, ket{rx(k)},
        TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  };

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto AB = binarize(rt(L"A", L"i_1", L"a_1") * rt(L"B", L"a_1", L"i_2"));
  auto ABc =
      binarize(conjugate(rt(L"A", L"i_1", L"a_1") * rt(L"B", L"a_1", L"i_2")));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(AB->hash_value() == ABc->hash_value());  // uniform {conj} hoists
  REQUIRE(ABc->canon_transform() == CanonTransform{.conj = true});
  REQUIRE_FALSE(AB->canon_transform().conj);

  // mixed marks salt
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto CD = binarize(rt(L"C", L"i_1", L"a_1") * rt(L"D", L"a_1", L"i_2"));
  auto CDc =
      binarize(rt(L"C", L"i_1", L"a_1") * conjugate(rt(L"D", L"a_1", L"i_2")));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(CD->hash_value() != CDc->hash_value());

  // a uniformly conjugated sum hoists the same way
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto S = binarize(rt(L"D", L"i_1", L"a_1") + rt(L"E", L"i_1", L"a_1"));
  auto Sc =
      binarize(conjugate(rt(L"D", L"i_1", L"a_1") + rt(L"E", L"i_1", L"a_1")));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(S->hash_value() == Sc->hash_value());
  REQUIRE(Sc->canon_transform().conj);

  SECTION("the adjoint of a contraction gets its own slot") {
    // over a complex basis the conjugate is the adjoint: every factor's
    // transform carries the bundle exchange, which respells the node's own
    // result, so the conjugation is not hoisted and the adjointed network
    // keeps its own slot
    auto X = ex<Tensor>(L"X", bra{L"i_1"}, ket{L"a_1"});
    auto Y = ex<Tensor>(L"Y", bra{L"a_1"}, ket{L"i_2"});
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto XY = binarize(X->clone() * Y->clone());
    auto XYc = binarize(conjugate(X->clone() * Y->clone()));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(XY->hash_value() != XYc->hash_value());
    REQUIRE_FALSE(XYc->canon_transform().conj);

    // the sum site reads the same rule: the adjoint of a sum salts
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto SXY = binarize(X->clone() + Y->clone());
    auto SXYc = binarize(conjugate(X->clone() + Y->clone()));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(SXY->hash_value() != SXYc->hash_value());
    REQUIRE_FALSE(SXYc->canon_transform().conj);
  }
}

TEST_CASE("conjugated_scalar_leaves", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Variable x{L"x"};
  Variable xs = x;
  xs.conjugate();
  REQUIRE(EvalExpr{x}.hash_value() == EvalExpr{xs}.hash_value());
  REQUIRE(EvalExpr{xs}.canon_transform().conj);
  REQUIRE(EvalExpr{x}.canon_transform().trivial());

  // conj(b^n) = conj(b)^n for integer n: the marker rides the transform
  auto pw = Power{ex<Variable>(L"y"), 2};
  auto pws = pw;
  pws.conjugate();
  REQUIRE(EvalExpr{pw}.hash_value() == EvalExpr{pws}.hash_value());
  REQUIRE(EvalExpr{pws}.canon_transform().conj);
  REQUIRE(EvalExpr{pw}.canon_transform().trivial());
}

// these cases drive binarize(ExprPtr) on purpose: the head layout they
// assert is exactly the one that overload derives from the source factors
SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
TEST_CASE("re_im_eval_nodes", "[EvalExpr][re-im]") {
  using namespace sequant;
  auto ctx_resetter = set_scoped_default_context(
      Context(get_default_context())
          .set_tensor_canonicalizer(
              L"", std::make_shared<DefaultTensorCanonicalizer>()));

  // The fold emission shape: Constant(2) * RealPart(s), s a
  // closed-contraction scalar network.
  auto g = deserialize(L"g{i_1;a_1}:N");
  auto t = deserialize(L"t{a_1;i_1}:N");
  auto s_expr = g->clone() * t->clone();

  auto two_re = ex<Constant>(2) * real_part(s_expr->clone());
  auto node = binarize(two_re);

  // locate the RealPart unary node in the scalar-wrapped root
  auto const& re = node.left()->op_type() ? node.left() : node.right();
  REQUIRE(re->op_type() == EvalOp::RealPart);
  REQUIRE(re->result_type() == ResultType::Scalar);
  REQUIRE_FALSE(re->is_primary());
  // Constant{1} sentinel right child; inner subtree on the left
  REQUIRE(re.right()->is_constant());
  REQUIRE(re.left()->is_product());

  // the inner subtree occupies the _same_ slot as an independent binarize of s
  auto inner_alone = binarize(s_expr->clone());
  REQUIRE(re.left()->hash_value() == inner_alone->hash_value());

  // Re, Im, and the bare inner all hash to distinct slots
  auto two_im = ex<Constant>(2) * imaginary_part(s_expr->clone());
  auto node_im = binarize(two_im);
  auto const& im = node_im.left()->op_type() ? node_im.left() : node_im.right();
  REQUIRE(im->op_type() == EvalOp::ImagPart);
  REQUIRE(re->hash_value() != im->hash_value());
  REQUIRE(re->hash_value() != inner_alone->hash_value());
  REQUIRE(im->hash_value() != inner_alone->hash_value());
}

TEST_CASE("tot_leaf_canonical_spelling_and_phase", "[eval_expr][tot][phase]") {
  // A tensor-of-tensors leaf is block-canonicalized in place like a flat
  // leaf: the stored spelling (what a leaf provider serves) and
  // canon_indices() carry the canonical slot order, the antisymmetric
  // reorder phase is the retrieval transform's byproduct, and the two
  // spellings share one slot. (HSeOH PNS sign regression, 2026-09-02: the
  // phase was recorded against an un-reordered spelling.)
  using namespace sequant;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  const Index i1(L"i↑_1"), i2(L"i↓_1");
  const Index a2(L"a↑_2", {i1, i2}), a3(L"a↑_3", {i1, i2});
  auto mk = [&](Index const& x, Index const& y) {
    return Tensor(L"t", bra{x, y}, ket{i1, i2}, Symmetry::Antisymm,
                  BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  EvalExpr e23{mk(a2, a3)}, e32{mk(a3, a2)};
  // both store the canonical spelling ...
  REQUIRE(e23.as_tensor().bra()[0].label() == e32.as_tensor().bra()[0].label());
  auto const& ci = e32.canon_indices();
  auto pos = [&](std::wstring_view l) {
    return std::find_if(ci.begin(), ci.end(),
                        [&](Index const& ix) { return ix.label() == l; }) -
           ci.begin();
  };
  REQUIRE(pos(e32.as_tensor().bra()[0].label()) <
          pos(e32.as_tensor().bra()[1].label()));
  // ... share one slot, and differ by the reorder phase only
  REQUIRE(e23.hash_value() == e32.hash_value());
  REQUIRE(e23.canon_transform().phase * e32.canon_transform().phase == -1);
  REQUIRE_FALSE(e23.canon_transform().conj);
  REQUIRE_FALSE(e32.canon_transform().conj);
  // a _mixed_-flavor bundle is stored space-major (up before down: the named
  // canonical order), whichever way it was written -- the raw canonical
  // vertex order sorts colors by their hash, which a provider cannot follow
  const Index b3(L"a↓_3", {i1, i2});
  EvalExpr eud{mk(a2, b3)}, edu{mk(b3, a2)};
  REQUIRE(eud.as_tensor().bra()[0].label() == L"a↑_2");
  REQUIRE(edu.as_tensor().bra()[0].label() == L"a↑_2");
  REQUIRE(eud.hash_value() == edu.hash_value());
  REQUIRE(eud.canon_transform().phase * edu.canon_transform().phase == -1);
}

TEST_CASE("leaf_reorder_phase_hoists_into_parents", "[eval_expr][tot][phase]") {
  // Products and sums that differ only by the slot order of an antisymmetric
  // ToT leaf: the leaf's reorder parity is a child _transform_ (its hash is
  // phase-blind), so the parents occupy _one_ slot. A product carries the
  // parity in its own transform (phases hoist multiplicatively); a sum
  // hoists a uniform parity and salts a mixed one -- otherwise one slot
  // would hold sign-different values for the two spellings.
  using namespace sequant;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  const Index i1(L"i↑_1"), i2(L"i↓_1");
  const Index a2(L"a↑_2", {i1, i2}), a3(L"a↑_3", {i1, i2});
  auto t = [&](Index const& x, Index const& y) {
    return ex<Tensor>(L"t", bra{x, y}, ket{i1, i2}, Symmetry::Antisymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  auto g = [&](std::wstring_view lbl) {
    return ex<Tensor>(lbl, bra{i1, i2}, ket{a2, a3}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  auto prod = [&](std::wstring_view lbl, Index const& x, Index const& y) {
    return ex<Product>(ExprPtrList{g(lbl), t(x, y)});
  };

  auto p23 = binarize(prod(L"g", a2, a3));
  auto p32 = binarize(prod(L"g", a3, a2));
  // the leaves share a slot and differ by the reorder parity ...
  REQUIRE(p23.right()->hash_value() == p32.right()->hash_value());
  REQUIRE(p23.right()->canon_phase() * p32.right()->canon_phase() == -1);
  // ... and so do the products, the parity hoisted into their transforms
  REQUIRE(p23->hash_value() == p32->hash_value());
  REQUIRE(p23->canon_phase() * p32->canon_phase() == -1);

  auto sum = [&](ExprPtr a, ExprPtr b) {
    return binarize(ex<Sum>(ExprPtrList{std::move(a), std::move(b)}));
  };
  auto s23 = sum(prod(L"g", a2, a3), prod(L"h", a2, a3));
  auto s32 = sum(prod(L"g", a3, a2), prod(L"h", a3, a2));
  auto s_mixed = sum(prod(L"g", a2, a3), prod(L"h", a3, a2));
  // uniform parity: one slot, the parity hoisted
  REQUIRE(s23->hash_value() == s32->hash_value());
  REQUIRE(s23->canon_phase() * s32->canon_phase() == -1);
  // mixed parity: not a whole-node transform -> its own slot
  REQUIRE(s_mixed->hash_value() != s23->hash_value());
  REQUIRE(s_mixed->canon_phase() == 1);
}

TEST_CASE("sum_slot_identity_covers_every_summand", "[eval_expr][sum]") {
  // a sum's slot depends on _all_ its summands: A + B and A + C are different
  // values; A + B and B + A are one value (the summand hash multiset is
  // order-blind) unless the leading summand's _layout_ differs
  using namespace sequant;
  auto sum_of = [](std::wstring_view a, std::wstring_view b) {
    return binarize(ex<Sum>(ExprPtrList{deserialize(a), deserialize(b)}));
  };
  auto const ab = sum_of(L"f{i_1;a_1}", L"g{i_1;a_1}");
  auto const ac = sum_of(L"f{i_1;a_1}", L"h{i_1;a_1}");
  auto const ba = sum_of(L"g{i_1;a_1}", L"f{i_1;a_1}");
  auto const ab2 = sum_of(L"f{i_2;a_2}", L"g{i_2;a_2}");
  REQUIRE(ab->hash_value() != ac->hash_value());
  // same summands, same layout, another order: one slot; a relabeled copy of
  // the same sum shares it too
  REQUIRE(ab->hash_value() == ba->hash_value());
  REQUIRE(ab->hash_value() == ab2->hash_value());
  // the result layout is the _first_ summand's: the same summands in another
  // order with a different leading layout are a different slot (the cached
  // array would be served in the wrong mode order otherwise)
  // (NonHermitian f keeps its written orientation, so the two leaves have
  // different layouts)
  auto const nonherm = [](std::wstring_view spec) {
    return deserialize<ExprPtr>(spec,
                                {.def_braket_symm = Hermiticity::NonHermitian});
  };
  auto const tf = binarize(
      ex<Sum>(ExprPtrList{nonherm(L"t{a_1;i_1}"), nonherm(L"f{i_1;a_1}")}));
  auto const ft = binarize(
      ex<Sum>(ExprPtrList{nonherm(L"f{i_1;a_1}"), nonherm(L"t{a_1;i_1}")}));
  REQUIRE(tf->canon_indices().front() != ft->canon_indices().front());
  REQUIRE(tf->hash_value() != ft->hash_value());
  // three summands: the _last_ one must count too
  auto const abc = binarize(ex<Sum>(ExprPtrList{deserialize(L"f{i_1;a_1}"),
                                                deserialize(L"g{i_1;a_1}"),
                                                deserialize(L"h{i_1;a_1}")}));
  auto const abd = binarize(ex<Sum>(ExprPtrList{deserialize(L"f{i_1;a_1}"),
                                                deserialize(L"g{i_1;a_1}"),
                                                deserialize(L"k{i_1;a_1}")}));
  REQUIRE(abc->hash_value() != abd->hash_value());
}

SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
TEST_CASE("denoted_expr_is_the_parent_network_spelling",
          "[eval_expr][denoted]") {
  // the denoted spelling is the leaf decoder's inverse: the {conj, swap}
  // channel is re-materialized as the '⁺' state of the stored array, a bare
  // {conj} as its '꙳'. The parent network sees the spelling as written, and
  // its states color the graph.
  using namespace sequant;

  // '⁺' over a NonHermitian core: stored bare + {conj, swap}; denoted is the
  // adjoint of the stored spelling
  Tensor w(L"F", bra{L"i_2"}, ket{L"a_1"});  // NonHermitian by default
  REQUIRE(w.adjoint() == 1);
  REQUIRE(w.adjointed());
  EvalExpr e{w};
  REQUIRE(e.canon_transform().conj);
  REQUIRE(e.canon_transform().braket_swap);
  REQUIRE_FALSE(e.as_tensor().adjointed());
  REQUIRE(e.denoted_expr()->as<Tensor>() == w);

  // a real-basis '꙳' kept by an indefinite parity: stored bare + {conj};
  // denoted is the stored spelling K-conjugated, slots in place
  Tensor s(L"F", bra{idx(L"i_1", Field::Real)}, ket{idx(L"a_1", Field::Real)},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  REQUIRE(s.kconjugate() == 1);
  REQUIRE(s.kconjugated());
  EvalExpr es{s};
  REQUIRE(es.canon_transform().conj);
  REQUIRE_FALSE(es.canon_transform().braket_swap);
  REQUIRE(es.denoted_expr()->as<Tensor>() == s);

  // a '꙳' over a complex basis is an array of its own: the state stays on
  // the stored spelling and there is nothing to re-materialize
  Tensor c(L"F", bra{L"i_1"}, ket{L"a_1"},
           TensorSymmetries{.conjugation_parity = ConjugationParity::None});
  REQUIRE(c.kconjugate() == 1);
  EvalExpr ec{c};
  REQUIRE(ec.canon_transform().trivial());
  REQUIRE(ec.denoted_expr()->as<Tensor>() == c);

  // an internal node inherits its tensor child's transform, but its
  // placeholder is built from the child's denoted spelling already, so
  // denoted_expr() hands it back untouched
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto const scaled = binarize(ex<Variable>(L"x") * ex<Tensor>(w));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(scaled->op_type().has_value());
  REQUIRE(scaled->is_tensor());
  REQUIRE(scaled->canon_transform().conj);
  REQUIRE(scaled->canon_transform().braket_swap);
  REQUIRE(scaled->denoted_expr()->as<Tensor>() == scaled->as_tensor());
}

TEST_CASE("tot_leaf_annotation_is_slot_faithful", "[eval_expr][tot]") {
  // A ToT leaf's canon_indices() is the layout the provider serves for the
  // _stored_ spelling: pure proto indices first (the pair labels of
  // C{mu;a<ij>}, laid out i,j,mu), then the plain slots in slot order, then
  // the proto-carrying slots. The occupied kets of t{a<ij>,b<ij>;j,i} are
  // both plain slots and proto indices; they take their slot position, since
  // the array's outer mode k pairs with inner mode k as a column -- read
  // i,j;a,b (the proto order) the same array would serve t^{ab}_{ij} for
  // t^{ab}_{ji} (nonrel PNO CSV-CCD exchange energy, 2026-09-11).
  using namespace sequant;
  const Index i1(L"i_1"), i2(L"i_2"), mu(L"a_1");
  const Index a(L"a_2", {i1, i2}), b(L"a_4", {i1, i2});
  auto t = [&](Index const& k0, Index const& k1) {
    return Tensor(L"t", bra{a, b}, ket{k0, k1}, Symmetry::Nonsymm,
                  BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
  };
  auto labels = [](EvalExpr const& e) {
    std::vector<std::wstring> out;
    for (auto const& ix : e.canon_indices()) out.emplace_back(ix.label());
    return out;
  };
  using L = std::vector<std::wstring>;
  REQUIRE(labels(EvalExpr{t(i1, i2)}) == L{L"i_1", L"i_2", L"a_2", L"a_4"});
  REQUIRE(labels(EvalExpr{t(i2, i1)}) == L{L"i_2", L"i_1", L"a_2", L"a_4"});
  REQUIRE(EvalExpr{t(i2, i1)}.indices_annot() == "i_2,i_1;a_2i_1i_2,a_4i_1i_2");
  // the two ket orders are different values (t^{ab}_{ji} vs t^{ab}_{ij});
  // their slot identity is label-blind (one provider array serves both), so
  // it is the annotation alone that tells the array's modes apart
  REQUIRE(EvalExpr{t(i1, i2)}.indices_annot() !=
          EvalExpr{t(i2, i1)}.indices_annot());
  // a CSV coefficient in either orientation: pair labels, expansion slot,
  // then the CSV slot
  Tensor const c_ket(L"C", bra{mu}, ket{a}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm);
  Tensor const c_bra(L"C", bra{a}, ket{mu}, Symmetry::Nonsymm,
                     BraKetSymmetry::Conjugate, ColumnSymmetry::Nonsymm);
  REQUIRE(labels(EvalExpr{c_ket}) == L{L"i_1", L"i_2", L"a_1", L"a_2"});
  REQUIRE(labels(EvalExpr{c_bra}) == L{L"i_1", L"i_2", L"a_1", L"a_2"});
  REQUIRE(EvalExpr{c_ket}.indices_annot() == "i_1,i_2,a_1;a_2i_1i_2");
}

// The symbolic spelling of an eval tree is the expression it denotes: a leaf
// stores the bare array with its states and its canonicalization sign on the
// transform, and to_expr puts both back.
TEST_CASE("to_expr_denotes_the_leaf", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));

  SECTION("an adjointed leaf") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    // the leaf stores the bare array under the adjoint channel ...
    REQUIRE_FALSE(node->as_tensor().adjointed());
    REQUIRE(node->canon_transform() ==
            CanonTransform{.conj = true, .braket_swap = true});
    // ... and denotes the spelling as written
    REQUIRE(*to_expr(node) == t);
  }

  SECTION("a K-conjugated leaf over a real basis") {
    auto rx = [](std::wstring_view l) { return idx(l, Field::Real); };
    auto r = ex<Tensor>(
        L"r", bra{rx(L"a_1")}, ket{rx(L"i_1")},
        TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    auto rc = conjugate(r);
    REQUIRE(rc->as<Tensor>().kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(rc);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE_FALSE(node->as_tensor().kconjugated());
    REQUIRE(node->canon_transform() == CanonTransform{.conj = true});
    REQUIRE(*to_expr(node) == *rc);
  }

  SECTION("a conjugated scalar leaf") {
    auto x = ex<Variable>(L"x");
    x->as<Variable>().conjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto node = binarize(x);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE_FALSE(node->as_variable().conjugated());
    REQUIRE(*to_expr(node) == *x);
  }

  SECTION("a leaf reordered with a sign") {
    // the block canonicalizer stores the mixed-space bra as i_3,a_1 with the
    // antisymmetric reorder's sign on the transform
    auto t = ex<Tensor>(L"t", bra{Index{L"a_1"}, Index{L"i_3"}},
                        ket{Index{L"i_1"}, Index{L"i_2"}}, Symmetry::Antisymm,
                        BraKetSymmetry::Nonsymm, ColumnSymmetry::Symm);
    auto v = ex<Tensor>(L"v", bra{L"i_1", L"i_2"}, ket{L"a_1", L"i_3"});
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto leaf = binarize(t);
    auto prod = binarize(t * v);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(leaf->canon_phase() == -1);
    auto const stored = leaf->as_tensor();
    REQUIRE(stored.bra()[0].label() == L"i_3");

    auto e = to_expr(leaf);
    REQUIRE(e->is<Product>());
    REQUIRE(e->as<Product>().scalar() == -1);
    REQUIRE(e->as<Product>().size() == 1);
    REQUIRE(*e->as<Product>().factor(0) == stored);

    auto pe = to_expr(prod);
    REQUIRE(pe->is<Product>());
    REQUIRE(pe->as<Product>().scalar() == -1);
    REQUIRE(pe->as<Product>().size() == 2);
    REQUIRE(*pe->as<Product>().factor(0) == stored);
    REQUIRE(*pe->as<Product>().factor(1) == *v);
  }
}

TEST_CASE("kramers_blind_erasure_helpers", "[eval_expr][kramers-blind]") {
  using namespace sequant;
  using namespace sequant::eval;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  auto parse = [](std::wstring_view s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Nonsymm,
                           .def_braket_symm = Hermiticity::NonHermitian});
  };
  // C(i↓_1,i↓_2,a_1; a↑_1<i↓_1 i↓_2>): outer pair slots blind, inner not
  // (↑ is the representative flavour, so ↓ labels are what gets erased)
  auto C = parse(L"C{i↓_1,i↓_2,a_1;a↑_1<i↓_1,i↓_2>}")->as<Tensor>();
  auto g = parse(L"g{i↓_1,a_2;a_1}")->as<Tensor>();
  KramersBlindness kb{.blind_slot =
                          [](Tensor const& t, std::size_t slot) {
                            return t.label() == L"C" && slot < 2;
                          },
                      .erase_space =
                          [](IndexSpace const& s) {
                            if (mbpt::to_spin(s.qns()) == mbpt::Spin::any)
                              return s;
                            return mbpt::make_spinalpha(Index(s, 1)).space();
                          }};
  SECTION("erasable set of a blind leaf") {
    auto er = erasable_indices(std::array{ExprPtr(ex<Tensor>(C))}, kb);
    REQUIRE(er.size() == 2);
    REQUIRE(er.count(Index(L"i↓_1")) == 1);
    REQUIRE(er.count(Index(L"i↓_2")) == 1);
  }
  SECTION("a non-blind occurrence pins the index") {
    auto er = erasable_indices(
        std::array{ExprPtr(ex<Tensor>(C)), ExprPtr(ex<Tensor>(g))}, kb);
    REQUIRE(er.size() == 1);
    REQUIRE(er.count(Index(L"i↓_2")) == 1);
  }
  SECTION("erased clone: spaces spin-free, protos rewritten, inner kept") {
    auto m = erasure_map(std::array{ExprPtr(ex<Tensor>(C))}, kb);
    REQUIRE(m.size() == 2);
    auto Ce = erase_indices(C, m);
    auto slots = Ce.const_slots() | ranges::to_vector;
    REQUIRE(slots[0].space() == Index(L"i↑_1").space());
    // placeholders are fresh temporaries minted in first-occurrence order
    REQUIRE(*slots[0].ordinal() >= Index::min_tmp_index());
    REQUIRE(*slots[1].ordinal() > *slots[0].ordinal());
    REQUIRE(slots[2].space() == Index(L"a_1").space());
    auto const& inner = slots[3];
    REQUIRE(inner.space() == Index(L"a↑_1").space());  // component kept
    REQUIRE(inner.proto_indices().size() == 2);
    REQUIRE(inner.proto_indices()[0].space() == Index(L"i↑_1").space());
    // the original is untouched
    REQUIRE((C.const_slots() | ranges::to_vector)[0].space() ==
            Index(L"i↓_1").space());
  }
  SECTION("empty hook erases nothing") {
    auto er = erasable_indices(std::array{ExprPtr(ex<Tensor>(C))},
                               KramersBlindness{});
    REQUIRE(er.empty());
  }
  SECTION("a blind composite slot nominates its protos (CSV spelling)") {
    // the CSV transform spells the projector C{a~; a<ij>}: the pair labels
    // occur only as protos of the PNS composite
    auto Cp = parse(L"C{a_1;a↑_1<i↓_1,i↓_2>}")->as<Tensor>();
    KramersBlindness kbp{
        .blind_slot =
            [](Tensor const& t, std::size_t slot) {
              auto const& ix = *(t.const_slots().begin() + slot);
              return t.label() == L"C" && ix.has_proto_indices();
            },
        .erase_space = kb.erase_space};
    auto er = erasable_indices(std::array{ExprPtr(ex<Tensor>(Cp))}, kbp);
    REQUIRE(er.size() == 2);
    auto m = erasure_map(std::array{ExprPtr(ex<Tensor>(Cp))}, kbp);
    REQUIRE(m.size() == 2);
    auto Ce = erase_indices(Cp, m);
    auto slots = Ce.const_slots() | ranges::to_vector;
    REQUIRE(slots[1].proto_indices()[0].space() == Index(L"i↑_1").space());
    REQUIRE(*slots[1].proto_indices()[0].ordinal() >= Index::min_tmp_index());
    REQUIRE(*slots[1].proto_indices()[1].ordinal() >
            *slots[1].proto_indices()[0].ordinal());
    // a non-blind composite (an amplitude) pins the protos
    auto t = parse(L"t{a↑_1<i↓_1,i↓_2>;a_2}")->as<Tensor>();
    auto er2 = erasable_indices(
        std::array{ExprPtr(ex<Tensor>(Cp)), ExprPtr(ex<Tensor>(t))}, kbp);
    REQUIRE(er2.empty());
    // a non-leaf factor is neutral
    auto er3 = erasable_indices(
        std::array{ExprPtr(ex<Tensor>(Cp)), ExprPtr(ex<Tensor>(t))}, kbp,
        container::svector<bool>{true, false});
    REQUIRE(er3.size() == 2);
  }
}

TEST_CASE("kramers_blind_leaf_identity", "[eval_expr][kramers-blind]") {
  using namespace sequant;
  using namespace sequant::eval;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  auto parse = [](std::wstring_view s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Nonsymm,
                           .def_braket_symm = Hermiticity::NonHermitian});
  };
  KramersBlindness kb{.blind_slot =
                          [](Tensor const& t, std::size_t slot) {
                            return t.label() == L"C" && slot < 2;
                          },
                      .erase_space =
                          [](IndexSpace const& s) {
                            if (mbpt::to_spin(s.qns()) == mbpt::Spin::any)
                              return s;
                            return mbpt::make_spinalpha(Index(s, 1)).space();
                          }};
  auto leaf = [&](std::wstring_view s, KramersBlindness const* k) {
    return EvalExpr(parse(s)->as<Tensor>(), k);
  };
  auto Cuu = leaf(L"C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", &kb);
  auto Cud = leaf(L"C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>}", &kb);
  auto Cdd = leaf(L"C{i↓_1,i↓_2,a_1;a↑_1<i↓_1,i↓_2>}", &kb);
  auto Cuu_dn = leaf(L"C{i↑_1,i↑_2,a_1;a↓_1<i↑_1,i↑_2>}", &kb);
  using Node = FullBinaryNode<EvalExpr>;
  TreeNodeEqualityComparator<Node> eq;
  // pair flavours are one identity
  REQUIRE(Cuu.hash_value() == Cud.hash_value());
  REQUIRE(Cuu.hash_value() == Cdd.hash_value());
  REQUIRE(eq(Node{Cuu}, Node{Cud}));
  REQUIRE(eq(Node{Cuu}, Node{Cdd}));
  // the PNS component stays value-distinctive
  REQUIRE(Cuu.hash_value() != Cuu_dn.hash_value());
  REQUIRE(!eq(Node{Cuu}, Node{Cuu_dn}));
  // labels / spelling untouched: the as-written flavours survive
  auto has_label = [](EvalExpr const& e, std::wstring_view lbl) {
    return ranges::any_of(
               e.canon_indices(),
               [&](Index const& i) { return i.full_label() == lbl; }) &&
           ranges::any_of(
               e.expr()->as<Tensor>().const_slots(),
               [&](Index const& i) { return i.full_label() == lbl; });
  };
  REQUIRE(has_label(Cud, L"i↓_1"));
  REQUIRE(has_label(Cud, L"i↑_2"));
  // without the hook nothing changes, and identities equal the hook-less ctor
  auto Cuu0 = leaf(L"C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", nullptr);
  auto Cud0 = leaf(L"C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>}", nullptr);
  REQUIRE(Cuu0.hash_value() != Cud0.hash_value());
  REQUIRE(Cuu0.hash_value() ==
          EvalExpr(parse(L"C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}")->as<Tensor>())
              .hash_value());
  // an inert hook (no blind slot) is the hook-less identity
  KramersBlindness inert{
      .blind_slot = [](Tensor const&, std::size_t) { return false; },
      .erase_space = kb.erase_space};
  REQUIRE(leaf(L"C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", &inert).hash_value() ==
          Cuu0.hash_value());
}

TEST_CASE("kramers_blind_product_identity", "[eval_expr][kramers-blind]") {
  using namespace sequant;
  using namespace sequant::eval;
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  auto parse = [](std::wstring_view s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Nonsymm,
                           .def_braket_symm = Hermiticity::NonHermitian});
  };
  BinarizationOptions opts{
      .kramers_blindness = {
          .blind_slot =
              [](Tensor const& t, std::size_t slot) {
                return t.label() == L"C" && slot < 2;
              },
          .erase_space =
              [](IndexSpace const& s) {
                if ((s.qns().to_int32() & mbpt::mask_v<mbpt::Spin>) == 0)
                  return s;
                if (mbpt::to_spin(s.qns()) == mbpt::Spin::any) return s;
                return mbpt::make_spinalpha(Index(s, 1)).space();
              }}};
  auto tree = [&](std::wstring_view head, std::wstring_view rhs,
                  BinarizationOptions const& o) {
    return binarize(ResultExpr{parse(head)->as<Tensor>(), parse(rhs)}, o);
  };
  using Node = FullBinaryNode<EvalExpr>;
  TreeNodeEqualityComparator<Node> eq;
  // PPL-like half projection (a_3 plays the aux index): one identity across
  // the pair flavours, spelled as the union residual blocks spell them
  auto Puu = tree(L"I{i↑_2,i↑_1,a_2,a_3;a↑_1<i↑_1,i↑_2>}",
                  L"g{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", opts);
  auto Pud = tree(L"I{i↑_2,i↓_1,a_2,a_3;a↑_2<i↑_2,i↓_1>}",
                  L"g{a_1,a_2,a_3} * C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>}", opts);
  auto Pdd = tree(L"I{i↓_2,i↓_1,a_2,a_3;a↑_1<i↓_1,i↓_2>}",
                  L"g{a_1,a_2,a_3} * C{i↓_1,i↓_2,a_1;a↑_1<i↓_1,i↓_2>}", opts);
  REQUIRE(Puu->hash_value() == Pud->hash_value());
  REQUIRE(Puu->hash_value() == Pdd->hash_value());
  REQUIRE(eq(Puu, Pud));
  REQUIRE(eq(Puu, Pdd));
  // the result keeps its as-written labels
  REQUIRE(ranges::any_of(Pud->canon_indices(), [](Index const& i) {
    return i.full_label() == L"i↓_1";
  }));
  // a flavoured non-blind leaf pins the index: g(i↑,..) vs g(i↓,..) differ
  auto Quu = tree(L"I{i↑_2,i↑_1,a_3;a↑_1<i↑_1,i↑_2>}",
                  L"g{a_1,i↑_1,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", opts);
  auto Qdd = tree(L"I{i↓_2,i↓_1,a_3;a↑_1<i↓_1,i↓_2>}",
                  L"g{a_1,i↓_1,a_3} * C{i↓_1,i↓_2,a_1;a↑_1<i↓_1,i↓_2>}", opts);
  REQUIRE(Quu->hash_value() != Qdd->hash_value());
  REQUIRE(!eq(Quu, Qdd));
  // the PNS component stays distinctive
  auto Puu_dn =
      tree(L"I{i↑_2,i↑_1,a_2,a_3;a↓_1<i↑_1,i↑_2>}",
           L"g{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↓_1<i↑_1,i↑_2>}", opts);
  REQUIRE(Puu->hash_value() != Puu_dn->hash_value());
  // hook off: nothing shared, identical to a hook-less binarize
  auto Puu0 = tree(L"I{i↑_2,i↑_1,a_2,a_3;a↑_1<i↑_1,i↑_2>}",
                   L"g{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}", {});
  auto Pud0 = tree(L"I{i↑_2,i↓_1,a_2,a_3;a↑_2<i↑_2,i↓_1>}",
                   L"g{a_1,a_2,a_3} * C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>}", {});
  REQUIRE(Puu0->hash_value() != Pud0->hash_value());
  REQUIRE(
      Puu0->hash_value() ==
      binarize(ResultExpr{
                   parse(L"I{i↑_2,i↑_1,a_2,a_3;a↑_1<i↑_1,i↑_2>}")->as<Tensor>(),
                   parse(L"g{a_1,a_2,a_3} * "
                         L"C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}")})
          ->hash_value());
  // the CSV spelling C{a~; a<ij>} (pair labels only as protos): blind on
  // the composite slot; f pins one pair label through a plain slot
  BinarizationOptions optsp{
      .kramers_blindness = {
          .blind_slot =
              [](Tensor const& t, std::size_t slot) {
                auto const& ix = *(t.const_slots().begin() + slot);
                return t.label() == L"C" && ix.has_proto_indices();
              },
          .erase_space =
              [](IndexSpace const& s) {
                if ((s.qns().to_int32() & mbpt::mask_v<mbpt::Spin>) == 0)
                  return s;
                if (mbpt::to_spin(s.qns()) == mbpt::Spin::any) return s;
                return mbpt::make_spinalpha(Index(s, 1)).space();
              }}};
  auto Cuu_p = tree(L"I{a_2,a_3;a↑_1<i↑_1,i↑_2>}",
                    L"g{a_1,a_2,a_3} * C{a_1;a↑_1<i↑_1,i↑_2>}", optsp);
  auto Cud_p = tree(L"I{a_2,a_3;a↑_2<i↑_2,i↓_1>}",
                    L"g{a_1,a_2,a_3} * C{a_1;a↑_2<i↑_2,i↓_1>}", optsp);
  auto Cdd_p = tree(L"I{a_2,a_3;a↑_1<i↓_1,i↓_2>}",
                    L"g{a_1,a_2,a_3} * C{a_1;a↑_1<i↓_1,i↓_2>}", optsp);
  REQUIRE(Cuu_p->hash_value() == Cud_p->hash_value());
  REQUIRE(Cuu_p->hash_value() == Cdd_p->hash_value());
  REQUIRE(eq(Cuu_p, Cud_p));
  auto Fuu_p =
      tree(L"I{i↑_1,a_3;a↑_1<i↑_1,i↑_2>}",
           L"g{a_1,a_2,a_3} * C{a_1;a↑_1<i↑_1,i↑_2>} * f{a_2;i↑_1}", optsp);
  auto Fdd_p =
      tree(L"I{i↓_1,a_3;a↑_1<i↓_1,i↓_2>}",
           L"g{a_1,a_2,a_3} * C{a_1;a↑_1<i↓_1,i↓_2>} * f{a_2;i↓_1}", optsp);
  REQUIRE(Fuu_p->hash_value() != Fdd_p->hash_value());
  // a sum of two blind products is one identity across the pair flavours
  auto Suu = tree(L"I{i↑_2,i↑_1,a_2,a_3;a↑_1<i↑_1,i↑_2>}",
                  L"g{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>} + "
                  L"f{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}",
                  opts);
  auto Sud = tree(L"I{i↑_2,i↓_1,a_2,a_3;a↑_2<i↑_2,i↓_1>}",
                  L"g{a_1,a_2,a_3} * C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>} + "
                  L"f{a_1,a_2,a_3} * C{i↑_2,i↓_1,a_1;a↑_2<i↑_2,i↓_1>}",
                  opts);
  REQUIRE(Suu->hash_value() == Sud->hash_value());
  REQUIRE(eq(Suu, Sud));
}

TEST_CASE("kramers_blind_guards", "[eval_expr][kramers-blind]") {
  using namespace sequant;
  using namespace sequant::eval;
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  auto parse = [](std::wstring_view s) {
    return deserialize(s, {.def_perm_symm = Symmetry::Nonsymm,
                           .def_braket_symm = Hermiticity::NonHermitian});
  };
  auto erase = [](IndexSpace const& s) {
    if ((s.qns().to_int32() & mbpt::mask_v<mbpt::Spin>) == 0) return s;
    if (mbpt::to_spin(s.qns()) == mbpt::Spin::any) return s;
    return mbpt::make_spinalpha(Index(s, 1)).space();
  };
  auto C = parse(L"C{i↓_1,i↓_2,a↓_1;a↑_1<i↓_1,i↓_2>}")->as<Tensor>();
  SECTION("a blind slot must be pure occupied") {
    KramersBlindness bad{
        .blind_slot = [](Tensor const&, std::size_t s) { return s == 2; },
        .erase_space = erase};
    REQUIRE_THROWS_AS(erasable_indices(std::array{ExprPtr(ex<Tensor>(C))}, bad),
                      std::invalid_argument);
  }
  SECTION("a blind composite slot must be indexed by pure-occupied protos") {
    auto D = parse(L"D{a_1;a↑_1<a↓_2>}")->as<Tensor>();
    KramersBlindness bad{
        .blind_slot = [](Tensor const&, std::size_t s) { return s == 1; },
        .erase_space = erase};
    REQUIRE_THROWS_AS(erasable_indices(std::array{ExprPtr(ex<Tensor>(D))}, bad),
                      std::invalid_argument);
    // with the pair labels also in plain slots, the plain slots govern: a
    // blind composite alone erases nothing here (its protos follow the
    // non-blind plain slots), and it is not a violation
    KramersBlindness ok{
        .blind_slot = [](Tensor const&, std::size_t s) { return s == 3; },
        .erase_space = erase};
    REQUIRE(erasable_indices(std::array{ExprPtr(ex<Tensor>(C))}, ok).empty());
  }
  SECTION("an index blind in one slot and not in another of one leaf") {
    auto D = parse(L"D{i↓_1,i↓_2;i↓_1}")->as<Tensor>();
    KramersBlindness bad{
        .blind_slot = [](Tensor const&, std::size_t s) { return s < 2; },
        .erase_space = erase};
    REQUIRE_THROWS_AS(erasable_indices(std::array{ExprPtr(ex<Tensor>(D))}, bad),
                      std::invalid_argument);
  }
  SECTION("inert hook == hook-less identity, leaf and tree") {
    BinarizationOptions inert{
        .kramers_blindness = {
            .blind_slot = [](Tensor const&, std::size_t) { return false; },
            .erase_space = erase}};
    auto e = parse(L"g{a_1,a_2,a_3} * C{i↑_1,i↑_2,a_1;a↑_1<i↑_1,i↑_2>}");
    auto head = parse(L"I{i↑_2,i↑_1,a_2,a_3;a↑_1<i↑_1,i↑_2>}")->as<Tensor>();
    auto t0 = binarize(ResultExpr{head, e}, {});
    auto t1 = binarize(ResultExpr{head, e}, inert);
    REQUIRE(t0->hash_value() == t1->hash_value());
    REQUIRE(!t1->has_identity_erasure());
    REQUIRE(t1->identity_tensor() == nullptr);
    using Node = FullBinaryNode<EvalExpr>;
    REQUIRE(TreeNodeEqualityComparator<Node>{}(t0, t1));
  }
}

#include <SeQuant/core/eval/eval_node.hpp>

TEST_CASE("kramers_flip_node", "[eval_expr][kramers-flip]") {
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  // EvalOp::KramersFlip: a unary wrapper (RealPart pattern: left = the
  // canonical ↑ node, right = the Constant(1) sentinel) denoting
  // phase · F(inner) with F the time-reversal flip over the given union modes
  using namespace sequant;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  const Index i1(L"i↑_1"), i2(L"i↑_2");
  const Index au(L"a↑_1", {i1, i2}), ad(L"a↓_1", {i1, i2});
  const Index a1(L"a_1"), a2(L"a_2"), a3(L"a_3");
  auto mk = [](std::wstring_view lbl, std::initializer_list<Index> b,
               std::initializer_list<Index> k) {
    return Tensor(lbl, bra(b), ket(k), Symmetry::Nonsymm,
                  BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm,
                  KramersSymmetry::TimeReversal);
  };
  // I{a_2,a_3;a↑_1<i↑_1,i↑_2>} = g{a_1;a_2,a_3} * C{a_1;a↑_1<i↑_1,i↑_2>}
  auto inner_expr = [&] {
    return ResultExpr{
        mk(L"I", {a2, a3}, {au}),
        ex<Product>(ExprPtrList{ex<Tensor>(mk(L"g", {a1}, {a2, a3})),
                                ex<Tensor>(mk(L"C", {a1}, {au}))})};
  };
  auto inner = binarize(inner_expr());
  Tensor const denoted = mk(L"I", {a2, a3}, {ad});
  using modes_t = container::svector<std::size_t>;

  auto wrap = make_kramers_flip_node(inner, modes_t{0, 1}, -1, denoted);
  REQUIRE(wrap->op_type() == EvalOp::KramersFlip);
  REQUIRE(wrap->result_type() == ResultType::Tensor);
  REQUIRE(wrap->kramers_flip_modes() == modes_t{0, 1});
  REQUIRE(wrap->kramers_flip_phase() == -1);
  REQUIRE(wrap.left()->hash_value() == inner->hash_value());
  REQUIRE(wrap.right()->is_constant());
  REQUIRE(wrap->hash_value() != inner->hash_value());
  // the wrapper denotes the flipped spelling with the child's layout
  REQUIRE(wrap->is_tensor());
  REQUIRE(wrap->as_tensor().label() == L"I");
  REQUIRE(wrap->canon_indices().size() == inner->canon_indices().size());
  REQUIRE(
      std::any_of(wrap->canon_indices().begin(), wrap->canon_indices().end(),
                  [&](Index const& ix) { return ix.space() == ad.space(); }));
  REQUIRE(
      std::none_of(wrap->canon_indices().begin(), wrap->canon_indices().end(),
                   [&](Index const& ix) { return ix.space() == au.space(); }));

  // same child, modes and phase => the same slot
  auto wrap2 = make_kramers_flip_node(binarize(inner_expr()), modes_t{0, 1}, -1,
                                      denoted);
  REQUIRE(wrap->hash_value() == wrap2->hash_value());
  using node_t = std::remove_cvref_t<decltype(wrap)>;
  TreeNodeEqualityComparator<node_t> same;
  REQUIRE(same(wrap, wrap2));

  // a different phase or mode set is a different value
  REQUIRE(
      make_kramers_flip_node(inner, modes_t{0, 1}, 1, denoted)->hash_value() !=
      wrap->hash_value());
  REQUIRE(
      make_kramers_flip_node(inner, modes_t{0}, -1, denoted)->hash_value() !=
      wrap->hash_value());

  // the linearized form spells the denoted (flipped) contraction
  auto lin = linearize_eval_node(wrap);
  REQUIRE(lin->is<Product>());
  bool has_down = false;
  lin->visit(
      [&](ExprPtr const& x) {
        if (!x->is<Tensor>()) return;
        for (auto const& ix : x->as<Tensor>().const_indices())
          if (ix.space() == ad.space()) has_down = true;
      },
      /*atoms_only=*/true);
  REQUIRE(has_down);
}

TEST_CASE("kramers_fold_intermediates", "[eval_expr][kramers-flip]") {
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  // BinarizationOptions::kramers_fold_intermediates: a down-majority
  // intermediate (all leaves time-reversal symmetric; the Kramers-blind pair
  // labels do not count) binarizes as a KramersFlip over its canonical (up)
  // partner's node, which the partner family builds: one contraction, one
  // O(size) flip
  using namespace sequant;
  using namespace sequant::eval;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  const Index i1(L"i↑_1"), i2(L"i↑_2");
  const Index au(L"a↑_1", {i1, i2}), ad(L"a↓_1", {i1, i2});
  const Index bu(L"a↑_2", {i1, i2}), bd(L"a↓_2", {i1, i2});
  const Index a1(L"a_1"), a2(L"a_2"), a3(L"a_3");
  auto mk = [](std::wstring_view lbl, std::initializer_list<Index> b,
               std::initializer_list<Index> k,
               KramersSymmetry ks = KramersSymmetry::TimeReversal) {
    return ex<Tensor>(lbl, bra(b), ket(k), Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm, ks);
  };
  // the pair labels are blind on the projector's composite slot (Phase 1)
  BinarizationOptions opts{
      .kramers_blindness = {
          .blind_slot =
              [](Tensor const& t, std::size_t slot) {
                return t.label() == L"C" && slot == 1;
              },
          .erase_space =
              [](IndexSpace const& s) {
                if (mbpt::to_spin(s.qns()) == mbpt::Spin::any) return s;
                return mbpt::make_spinalpha(Index(s, 1)).space();
              }},
      .kramers_fold_intermediates = true};
  using modes_t = container::svector<std::size_t>;

  // I{a_2,a_3;a<i↑_1,i↑_2>} = g{a_1;a_2,a_3} * C{a_1;a<i↑_1,i↑_2>}
  auto I = [&](Index const& comp) {
    return Tensor(L"I", bra{a2, a3}, ket{comp});
  };
  auto rhs = [&](Index const& comp,
                 KramersSymmetry ks = KramersSymmetry::TimeReversal) {
    return ex<Product>(
        ExprPtrList{mk(L"g", {a1}, {a2, a3}, ks), mk(L"C", {a1}, {comp}, ks)});
  };
  // a ResultExpr root keeps the head's spelling: it never folds as a whole
  // (its factors may); the fold decision is exercised through the internal
  // entry a root Product's factors and a Sum's summands' factors go through
  REQUIRE(binarize(ResultExpr{I(ad), rhs(ad)}, opts)->op_type() ==
          EvalOp::Product);
  auto bin = [](ExprPtr const& e, IndexSet const& ext,
                BinarizationOptions const& o) {
    std::size_t counter = 0;
    return impl::binarize(e, ext, o, counter);
  };
  auto up = bin(rhs(au), IndexSet{au, a2, a3}, opts);
  auto dn = bin(rhs(ad), IndexSet{ad, a2, a3}, opts);
  REQUIRE(up->op_type() == EvalOp::Product);
  REQUIRE(dn->op_type() == EvalOp::KramersFlip);
  REQUIRE(dn.left()->hash_value() == up->hash_value());  // one shared family
  REQUIRE(dn->kramers_flip_phase() == -1);               // one down composite
  // the union legs a_2, a_3 as OUTER mode positions: the pair labels
  // i↑_1, i↑_2 are outer (CSV pair) modes too, the composite is inner
  REQUIRE(dn->kramers_flip_modes() == modes_t{2, 3});
  REQUIRE(dn->as_tensor().label() == L"I");
  REQUIRE(
      std::any_of(dn->canon_indices().begin(), dn->canon_indices().end(),
                  [&](Index const& ix) { return ix.space() == ad.space(); }));

  // no fold without the option, or with a leaf that is not time-reversal
  // symmetric
  auto nofold = opts;
  nofold.kramers_fold_intermediates = false;
  REQUIRE(bin(rhs(ad), IndexSet{ad, a2, a3}, nofold)->op_type() ==
          EvalOp::Product);
  REQUIRE(bin(rhs(ad, KramersSymmetry::Nonsymm), IndexSet{ad, a2, a3}, opts)
              ->op_type() == EvalOp::Product);

  // a tie (one up, one down composite) resolves by the flavour string in the
  // flavour-blind canonical order of the externals: (down_1, up_2) folds onto
  // (up_1, down_2), never both ways
  auto rhs2 = [&](Index const& x, Index const& y) {
    return ex<Product>(ExprPtrList{mk(L"g", {a1, a2}, {a3}),
                                   mk(L"C", {a1}, {x}), mk(L"C", {a2}, {y})});
  };
  auto t_ud = bin(rhs2(au, bd), IndexSet{au, bd, a3}, opts);
  auto t_du = bin(rhs2(ad, bu), IndexSet{ad, bu, a3}, opts);
  REQUIRE(t_ud->op_type() == EvalOp::Product);
  REQUIRE(t_du->op_type() == EvalOp::KramersFlip);
  REQUIRE(t_du.left()->hash_value() == t_ud->hash_value());
  REQUIRE(t_du->kramers_flip_phase() == -1);
  REQUIRE(t_du->kramers_flip_modes() ==
          modes_t{2});  // a_3 after the pair modes

  // a sum of down-majority terms folds as one node onto the up sum
  auto sum_of = [&](Index const& comp) {
    return ex<Sum>(ExprPtrList{
        rhs(comp), ex<Product>(ExprPtrList{mk(L"h", {a1}, {a2, a3}),
                                           mk(L"C", {a1}, {comp})})});
  };
  auto s_up = bin(sum_of(au), IndexSet{au, a2, a3}, opts);
  auto s_dn = bin(sum_of(ad), IndexSet{ad, a2, a3}, opts);
  REQUIRE(s_up->op_type() == EvalOp::Sum);
  REQUIRE(s_dn->op_type() == EvalOp::KramersFlip);
  REQUIRE(s_dn.left()->hash_value() == s_up->hash_value());
  REQUIRE(s_dn.left()->op_type() == EvalOp::Sum);

  // the flip is DEEP: the flavoured pair labels inside an unflavoured
  // (Kramers-union) composite follow their plain-slot occurrences, so every
  // leaf of the canonical partner is spelled consistently (an amplitude leaf
  // t{a_2<i↓_1,i_2>; i↓_1, i_2} becomes t{a_2<i↑_1,i_2>; i↑_1, i_2})
  // (X, not the blind projector C: its bra composite would carry the pair
  // labels in a non-blind slot, which the blindness guard rejects)
  const Index i1d(L"i↓_1"), iu2(L"i_2");
  const Index au3(L"a_3", {i1d, iu2}), cd(L"a↓_1", {i1d, iu2});
  auto deep = bin(ex<Product>(ExprPtrList{mk(L"X", {au3}, {cd}),
                                          mk(L"t", {au3}, {i1d, iu2})}),
                  IndexSet{i1d, iu2, cd}, opts);
  REQUIRE(deep->op_type() == EvalOp::KramersFlip);
  REQUIRE(deep->kramers_flip_phase() == 1);  // two down externals
  auto no_down = [&](auto const& node, auto& self) -> bool {
    if (node.leaf()) {
      if (!node->is_tensor()) return true;
      for (auto const& ix : node->as_tensor().const_indices()) {
        if (ix.space() == i1d.space()) return false;
        for (auto const& p : ix.proto_indices())
          if (p.space() == i1d.space()) return false;
      }
      return true;
    }
    return self(node.left(), self) && self(node.right(), self);
  };
  REQUIRE(no_down(deep.left(), no_down));

  // a Sum's direct summands keep the Sum's labels: a summand that would fold
  // on its own (its pair labels are blind inside it) does not when another
  // summand pins them and the Sum as a whole is canonical; nested factors
  // still may
  auto blind_occ = opts;
  blind_occ.kramers_blindness.blind_slot = [isr](Tensor const& t,
                                                 std::size_t slot) {
    if (t.label() != L"C") return false;
    auto const& ix = *(t.const_slots().begin() + slot);
    return ix.has_proto_indices() || isr->is_pure_occupied(ix.space());
  };
  const Index cdu(L"a↓_1", {i1, i2});  // down column, up pair labels
  auto p_pinned =
      ex<Product>(ExprPtrList{mk(L"g", {a3}, {i1, i2}), mk(L"C", {a3}, {cdu})});
  auto p_blind = ex<Product>(
      ExprPtrList{mk(L"h", {a3}, {}), mk(L"C", {i1, i2, a3}, {cdu})});
  auto head = Tensor(L"I", bra{i1, i2}, ket{cdu});
  // alone, the blind summand is down-majority (its pair labels do not count)
  REQUIRE(bin(p_blind, IndexSet{i1, i2, cdu}, blind_occ)->op_type() ==
          EvalOp::KramersFlip);
  // in the Sum the pinned pair labels make the whole canonical (1 down, 2 up)
  auto s_mixed = binarize(
      ResultExpr{head, ex<Sum>(ExprPtrList{p_pinned, p_blind})}, blind_occ);
  REQUIRE(s_mixed->op_type() == EvalOp::Sum);
  auto no_flip = [&](auto const& node, auto& self) -> bool {
    if (node->op_type() == EvalOp::KramersFlip) return false;
    if (node.leaf()) return true;
    return self(node.left(), self) && self(node.right(), self);
  };
  REQUIRE(no_flip(s_mixed, no_flip));

  // a product whose factor folds contracts THROUGH the factor's denoted
  // (flipped-flavour) tensor: the parent's network holds the wrapper's
  // tensor, its annotations are consistent with the wrapper's labels, and
  // the parent's result is the head
  auto bracket = ex<Product>(ExprPtrList{mk(L"g", {a1}, {a2, a3}),
                                         mk(L"C", {a1}, {ad})});  // I{a2,a3;a↓}
  const Index a4(L"a_4");
  auto outer = ex<Product>(ExprPtrList{
      bracket, mk(L"Y", {a2, a3}, {a4})});  // J{a4;a↓<i1,i2>} = I * Y
  auto pj = bin(outer, IndexSet{a4, ad}, opts);
  REQUIRE(pj->op_type() == EvalOp::KramersFlip);  // the whole (1 down) folds
  // ... but not as a ResultExpr root
  REQUIRE(binarize(ResultExpr{Tensor(L"J", bra{a4}, ket{ad}), outer}, opts)
              ->op_type() == EvalOp::Product);
  // the same product where the head is up: the bracket alone folds, the
  // outer product stays a Product whose left child is the wrapper
  auto bracket_u =
      ex<Product>(ExprPtrList{mk(L"g", {a1}, {a2, a3}), mk(L"C", {a1}, {ad})});
  auto outer_u = ex<Product>(
      Product{1, ExprPtrList{bracket_u, mk(L"X", {i1, i2, ad}, {a2, a3, au})},
              Product::Flatten::No});  // K{i1,i2;a↑}, the bracket kept
  auto pk =
      binarize(ResultExpr{Tensor(L"K", bra{i1, i2}, ket{au}), outer_u}, opts);
  REQUIRE(pk->op_type() == EvalOp::Product);
  bool saw_wrapper = false;
  for (auto const& ch : {pk.left(), pk.right()}) {
    if (ch->op_type() == EvalOp::KramersFlip) {
      saw_wrapper = true;
      // the wrapper's labels are the bracket's as written (a↓ column)
      REQUIRE(std::any_of(
          ch->canon_indices().begin(), ch->canon_indices().end(),
          [&](Index const& ix) { return ix.space() == ad.space(); }));
    }
  }
  REQUIRE(saw_wrapper);
  // every index of the product's children is either shared or in the result
  auto labels = [](auto const& n) {
    container::set<std::wstring> out;
    for (auto const& ix : n->canon_indices())
      out.insert(std::wstring(ix.full_label()));
    return out;
  };
  auto const L = labels(pk.left()), R = labels(pk.right()), T = labels(pk);
  for (auto const& l : R) REQUIRE((L.contains(l) || T.contains(l)));
  for (auto const& l : L) REQUIRE((R.contains(l) || T.contains(l)));
}

TEST_CASE("kramers_fold_batch_annotations",
          "[eval_expr][kramers-flip][batched-here]") {
  auto lax_braket = sequant::tests::scoped_lax_braket_pairing();
  // The optimizer's per-node batch annotations
  // (BinarizationOptions::node_batch_axes) name the term's own indices. A
  // folded family is binarized from its flipped spelling, so the entries its
  // nodes consume must be flipped with it; the KramersFlip wrapper -- the
  // family's root in the original spelling -- takes the root's External
  // slicing and loop opens in that spelling, while a Contracted open stays at
  // its reduction node inside the flip.
  using namespace sequant;
  using namespace sequant::eval;
  auto isr = mbpt::make_min_sr_spaces(mbpt::SpinConvention::None);
  mbpt::add_fermi_spin(*isr);
  Context ctx = get_default_context();
  ctx.set(isr);
  auto resetter = set_scoped_default_context(ctx);
  const Index i1(L"i↑_1"), i2(L"i↑_2");
  const Index au(L"a↑_1", {i1, i2}), ad(L"a↓_1", {i1, i2});
  const Index a1(L"a_1"), a2(L"a_2"), a3(L"a_3");
  auto mk = [](std::wstring_view lbl, std::initializer_list<Index> b,
               std::initializer_list<Index> k) {
    return ex<Tensor>(lbl, bra(b), ket(k), Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm,
                      KramersSymmetry::TimeReversal);
  };
  BinarizationOptions opts{
      .kramers_blindness = {
          .blind_slot =
              [](Tensor const& t, std::size_t slot) {
                return t.label() == L"C" && slot == 1;
              },
          .erase_space =
              [](IndexSpace const& s) {
                if (mbpt::to_spin(s.qns()) == mbpt::Spin::any) return s;
                return mbpt::make_spinalpha(Index(s, 1)).space();
              }},
      .kramers_fold_intermediates = true};
  // one DP node, I{a_2,a_3;a<i,i>} = g{a_1;a_2,a_3} * C{a_1;a<i,i>}: sliced
  // on the pair label i↑_1 (External, opened here) and on a_1 (Contracted:
  // its reduction is here)
  using axes_t = container::svector<std::pair<Index, BatchModeType>>;
  NodeBatchAnnotation ann;
  ann.axes =
      axes_t{{i1, BatchModeType::External}, {a1, BatchModeType::Contracted}};
  ann.opened_here = ann.axes;
  ann.effective_count = 3;
  ann.order_aware = true;
  opts.node_batch_axes = {ann};
  auto rhs = [&](Index const& comp) {
    return ex<Product>(
        ExprPtrList{mk(L"g", {a1}, {a2, a3}), mk(L"C", {a1}, {comp})});
  };
  auto bin = [](ExprPtr const& e, IndexSet const& ext,
                BinarizationOptions const& o) {
    std::size_t counter = 0;
    auto node = impl::binarize(e, ext, o, counter);
    REQUIRE(counter == o.node_batch_axes.size());
    return node;
  };
  // the canonical (up) family keeps its annotation as given
  auto up = bin(rhs(au), IndexSet{au, a2, a3}, opts);
  REQUIRE(up->op_type() == EvalOp::Product);
  CHECK(up->node_slice_mask() == ann.axes);
  CHECK(up->batch_loops_opened_here() == ann.opened_here);

  // the down family: a KramersFlip over the up family's node
  auto dn = bin(rhs(ad), IndexSet{ad, a2, a3}, opts);
  REQUIRE(dn->op_type() == EvalOp::KramersFlip);
  REQUIRE(dn.left()->hash_value() == up->hash_value());
  Index const i1_flipped = kramers_flipped_deep(i1, *isr);
  REQUIRE(i1_flipped != i1);
  // the inner node reads its annotation in the flipped spelling: sliced on
  // the flipped pair label, on the unflavoured a_1 as is; the External loop
  // is no longer opened here, the Contracted one still is
  CHECK(dn.left()->node_slice_mask() ==
        axes_t{{i1_flipped, BatchModeType::External},
               {a1, BatchModeType::Contracted}});
  CHECK(dn.left()->batch_loops_opened_here() ==
        axes_t{{a1, BatchModeType::Contracted}});
  CHECK(dn.left()->batch_effective_count() == 3);
  CHECK(dn.left()->batch_order_aware());
  // the wrapper is the family root in the original spelling: sliced on the
  // External axis, which it opens
  CHECK(dn->node_slice_mask() == axes_t{{i1, BatchModeType::External}});
  CHECK(dn->batch_loops_opened_here() == axes_t{{i1, BatchModeType::External}});
  CHECK(dn->batch_effective_count() == 3);
  CHECK(dn->batch_order_aware());
}

TEST_CASE("binarize-broadcast-sum", "[EvalExpr][broadcast]") {
  using namespace sequant;
  using Node = FullBinaryNode<EvalExpr>;

  // the summands of the Sum chain at root, in order
  auto summands_of = [](Node const& root) {
    container::svector<Node const*> result;
    Node const* n = &root;
    while (!n->leaf() && (*n)->op_type() == EvalOp::Sum) {
      result.push_back(&n->right());
      n = &n->left();
    }
    result.push_back(n);
    std::ranges::reverse(result);
    return result;
  };
  auto indices_of = [](ResultExpr const& res) {
    return res.indices() | ranges::to<EvalExpr::index_vector>;
  };
  // n is a Broadcast node over the head of res
  auto broadcasts_over = [&indices_of](Node const& n, ResultExpr const& res) {
    if (n.leaf() || !n->is_broadcast()) return false;
    REQUIRE(n->op_type() == EvalOp::Broadcast);
    REQUIRE(n.right()->is_constant());
    REQUIRE(n->canon_indices() == indices_of(res));
    auto const& t = n->as_tensor();
    auto const head = res.result_as_tensor();
    REQUIRE(t.symmetry() == Symmetry::Nonsymm);
    REQUIRE(t.braket_symmetry() == BraKetSymmetry::Nonsymm);
    return std::ranges::equal(t.bra(), head.bra()) &&
           std::ranges::equal(t.ket(), head.ket()) &&
           std::ranges::equal(t.aux(), head.aux());
  };

  const ResultExpr d4 = deserialize<ResultExpr>(
      L"D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}");

  SECTION("a Sum chain over Broadcast nodes of the summands") {
    auto const tree = binarize(d4);
    REQUIRE(tree->op_type() == EvalOp::Sum);
    REQUIRE(tree->as_tensor() == d4.result_as_tensor());
    auto const summands = summands_of(tree);
    REQUIRE(summands.size() == 4);
    auto const& eps = d4.expression()->as<Sum>().summands();
    for (std::size_t s = 0; s != 4; ++s) {
      CAPTURE(s);
      Node const& b = *summands[s];
      REQUIRE(broadcasts_over(b, d4));
      Node const& operand = b.left();
      if (s < 2) {
        // ε{;a_k}
        REQUIRE(operand.leaf());
        REQUIRE(*operand->expr() == *eps[s]);
      } else {
        // −1·ε{i_k;}
        REQUIRE(operand->op_type() == EvalOp::Product);
        REQUIRE(operand.left().leaf());
        REQUIRE(*operand.left()->expr() == *eps[s]->as<Product>().factor(0));
        REQUIRE(operand.right()->is_constant());
        REQUIRE(operand.right()->as_constant().value() == -1);
      }
      REQUIRE(b->canon_phase() == operand->canon_phase());
    }
  }

  SECTION("to_expr spells the broadcast summands as written") {
    auto const tree = binarize(d4);
    REQUIRE(*to_expr(tree) == *d4.expression());
  }

  SECTION("bra-ket-swapped supply counts as present") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"D{i_1;a_1} = ε{a_1;} - ε{i_1;}");
    auto const summands_tree = binarize(res);
    auto const summands = summands_of(summands_tree);
    REQUIRE(summands.size() == 2);
    for (Node const* b : summands) REQUIRE(broadcasts_over(*b, res));
    // ε{a_1;} lacks i_1 only, ε{i_1;} lacks a_1 only
    REQUIRE(summands[0]->left()->canon_indices() ==
            EvalExpr::index_vector{Index{L"a_1"}});
    REQUIRE(summands[1]->left()->canon_indices() ==
            EvalExpr::index_vector{Index{L"i_1"}});
  }

  SECTION("scalar summands") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"D{i_1;a_1} = ε{;a_1} - ε{i_1;} - E");
    auto const summands_tree = binarize(res);
    auto const summands = summands_of(summands_tree);
    REQUIRE(summands.size() == 3);
    Node const& e = *summands[2];
    REQUIRE(broadcasts_over(e, res));
    // −1·E, a scalar
    REQUIRE(e.left()->is_scalar());
    REQUIRE(e.left()->op_type() == EvalOp::Product);
    REQUIRE(e.left().left()->is_variable());
    REQUIRE(e.left().left()->as_variable() == Variable(L"E"));

    // a Constant summand
    const ResultExpr c = deserialize<ResultExpr>(L"R{a_1} = 2 + t{a_1}");
    auto const csummands_tree = binarize(c);
    auto const csummands = summands_of(csummands_tree);
    REQUIRE(csummands.size() == 2);
    auto const n_bcast = std::ranges::count_if(
        csummands, [&](Node const* n) { return broadcasts_over(*n, c); });
    REQUIRE(n_bcast == 1);
  }

  SECTION("a summand supplying every head index is not wrapped") {
    const ResultExpr complete = deserialize<ResultExpr>(
        L"R{i_1;a_1} = f{i_1;a_1} + g{i_1;i_2} t{i_2;a_1}");
    auto const complete_tree = binarize(complete);
    REQUIRE(std::ranges::none_of(summands_of(complete_tree), [](Node const* n) {
      return !n->leaf() && (*n)->is_broadcast();
    }));
    // a head index used twice (a Hadamard product) is supplied
    const ResultExpr hadamard = deserialize<ResultExpr>(
        L"R{i_1;a_1} = f{i_1;a_1} + M{i_1;a_1} D{a_1;i_1}");
    REQUIRE(is_valid(hadamard));
    auto const hadamard_tree = binarize(hadamard);
    REQUIRE(std::ranges::none_of(summands_of(hadamard_tree), [](Node const* n) {
      return !n->leaf() && (*n)->is_broadcast();
    }));
  }

  SECTION("a non-Sum rhs is not wrapped") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"R{i_1;a_1} = 2 g{i_1;i_2} t{i_2;a_1}");
    auto const tree = binarize(res);
    REQUIRE(tree->op_type() == EvalOp::Product);
    tree.visit(
        [](Node const& n) { REQUIRE((n.leaf() || !n->is_broadcast())); });
  }

  SECTION("a non-Sum right-hand side cannot broadcast") {
    for (std::wstring_view spec :
         {L"D{i_1;a_1} = ε{;a_1}", L"D{i_1;a_1} = 2 ε{;a_1}", L"D{i_1;a_1} = E",
          L"D{i_1;a_1} = ε{;a_1} ε{;a_1}"}) {
      CAPTURE(toUtf8(spec));
      const ResultExpr res = deserialize<ResultExpr>(spec);
      REQUIRE_FALSE(is_valid(res));
      REQUIRE_THROWS_WITH(
          binarize(res),
          Catch::Matchers::ContainsSubstring("lacks head indices i_1") &&
              Catch::Matchers::ContainsSubstring(
                  "broadcast is only defined for a Sum under a head"));
    }
  }

  SECTION("a broadcast summand carrying a non-head index") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"D{i_1;a_1} = ε{;a_1} + f{i_1,i_2;i_2}");
    REQUIRE_THROWS_WITH(
        binarize(res),
        Catch::Matchers::ContainsSubstring("carries index i_2") &&
            Catch::Matchers::ContainsSubstring(
                "that is not an index of the result"));
  }

  SECTION("a nested Sum broadcasts at the top level only") {
    const ResultExpr res = deserialize<ResultExpr>(
        L"R{i_1;a_1} = (f{i_1;} + g{i_1;}) c + h{i_1;a_1}");
    auto const summands_tree = binarize(res);
    auto const summands = summands_of(summands_tree);
    REQUIRE(summands.size() == 2);
    REQUIRE(broadcasts_over(*summands[0], res));
    REQUIRE(summands[0]->left()->canon_indices() ==
            EvalExpr::index_vector{Index{L"i_1"}});
    REQUIRE(summands[1]->leaf());
  }

  SECTION("a symmetrizer supplies its indices") {
    // the binarized summand carries the symmetrizer's indices, so a Broadcast
    // node would add no mode
    auto const no_broadcast = [](Node const& root) {
      bool found = false;
      root.visit([&](Node const& n) {
        if (!n.leaf() && n->is_broadcast()) found = true;
      });
      return !found;
    };
    const ResultExpr top = deserialize<ResultExpr>(
        L"R{a1,a2;i1,i2} = Â{i1,i2;a1,a2} t{a1;i1} ε{a2;} + g{a1,a2;i1,i2}");
    REQUIRE(no_broadcast(binarize(top)));
    const ResultExpr nested = deserialize<ResultExpr>(
        L"R{a1,a2;i1,i2} = c (Â{i1,i2;a1,a2} t{a1;i1} ε{a2;} + "
        L"Â{i1,i2;a1,a2} t{a2;i1} ε{a1;}) + g{a1,a2;i1,i2}");
    REQUIRE(no_broadcast(binarize(nested)));
  }

  SECTION("aux index cannot be broadcast") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"R{i_1;;a_2} = f{i_1;;a_2} + g{i_1;}");
    REQUIRE_FALSE(is_valid(res));
    REQUIRE_THROWS_AS(binarize(res), Exception);

    // a summand that carries the aux index may lack a bra/ket head index
    const ResultExpr ok =
        deserialize<ResultExpr>(L"R{i_1,i_2;;a_2} = f{i_1;;a_2} + g{i_2;;a_2}");
    REQUIRE(is_valid(ok));
    auto const summands_tree = binarize(ok);
    auto const summands = summands_of(summands_tree);
    REQUIRE(summands.size() == 2);
    for (Node const* b : summands) REQUIRE(broadcasts_over(*b, ok));
  }

  SECTION("hash: the operand and its placement in the result") {
    auto const tree = binarize(d4);
    auto const summands = summands_of(tree);
    // ε{;a_1} and ε{;a_2} are one value, broadcast to different modes
    REQUIRE(summands[0]->left()->hash_value() ==
            summands[1]->left()->hash_value());
    REQUIRE((*summands[0])->hash_value() != (*summands[1])->hash_value());

    // the same operand broadcast alike: one value
    auto const first_bcast = [&summands_of](ResultExpr const& res) {
      auto const tree = binarize(res);
      return (*summands_of(tree).front())->hash_value();
    };
    auto const d = deserialize<ResultExpr>(L"D{i_1;a_1} = ε{;a_1} - ε{i_1;}");
    auto const f =
        deserialize<ResultExpr>(L"F{i_1;a_1} = ε{;a_1} + t{i_1;a_1}");
    REQUIRE(first_bcast(d) == first_bcast(f));
    // ... over a different set of modes: different values
    auto const g =
        deserialize<ResultExpr>(L"G{i_1,i_2;a_1} = ε{;a_1} + t{i_1,i_2;a_1}");
    REQUIRE(first_bcast(d) != first_bcast(g));
  }
}

TEST_CASE("binarize-system", "[EvalExpr][inline]") {
  // a real-basis registry: the tensors here pin BraKetSymmetry::Symm, which
  // over the (default) complex basis no trait derives
  [[maybe_unused]] auto const real_basis_scope = scoped_real_field_context();
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  using Catch::Matchers::ContainsSubstring;
  using Node = FullBinaryNode<EvalExpr>;
  using System = container::svector<ResultExpr>;
  using Flags = container::svector<bool>;
  auto _ = set_scoped_default_context(lax_braket_context());

  // the nodes of the subtree at root satisfying pred, in pre-order
  // (FullBinaryNode::visit walks a whole tree, not a subtree)
  auto find_all = [](Node const& root, auto pred) {
    container::svector<Node const*> result;
    auto descend = [&](auto& self, Node const& n) -> void {
      if (pred(n)) result.push_back(&n);
      if (n.leaf()) return;
      self(self, n.left());
      self(self, n.right());
    };
    descend(descend, root);
    return result;
  };
  auto find = [&find_all](Node const& root, auto pred) -> Node const* {
    auto const all = find_all(root, pred);
    return all.empty() ? nullptr : all.front();
  };
  auto labelled = [](std::wstring_view label, bool leaf) {
    return [label, leaf](Node const& n) {
      return n.leaf() == leaf && n->is_tensor() &&
             n->as_tensor().label() == label;
    };
  };
  auto variable = [](std::wstring_view label, bool leaf) {
    return [label, leaf](Node const& n) {
      return n.leaf() == leaf && n->expr()->is<Variable>() &&
             n->expr()->as<Variable>().label() == label;
    };
  };
  auto leaf_labels = [](Node const& root) {
    std::multiset<std::string> result;
    root.visit_leaf([&](Node const& n) { result.insert(n->label()); });
    return result;
  };
  // x equals tensor y, or y with bra and ket swapped (a leaf of a bra-ket
  // symmetric tensor may hold either)
  auto same_tensor = [](ExprPtr const& x, ExprPtr const& y) {
    auto const& t = y->as<Tensor>();
    return *x == *y ||
           *x == Tensor(t.label(), bra(t.ket()), ket(t.bra()), aux(t.aux()),
                        TensorSymmetries{.perm = t.symmetry(),
                                         .braket = t.braket_symmetry(),
                                         .column = t.column_symmetry()});
  };
  // x is m ⊙ pow(d,-1)
  auto is_hadamard_power = [&same_tensor](ExprPtr const& x, ExprPtr const& m,
                                          ExprPtr const& d) {
    if (!x->is<Product>() || x->as<Product>().factors().size() != 2)
      return false;
    auto const& f = x->as<Product>().factors();
    return same_tensor(f[0], m) && f[1]->is<TensorPower>() &&
           f[1]->as<TensorPower>().exponent() == rational{-1} &&
           same_tensor(f[1]->as<TensorPower>().base(), d);
  };
  auto antisymm_ell = [](auto b, auto k) {
    return ex<Tensor>(relabelled(ell_definition().result_as_tensor(),
                                 std::move(b), std::move(k)));
  };

  SECTION("a definition becomes one factor of its consumer") {
    System const sys{ell_definition(), consumer_E()};
    auto const res = binarize(sys);
    REQUIRE(res.trees.size() == 2);
    REQUIRE(res.consumed == Flags{true, false});
    REQUIRE(res.roots().size() == 1);
    Node const& e = res.trees[1];
    REQUIRE(&res.roots()[0].get() == &e);
    REQUIRE(e->op_type() == EvalOp::Product);

    Node const& ell = e.right();
    REQUIRE_FALSE(ell.leaf());
    REQUIRE(ell->as_tensor() == relabelled(ell_definition().result_as_tensor(),
                                           bra{L"i_1", L"i_2"},
                                           ket{L"a_1", L"a_2"}));
    REQUIRE(find(ell, labelled(L"M", true)));
    REQUIRE(find(ell, labelled(L"D", true)));
    REQUIRE(
        find(ell, [](Node const& n) { return n->op_type() == EvalOp::Power; }));

    // the leaves and the Hadamard subtree of the Hadamard term; the root
    // contracts L with the head ℓ instead of with M and pow(D,-1)
    auto const ref = binarize(hadamard_power_term());
    REQUIRE(ell->hash_value() == ref.right()->hash_value());
    REQUIRE(leaf_labels(e) == leaf_labels(ref));
    REQUIRE(ref.right()->as_tensor().label() != L"ℓ");

    // the inlined root alone is an equation boundary
    REQUIRE(ell->is_equation_boundary());
    REQUIRE_FALSE(e->is_equation_boundary());
    REQUIRE(find_all(ell.left(), [](Node const& n) {
              return n->is_equation_boundary();
            }).empty());
    REQUIRE(find_all(ell.right(), [](Node const& n) {
              return n->is_equation_boundary();
            }).empty());
  }

  SECTION("definitions inline recursively, in any order") {
    System const sys{D_broadcast(), ell_definition(), consumer_E()};
    auto const fwd = binarize(sys);
    REQUIRE(fwd.consumed == Flags{true, true, false});
    Node const& e = fwd.trees[2];
    REQUIRE_FALSE(find(e, labelled(L"D", true)));
    Node const* d = find(e, labelled(L"D", false));
    REQUIRE(d);
    REQUIRE((*d)->op_type() == EvalOp::Sum);
    REQUIRE((*d)->as_tensor() == power_base()->as<Tensor>());
    REQUIRE(find_all(*d, [](Node const& n) {
              return !n.leaf() && n->is_broadcast();
            }).size() == 4);

    // to_expr keeps the inlined equations below a node opaque
    auto product_of = [](ExprPtr const& l, ExprPtr const& r) {
      Product p;
      p.append(1, l, Product::Flatten::No);
      p.append(1, r, Product::Flatten::No);
      return ex<Product>(std::move(p));
    };
    REQUIRE(*to_expr(e) ==
            *product_of(contractor(), antisymm_ell(bra{L"i_1", L"i_2"},
                                                   ket{L"a_1", L"a_2"})));
    REQUIRE(is_hadamard_power(
        to_expr(e.right()),
        antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}),
        power_base()));
    Node const* power =
        find(e, [](Node const& n) { return n->op_type() == EvalOp::Power; });
    REQUIRE(power);
    REQUIRE(*to_expr(*power) == *pow(power_base(), rational{-1}));
    // a boundary node itself expands to its rhs
    REQUIRE(to_expr(*d)->is<Sum>());

    System const rsys{consumer_E(), ell_definition(), D_broadcast()};
    auto const rev = binarize(rsys);
    REQUIRE(rev.consumed == Flags{false, true, true});
    for (std::size_t k = 0; k != 3; ++k) {
      CAPTURE(k);
      Node const& f = fwd.trees[k];
      Node const& r = rev.trees[2 - k];
      REQUIRE(f->hash_value() == r->hash_value());
      REQUIRE(f.size() == r.size());
      REQUIRE(leaf_labels(f) == leaf_labels(r));
    }
  }

  SECTION("a leaf of another shape stays a leaf") {
    auto const rank = ResultExpr(
        Tensor(L"R", bra{L"i_1"}, ket{L"a_1"}),
        ex<Constant>(2) * ex<Tensor>(L"ℓ", bra{L"i_1"}, ket{L"a_1"}));
    // spaces that match neither as written nor with bra and ket swapped
    auto const spaces =
        ResultExpr(Tensor(L"R", bra{L"i_1", L"a_2"}, ket{L"a_1", L"i_2"}),
                   ex<Constant>(2) *
                       antisymm_ell(bra{L"i_1", L"a_2"}, ket{L"a_1", L"i_2"}));
    for (auto const& consumer : {rank, spaces}) {
      System const sys{ell_definition(), consumer};
      auto const res = binarize(sys);
      REQUIRE(res.consumed == Flags{false, false});
      REQUIRE(find(res.trees[1], labelled(L"ℓ", true)));
      REQUIRE_FALSE(find(res.trees[1], labelled(L"ℓ", false)));
    }
  }

  SECTION("a bra-ket-swapped use of a head") {
    // E{i_2;i_3} = L{a_1,a_2;i_3,i_1} ℓ{a_1,a_2;i_2,i_1}, ℓ as canonicalize()
    // spells it: bra and ket swapped, which names the bra-ket symmetric ℓ
    auto const L = antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_3", L"i_1"});
    auto const swapped = ResultExpr(
        head(L"E"), L * antisymm_ell(bra{L"a_1", L"a_2"}, ket{L"i_2", L"i_1"}));
    System const sys{ell_definition(), swapped};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    Node const* ell = find(res.trees[1], labelled(L"ℓ", false));
    REQUIRE(ell);
    // the head's bra takes the leaf's ket and vice versa
    REQUIRE((*ell)->as_tensor() ==
            relabelled(ell_definition().result_as_tensor(), bra{L"i_2", L"i_1"},
                       ket{L"a_1", L"a_2"}));
    Node const* m = find(*ell, labelled(L"M", true));
    REQUIRE(m);
    REQUIRE(same_tensor((*m)->expr(), antisymm(L"M", bra{L"i_2", L"i_1"},
                                               ket{L"a_1", L"a_2"})));
    // the same tree as the use spelled like the head
    auto const plain = ResultExpr(
        head(L"E"), L * antisymm_ell(bra{L"i_2", L"i_1"}, ket{L"a_1", L"a_2"}));
    System const psys{ell_definition(), plain};
    REQUIRE(binarize(psys).trees[1]->hash_value() ==
            res.trees[1]->hash_value());

    // the canonicalized consumer_E() inlines too
    auto const canon = ResultExpr(
        head(L"E"), canonicalize(consumer_E().expression()->clone()));
    System const csys{ell_definition(), canon};
    auto const cres = binarize(csys);
    REQUIRE(cres.consumed == Flags{true, false});
    REQUIRE(find(cres.trees[1], labelled(L"ℓ", false)));

    // K{i_1;a_1} = t{i_1;a_1}, read as K{a_1;i_1}
    auto k_system = [](BraKetSymmetry bk) {
      TensorSymmetries const syms{.perm = Symmetry::Nonsymm, .braket = bk};
      return System{
          ResultExpr(Tensor(L"K", bra{L"i_1"}, ket{L"a_1"}, syms),
                     ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"})),
          ResultExpr(Tensor(L"R", bra{L"i_1"}, ket{L"a_1"}),
                     ex<Constant>(2) *
                         ex<Tensor>(L"K", bra{L"a_1"}, ket{L"i_1"}, syms))};
    };
    // a conjugate-symmetric head read swapped is its conjugate
    REQUIRE_THROWS_WITH(
        binarize(k_system(BraKetSymmetry::Conjugate)),
        ContainsSubstring("inline: bra-ket-swapped use of the "
                          "conjugate-symmetric head K is not supported"));
    // a head without bra-ket symmetry read swapped is another object
    auto const nonsymm = binarize(k_system(BraKetSymmetry::Nonsymm));
    REQUIRE(nonsymm.consumed == Flags{false, false});
    REQUIRE(find(nonsymm.trees[1], labelled(L"K", true)));

    // a head matching as written wins over one matching swapped: with the
    // conjugate X{a_1;i_1} = s{a_1;i_1} defined next to X{i_1;a_1} =
    // t{i_1;a_1}, X{a_1;i_1} reads the former, in either order
    auto x_pair = [](BraKetSymmetry bk) {
      TensorSymmetries const syms{.perm = Symmetry::Nonsymm, .braket = bk};
      return std::array{ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}, syms),
                                   ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"})),
                        ResultExpr(Tensor(L"X", bra{L"a_1"}, ket{L"i_1"}, syms),
                                   ex<Tensor>(L"s", bra{L"a_1"}, ket{L"i_1"}))};
    };
    TensorSymmetries const conj{.perm = Symmetry::Nonsymm,
                                .braket = BraKetSymmetry::Conjugate};
    auto const r = ResultExpr(
        Tensor(L"R", bra{L"a_1"}, ket{L"i_1"}),
        ex<Constant>(2) * ex<Tensor>(L"X", bra{L"a_1"}, ket{L"i_1"}, conj));
    auto const [x, x_conj] = x_pair(BraKetSymmetry::Conjugate);
    for (auto const& [sys, flags] :
         {std::pair{System{x, x_conj, r}, Flags{false, true, false}},
          std::pair{System{x_conj, x, r}, Flags{true, false, false}}}) {
      std::optional<InlinedSystem<EvalExpr>> cres;
      REQUIRE_NOTHROW(cres.emplace(binarize(sys)));
      REQUIRE(cres->consumed == flags);
      Node const* xn = find(cres->trees[2], labelled(L"X", false));
      REQUIRE(xn);
      REQUIRE((*xn)->as_tensor() == x_conj.result_as_tensor());
      REQUIRE(find(*xn, labelled(L"s", true)));
      REQUIRE_FALSE(find(*xn, labelled(L"t", true)));
    }
    // a bra-ket symmetric head and its swapped spelling are one tensor
    auto const symm_pair = x_pair(BraKetSymmetry::Symm);
    REQUIRE_THROWS_WITH(binarize(System{symm_pair[0], symm_pair[1]}),
                        ContainsSubstring("inline: duplicate definition of X"));
  }

  SECTION("a matching leaf with other symmetries throws") {
    auto const e =
        ResultExpr(head(L"E"),
                   contractor() * ex<Tensor>(L"ℓ", bra{L"i_1", L"i_2"},
                                             ket{L"a_1", L"a_2"}, symmetric()));
    System const sys{ell_definition(), e};
    REQUIRE_THROWS_WITH(binarize(sys),
                        ContainsSubstring("declares different symmetries"));
  }

  SECTION("definition dummies stay distinct from the consumer's indices") {
    auto nonsymm = [](std::wstring_view label, auto b, auto k) {
      return ex<Tensor>(label, std::move(b), std::move(k),
                        TensorSymmetries{.perm = Symmetry::Nonsymm});
    };
    // λ{i_1,i_2;a_1,a_2} = M{i_1,i_2;a_3,a_2} X{a_3;a_1}
    auto const lambda = ResultExpr(
        nonsymm(L"λ", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"})->as<Tensor>(),
        nonsymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_3", L"a_2"}) *
            nonsymm(L"X", bra{L"a_3"}, ket{L"a_1"}));
    // E{i_2;i_3} = L{a_3,a_4;i_1,i_3} λ{i_1,i_2;a_3,a_4}
    auto const e = ResultExpr(
        head(L"E"),
        antisymm(L"L", bra{L"a_3", L"a_4"}, ket{L"i_1", L"i_3"}) *
            ex<Tensor>(relabelled(lambda.result_as_tensor(),
                                  bra{L"i_1", L"i_2"}, ket{L"a_3", L"a_4"})));
    System const sys{lambda, e};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    Node const& root = res.trees[1];
    Node const* l = find(root, labelled(L"λ", false));
    REQUIRE(l);
    Node const* m = find(*l, labelled(L"M", true));
    Node const* x = find(*l, labelled(L"X", true));
    REQUIRE(m);
    REQUIRE(x);
    Tensor const& mt = (*m)->as_tensor();
    Tensor const& xt = (*x)->as_tensor();
    REQUIRE(mt.bra()[0] == Index{L"i_1"});
    REQUIRE(mt.bra()[1] == Index{L"i_2"});
    REQUIRE(mt.ket()[1] == Index{L"a_4"});
    REQUIRE(xt.ket()[0] == Index{L"a_3"});
    // the dummy is contracted between M and X and is no other index of the
    // system
    Index const& dummy = xt.bra()[0];
    REQUIRE(mt.ket()[0] == dummy);
    for (auto label : {L"i_1", L"i_2", L"i_3", L"a_1", L"a_2", L"a_3", L"a_4"})
      REQUIRE(dummy != Index{label});
  }

  SECTION("a proto-indexed dummy keeps its space and follows the head") {
    Index const a2i1(L"a_2", {L"i_1"});
    // X{i_1;a_1} = S{a_1;a_2<i_1>} t{i_1;a_2<i_1>}
    auto const x = ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}),
                              ex<Tensor>(L"S", bra{L"a_1"}, ket{a2i1}) *
                                  ex<Tensor>(L"t", bra{L"i_1"}, ket{a2i1}));
    // F{i_3;a_1} = 2 X{i_3;a_1}
    auto const f = ResultExpr(
        Tensor(L"F", bra{L"i_3"}, ket{L"a_1"}),
        ex<Constant>(2) * ex<Tensor>(L"X", bra{L"i_3"}, ket{L"a_1"}));
    System const sys{x, f};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    Node const* xn = find(res.trees[1], labelled(L"X", false));
    REQUIRE(xn);
    Node const* sn = find(*xn, labelled(L"S", true));
    Node const* tn = find(*xn, labelled(L"t", true));
    REQUIRE(sn);
    REQUIRE(tn);
    Index const& dummy = (*sn)->as_tensor().ket()[0];
    REQUIRE((*tn)->as_tensor().ket()[0] == dummy);
    REQUIRE((*tn)->as_tensor().bra()[0] == Index{L"i_3"});
    REQUIRE(dummy.label() != a2i1.label());
    REQUIRE(dummy.space() == a2i1.space());
    REQUIRE(dummy.proto_indices().size() == 1);
    REQUIRE(dummy.proto_indices()[0] == Index{L"i_3"});
  }

  SECTION("a proto-indexed slot matches only a proto-indexed slot") {
    Index const a1i1(L"a_1", {L"i_1"});
    Index const a3i2(L"a_3", {L"i_2"});
    // X{i_1;a_1} = t{i_1;a_1}: a CSV use X{i_2;a_3<i_2>} is another tensor
    auto const plain_x = ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}),
                                    ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    // X{i_1;a_1<i_1>} = t{i_1;a_1<i_1>}
    auto const csv_x = ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{a1i1}),
                                  ex<Tensor>(L"t", bra{L"i_1"}, ket{a1i1}));
    // R{i_2;a_3<i_2>} = 2 X{i_2;a_3<i_2>}
    auto const r =
        ResultExpr(Tensor(L"R", bra{L"i_2"}, ket{a3i2}),
                   ex<Constant>(2) * ex<Tensor>(L"X", bra{L"i_2"}, ket{a3i2}));
    System const mismatch{plain_x, r};
    auto const res = binarize(mismatch);
    REQUIRE(res.consumed == Flags{false, false});
    REQUIRE(find(res.trees[1], labelled(L"X", true)));

    System const match{csv_x, r};
    auto const mres = binarize(match);
    REQUIRE(mres.consumed == Flags{true, false});
    Node const* t = find(mres.trees[1], labelled(L"t", true));
    REQUIRE(t);
    REQUIRE((*t)->as_tensor().ket()[0] == a3i2);
  }

  SECTION("two uses of one definition share a hash") {
    // G{i_2;i_3} = L{a_1,a_2;i_1,i_3} ℓ{i_1,i_2;a_1,a_2}
    //            + K{a_1,a_2;i_1,i_2} ℓ{i_1,i_3;a_1,a_2}
    auto const g = ResultExpr(
        head(L"G"),
        contractor() * antisymm_ell(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) +
            antisymm(L"K", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
                antisymm_ell(bra{L"i_1", L"i_3"}, ket{L"a_1", L"a_2"}));
    System const sys{ell_definition(), g};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    auto const uses = find_all(res.trees[1], labelled(L"ℓ", false));
    REQUIRE(uses.size() == 2);
    REQUIRE((*uses[0])->hash_value() == (*uses[1])->hash_value());
    REQUIRE((*uses[0])->as_tensor().bra()[1] !=
            (*uses[1])->as_tensor().bra()[1]);
    // each use expands to the definition's rhs over its own indices
    for (auto const* use : uses) {
      Tensor const& head = (*use)->as_tensor();
      REQUIRE(is_hadamard_power(
          to_expr(*use), antisymm(L"M", bra(head.bra()), ket(head.ket())),
          ex<Tensor>(L"D", bra(head.bra()), ket(head.ket()), symmetric())));
    }
  }

  SECTION("a Sum chain keeps an inlined Sum opaque") {
    auto t = [](std::wstring_view label) {
      return ex<Tensor>(label, bra{L"i_1"}, ket{L"a_1"});
    };
    // X{i_1;a_1} = f{i_1;a_1} + g{i_1;a_1}
    auto const x =
        ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}), t(L"f") + t(L"g"));
    // R{i_1;a_1} = X{i_1;a_1} + h{i_1;a_1} + k{i_1;a_1}, X at either end of
    // the chain
    for (bool x_first : {true, false}) {
      auto const rhs = x_first
                           ? ex<Sum>(ExprPtrList{t(L"X"), t(L"h"), t(L"k")})
                           : ex<Sum>(ExprPtrList{t(L"h"), t(L"k"), t(L"X")});
      auto const r = ResultExpr(Tensor(L"R", bra{L"i_1"}, ket{L"a_1"}), rhs);
      // the definition listed before and after its consumer
      for (bool def_first : {true, false}) {
        CAPTURE(x_first, def_first);
        System const sys = def_first ? System{x, r} : System{r, x};
        auto const res = binarize(sys);
        REQUIRE(res.roots().size() == 1);
        auto const expr = to_expr(res.roots()[0].get());
        REQUIRE(expr->is<Sum>());
        auto const& summands = expr->as<Sum>().summands();
        REQUIRE(summands.size() == 3);
        for (auto label : {L"X", L"h", L"k"})
          REQUIRE(std::ranges::count_if(summands, [&](ExprPtr const& s) {
                    return *s == *t(label);
                  }) == 1);
      }
    }
  }

  SECTION("a cyclic definition throws") {
    auto const x = ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}),
                              ex<Tensor>(L"Y", bra{L"i_1"}, ket{L"a_1"}));
    auto const y = ResultExpr(
        Tensor(L"Y", bra{L"i_1"}, ket{L"a_1"}),
        ex<Constant>(2) * ex<Tensor>(L"X", bra{L"i_1"}, ket{L"a_1"}));
    System const cycle{x, y};
    REQUIRE_THROWS_WITH(binarize(cycle),
                        ContainsSubstring("inline: cyclic definition of"));

    auto const self = ResultExpr(
        Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}),
        ex<Tensor>(L"X", bra{L"i_1"}, ket{L"a_1"}) * ex<Variable>(L"c"));
    System const selfref{self};
    REQUIRE_THROWS_WITH(binarize(selfref),
                        ContainsSubstring("inline: cyclic definition of"));
  }

  SECTION("scalar heads inline into Variable leaves") {
    // E = L{a_1,a_2;i_1,i_2} ℓ{i_1,i_2;a_1,a_2}
    auto const e =
        ResultExpr(Variable(L"E"),
                   antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
                       antisymm_ell(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}));
    // F{i_1;a_1} = E t{i_1;a_1}
    auto const f = ResultExpr(
        Tensor(L"F", bra{L"i_1"}, ket{L"a_1"}),
        ex<Variable>(L"E") * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const sys{e, f};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    REQUIRE_FALSE(find(res.trees[1], variable(L"E", true)));
    Node const* inlined = find(res.trees[1], variable(L"E", false));
    REQUIRE(inlined);
    REQUIRE((*inlined)->op_type() == EvalOp::Product);
    REQUIRE(find(*inlined, labelled(L"L", true)));

    auto const self =
        ResultExpr(Variable(L"E"), ex<Constant>(2) * ex<Variable>(L"E"));
    System const selfref{self};
    REQUIRE_THROWS_WITH(binarize(selfref),
                        ContainsSubstring("inline: cyclic definition of"));

    // a conjugated use of a scalar head is not inlined silently
    auto e_conj = ex<Variable>(L"E");
    e_conj->as<Variable>().conjugate();
    auto const fc =
        ResultExpr(Tensor(L"F", bra{L"i_1"}, ket{L"a_1"}),
                   e_conj * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const conj{e, fc};
    REQUIRE_THROWS_WITH(
        binarize(conj),
        ContainsSubstring("inline: conjugated use of head E is not supported"));

    // a scalar head inside a scalar power is neither a leaf nor inlined
    auto const gp = ResultExpr(Tensor(L"G", bra{L"i_1"}, ket{L"a_1"}),
                               ex<Power>(ex<Variable>(L"E"), 2) *
                                   ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const powered{e, gp};
    REQUIRE_THROWS_WITH(binarize(powered),
                        ContainsSubstring("inline: head E cannot be used "
                                          "inside a scalar power"));

    // c has no definition
    auto const g = ResultExpr(
        Tensor(L"G", bra{L"i_1"}, ket{L"a_1"}),
        ex<Variable>(L"c") * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const undefined{e, g};
    auto const ures = binarize(undefined);
    REQUIRE(ures.consumed == Flags{false, false});
    REQUIRE(find(ures.trees[1], variable(L"c", true)));
  }

  SECTION("a scalar head feeds a tensor-tensor contraction") {
    // E = L{a_1,a_2;i_1,i_2} ℓ{i_1,i_2;a_1,a_2}
    auto const e =
        ResultExpr(Variable(L"E"),
                   antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
                       antisymm_ell(bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}));
    // H{i_1;a_1} = E t{i_1;a_2} u{a_2;a_1}
    auto const h = ResultExpr(Tensor(L"H", bra{L"i_1"}, ket{L"a_1"}),
                              ex<Variable>(L"E") *
                                  ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_2"}) *
                                  ex<Tensor>(L"u", bra{L"a_2"}, ket{L"a_1"}));
    System const sys{e, h};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false});
    Node const& root = res.trees[1];
    REQUIRE(root->op_type() == EvalOp::Product);
    REQUIRE(contracted_indices(root) == Index::index_vector{Index{L"a_2"}});
    // the inlined E is a boundary in the operand that carries a_2
    Node const* inlined = find(root, variable(L"E", false));
    REQUIRE(inlined);
    REQUIRE((*inlined)->is_equation_boundary());
    REQUIRE((*inlined)->op_type() == EvalOp::Product);
    Node const& operand =
        find(root.left(), variable(L"E", false)) ? root.left() : root.right();
    REQUIRE(find(operand, labelled(L"t", true)));
  }

  SECTION("a scalar Power in a definition is inlined") {
    auto const c2 = ex<Power>(ex<Variable>(L"c"), 2);
    // X{i_1;a_1} = c^2 t{i_1;a_1}
    auto const x = ResultExpr(Tensor(L"X", bra{L"i_1"}, ket{L"a_1"}),
                              c2 * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    // F{i_2;a_2} = 2 X{i_2;a_2}
    auto const f = ResultExpr(
        Tensor(L"F", bra{L"i_2"}, ket{L"a_2"}),
        ex<Constant>(2) * ex<Tensor>(L"X", bra{L"i_2"}, ket{L"a_2"}));
    // E = c^-1, G{i_1;a_1} = E t{i_1;a_1}
    auto const e =
        ResultExpr(Variable(L"E"), ex<Power>(ex<Variable>(L"c"), -1));
    auto const g = ResultExpr(
        Tensor(L"G", bra{L"i_1"}, ket{L"a_1"}),
        ex<Variable>(L"E") * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const sys{x, f, e, g};
    auto const res = binarize(sys);
    REQUIRE(res.consumed == Flags{true, false, true, false});
    Node const* xn = find(res.trees[1], labelled(L"X", false));
    REQUIRE(xn);
    Node const* t = find(*xn, labelled(L"t", true));
    REQUIRE(t);
    REQUIRE((*t)->as_tensor().bra()[0] == Index{L"i_2"});
    REQUIRE((*t)->as_tensor().ket()[0] == Index{L"a_2"});
    REQUIRE(find(*xn, [](Node const& n) { return n->is_power(); }));
    Node const* en = find(res.trees[3], variable(L"E", false));
    REQUIRE(en);
    REQUIRE(find(*en, [](Node const& n) { return n->is_power(); }));
  }

  SECTION("batch-axis stamping is applied to the returned trees") {
    System const sys{ell_definition(), consumer_E()};
    // one annotation list cannot describe several equations
    auto const opts = BinarizationOptions{
        .node_batch_axes = {NodeBatchAnnotation{
            .axes = {{batched_occ(), BatchModeType::Contracted}}}}};
    REQUIRE_THROWS_WITH(
        binarize(sys, opts),
        ContainsSubstring("node_batch_axes cannot describe several equations"));

    auto res = binarize(sys);
    Node& e = res.trees[1];
    stamp_occ_loop(e);
    REQUIRE(e->node_slice_mask().size() == 1);
    REQUIRE(e->batch_loops_opened_here().size() == 1);
    REQUIRE(find_all(e.right(), [](Node const& n) {
              return !n->node_slice_mask().empty();
            }).empty());
  }
}

TEST_CASE("binarize-system-strict-braket", "[EvalExpr][inline]") {
  // a real-basis registry: the tensors here pin BraKetSymmetry::Symm, which
  // over the (default) complex basis no trait derives
  [[maybe_unused]] auto const real_basis_scope = scoped_real_field_context();
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  using Node = FullBinaryNode<EvalExpr>;
  // E's factors are L and the head ℓ, not the Hadamard chain L, M, pow(D,-1)
  // the strict bra-ket rule rejects
  container::svector<ResultExpr> const sys{ell_definition(), consumer_E()};
  std::optional<InlinedSystem<EvalExpr>> strict;
  REQUIRE_NOTHROW(strict.emplace(binarize(sys)));
  REQUIRE_THROWS_AS(binarize(hadamard_power_term()), Exception);

  auto lax = [&] {
    auto _ = set_scoped_default_context(lax_braket_context());
    return binarize(sys);
  }();
  Node const& e = strict->trees[1];
  Node const& e_lax = lax.trees[1];
  REQUIRE(e->hash_value() == e_lax->hash_value());
  REQUIRE(e.size() == e_lax.size());
  REQUIRE(e.right()->is_equation_boundary());
  REQUIRE(e.right()->hash_value() == e_lax.right()->hash_value());
}
