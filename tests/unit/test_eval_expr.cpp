#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"
#include "tensor_power_fixtures.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <algorithm>
#include <array>
#include <initializer_list>
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

TEST_CASE("eval_expr", "[EvalExpr]") {
  using namespace std::string_literals;
  using sequant::EvalExpr;
  using namespace sequant;

  SECTION("Constructors") {
    auto t1 = parse_tensor(L"t_{i1, i2}^{a1, a2}");

    REQUIRE_NOTHROW(EvalExpr{t1});

    auto p1 = deserialize(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};

    REQUIRE_NOTHROW(EvalExpr{Variable{L"λ"}});

    REQUIRE_NOTHROW(EvalExpr{Constant{1}});

    REQUIRE_NOTHROW(
        EvalExpr{Power(ex<Constant>(rational{1, 2}), rational{1, 2})});
    REQUIRE_NOTHROW(EvalExpr{Power(ex<Variable>(L"x"), rational{3, 1})});
  }

  SECTION("EvalExpr::EvalOp types") {
    auto t1 = parse_tensor(L"t_{i1, i2}^{a1, a2}");

    auto x1 = EvalExpr(t1);

    REQUIRE(!x1.op_type());

    auto p1 = deserialize(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};

    auto x2 = EvalExpr(deserialize(L"1/2")->as<Constant>());
    REQUIRE(!x2.op_type());

    REQUIRE(!EvalExpr{Variable{L"λ"}}.op_type());

    REQUIRE(!EvalExpr{Power(ex<Constant>(rational{1, 2}), rational{1, 2})}
                 .op_type());
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

  SECTION(
      "scalar * tensor product node inherits the tensor operand canon_phase") {
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

  SECTION("Adjoint op") {
    // A Nonsymm-braket tensor's adjoint() relabels its label with U+207A '⁺'
    // (Tensor::adjoint() at expressions/tensor.cpp:25-41). When that leaf
    // reaches binarize, the resulting eval-tree should expose the adjoint as
    // a first-class IR op (EvalOp::Adjoint) holding the bare-label tensor as
    // its single operand — not as a leaf with the marker still in the label.
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    REQUIRE(t.label() == L"t");
    Tensor t_adj = t;
    t_adj.adjoint();
    REQUIRE(t_adj.label() == L"t⁺");
    REQUIRE(t_adj.bra().at(0).label() == L"i_1");
    REQUIRE(t_adj.ket().at(0).label() == L"a_1");

    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // The tree should NOT be a leaf — the adjoint marker on the leaf label is
    // a structural property the IR should surface as a unary op.
    REQUIRE_FALSE(tree.leaf());
    REQUIRE(tree->op_type() == EvalOp::Adjoint);
    REQUIRE(tree->is_tensor());  // adjoint of a tensor is still a tensor

    // Left child: the bare-label leaf with the marker stripped and bra/ket
    // swapped back to the operand's natural orientation.
    REQUIRE(tree.left().leaf());
    REQUIRE(tree.left()->is_tensor());
    REQUIRE(tree.left()->as_tensor().label() == L"t");
    REQUIRE(tree.left()->as_tensor().bra().at(0).label() == L"a_1");
    REQUIRE(tree.left()->as_tensor().ket().at(0).label() == L"i_1");

    // Right child: a Constant(1) sentinel — present so the full-binary-tree
    // invariant holds (every non-leaf has both children); ignored by
    // evaluate's Adjoint-case dispatch.
    REQUIRE(tree.right().leaf());
    REQUIRE(tree.right()->is_constant());
    REQUIRE(tree.right()->as_constant().value() == Constant::scalar_type{1});

    // The Adjoint node's canon_indices are the *original* (marker-bearing)
    // tensor's index order — that's what a downstream Sum/Product parent
    // sees as the operand's slot order.
    auto const& adj_canon = tree->canon_indices();
    REQUIRE(adj_canon.size() == 2);
    REQUIRE(adj_canon[0].label() == L"i_1");
    REQUIRE(adj_canon[1].label() == L"a_1");

    // The bare-leaf operand's canon_indices are the swapped order.
    auto const& bare_canon = tree.left()->canon_indices();
    REQUIRE(bare_canon.size() == 2);
    REQUIRE(bare_canon[0].label() == L"a_1");
    REQUIRE(bare_canon[1].label() == L"i_1");

    // Hash uniqueness: Adjoint(t) and the bare t leaf must have distinct
    // hashes (else cache collisions). And Adjoint(Adjoint(t)) ≡ t.
    auto bare_leaf = ex<Tensor>(t);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto bare_tree = binarize(bare_leaf);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(bare_tree.leaf());
    REQUIRE(bare_tree->hash_value() != tree->hash_value());

    // Hermitian (BraKetSymmetry::Conjugate/Symm) tensors don't carry the
    // marker — Tensor::adjoint() guards the relabel on Nonsymm only — so
    // binarize gives a plain leaf (no Adjoint op).
    Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
             BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    Tensor g_adj = g;
    g_adj.adjoint();
    REQUIRE(g_adj.label() == L"g");  // no marker added for Conjugate
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree = binarize(ex<Tensor>(g_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree.leaf());  // no Adjoint wrapper
  }

  SECTION("Adjoint op in a binarized term") {
    // Regression: a tensor leaf can carry the adjoint marker U+207A '⁺' in its
    // label without having been produced by Tensor::adjoint() — e.g. when built
    // directly from a label string that already ends in '⁺'. Adjointness is
    // tracked solely by this label marker, so such a leaf is a valid adjoint.
    //
    // binarize() keys off the '⁺' label to surface the marker-bearing leaf as
    // EvalOp::Adjoint (see the "Adjoint op" section above), stripping the
    // marker via Tensor::adjoint(). That path must tolerate a leaf that never
    // went through Tensor::adjoint().
    auto expr = deserialize(L"1/2 g{i_1,i_2;a_1,a_2} t⁺{a_1;i_1} t{a_2;i_2}",
                            {.def_braket_symm = BraKetSymmetry::Nonsymm});
    REQUIRE(expr->is<Product>());

    bool has_marker_leaf = false;
    for (auto const& factor : expr->as<Product>().factors())
      has_marker_leaf |=
          factor->is<Tensor>() && factor->as<Tensor>().label() == L"t⁺";
    REQUIRE(has_marker_leaf);

    // binarize() must not throw on this term:
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    REQUIRE_NOTHROW(binarize(expr));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }

  SECTION("Power op") {
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

    // power of an adjoint-marked base: Power over Adjoint over the bare leaf
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor t_adj = t;
    t_adj.adjoint();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto pa = binarize(pow(ex<Tensor>(t_adj), rational{2}));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(pa->op_type() == EvalOp::Power);
    REQUIRE(pa.left()->op_type() == EvalOp::Adjoint);
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

    // without external indices: i1,i2,a3 are all contracted
    // result has only {a1,a2}
    {
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

    // with external={i1,i2}: only a3 is contracted
    // result has {a1,a2,i1,i2} with i1,i2 in aux
    {
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

#if 0
  SECTION("Symmetry of sum") {
    auto tensor = [](Symmetry s) {
      return deserialize(L"I_{i1,i2}^{a1,a2}", s)->as<Tensor>();
    };

    auto symmetry = [](const EvalExpr& x) {
      return x.expr()->as<Tensor>().symmetry();
    };

    auto imed = [](const Tensor& t1, const Tensor& t2) {
      return result_expr(EvalExpr{t1}, EvalExpr{t2}, EvalOp::Sum);
    };

    const auto t1 = tensor(Symmetry::Antisymm);
    const auto t2 = tensor(Symmetry::Antisymm);

    const auto t3 = tensor(Symmetry::Symm);
    const auto t4 = tensor(Symmetry::Symm);

    const auto t5 = tensor(Symmetry::Nonsymm);
    const auto t6 = tensor(Symmetry::Nonsymm);

    // sum of two antisymm tensors.
    REQUIRE(symmetry(imed(t1, t2)) == Symmetry::Antisymm);

    // sum of one antisymm and one symmetric tensors
    REQUIRE(symmetry(imed(t1, t3)) == Symmetry::Symm);

    // sum of two symmetric tensors
    REQUIRE(symmetry(imed(t3, t4)) == Symmetry::Symm);

    // sum of an antisymmetric and a nonsymmetric tensors
    REQUIRE(symmetry(imed(t1, t5)) == Symmetry::Nonsymm);

    // sum of one symmetric and one nonsymmetric tensors
    REQUIRE(symmetry(imed(t3, t5)) == Symmetry::Nonsymm);

    // sum of two nonsymmetric tensors
    REQUIRE(symmetry(imed(t5, t6)) == Symmetry::Nonsymm);
  }
#endif

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

TEST_CASE("lower_broadcast_sum", "[EvalExpr][broadcast]") {
  using namespace sequant;

  auto one = [](std::wstring_view label, SlotType slot) {
    return make_ones(Index{label}, slot)->as<Tensor>();
  };
  // the all-ones factors of a lowered summand, in order
  auto ones_of = [](const ExprPtr& summand) {
    container::svector<Tensor> result;
    for (const ExprPtr& factor : summand->as<Product>().factors()) {
      if (factor->is<Tensor>() && is_ones(factor->as<Tensor>())) {
        result.push_back(factor->as<Tensor>());
      }
    }
    return result;
  };
  using Ones = container::svector<Tensor>;

  const ResultExpr d4 = deserialize<ResultExpr>(
      L"D{i_1,i_2;a_1,a_2} = ε{;a_1} + ε{;a_2} - ε{i_1;} - ε{i_2;}");

  SECTION("outer products with ones leaves in head order and slots") {
    const ResultExpr lowered = lower_broadcast_sum(d4);
    REQUIRE(lowered.result_as_tensor() == d4.result_as_tensor());
    REQUIRE(lowered.expression()->is<Sum>());
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 4);

    const auto& eps = d4.expression()->as<Sum>().summands();
    const auto i1 = one(L"i_1", SlotType::Bra);
    const auto i2 = one(L"i_2", SlotType::Bra);
    const auto a1 = one(L"a_1", SlotType::Ket);
    const auto a2 = one(L"a_2", SlotType::Ket);
    const std::array<Ones, 4> expected{Ones{i1, i2, a2}, Ones{i1, i2, a1},
                                       Ones{i2, a1, a2}, Ones{i1, a1, a2}};
    for (std::size_t s = 0; s != 4; ++s) {
      CAPTURE(s);
      const ExprPtr& summand = sum.summand(s);
      REQUIRE(summand->is<Product>());
      const Product& prod = summand->as<Product>();
      REQUIRE(prod.size() == 4);
      // the ε factor of the input summand comes first
      const ExprPtr& eps_in =
          eps[s]->is<Product>() ? eps[s]->as<Product>().factor(0) : eps[s];
      REQUIRE(*prod.factor(0) == *eps_in);
      REQUIRE(ones_of(summand) == expected[s]);
      // −ε{i_1;} and −ε{i_2;} keep their sign
      REQUIRE(prod.scalar() == (s < 2 ? 1 : -1));
    }
  }

  SECTION("idempotent") {
    const ResultExpr lowered = lower_broadcast_sum(d4);
    REQUIRE(*lower_broadcast_sum(lowered).expression() ==
            *lowered.expression());

    // a Sum whose summands all cover the head is unchanged
    const ResultExpr complete = deserialize<ResultExpr>(
        L"R{i_1;a_1} = f{i_1;a_1} + g{i_1;i_2} t{i_2;a_1}");
    REQUIRE(*lower_broadcast_sum(complete).expression() ==
            *complete.expression());
  }

  SECTION("non-Sum unchanged") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"R{i_1;a_1} = 2 g{i_1;i_2} t{i_2;a_1}");
    const ResultExpr lowered = lower_broadcast_sum(res);
    REQUIRE(*lowered.expression() == *res.expression());
    REQUIRE(lowered.result_as_tensor() == res.result_as_tensor());
  }

  SECTION("bra-ket-swapped supply counts as present") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"D{i_1;a_1} = ε{a_1;} - ε{i_1;}");
    const ResultExpr lowered = lower_broadcast_sum(res);
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 2);
    REQUIRE(ones_of(sum.summand(0)) == Ones{one(L"i_1", SlotType::Bra)});
    REQUIRE(ones_of(sum.summand(1)) == Ones{one(L"a_1", SlotType::Ket)});
    REQUIRE(sum.summand(1)->as<Product>().scalar() == -1);
  }

  SECTION("scalar summands") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"D{i_1;a_1} = ε{;a_1} - ε{i_1;} - E");
    const ResultExpr lowered = lower_broadcast_sum(res);
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 3);
    const Product& e = sum.summand(2)->as<Product>();
    REQUIRE(e.size() == 3);
    REQUIRE(*e.factor(0) == Variable(L"E"));
    REQUIRE(ones_of(sum.summand(2)) ==
            Ones{one(L"i_1", SlotType::Bra), one(L"a_1", SlotType::Ket)});
    REQUIRE(e.scalar() == -1);

    // a Constant summand becomes the scaled ones leaves
    const ResultExpr c = deserialize<ResultExpr>(L"R{a_1} = 2 + t{a_1}");
    const ResultExpr clowered = lower_broadcast_sum(c);
    const Sum& csum = clowered.expression()->as<Sum>();
    REQUIRE(csum.size() == 2);
    const auto cit = std::ranges::find_if(csum.summands(), [](const auto& s) {
      return s->template is<Product>() &&
             s->template as<Product>().scalar() == 2;
    });
    REQUIRE(cit != csum.summands().end());
    REQUIRE((*cit)->as<Product>().size() == 1);
    REQUIRE(ones_of(*cit) == Ones{one(L"a_1", SlotType::Bra)});
  }

  SECTION("a nested Sum broadcasts at the top level only") {
    const ResultExpr res = deserialize<ResultExpr>(
        L"R{i_1;a_1} = (f{i_1;} + g{i_1;}) c + h{i_1;a_1}");
    const ResultExpr lowered = lower_broadcast_sum(res);
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 2);
    REQUIRE(*sum.summand(1) == *res.expression()->as<Sum>().summand(1));
    const Product& prod = sum.summand(0)->as<Product>();
    REQUIRE(prod.size() == 3);
    REQUIRE(*prod.factor(0) == *deserialize(L"f{i_1;} + g{i_1;}"));
    REQUIRE(*prod.factor(1) == Variable(L"c"));
    REQUIRE(ones_of(sum.summand(0)) == Ones{one(L"a_1", SlotType::Ket)});
  }

  SECTION("a head index used twice in a summand is supplied") {
    const ResultExpr res = deserialize<ResultExpr>(
        L"R{i_1;a_1} = f{i_1;a_1} + M{i_1;a_1} D{a_1;i_1}");
    REQUIRE(is_valid(res));
    const ResultExpr lowered = lower_broadcast_sum(res);
    REQUIRE(*lowered.expression() == *res.expression());
    REQUIRE_NOTHROW(binarize(res));
  }

  SECTION("a symmetrizer in a nested Sum supplies no index") {
    const ResultExpr res = deserialize<ResultExpr>(
        L"R{a1,a2;i1,i2} = c (Â{i1,i2;a1,a2} t{a1;i1} ε{a2;} + "
        L"Â{i1,i2;a1,a2} t{a2;i1} ε{a1;}) + g{a1,a2;i1,i2}");
    const ResultExpr lowered = lower_broadcast_sum(res);
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 2);
    REQUIRE(ones_of(sum.summand(0)) == Ones{one(L"i_2", SlotType::Ket)});
    REQUIRE(*sum.summand(1) == *res.expression()->as<Sum>().summand(1));
  }

  SECTION("aux index cannot be broadcast") {
    const ResultExpr res =
        deserialize<ResultExpr>(L"R{i_1;;a_2} = f{i_1;;a_2} + g{i_1;}");
    REQUIRE_FALSE(is_valid(res));
    REQUIRE_THROWS_AS(lower_broadcast_sum(res), Exception);

    // a summand that carries the aux index may lack a bra/ket head index
    const ResultExpr ok =
        deserialize<ResultExpr>(L"R{i_1,i_2;;a_2} = f{i_1;;a_2} + g{i_2;;a_2}");
    REQUIRE(is_valid(ok));
    const ResultExpr lowered = lower_broadcast_sum(ok);
    const Sum& sum = lowered.expression()->as<Sum>();
    REQUIRE(sum.size() == 2);
    REQUIRE(ones_of(sum.summand(0)) == Ones{one(L"i_2", SlotType::Bra)});
    REQUIRE(ones_of(sum.summand(1)) == Ones{one(L"i_1", SlotType::Bra)});
  }

  SECTION("binarize") {
    std::optional<FullBinaryNode<EvalExpr>> tree;
    REQUIRE_NOTHROW(tree.emplace(binarize(d4)));
    REQUIRE((*tree)->op_type() == EvalOp::Sum);
    REQUIRE((*tree)->as_tensor() == d4.result_as_tensor());
    std::size_t n_eps = 0, n_ones = 0;
    tree->visit_leaf([&](const auto& node) {
      // the −1 prefactors are scalar leaves
      if (!node->is_tensor()) {
        REQUIRE(node->is_constant());
        return;
      }
      if (is_ones(node->as_tensor()))
        ++n_ones;
      else if (node->as_tensor().label() == L"ε")
        ++n_eps;
    });
    REQUIRE(n_eps == 4);
    REQUIRE(n_ones == 12);
  }
}

TEST_CASE("binarize a system of ResultExprs", "[EvalExpr][inline]") {
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
              return n.leaf() && n->is_tensor() && is_ones(n->as_tensor());
            }).size() == 12);

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
    auto const spaces =
        ResultExpr(Tensor(L"R", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}),
                   ex<Constant>(2) *
                       antisymm_ell(bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}));
    for (auto const& consumer : {rank, spaces}) {
      System const sys{ell_definition(), consumer};
      auto const res = binarize(sys);
      REQUIRE(res.consumed == Flags{false, false});
      REQUIRE(find(res.trees[1], labelled(L"ℓ", true)));
      REQUIRE_FALSE(find(res.trees[1], labelled(L"ℓ", false)));
    }
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

    // c has no definition
    auto const g = ResultExpr(
        Tensor(L"G", bra{L"i_1"}, ket{L"a_1"}),
        ex<Variable>(L"c") * ex<Tensor>(L"t", bra{L"i_1"}, ket{L"a_1"}));
    System const undefined{e, g};
    auto const ures = binarize(undefined);
    REQUIRE(ures.consumed == Flags{false, false});
    REQUIRE(find(ures.trees[1], variable(L"c", true)));
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

TEST_CASE("binarize a system under the strict bra-ket context",
          "[EvalExpr][inline]") {
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
