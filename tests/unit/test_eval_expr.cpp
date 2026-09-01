#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <initializer_list>
#include <memory>
#include <set>
#include <string>
#include <string_view>

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
    REQUIRE(root_expr.as<Tensor>() == Tensor(L"Amplitude",
                                             bra(IndexList{L"i_1"}),
                                             ket(IndexList{L"a_1"})));
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

    // Hermitian (BraKetSymmetry::Conjugate/Symm) tensors don't get the '⁺'
    // label decoration — Tensor::adjoint() guards the relabel on Nonsymm
    // only. At the eval boundary a flat Conjugate leaf is NOT folded onto
    // one orientation (that is the lazy-conj eval follow-up): both
    // orientations binarize to plain unmarked leaves in their as-written
    // spelling. An already-starred spelling, however, is served through an
    // EvalOp::Adjoint node over its unmarked VALUE-orientation operand.
    Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
             BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    Tensor g_adj = g;
    g_adj.adjoint();
    REQUIRE(g_adj.label() == L"g");  // no label marker added for Conjugate
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree = binarize(ex<Tensor>(g_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree.leaf());
    REQUIRE_FALSE(g_tree->as_tensor().conjugated());
    REQUIRE(g_tree->as_tensor().bra().at(0).label() == L"p_3");  // as written
    auto g_tree2 = binarize(ex<Tensor>(g));
    REQUIRE(g_tree2.leaf());
    REQUIRE_FALSE(g_tree2->as_tensor().conjugated());
    // a starred spelling is served via Adjoint over the value orientation
    Tensor g_star = g;
    g_star.conjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree3 = binarize(ex<Tensor>(g_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE_FALSE(g_tree3.leaf());
    REQUIRE(g_tree3->op_type() == EvalOp::Adjoint);
    REQUIRE(g_tree3.left().leaf());
    REQUIRE_FALSE(g_tree3.left()->as_tensor().conjugated());
    // bare operand = value orientation of g^*: bra/ket swapped back
    REQUIRE(g_tree3.left()->as_tensor().bra().at(0).label() == L"p_3");
    REQUIRE(g_tree3.right()->is_constant());  // sentinel
    REQUIRE(g_tree3->hash_value() != g_tree3.left()->hash_value());
  }

  SECTION("starred non-Conjugate leaves") {
    // The conjugation marker is first-class and can land on leaves whose
    // braket symmetry is not Conjugate. Symm: conj is the identity in value,
    // so the marker just drops and the leaf is served unmarked. Nonsymm:
    // conj(t) has no slot spelling and no Adjoint-served equivalent, so
    // binarize must refuse loudly instead of serving the wrong value
    // (lazy-conj eval is the follow-up).
    Tensor s(L"s", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Symm, ColumnSymmetry::Nonsymm);
    Tensor s_star = s;
    s_star.conjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto s_tree = binarize(ex<Tensor>(s_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(s_tree.leaf());
    REQUIRE_FALSE(s_tree->as_tensor().conjugated());

    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor t_star = t;
    t_star.conjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    REQUIRE_THROWS_AS(binarize(ex<Tensor>(t_star)), std::logic_error);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    // a directly-constructed EvalExpr must not alias t and t* onto one cache
    // slot: for Nonsymm they are value-distinct
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_star}.hash_value());

    // '⁺'-labeled AND marked (the symbolic transpose conj(adjoint(t))): the
    // marker refusal must fire before the '⁺' label channel, else the label
    // channel builds Adjoint over a still-marked bare leaf and the marker is
    // silently ignored
    Tensor t_adj_star = t;
    t_adj_star.adjoint();    // '⁺' label + slot swap
    t_adj_star.conjugate();  // marker on top: t^T, not evaluable
    REQUIRE(t_adj_star.label() == L"t⁺");
    REQUIRE(t_adj_star.conjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    REQUIRE_THROWS_AS(binarize(ex<Tensor>(t_adj_star)), std::logic_error);
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
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

// At the eval boundary a BraKetSymmetry::Conjugate leaf is
// orientation-sensitive: neither the flat (block-canonicalization) branch
// nor the ToT (canonicalize_slots) branch folds bra<->ket orientations --
// leaves keep their as-written spelling, and only an already-starred
// spelling carries a marker (served by binarize through an EvalOp::Adjoint
// node). The conjugation marker still colors the ToT graph, so T and T*
// stay distinct. Folding fresh leaves onto one orientation-shared slot is
// the lazy-conj eval follow-up.
TEST_CASE("conjugate eval fold", "[eval_expr][conjugate-fold]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));

  // A proto-indexed (Tensor-of-Tensor) Conjugate leaf and its bra<->ket swap.
  // These evaluate to complex conjugates of each other. Proto indices route
  // the leaf ctor through canonicalize_slots; a flat tensor takes the
  // block-canonicalization branch instead. Neither path folds at the eval
  // boundary.
  auto C = deserialize(L"C{a_1<i_1>;i_2}:N-C-S")->as<Tensor>();
  REQUIRE(ranges::any_of(C.const_indices(), &Index::has_proto_indices));
  auto C_swap = C;
  C_swap.adjoint();  // swaps bra<->ket; no '⁺' marker for Conjugate
  REQUIRE(C_swap.label() == L"C");

  auto is_conj_leaf = [](EvalExpr const& e) {
    return e.expr()->is<Tensor>() && e.expr()->as<Tensor>().conjugated();
  };

  SECTION("leaf identity: orientations are distinct at eval") {
    EvalExpr a{C};
    EvalExpr b{C_swap};
    // no fold: distinct spellings, distinct slots, no markers
    REQUIRE(a.hash_value() != b.hash_value());
    REQUIRE_FALSE(is_conj_leaf(a));
    REQUIRE_FALSE(is_conj_leaf(b));
    // a starred ToT spelling keeps its marker, and the marker colors the
    // graph: C and C* stay distinct
    auto C_star = C;
    C_star.conjugate();
    EvalExpr s{C_star};
    REQUIRE(is_conj_leaf(s));
    REQUIRE(s.hash_value() != a.hash_value());
  }

  SECTION("flat (block-canon) Conjugate leaf does NOT fold at eval") {
    // A flat (protoindex-free) Conjugate leaf takes the block-canonicalization
    // branch, not canonicalize_slots. There the conjugate-braket fold is
    // DISABLED: leaves keep their as-written orientation so leaf yielders and
    // evaluators need no conjugation awareness -- folding flat leaves onto
    // one orientation-shared slot is the lazy-conj eval follow-up. A marked
    // spelling keeps its marker, shares the slot of its unstarred spelling
    // (leaf-hash invariant), and binarize serves it through EvalOp::Adjoint.
    auto F = deserialize(L"C{a_1;i_1}:N-C-S")->as<Tensor>();
    REQUIRE_FALSE(ranges::any_of(F.const_indices(), &Index::has_proto_indices));
    auto F_swap = F;
    F_swap.adjoint();
    REQUIRE(F_swap.label() == L"C");

    EvalExpr fa{F};
    EvalExpr fb{F_swap};
    // no fold: both orientations stay unmarked, in their own spelling
    REQUIRE_FALSE(is_conj_leaf(fa));
    REQUIRE_FALSE(is_conj_leaf(fb));
    // a starred spelling hashes onto its unstarred spelling's slot
    auto F_star = F;
    F_star.conjugate();
    EvalExpr fs{F_star};
    REQUIRE(is_conj_leaf(fs));
    REQUIRE(fs.hash_value() == fa.hash_value());
  }
}

TEST_CASE("eval_expr_conjugation_marker_identity",
          "[EvalExpr][conjugate-fold]") {
  // C(a_1;p) C*(a_2;p) and C*(a_1;p) C(a_2;p) are one tensor S up to the
  // named relabeling a_1 <-> a_2 (its slots are "index of the unconjugated
  // factor" and "index of the conjugated factor"), so they may share one
  // eval-node hash and cache slot -- but only if their canonical layouts put
  // the unconjugated factor's index in the same slot. With the conjugation
  // marker invisible to the graph coloring, the two spellings were the same
  // colored graph with an automorphism exchanging the factors, and bliss
  // pinned the slot order by the index labels alone: the same buffer was then
  // read as S by one occurrence and as S^T* by the other (identical only for
  // real C). Regression for the Kramers-restricted CSV inter-pair overlap.
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  auto C = [](std::wstring_view ext) {
    return ex<Tensor>(L"C", bra{Index{ext}}, ket{Index{L"p_1"}},
                      Symmetry::Nonsymm, BraKetSymmetry::Conjugate,
                      ColumnSymmetry::Symm);
  };
  auto Cstar = [&C](std::wstring_view ext) {
    auto t = C(ext);
    t->as<Tensor>().conjugate();
    return t;
  };

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto A = binarize(C(L"a_1") * Cstar(L"a_2"));  // S(a_1, a_2)
  auto B = binarize(Cstar(L"a_1") * C(L"a_2"));  // S(a_2, a_1)
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

  // the genuinely different tensor C(a_1) C(a_2) (no conjugation) is kept apart
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto D = binarize(C(L"a_1") * C(L"a_2"));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(D->hash_value() != A->hash_value());
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
  // structural salt: conj/swap enter, phase does NOT (hoistable)
  REQUIRE(CanonTransform{.phase = -1}.structural_salt() ==
          CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != CanonTransform{}.structural_salt());
  REQUIRE(c.structural_salt() != s.structural_salt());
  REQUIRE(c.structural_salt() != cs.structural_salt());
}
