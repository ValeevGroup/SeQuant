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

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <optional>
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

namespace {
/// an Index whose space carries the given field (the default field is Complex)
sequant::Index idx(std::wstring_view label, sequant::Field field) {
  sequant::Index i(label);
  sequant::IndexSpace sp = i.space();
  sp.field(field);
  return sequant::Index(label, sp);
}
}  // namespace

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
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(res)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr.is<Variable>());
    REQUIRE(root_expr.as<Variable>().label() == L"E");

    // The binarized tree shall respect the indexing of the ResultExpr
    // (the result's `S` braket letter is derivable only over a real basis)
    auto real_basis = sequant::tests::scoped_real_basis();
    res = deserialize<ResultExpr>(
        L"Result{a2;i2}:A-S-S = g{i1,i2;a1,a2} t{a1;i1}");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(res)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr.is<Tensor>());
    REQUIRE(root_expr.as<Tensor>() ==
            Tensor(L"Result", bra(IndexList{L"a_2"}), ket(IndexList{L"i_2"}),
                   Symmetry::Antisymm, BraKetSymmetry::Symm,
                   ColumnSymmetry::Symm));

    // continued ->  check that changing indexing in result changes indexing in
    // tree
    res = deserialize<ResultExpr>(
        L"Result{i2;a2}:A-S-S = g{i1,i2;a1,a2} t{a1;i1}");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(res)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr.is<Tensor>());
    REQUIRE(root_expr.as<Tensor>() ==
            Tensor(L"Result", bra(IndexList{L"i_2"}), ket(IndexList{L"a_2"}),
                   Symmetry::Antisymm, BraKetSymmetry::Symm,
                   ColumnSymmetry::Symm));

    // The name-respecting property shall also hold for terminals
    res = deserialize<ResultExpr>(L"Other = Var");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(res)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(root_expr.is<Variable>());
    REQUIRE(root_expr.as<Variable>().label() == L"Other");

    res = deserialize<ResultExpr>(L"Amplitude{i1;a1} = t{a1;i1}");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    root_expr = binarize(res)->expr();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
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
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
      auto root = binarize(res);
      SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
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
    // A Nonsymm-braket tensor's adjoint() sets the adjointed state (printed
    // as the '⁺' mark). When that leaf
    // reaches binarize, the resulting eval-tree should expose the adjoint as
    // a first-class IR op (EvalOp::Adjoint) holding the bare-label tensor as
    // its single operand — not as a leaf with the marker still in the label.
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    REQUIRE(t.label() == L"t");
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);
    REQUIRE(t_adj.label() == L"t");
    REQUIRE(t_adj.adjointed());
    REQUIRE(t_adj.decorated_label() == L"t⁺");
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

    // A Hermitian tensor's adjoint is a bra/ket exchange: the '⁺' state
    // normalizes away, so both orientations binarize to plain leaves in
    // their as-written spelling (a leaf keeps its orientation; the eval
    // boundary never folds one orientation onto the other).
    Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
             BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    Tensor g_adj = g;
    REQUIRE(g_adj.adjoint() == 1);
    REQUIRE_FALSE(g_adj.adjointed());
    REQUIRE(g_adj.label() == L"g");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree = binarize(ex<Tensor>(g_adj));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree.leaf());
    REQUIRE_FALSE(g_tree->as_tensor().kconjugated());
    REQUIRE(g_tree->as_tensor().bra().at(0).label() == L"p_3");  // as written
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree2 = binarize(ex<Tensor>(g));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree2.leaf());
    REQUIRE_FALSE(g_tree2->as_tensor().kconjugated());
    REQUIRE(g_tree2->as_tensor().bra().at(0).label() == L"p_1");

    // a '꙳' on an even-parity tensor normalizes away with the slots in
    // place: the same plain leaf as g
    Tensor g_star = g;
    REQUIRE(g_star.kconjugate() == 1);
    REQUIRE_FALSE(g_star.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree3 = binarize(ex<Tensor>(g_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree3.leaf());
    REQUIRE(g_tree3->as_tensor().bra().at(0).label() == L"p_1");
    REQUIRE(g_tree3->hash_value() == g_tree2->hash_value());
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

    // nor t and t⁺, nor t and t⁺꙳
    Tensor t_adj = t;
    REQUIRE(t_adj.adjoint() == 1);
    REQUIRE(t_adj.adjointed());
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_adj}.hash_value());
    Tensor t_adj_star = t;
    REQUIRE(t_adj_star.adjoint() == 1);
    REQUIRE(t_adj_star.kconjugate() == 1);
    REQUIRE(t_adj_star.adjointed());
    REQUIRE(t_adj_star.kconjugated());
    REQUIRE(EvalExpr{t}.hash_value() != EvalExpr{t_adj_star}.hash_value());
    REQUIRE(EvalExpr{t_adj}.hash_value() != EvalExpr{t_adj_star}.hash_value());

    // both states: Adjoint over the '꙳' leaf, which is its own array
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto both = binarize(ex<Tensor>(t_adj_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(both->op_type() == EvalOp::Adjoint);
    REQUIRE(both.left().leaf());
    REQUIRE_FALSE(both.left()->as_tensor().adjointed());
    REQUIRE(both.left()->as_tensor().kconjugated());
    REQUIRE(both.left()->as_tensor().bra().at(0).label() == L"a_1");
    REQUIRE(both.left()->hash_value() == t_star_tree->hash_value());
  }

  SECTION("Adjoint op in a binarized term") {
    // Regression: a tensor leaf can carry the Adjoint modifier without having
    // been produced by Tensor::adjoint() — e.g. when built from a label
    // string ending in '⁺' (the constructor adopts the mark into the bits).
    //
    // binarize() reads the states to surface such a leaf as
    // EvalOp::Adjoint (see the "Adjoint op" section above), stripping the
    // marker via Tensor::adjoint(). That path must tolerate a leaf that never
    // went through Tensor::adjoint().
    auto expr = deserialize(L"1/2 g{i_1,i_2;a_1,a_2} t⁺{a_1;i_1} t{a_2;i_2}",
                            {.def_braket_symm = BraKetSymmetry::Nonsymm});
    REQUIRE(expr->is<Product>());

    bool has_marker_leaf = false;
    for (auto const& factor : expr->as<Product>().factors())
      has_marker_leaf |=
          factor->is<Tensor>() && factor->as<Tensor>().adjointed();
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
// nor the ToT (canonicalize_slots) branch exchanges a Conjugate tensor's
// bundles, so the two orientations are distinct leaves in their as-written
// spelling. A kept '꙳' state (parity None) is part of the leaf's identity:
// the state byte enters the flat leaf's hash and colours the ToT graph, so
// C and C꙳ never share a cache slot.
TEST_CASE("conjugate eval fold", "[eval_expr][conjugate-fold]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));

  // A proto-indexed (Tensor-of-Tensor) Conjugate leaf and its bra<->ket swap.
  // These evaluate to complex conjugates of each other. Proto indices route
  // the leaf ctor through canonicalize_slots; a flat tensor takes the
  // block-canonicalization branch instead. Neither path exchanges the
  // bundles. The parity is None so that a '꙳' put on the tensor stays.
  auto C = deserialize(L"C{a_1<i_1>;i_2}:N-C-S-N")->as<Tensor>();
  REQUIRE(ranges::any_of(C.const_indices(), &Index::has_proto_indices));
  auto C_swap = C;
  // swaps bra<->ket; a Hermitian tensor's '⁺' normalizes away
  REQUIRE(C_swap.adjoint() == 1);
  REQUIRE_FALSE(C_swap.adjointed());
  REQUIRE(C_swap.label() == L"C");

  auto is_conj_leaf = [](EvalExpr const& e) {
    return e.expr()->is<Tensor>() && e.expr()->as<Tensor>().kconjugated();
  };

  SECTION("leaf identity: orientations are distinct at eval") {
    EvalExpr a{C};
    EvalExpr b{C_swap};
    // distinct spellings, distinct slots, no states
    REQUIRE(a.hash_value() != b.hash_value());
    REQUIRE_FALSE(is_conj_leaf(a));
    REQUIRE_FALSE(is_conj_leaf(b));
    // a K-conjugated ToT spelling keeps its state, and the state colours the
    // graph: C and C꙳ stay distinct
    auto C_star = C;
    REQUIRE(C_star.kconjugate() == 1);
    EvalExpr s{C_star};
    REQUIRE(is_conj_leaf(s));
    REQUIRE(s.hash_value() != a.hash_value());
  }

  SECTION("flat (block-canon) Conjugate leaf keeps its orientation") {
    // A flat (protoindex-free) Conjugate leaf takes the block-canonicalization
    // branch, not canonicalize_slots. The canonicalizer never exchanges a
    // Conjugate tensor's bundles: leaves keep their as-written orientation,
    // so leaf yielders need no conjugation awareness. A K-conjugated
    // spelling keeps its state, and its leaf hash differs from its unmarked
    // twin's (the state is part of the leaf's value identity, see
    // hash_terminal_tensor).
    auto F = deserialize(L"C{a_1;i_1}:N-C-S-N")->as<Tensor>();
    REQUIRE_FALSE(ranges::any_of(F.const_indices(), &Index::has_proto_indices));
    auto F_swap = F;
    REQUIRE(F_swap.adjoint() == 1);
    REQUIRE_FALSE(F_swap.adjointed());
    REQUIRE(F_swap.label() == L"C");

    EvalExpr fa{F};
    EvalExpr fb{F_swap};
    // both orientations stay unmarked, in their own spelling
    REQUIRE_FALSE(is_conj_leaf(fa));
    REQUIRE_FALSE(is_conj_leaf(fb));
    REQUIRE(fa.hash_value() != fb.hash_value());
    auto F_star = F;
    REQUIRE(F_star.kconjugate() == 1);
    EvalExpr fs{F_star};
    REQUIRE(is_conj_leaf(fs));
    REQUIRE(fs.hash_value() != fa.hash_value());
  }
}

TEST_CASE("leaf lowering of the core states", "[eval_expr]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  SECTION("adjointed leaf: Adjoint over the bare leaf") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree->op_type() == EvalOp::Adjoint);
    REQUIRE(tree.left().leaf());
    REQUIRE_FALSE(tree.left()->as_tensor().adjointed());
    REQUIRE(tree.left()->as_tensor().bra()[0].label() == L"a_1");
  }
  SECTION("adjointed leaf rebuilt onto a real basis: still Adjoint") {
    // index-mutating APIs do not re-normalize, so a '⁺' can stand over a
    // real basis; the lowering is value-correct in every basis and must not
    // assert
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    container::map<Index, Index> to_real{
        {Index{L"a_1"}, idx(L"a_1", Field::Real)},
        {Index{L"i_1"}, idx(L"i_1", Field::Real)}};
    REQUIRE(t.transform_indices(to_real));
    t.reset_tags();  // transform_indices tags the replaced indices
    REQUIRE(t.base_field() == Field::Real);
    REQUIRE(t.adjointed());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree->op_type() == EvalOp::Adjoint);
    REQUIRE(tree.left().leaf());
    REQUIRE_FALSE(tree.left()->as_tensor().adjointed());
    REQUIRE_FALSE(tree.left()->as_tensor().kconjugated());
    REQUIRE(tree.left()->as_tensor().bra()[0].label() == L"a_1");
    REQUIRE(tree->canon_indices() != tree.left()->canon_indices());
  }
  SECTION(
      "K-conjugated leaf over a real basis: Adjoint with an identity "
      "layout") {
    Index a = idx(L"a_1", Field::Real), i = idx(L"i_1", Field::Real);
    Tensor t(L"t", bra{a}, ket{i},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.kconjugate() == 1);
    REQUIRE(t.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree->op_type() == EvalOp::Adjoint);
    REQUIRE(tree.left().leaf());
    REQUIRE_FALSE(tree.left()->as_tensor().kconjugated());
    REQUIRE(tree->canon_indices() == tree.left()->canon_indices());
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
    REQUIRE(b->as_tensor().kconjugated());
    REQUIRE(a->hash_value() != b->hash_value());
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

// The cases below build eval trees straight from expressions: the head layout
// is irrelevant to what they check (slot identity, transforms, phases), so
// the deprecated binarize(ExprPtr) is used on purpose.
SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
TEST_CASE("eval_expr_conjugation_marker_identity",
          "[EvalExpr][conjugate-fold]") {
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
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

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
  // '⁺' stays, and C{a_1;p} C⁺{p;a_2} (an Adjoint node inside the product)
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
  REQUIRE(D.right()->op_type() == EvalOp::Adjoint);
  REQUIRE(D->hash_value() != E->hash_value());
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

SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

// Task 5 (multiroot-single-dag-eval): binarize(Sum const&, ...)'s make_sum
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
