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

TEST_CASE("eval_expr", "[EvalExpr]") {
  using namespace std::string_literals;
  using sequant::EvalExpr;
  using namespace sequant;
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

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

TEST_CASE("eval_expr_outer_proto_position",
          "[EvalExpr][outer-proto-position]") {
  using namespace sequant;

  // A proto-only spectator node: both legs are PNO-like composites
  // a_1<i_1,i_2> and a_2<i_1,i_2> that share the occ protos i_1, i_2. The occ
  // indices appear ONLY as protoindices; there is no top-level occ slot. Yet
  // the tensor-of-tensor array's outer modes are exactly these distinct protos,
  // in canonical order, which is what slice_mode()/write_into_slice() address.
  //
  // Build the node directly via the low-level EvalExpr ctor with canon_indices
  // = the two composite legs, mirroring a real contraction intermediate (e.g.
  // the PPL residual W) whose proto-only externals are dropped from the result
  // slots by binarize (target_indices, eval_expr.cpp). The leaf ctor is NOT
  // used here: its per-tensor slot canonicalization promotes the shared protos
  // to top-level named slots, which is not the proto-only shape under test.
  auto const tnsr = parse_tensor(L"W{a_1<i_1,i_2>;a_2<i_1,i_2>}",
                                 {.def_perm_symm = Symmetry::Nonsymm});
  EvalExpr::index_vector composite_legs;
  for (auto const& ix : tnsr.const_indices())
    if (ix.has_proto_indices()) composite_legs.push_back(ix);
  EvalExpr const node{EvalOp::Product,
                      ResultType::Tensor,
                      tnsr.clone(),
                      composite_legs,  //
                      /*phase=*/1,     //
                      /*hash=*/0,      //
                      /*connectivity=*/nullptr};

  // The node is a tensor-of-tensor whose canonical legs are all composite (they
  // each carry proto indices) -- so there is no top-level plain occ slot.
  REQUIRE(node.tot());
  auto const& legs = node.canon_indices();
  REQUIRE(legs.size() == 2);
  REQUIRE(std::all_of(legs.begin(), legs.end(),
                      [](Index const& ix) { return ix.has_proto_indices(); }));

  // The occ protos, in the order they appear on a leg == their canonical outer
  // order (both legs share the same, canonicalized, proto list).
  Index const i1 = legs.front().proto_indices().front();
  Index const i2 = legs.front().proto_indices().back();
  REQUIRE(i1 != i2);

  // Distinct outer protos in canonical order: {i1, i2} -> positions 0, 1.
  auto const p1 = outer_proto_position(node, i1);
  REQUIRE(p1.has_value());
  REQUIRE(*p1 == 0);

  auto const p2 = outer_proto_position(node, i2);
  REQUIRE(p2.has_value());
  REQUIRE(*p2 == 1);

  // An index that is not one of the node's outer protos -> nullopt.
  Index const absent{L"i_5"};
  REQUIRE_FALSE(outer_proto_position(node, absent).has_value());

  // The occ proto has no top-level slot: it does not appear as a top-level
  // entry of canon_indices() (the proto-only property that makes index_position
  // -- which scans top-level slots only -- miss it).
  REQUIRE(std::none_of(legs.begin(), legs.end(),
                       [&](Index const& ix) { return ix == i1; }));
}
