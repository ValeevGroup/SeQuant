#include <SeQuant/core/expressions/complex.hpp>
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
    // A Nonsymm-braket tensor's adjoint() sets the '⁺' (adjointed) state.
    // Such a leaf binarizes to a plain LEAF holding the bare spelling on the
    // bare tensor's slot: the adjointness rides in the leaf's CanonTransform
    // ({conj, braket_swap}) and is applied on retrieval.
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
    auto bare_tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

    REQUIRE(tree.leaf());
    REQUIRE(tree->is_tensor());
    REQUIRE(tree->as_tensor().label() == L"t");  // bare spelling stored
    REQUIRE(tree->as_tensor().bra().at(0).label() == L"a_1");
    REQUIRE(tree->as_tensor().ket().at(0).label() == L"i_1");
    REQUIRE(tree->canon_transform().conj);
    REQUIRE(tree->canon_transform().braket_swap);

    // slot identity: the adjoint shares the bare tensor's cache slot but is
    // structurally distinct via the transform salt
    REQUIRE(bare_tree.leaf());
    REQUIRE(bare_tree->hash_value() == tree->hash_value());
    REQUIRE(bare_tree->canon_transform().trivial());

    // A Hermitian (BraKetSymmetry::Conjugate) tensor never keeps a '⁺'
    // state; both plain orientations land on ONE canonical slot, the
    // non-canonical spelling carrying the fold map {conj, braket_swap}
    // (the adjoint is the identity on a Hermitian value), and a '꙳'
    // spelling composes a pure {conj} on top.
    Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
             BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
    Tensor g_adj = g;
    REQUIRE(g_adj.adjoint() == 1);
    REQUIRE_FALSE(g_adj.adjointed());
    REQUIRE(g_adj.label() == L"g");
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree = binarize(ex<Tensor>(g_adj));
    auto g_tree2 = binarize(ex<Tensor>(g));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree.leaf());
    REQUIRE(g_tree2.leaf());
    REQUIRE_FALSE(g_tree->as_tensor().kconjugated());  // states never stored
    REQUIRE_FALSE(g_tree2->as_tensor().kconjugated());
    REQUIRE(g_tree->hash_value() == g_tree2->hash_value());  // one slot
    // exactly one of the two spellings is non-canonical: it carries the fold
    REQUIRE(g_tree->canon_transform().trivial() !=
            g_tree2->canon_transform().trivial());

    // '꙳' spelling: same slot, the state composed as a pure conj bit
    Tensor g_star = g;
    REQUIRE(g_star.kconjugate() == 1);
    REQUIRE_FALSE(g_star.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto g_tree3 = binarize(ex<Tensor>(g_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(g_tree3.leaf());
    REQUIRE_FALSE(g_tree3->as_tensor().kconjugated());
    REQUIRE(g_tree3->hash_value() == g_tree2->hash_value());
    REQUIRE(
        compose(g_tree3->canon_transform(), g_tree2->canon_transform()).conj);
  }

  SECTION("starred non-Conjugate leaves") {
    // The '꙳' state is first-class and can land on leaves whose braket
    // symmetry is not Conjugate. Symm: conj is the identity in value, so the
    // state just drops. Nonsymm: conj(t) is value-distinct and is served
    // from t's slot through a {conj} transform.
    Tensor s(L"s", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Symm, ColumnSymmetry::Nonsymm);
    Tensor s_star = s;
    [[maybe_unused]] auto s_sign = s_star.kconjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto s_tree = binarize(ex<Tensor>(s_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(s_tree.leaf());
    REQUIRE_FALSE(s_tree->as_tensor().kconjugated());
    REQUIRE_FALSE(s_tree->canon_transform().conj);  // Symm state dropped

    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
             BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
    Tensor t_star = t;
    [[maybe_unused]] auto t_sign = t_star.kconjugate();
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto t_tree = binarize(ex<Tensor>(t_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(t_tree.leaf());
    REQUIRE(t_tree->canon_transform().conj);
    REQUIRE_FALSE(t_tree->canon_transform().braket_swap);

    // slot identity: one slot for t and t꙳, the conj in the transform
    REQUIRE(EvalExpr{t}.hash_value() == EvalExpr{t_star}.hash_value());
    REQUIRE(EvalExpr{t_star}.canon_transform().conj);

    // '⁺' AND '꙳': conj(adjoint(t)) is the symbolic transpose t^T -- the two
    // channels compose to a pure {braket_swap}
    Tensor t_adj_star = t;
    [[maybe_unused]] auto adj_sign = t_adj_star.adjoint();
    [[maybe_unused]] auto star_sign = t_adj_star.kconjugate();
    REQUIRE(t_adj_star.decorated_label() == L"t⁺꙳");
    REQUIRE(t_adj_star.adjointed());
    REQUIRE(t_adj_star.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tt = binarize(ex<Tensor>(t_adj_star));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tt.leaf());
    REQUIRE(tt->hash_value() == EvalExpr{t}.hash_value());
    REQUIRE_FALSE(tt->canon_transform().conj);
    REQUIRE(tt->canon_transform().braket_swap);
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
  SECTION("Adjoint op in a binarized term") {
    // Regression: a tensor leaf can carry the Adjoint modifier without having
    // been produced by Tensor::adjoint() — e.g. when built from a label
    // string ending in '⁺' (the constructor adopts the mark into the bits).
    //
    // The leaf ctor keys off the '⁺' state to strip it into the bare
    // spelling + a {conj, braket_swap} transform (see the "Adjoint op"
    // section above). That path must tolerate a leaf that never went
    // through Tensor::adjoint().
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

  SECTION("leaf identity: orientations share one slot via the transform") {
    EvalExpr a{C};
    EvalExpr b{C_swap};
    // fold ON: both orientations land on one canonical slot; the
    // non-canonical spelling carries the fold map, and the '꙳' state never
    // survives on the stored leaf
    REQUIRE(a.hash_value() == b.hash_value());
    REQUIRE_FALSE(is_conj_leaf(a));
    REQUIRE_FALSE(is_conj_leaf(b));
    REQUIRE(a.canon_transform().trivial() != b.canon_transform().trivial());
    // a '꙳' ToT spelling: same slot, transforms differ by exactly conj
    auto C_star = C;
    REQUIRE(C_star.kconjugate() == 1);
    EvalExpr s{C_star};
    REQUIRE_FALSE(is_conj_leaf(s));
    REQUIRE(s.hash_value() == a.hash_value());
    REQUIRE(compose(s.canon_transform(), a.canon_transform()).conj);
  }

  SECTION("flat (block-canon) Conjugate leaf FOLDS at eval") {
    // A flat (protoindex-free) Conjugate leaf takes the block-canonicalize
    // branch WITH the fold: both plain orientations land on one canonical
    // slot (the non-canonical one carrying the {conj,swap} fold map), and a
    // '꙳' spelling composes a pure {conj} bit on top -- all served on
    // retrieval, no state ever stored on the leaf.
    auto F = deserialize(L"C{a_1;i_1}:N-C-S-N")->as<Tensor>();
    REQUIRE_FALSE(ranges::any_of(F.const_indices(), &Index::has_proto_indices));
    auto F_swap = F;
    REQUIRE(F_swap.adjoint() == 1);
    REQUIRE_FALSE(F_swap.adjointed());
    REQUIRE(F_swap.label() == L"C");

    EvalExpr fa{F};
    EvalExpr fb{F_swap};
    REQUIRE_FALSE(is_conj_leaf(fa));
    REQUIRE_FALSE(is_conj_leaf(fb));
    REQUIRE(fa.hash_value() == fb.hash_value());  // one slot
    REQUIRE(fa.canon_transform().trivial() != fb.canon_transform().trivial());
    // starred spelling: same slot, transforms differ by exactly conj
    auto F_star = F;
    REQUIRE(F_star.kconjugate() == 1);
    EvalExpr fs{F_star};
    REQUIRE_FALSE(is_conj_leaf(fs));
    REQUIRE(fs.hash_value() == fa.hash_value());
    REQUIRE(compose(fs.canon_transform(), fa.canon_transform()).conj);
  }
}

TEST_CASE("leaf lowering of the core states", "[eval_expr]") {
  using namespace sequant;
  auto ctx = set_scoped_default_context(
      Context{get_default_context()}.set(AssertStrictBraKetSymmetry::No));
  SECTION("adjointed leaf: the bare leaf plus a transform") {
    Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"});
    REQUIRE(t.adjoint() == 1);
    REQUIRE(t.adjointed());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree.leaf());
    REQUIRE_FALSE(tree->as_tensor().adjointed());
    REQUIRE(tree->as_tensor().bra()[0].label() == L"a_1");
  }
  SECTION("adjointed leaf moved onto a real basis: it arrives normalized") {
    // transform_indices re-normalizes the states over the new slots' field, so
    // the lowering never sees a '⁺' over a real basis: the coset rule has
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
    REQUIRE(even_tree->as_tensor().bra()[0].label() == L"a_1");

    // an indefinite parity keeps it: the conjugation rides the transform
    Tensor none = moved_onto_real(ConjugationParity::None);
    REQUIRE(none.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto none_tree = binarize(ex<Tensor>(none));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(none_tree.leaf());
    REQUIRE_FALSE(none_tree->as_tensor().adjointed());
    REQUIRE_FALSE(none_tree->as_tensor().kconjugated());
    REQUIRE(none_tree->as_tensor().bra()[0].label() == L"a_1");
  }
  SECTION("K-conjugated leaf over a real basis: a conjugating transform") {
    Index a = idx(L"a_1", Field::Real), i = idx(L"i_1", Field::Real);
    Tensor t(L"t", bra{a}, ket{i},
             TensorSymmetries{.conjugation_parity = ConjugationParity::None});
    REQUIRE(t.kconjugate() == 1);
    REQUIRE(t.kconjugated());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    auto tree = binarize(ex<Tensor>(t));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    REQUIRE(tree.leaf());
    REQUIRE_FALSE(tree->as_tensor().kconjugated());
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
  REQUIRE(D.right().leaf());
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

TEST_CASE("eval_expr_carries_canon_transform", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  EvalExpr ee{t};
  REQUIRE(ee.canon_transform().trivial());
  REQUIRE(ee.canon_phase() == ee.canon_transform().phase);  // compat accessor
}

TEST_CASE("leaf_slot_identity_is_canonical_spelling",
          "[EvalExpr][conj-transform]") {
  using namespace sequant;
  Tensor t(L"t", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor ts = t;
  REQUIRE(ts.kconjugate() == 1);
  // t and t^* SHARE one slot; the conj rides in the transform
  REQUIRE(EvalExpr{t}.hash_value() == EvalExpr{ts}.hash_value());
  REQUIRE(EvalExpr{ts}.canon_transform().conj);
  REQUIRE_FALSE(EvalExpr{t}.canon_transform().conj);
}

TEST_CASE("leaf_transform_channels", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  // Conjugate: folded (starred+swapped) spelling -> canonical slot +
  // {conj,swap}
  Tensor g(L"g", bra{L"p_1", L"p_2"}, ket{L"p_3", L"p_4"}, Symmetry::Nonsymm,
           BraKetSymmetry::Conjugate, ColumnSymmetry::Symm);
  Tensor g_folded = g;
  REQUIRE(g_folded.kconjugate() == 1);
  REQUIRE(g_folded.adjoint() == 1);  // pure swap for Conjugate
  EvalExpr eg{g}, egf{g_folded};
  REQUIRE(eg.hash_value() == egf.hash_value());  // one slot
  // the '꙳' contributes {conj}, the fold contributes {conj,swap}: the folded
  // (starred+swapped) spelling's net map differs from the plain spelling's
  // by exactly {braket_swap} -- on Hermitian values the pure transpose,
  // which is precisely the fold identity's conjugation
  auto const delta = compose(eg.canon_transform(), egf.canon_transform());
  REQUIRE(delta.braket_swap);
  REQUIRE_FALSE(delta.conj);

  // '⁺' Nonsymm adjoint: label stripped, {conj, swap}
  Tensor t(L"t", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  Tensor t_adj = t;
  REQUIRE(t_adj.adjoint() == 1);
  EvalExpr et{t}, eta{t_adj};
  REQUIRE(et.hash_value() == eta.hash_value());
  REQUIRE(eta.canon_transform().conj);
  REQUIRE(eta.canon_transform().braket_swap);
  REQUIRE(eta.as_tensor().decorated_label() == L"t");  // bare spelling

  // '⁺' + '꙳' = pure transpose {swap}
  Tensor t_adj_star = t_adj;
  REQUIRE(t_adj_star.kconjugate() == 1);
  EvalExpr etas{t_adj_star};
  REQUIRE(etas.hash_value() == et.hash_value());
  REQUIRE_FALSE(etas.canon_transform().conj);
  REQUIRE(etas.canon_transform().braket_swap);

  // Symm '꙳': dropped
  Tensor s(L"s", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
           BraKetSymmetry::Symm, ColumnSymmetry::Symm);
  Tensor s_star = s;
  REQUIRE(s_star.kconjugate() == 1);
  REQUIRE(EvalExpr{s_star}.canon_transform() == EvalExpr{s}.canon_transform());

  // ToT: the same channels via canonicalize_slots' conjugated_tensors report
  auto ct = deserialize(L"C{a_1<i_1>;i_1}:N-C-S")->as<Tensor>();
  REQUIRE(ranges::any_of(ct.const_indices(), &Index::has_proto_indices));
  Tensor ct_star = ct;
  REQUIRE(ct_star.kconjugate() == 1);
  EvalExpr ec{ct}, ecs{ct_star};
  REQUIRE(ec.hash_value() == ecs.hash_value());  // one slot
  REQUIRE(compose(ec.canon_transform(), ecs.canon_transform()).conj);
}

TEST_CASE("conj_hoisting_structural_identity", "[EvalExpr][conj-transform]") {
  using namespace sequant;
  auto A = ex<Tensor>(L"A", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  auto B = ex<Tensor>(L"B", bra{L"a_1"}, ket{L"i_1"}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);

  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto AB = binarize(A->clone() * B->clone());
  auto ABc = binarize(conjugate(A->clone() * B->clone()));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  // uniform conj HOISTS: one slot, root transform conj
  REQUIRE(AB->hash_value() == ABc->hash_value());
  REQUIRE(ABc->canon_transform().conj);
  REQUIRE_FALSE(AB->canon_transform().conj);

  // mixed conj SALTS: C·D^* keeps its own identity vs C·D
  auto Cx = ex<Tensor>(L"C", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                       BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  auto Dx = ex<Tensor>(L"D", bra{L"a_1"}, ket{L"i_2"}, Symmetry::Nonsymm,
                       BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto CD = binarize(Cx->clone() * Dx->clone());
  auto CDc = binarize(Cx->clone() * conjugate(Dx->clone()));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(CD->hash_value() != CDc->hash_value());

  // sum level: a uniformly conjugated SUM of products (the
  // \mathcal{T}-partner shape) hoists onto the unconjugated sum's slot with a
  // conj transform
  auto D = ex<Tensor>(L"D", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  auto E = ex<Tensor>(L"E", bra{L"i_1"}, ket{L"a_1"}, Symmetry::Nonsymm,
                      BraKetSymmetry::Nonsymm, ColumnSymmetry::Nonsymm);
  auto sum = D->clone() + E->clone();
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto S = binarize(sum);
  auto Sc = binarize(conjugate(sum));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(S->hash_value() == Sc->hash_value());
  REQUIRE(Sc->canon_transform().conj);
  REQUIRE_FALSE(S->canon_transform().conj);
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
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

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

  // the inner subtree occupies the SAME slot as an independent binarize of s
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
  // A tensor-of-tensors leaf is block-canonicalized IN PLACE like a flat
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
  // a MIXED-flavor bundle is stored space-major (up before down: the named
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
  // ToT leaf: the leaf's reorder parity is a child TRANSFORM (its hash is
  // phase-blind), so the parents occupy ONE slot. A product carries the
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
  // a sum's slot depends on ALL its summands: A + B and A + C are different
  // values; A + B and B + A are one value (the summand hash multiset is
  // order-blind) unless the leading summand's LAYOUT differs
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
  // the result layout is the FIRST summand's: the same summands in another
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
  // three summands: the LAST one must count too
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
  // the denoted spelling re-materializes the transform syntactically: the
  // conj bit becomes the marker (it colors the parent's graph), a swap is
  // swapped back
  using namespace sequant;
  // flat Hermitian leaf written in the non-canonical orientation: stored
  // swapped + {conj, swap}; denoted = swapped back AND marked
  auto const w = deserialize(L"F{i_2;i_1}:N-C-S")->as<Tensor>();
  EvalExpr e{w};
  if (e.canon_transform().braket_swap) {
    REQUIRE(e.canon_transform().conj);
    auto w_marked = w;
    REQUIRE(w_marked.kconjugate() == 1);
    REQUIRE(e.denoted_expr()->as<Tensor>() == w_marked);
  }
  // the marked spelling of the canonical orientation: stored unmarked +
  // {conj}; denoted = as written, marked
  auto s = deserialize(L"F{i_1;i_2}:N-C-S")->as<Tensor>();
  REQUIRE(s.kconjugate() == 1);
  EvalExpr es{s};
  REQUIRE(es.canon_transform().conj);
  REQUIRE(es.denoted_expr()->as<Tensor>() == s);
}

TEST_CASE("tot_leaf_annotation_is_slot_faithful", "[eval_expr][tot]") {
  // A ToT leaf's canon_indices() is the layout the provider serves for the
  // STORED spelling: pure proto indices first (the pair labels of
  // C{mu;a<ij>}, laid out i,j,mu), then the plain slots IN SLOT ORDER, then
  // the proto-carrying slots. The occupied kets of t{a<ij>,b<ij>;j,i} are
  // both plain slots and proto indices; they take their SLOT position, since
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
