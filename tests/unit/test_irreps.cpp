#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/irreps/algebra.hpp>
#include <SeQuant/core/irreps/any.hpp>
#include <SeQuant/core/irreps/eval_tree.hpp>
#include <SeQuant/core/irreps/propagation.hpp>
#include <SeQuant/core/irreps/table_algebra.hpp>
#include <SeQuant/core/irreps/z2_algebra.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/expr_matcher.hpp>

#include <compare>
#include <format>
#include <map>
#include <string>
#include <vector>

namespace {

/// IRREPs of the cyclic group Z3 (which has complex IRREPs)
struct Z3Irrep {
  int value = 0;

  auto operator<=>(const Z3Irrep &) const = default;
};

struct Z3Algebra {
  using irrep_type = Z3Irrep;

  irrep_type totally_symmetric() const { return {0}; }
  sequant::IrrepSet<irrep_type> product(const irrep_type &lhs,
                                        const irrep_type &rhs) const {
    return {Z3Irrep{(lhs.value + rhs.value) % 3}};
  }
  irrep_type conjugate(const irrep_type &irrep) const {
    return {(3 - irrep.value) % 3};
  }
};

}  // namespace

template <>
struct std::formatter<Z3Irrep> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
  auto format(const Z3Irrep &irrep, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "G{}", irrep.value);
  }
};

static_assert(sequant::IrrepAlgebra<Z3Algebra>);

namespace custom {

/// IRREP type stringified via an overload of irrep_to_string (found via ADL)
struct NamedIrrep {
  std::string name;

  auto operator<=>(const NamedIrrep &) const = default;
};

std::string irrep_to_string(const NamedIrrep &irrep) {
  return "<" + irrep.name + ">";
}

struct NamedAlgebra {
  using irrep_type = NamedIrrep;

  irrep_type totally_symmetric() const { return {"A"}; }
  sequant::IrrepSet<irrep_type> product(const irrep_type &lhs,
                                        const irrep_type &rhs) const {
    return {lhs == totally_symmetric() ? rhs : lhs};
  }
  irrep_type conjugate(const irrep_type &irrep) const { return irrep; }
};

}  // namespace custom

static_assert(sequant::IrrepAlgebra<custom::NamedAlgebra>);

TEST_CASE("irreps", "[irreps]") {
  using namespace sequant;

  SECTION("IrrepSet") {
    using Set = IrrepSet<std::string>;

    Set set{"b", "a", "b"};
    REQUIRE(set.size() == 2);
    REQUIRE(to_string(set) == "{a,b}");
    REQUIRE(set.contains("a"));
    REQUIRE(!set.contains("c"));
    REQUIRE(!set.empty());
    REQUIRE(!set.is_any());

    REQUIRE(Set{}.empty());
    REQUIRE(to_string(Set{}) == "{}");

    Set any = Set::any();
    REQUIRE(any.is_any());
    REQUIRE(!any.empty());
    REQUIRE(any.contains("whatever"));
    REQUIRE(to_string(any) == "*");

    Set united = set;
    united.unite(Set{"c"});
    REQUIRE(united == Set{"a", "b", "c"});
    united.unite(any);
    REQUIRE(united == any);
    REQUIRE(Set{"a"}.unite(Set{}) == Set{"a"});

    REQUIRE(Set{"a"} < Set{"b"});
    REQUIRE(Set{"z"} < any);
  }

  SECTION("Z2Algebra (symbolic)") {
    Z2Algebra<> alg;
    using Irrep = Z2Algebra<>::irrep_type;
    using Set = IrrepSet<Irrep>;

    const Irrep gx = make_irrep("Gx");
    const Irrep gy = make_irrep("Gy");

    REQUIRE(irrep_to_string(alg.totally_symmetric()) == "1");
    REQUIRE(irrep_to_string(gx) == "Gx");

    REQUIRE(alg.product(gx, gy) == Set{make_irrep("Gx*Gy")});
    REQUIRE(alg.product(gy, gx) == alg.product(gx, gy));
    REQUIRE(alg.product(gx, gx) == Set{alg.totally_symmetric()});
    REQUIRE(alg.product(gx, alg.totally_symmetric()) == Set{gx});
    REQUIRE(irrep_to_string(make_irrep("Gy*Gx")) == "Gx*Gy");
    REQUIRE(make_irrep("Gx*Gy*Gx") == gy);

    // Associativity
    const Irrep gz = make_irrep("Gz");
    REQUIRE(product(alg, alg.product(gx, gy), Set{gz}) ==
            product(alg, Set{gx}, alg.product(gy, gz)));

    REQUIRE(alg.conjugate(gx) == gx);

    // Direct construction canonicalizes the symbols
    REQUIRE(Irrep(0, {"Gy", "Gx"}) == make_irrep("Gx*Gy"));
    REQUIRE(Irrep(0, {"Gy", "Gx", "Gy"}) == gx);
    REQUIRE(Irrep(0, {"Gx", "Gx"}) == alg.totally_symmetric());
    REQUIRE(alg.product(Irrep(0, {"Gy", "Gx"}), gx) == Set{gy});
    REQUIRE_THROWS_AS(Irrep(1), Exception);
    REQUIRE_THROWS_AS(Irrep(0, {""}), Exception);
    REQUIRE_THROWS_AS(Irrep(0, {"Gx*Gy"}), Exception);
    REQUIRE_THROWS_AS(Irrep(0, {"1"}), Exception);
    REQUIRE_THROWS_AS(Z2Irrep<point_groups::D2>(0, {"B1"}), Exception);
  }

  SECTION("Z2Algebra (relations)") {
    using Irrep = Z2Algebra<>::irrep_type;
    using Set = IrrepSet<Irrep>;
    const Irrep my1 = make_irrep("My1");
    const Irrep my2 = make_irrep("My2");
    const Irrep my3 = make_irrep("My3");
    const Irrep x = make_irrep("X");

    Z2Algebra<> alg({{make_irrep("My1*My2"), my3}});

    REQUIRE(alg.product(my1, my2) == Set{my3});
    REQUIRE(alg.product(my2, my1) == Set{my3});
    REQUIRE(alg.product(my1, my3) == Set{my2});
    REQUIRE(alg.product(my3, my2) == Set{my1});
    REQUIRE(alg.product(my3, my3) == Set{alg.totally_symmetric()});
    REQUIRE(alg.reduce(make_irrep("My1*My2*My3")) == alg.totally_symmetric());
    // Unrelated symbols remain independent
    REQUIRE(alg.product(my1, x) == Set{make_irrep("My1*X")});
    REQUIRE(alg.product(alg.reduce(make_irrep("My1*X")), my2) ==
            Set{make_irrep("My3*X")});

    // Reduction is compatible with multiplication and all members of the same
    // equivalence class reduce to the same representation
    {
      const std::vector<Irrep> elements = {alg.totally_symmetric(),
                                           my1,
                                           my2,
                                           my3,
                                           x,
                                           make_irrep("My1*My2"),
                                           make_irrep("My1*My3"),
                                           make_irrep("My2*My3*X"),
                                           make_irrep("My1*My2*My3*X")};
      for (const Irrep &a : elements) {
        REQUIRE(alg.reduce(alg.reduce(a)) == alg.reduce(a));
        for (const Irrep &b : elements) {
          const Irrep plain = *Z2Algebra<>{}.product(a, b).begin();
          REQUIRE(alg.product(a, b) ==
                  alg.product(alg.reduce(a), alg.reduce(b)));
          REQUIRE(alg.product(a, b) == Set{alg.reduce(plain)});
        }
        // Multiplying with a relation doesn't change the class
        REQUIRE(
            alg.reduce(
                *Z2Algebra<>{}.product(a, make_irrep("My1*My2*My3")).begin()) ==
            alg.reduce(a));
      }
    }

    // Too many independent relations
    {
      std::vector<std::pair<Irrep, Irrep>> relations;
      for (std::size_t i = 0; i <= Z2Algebra<>::max_relations; ++i) {
        relations.emplace_back(make_irrep("S" + std::to_string(i)),
                               make_irrep("T" + std::to_string(i)));
      }
      REQUIRE_THROWS_AS(Z2Algebra<>(relations), Exception);
      relations.pop_back();
      REQUIRE_NOTHROW(Z2Algebra<>(relations));
    }

    // Without relations, products remain unevaluated
    REQUIRE(Z2Algebra<>{}.product(my1, my2) == Set{make_irrep("My1*My2")});

    // Redundant relations are fine
    Z2Algebra<> redundant({{make_irrep("My1*My2"), my3},
                           {make_irrep("My2*My3"), my1},
                           {my1, my1}});
    REQUIRE(redundant.product(my1, my2) == Set{my3});

    // Relations involving concrete IRREPs
    using point_groups::D2;
    using D2Set = IrrepSet<Z2Irrep<D2>>;
    Z2Algebra<D2> d2({{make_irrep<D2>("Gp*Gq"), make_irrep<D2>("B1")}});
    REQUIRE(d2.product(make_irrep<D2>("Gp"), make_irrep<D2>("Gq")) ==
            D2Set{make_irrep<D2>("B1")});
    REQUIRE(d2.product(make_irrep<D2>("Gp"), make_irrep<D2>("B1")) ==
            D2Set{make_irrep<D2>("Gq")});
    REQUIRE(irrep_to_string(
                *d2.product(make_irrep<D2>("Gp*B2"), make_irrep<D2>("Gq"))
                     .begin()) == "B3");

    // Contradictions
    REQUIRE_THROWS_AS(
        Z2Algebra<D2>({{make_irrep<D2>("B1"), make_irrep<D2>("B2")}}),
        Exception);
    REQUIRE_THROWS_AS(
        Z2Algebra<D2>({{make_irrep<D2>("Gp"), make_irrep<D2>("B1")},
                       {make_irrep<D2>("Gp"), make_irrep<D2>("B2")}}),
        Exception);

    // Resolved leaf IRREPs are brought into canonical form
    IrrepContext<Z2Algebra<>> ctx(
        alg, resolve_by_label<Irrep>(
                 {{L"t", Set{make_irrep("My1*My2")}}, {L"v", Set{my1}}}));
    REQUIRE(irreps_of(*deserialize(L"t{a1;i1}"), ctx) == Set{my3});
    REQUIRE(irreps_of(*deserialize(L"t{a1;i1} v{i1;a2}"), ctx) == Set{my2});
  }

  SECTION("Z2Algebra (concrete + symbolic)") {
    using point_groups::D2;
    using point_groups::D2h;

    Z2Algebra<D2> d2;
    using Set = IrrepSet<Z2Irrep<D2>>;

    REQUIRE(irrep_to_string(d2.totally_symmetric()) == "A");
    REQUIRE(d2.product(make_irrep<D2>("B1"), make_irrep<D2>("B2")) ==
            Set{make_irrep<D2>("B3")});
    REQUIRE(irrep_to_string(make_irrep<D2>("B1*Gx*B2")) == "B3*Gx");
    REQUIRE(irrep_to_string(make_irrep<D2>("B1*Gx*B1")) == "Gx");
    REQUIRE(irrep_to_string(make_irrep<D2>("B1*Gx*B1*Gx")) == "A");

    Z2Algebra<D2h> d2h;
    auto mul = [&](std::string_view lhs, std::string_view rhs) {
      return irrep_to_string(
          *d2h.product(make_irrep<D2h>(lhs), make_irrep<D2h>(rhs)).begin());
    };
    REQUIRE(mul("B3u", "B2u") == "B1g");
    REQUIRE(mul("B1g", "B1u") == "Au");
    REQUIRE(mul("B3u", "B3g") == "Au");
    REQUIRE(mul("B2g", "B1g") == "B3g");
    REQUIRE(mul("Au", "Au") == "Ag");

    Z2Algebra<point_groups::C2v> c2v;
    REQUIRE(irrep_to_string(*c2v.product(make_irrep<point_groups::C2v>("B1"),
                                         make_irrep<point_groups::C2v>("B2"))
                                 .begin()) == "A2");

    // Sets
    Set b1_or_b2{make_irrep<D2>("B1"), make_irrep<D2>("B2")};
    REQUIRE(to_string(b1_or_b2) == "{B2,B1}");
    REQUIRE(product(d2, b1_or_b2, Set{make_irrep<D2>("B3")}) == b1_or_b2);
    REQUIRE(product(d2, b1_or_b2, Set{}).empty());
    REQUIRE(product(d2, b1_or_b2, Set::any()).is_any());
    REQUIRE(product(d2, Set{}, Set::any()).empty());
  }

  SECTION("TableAlgebra") {
    using Set = IrrepSet<std::string>;

    // C3v
    TableAlgebra c3v("A1", {{{"A2", "A2"}, Set{"A1"}},
                            {{"A2", "E"}, Set{"E"}},
                            {{"E", "E"}, Set{"A1", "A2", "E"}}});

    REQUIRE(c3v.totally_symmetric() == "A1");
    REQUIRE(c3v.product("E", "E") == Set{"A1", "A2", "E"});
    REQUIRE(c3v.product("E", "A2") == Set{"E"});
    REQUIRE(c3v.product("A1", "E") == Set{"E"});
    REQUIRE(c3v.product("A2", "A1") == Set{"A2"});
    REQUIRE(product(c3v, Set{"A2", "E"}, Set{"E"}) == Set{"A1", "A2", "E"});
    REQUIRE(c3v.conjugate("E") == "E");

    // Missing entry
    TableAlgebra incomplete("A1", {{{"A2", "A2"}, Set{"A1"}}});
    REQUIRE_THROWS_AS(incomplete.product("A2", "E"), Exception);

    // Inconsistent tables
    REQUIRE_THROWS_AS((TableAlgebra("A1", {{{"A1", "E"}, Set{"A2"}}})),
                      Exception);
    REQUIRE_THROWS_AS((TableAlgebra("A1", {{{"A2", "E"}, Set{"E"}},
                                           {{"E", "A2"}, Set{"A1"}}})),
                      Exception);
    REQUIRE_THROWS_AS((TableAlgebra("A1", {{{"A2", "E"}, Set::any()}})),
                      Exception);

    // Empty products are invalid
    REQUIRE_THROWS_AS((TableAlgebra("A1", {{{"A2", "E"}, Set{}}})), Exception);

    // Conjugation map (only one direction of a conjugate pair is needed)
    TableAlgebra complex("A", {{{"E1", "E2"}, Set{"A"}}}, {{"E1", "E2"}});
    REQUIRE(complex.conjugate("E1") == "E2");
    REQUIRE(complex.conjugate("E2") == "E1");
    REQUIRE(complex.conjugate("A") == "A");
    REQUIRE(conjugate(complex, Set{"E1", "A"}) == Set{"A", "E2"});
    REQUIRE_NOTHROW(TableAlgebra("A", {}, {{"E1", "E2"}, {"E2", "E1"}}));
    REQUIRE_THROWS_AS((TableAlgebra("A", {}, {{"A", "E1"}})), Exception);
    REQUIRE_THROWS_AS((TableAlgebra("A", {}, {{"E1", "E2"}, {"E2", "E3"}})),
                      Exception);
  }

  SECTION("custom algebra") {
    Z3Algebra alg;
    using Set = IrrepSet<Z3Irrep>;

    REQUIRE(irrep_to_string(Z3Irrep{2}) == "G2");

    // Stringification via an ADL-found overload
    REQUIRE(to_string(IrrepSet<custom::NamedIrrep>{{"X"}, {"Y"}}) ==
            "{<X>,<Y>}");
    REQUIRE(AnyIrrep(custom::NamedIrrep{"X"}).to_string() == "<X>");
    REQUIRE(to_string(product(alg, Set{Z3Irrep{1}}, Set{Z3Irrep{1}})) ==
            "{G2}");
    REQUIRE(conjugate(alg, Set{Z3Irrep{1}, Z3Irrep{0}}) ==
            Set{Z3Irrep{0}, Z3Irrep{2}});
  }

  SECTION("resolvers") {
    using Irrep = Z2Irrep<point_groups::D2>;
    using Set = IrrepSet<Irrep>;
    const Set b1{make_irrep<point_groups::D2>("B1")};
    const Set b2{make_irrep<point_groups::D2>("B2")};

    const Tensor t_ov = deserialize(L"t{i1;a1}")->as<Tensor>();
    const Tensor t_vo = deserialize(L"t{a1;i1}")->as<Tensor>();
    const Tensor f = deserialize(L"f{i1;i2}")->as<Tensor>();

    auto by_label = resolve_by_label<Irrep>({{L"t", b1}});
    REQUIRE(by_label(t_ov) == b1);
    REQUIRE(by_label(f).is_any());

    std::map<Tensor, Set, TensorBlockLessThanComparator> blocks;
    blocks.emplace(deserialize(L"t{i2;a3}")->as<Tensor>(), b1);
    blocks.emplace(t_vo, b2);
    auto by_block = resolve_by_block<Irrep>(std::move(blocks));
    REQUIRE(by_block(t_ov) == b1);
    REQUIRE(by_block(t_vo) == b2);
    REQUIRE(by_block(f).is_any());

    auto by_pattern = resolve_by_pattern<Irrep>(
        {{ExprMatcher(deserialize(L"t{i1;a1}")->as<Tensor>(),
                      {.tensor_cmp = TensorComparison::Block}),
          b2}});
    REQUIRE(by_pattern(deserialize(L"t{i3;a2}")->as<Tensor>()) == b2);
    REQUIRE(by_pattern(t_vo).is_any());
  }

  SECTION("propagation") {
    using point_groups::D2;
    using Irrep = Z2Irrep<D2>;
    using Set = IrrepSet<Irrep>;
    auto irrep = [](std::string_view name) { return make_irrep<D2>(name); };

    IrrepContext<Z2Algebra<D2>> ctx(
        resolve_by_label<Irrep>({{L"t", Set{irrep("B1")}},
                                 {L"v", Set{irrep("B2")}},
                                 {L"f", Set{irrep("A")}},
                                 {L"m", Set{irrep("B1"), irrep("B2")}},
                                 {L"p", Set{irrep("Gp")}}}));

    auto irreps = [&](std::wstring_view spec) {
      return irreps_of(*deserialize(spec), ctx);
    };

    REQUIRE(irreps(L"t{a1;i1}") == Set{irrep("B1")});
    REQUIRE(irreps(L"t{a1;i1} v{i1;a2}") == Set{irrep("B3")});
    // Element-wise power of a tensor spans the IRREPs of its base
    ExprPtr m = deserialize(L"m{a1;i1}");
    REQUIRE(irreps_of(*pow(m, rational{-1}), ctx) == irreps_of(*m, ctx));
    REQUIRE(irreps_of(*m, ctx) == Set{irrep("B1"), irrep("B2")});
    REQUIRE(irreps(L"1/2 * t{a1;i1} v{i1;a2} f{a2;a3}") == Set{irrep("B3")});
    REQUIRE(irreps(L"t{a1;i1} + v{a1;i1}") == Set{irrep("B1"), irrep("B2")});
    REQUIRE(irreps(L"t{a1;i1} m{i1;a2}") == Set{irrep("A"), irrep("B3")});
    REQUIRE(irreps(L"t{a1;i1} p{i1;a2}") == Set{irrep("B1*Gp")});
    REQUIRE(irreps(L"p{a1;i1} p{i1;a2}") == Set{irrep("A")});
    // Unknown tensors span any IRREP
    REQUIRE(irreps(L"t{a1;i1} x{i1;a2}").is_any());
    // Scalars are totally symmetric
    REQUIRE(irreps(L"2 * c * t{a1;i1}") == Set{irrep("B1")});
    // Reserved tensors are totally symmetric
    REQUIRE(irreps(L"t{a1;i1} δ{i1;i2}") == Set{irrep("B1")});
    // ... but an all-ones tensor is nonzero in every IRREP block
    REQUIRE(irreps_of(*make_ones(Index(L"a_1"), SlotType::Bra), ctx).is_any());
    REQUIRE(irreps_of(*make_ones(Index(L"i_1"), SlotType::Ket), ctx).is_any());
    // Vanishing expressions span no IRREP
    REQUIRE(irreps_of(Constant(0), ctx).empty());
    REQUIRE(irreps_of(Sum{}, ctx).empty());
    REQUIRE(irreps_of(Product(0, {deserialize(L"t{a1;i1}")}), ctx).empty());
    // Powers (of scalars) are totally symmetric
    REQUIRE(irreps(L"c^(2) * t{a1;i1}") == Set{irrep("B1")});

    // Equality of propagation on plain expressions and on evaluation trees
    for (std::wstring_view spec : {L"t{a1;i1} v{i1;a2} f{a2;a3} m{a3;i2}",
                                   L"t{a1;i1} v{i1;i2} + p{a1;a2} f{a2;i2}",
                                   L"t{a1;i1} m{i1;i2} + v{a1;i3} f{i3;i2}",
                                   // a broadcast summand
                                   L"t{a1;i2} + m{a1;}"}) {
      ResultExpr res(deserialize(L"R{a1;i2}")->as<Tensor>(), deserialize(spec));
      auto tree = binarize(res);

      REQUIRE(irreps_of(tree, ctx) == irreps_of(res, ctx));
    }

    // Callback is called for every node
    {
      ResultExpr res(deserialize(L"R{a1;i2}")->as<Tensor>(),
                     deserialize(L"t{a1;i1} v{i1;i2} + p{a1;a2} f{a2;i2}"));
      auto tree = binarize(res);
      std::size_t n_nodes = 0;
      irreps_of(tree, ctx, [&](const auto &, const Set &) { n_nodes++; });
      REQUIRE(n_nodes == tree.size());
    }
  }

  SECTION("adjoint") {
    Z3Algebra alg;
    using Set = IrrepSet<Z3Irrep>;
    IrrepContext<Z3Algebra> ctx(
        alg, resolve_by_label<Z3Irrep>(
                 {{L"t", Set{Z3Irrep{1}}}, {L"v", Set{Z3Irrep{1}}}}));

    const io::serialization::DeserializationOptions nonsymm{
        .def_perm_symm = Symmetry::Nonsymm,
        .def_braket_symm = BraKetSymmetry::Nonsymm,
        .def_col_symm = ColumnSymmetry::Nonsymm};

    ExprPtr expr = deserialize(L"t{a1;i1} v{i1;a2}", nonsymm);
    REQUIRE(irreps_of(*expr, ctx) == Set{Z3Irrep{2}});

    // Adjoint of t -> conjugate IRREP
    ExprPtr adjointed = deserialize(L"t{i1;a1}", nonsymm);
    adjointed->adjoint();
    REQUIRE(adjointed->as<Tensor>().label() != L"t");
    ExprPtr prod = adjointed * deserialize(L"v{i1;a2}", nonsymm);
    REQUIRE(irreps_of(*prod, ctx) == Set{Z3Irrep{0}});

    ResultExpr res(deserialize(L"R{a2;a1}")->as<Tensor>(), prod);
    auto tree = binarize(res);
    REQUIRE(irreps_of(tree, ctx) == Set{Z3Irrep{0}});
  }

  SECTION("type erasure") {
    using point_groups::D2;
    using Irrep = Z2Irrep<D2>;
    using Set = IrrepSet<Irrep>;
    auto irrep = [](std::string_view name) { return make_irrep<D2>(name); };

    IrrepContext<Z2Algebra<D2>> typed(resolve_by_label<Irrep>(
        {{L"t", Set{irrep("B1")}}, {L"m", Set{irrep("B1"), irrep("B2")}}}));
    AnyIrrepContext erased(typed);

    for (std::wstring_view spec :
         {L"t{a1;i1}", L"t{a1;i1} m{i1;a2}", L"t{a1;i1} + m{a1;i1}",
          L"t{a1;i1} x{i1;a2}"}) {
      ExprPtr expr = deserialize(spec);
      AnyIrrepSet result = irreps_of(*expr, erased);
      const Set *typed_result = result.as<Irrep>();
      REQUIRE(typed_result != nullptr);
      REQUIRE(*typed_result == irreps_of(*expr, typed));
      REQUIRE(result.to_string() == to_string(*typed_result));
    }

    // Deep copy
    AnyIrrepContext copy(erased);
    REQUIRE(copy.as<Z2Algebra<D2>>() != nullptr);
    REQUIRE(copy.as<Z2Algebra<D2>>() != erased.as<Z2Algebra<D2>>());
    REQUIRE(copy.as<Z2Algebra<>>() == nullptr);

    // Untyped empty sets
    AnyIrrepSet empty;
    REQUIRE(empty.empty());
    REQUIRE(empty == AnyIrrepSet(Set{}));
    REQUIRE(empty != AnyIrrepSet(Set{irrep("B1")}));
    empty.unite(AnyIrrepSet(Set{irrep("B1")}));
    REQUIRE(empty == AnyIrrepSet(Set{irrep("B1")}));

    // Type mismatch
    AnyIrrepSet other_type(IrrepSet<std::string>{"A"});
    REQUIRE(other_type != empty);
    REQUIRE_THROWS_AS(empty.unite(other_type), Exception);
    REQUIRE(AnyIrrep(std::string("B1")) != AnyIrrep(irrep("B1")));
    REQUIRE(AnyIrrep(irrep("B1")) == AnyIrrep(irrep("B1")));

    // Empty sets are equal regardless of their type (transitivity)
    REQUIRE(AnyIrrepSet(Set{}) == AnyIrrepSet());
    REQUIRE(AnyIrrepSet() == AnyIrrepSet(IrrepSet<std::string>{}));
    REQUIRE(AnyIrrepSet(Set{}) == AnyIrrepSet(IrrepSet<std::string>{}));
    REQUIRE(AnyIrrepSet(Set::any()) != AnyIrrepSet(Set{}));

    // Iteration
    std::vector<std::string> names;
    AnyIrrepSet(Set{irrep("B1"), irrep("B2*Gx")})
        .for_each([&](const AnyIrrep &current) {
          names.push_back(current.to_string());
          REQUIRE(current.as<Irrep>() != nullptr);
        });
    REQUIRE(names == std::vector<std::string>{"B2*Gx", "B1"});
  }
}
