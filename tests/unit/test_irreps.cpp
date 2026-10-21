#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/irreps/algebra.hpp>
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

    // Conjugation map
    TableAlgebra complex("A", {{{"E1", "E2"}, Set{"A"}}}, {{"E1", "E2"}});
    REQUIRE(complex.conjugate("E1") == "E2");
    REQUIRE(complex.conjugate("A") == "A");
  }

  SECTION("custom algebra") {
    Z3Algebra alg;
    using Set = IrrepSet<Z3Irrep>;

    REQUIRE(irrep_to_string(Z3Irrep{2}) == "G2");
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

    // Equality of propagation on plain expressions and on evaluation trees
    for (std::wstring_view spec : {L"t{a1;i1} v{i1;a2} f{a2;a3} m{a3;i2}",
                                   L"t{a1;i1} v{i1;i2} + p{a1;a2} f{a2;i2}",
                                   L"t{a1;i1} m{i1;i2} + v{a1;i3} f{i3;i2}"}) {
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
}
