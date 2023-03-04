#include "catch.hpp"

#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/eval_seq.hpp>
#include <SeQuant/core/parse_expr.hpp>

// validates if x is constructible from tspec using parse_expr
auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
  return x.to_latex() ==
         sequant::parse_expr(tspec, sequant::Symmetry::antisymm)->to_latex();
};

TEST_CASE("TEST EVAL_NODE", "[EvalNode]") {
  using namespace sequant;

  auto parse_expr_antisymm = [](auto const& xpr) {
    return parse_expr(xpr, Symmetry::antisymm);
  };

  SECTION("product") {
    // 1/16 * (A * B) * C
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");

    auto node1 = to_eval_node(p1);

    REQUIRE(node1->scalar() == Constant{1.0 / 16});

    REQUIRE(validate_tensor(node1->tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node1.left()->tensor(), L"I_{a1,a2}^{a3,a4}"));

    REQUIRE(validate_tensor(node1.right()->tensor(), L"t_{a3,a4}^{i1,i2}"));

    REQUIRE(
        validate_tensor(node1.left().left()->tensor(), L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(
        validate_tensor(node1.left().right()->tensor(), L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(node1.right().leaf());
    REQUIRE(node1.left().left().leaf());
    REQUIRE(node1.left().right().leaf());
    REQUIRE(node1->op() == EvalOp::Prod);
    REQUIRE(node1.left()->op() == EvalOp::Prod);

    // 1/16 * A * (B * C)
    auto node2p = Product{p1->as<Product>().scalar(), {}};
    node2p.append(p1->at(0));
    node2p.append(ex<Product>(Product{p1->at(1), p1->at(2)}));

    auto const node2 = to_eval_node(ex<Product>(node2p));

    REQUIRE(node2->scalar() == Constant{1.0 / 16});

    REQUIRE(validate_tensor(node2->tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(validate_tensor(node2.left()->tensor(), L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node2.right()->tensor(),
                            L"I_{a3,a4,a1,a2}^{i1,i2,i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().left()->tensor(), L"t_{a1,a2}^{i3,i4}"));

    REQUIRE(
        validate_tensor(node2.right().right()->tensor(), L"t_{a3,a4}^{i1,i2}"));

    REQUIRE(node2.left().leaf());
    REQUIRE(node2.right().right().leaf());
    REQUIRE(node2.right().left().leaf());
    REQUIRE(node2->op() == EvalOp::Prod);
    REQUIRE(node2.right()->op() == EvalOp::Prod);
  }

  SECTION("sum") {
    auto const sum1 = parse_expr_antisymm(
        L"X^{i1,i2}_{a1,a2} "
        L"+ Y^{i1, i2}_{a1,a2}"
        L"+ g_{i3,a1}^{i1,i2} * t_{a2}^{i3}");

    auto const node1 = to_eval_node(sum1);
    REQUIRE(node1->op() == EvalOp::Sum);
    REQUIRE(node1.left()->op() == EvalOp::Sum);
    REQUIRE(validate_tensor(node1.left()->tensor(), L"I^{i1,i2}_{a1,a2}"));
    REQUIRE(
        validate_tensor(node1.left().left()->tensor(), L"X^{i1,i2}_{a1,a2}"));
    REQUIRE(
        validate_tensor(node1.left().right()->tensor(), L"Y^{i1,i2}_{a1,a2}"));

    REQUIRE(node1.right()->op() == EvalOp::Prod);
    REQUIRE((validate_tensor(node1.right()->tensor(), L"I_{a2,a1}^{i1,i2}") ||
             validate_tensor(node1.right()->tensor(), L"I_{a1,a2}^{i2,i1}")));
    REQUIRE(
        validate_tensor(node1.right().left()->tensor(), L"g_{i3,a1}^{i1,i2}"));
    REQUIRE(validate_tensor(node1.right().right()->tensor(), L"t_{a2}^{i3}"));
  }

  SECTION("to_expr") {
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");

    auto p1_after = ex<Product>(1. / 16, ExprPtrList{});
    p1_after->as<Product>().append(
        ex<Product>(ExprPtrList{p1->at(0)->clone(), p1->at(1)->clone()}));
    p1_after->as<Product>().append(p1->at(2)->clone());

    auto const n1 = to_eval_node(p1);

    REQUIRE(to_expr(n1)->to_latex() == p1_after->to_latex());

    auto const p2 = parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}");

    auto const n2 = to_eval_node(p2);
    REQUIRE(to_expr(n2)->to_latex() == p2->to_latex());
  }

  SECTION("linearize_eval_node") {
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");
    REQUIRE(linearize_eval_node(to_eval_node(p1))->to_latex() ==
            p1->to_latex());

    auto const p2 = parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}");
    REQUIRE(linearize_eval_node(to_eval_node(p2))->to_latex() ==
            parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}")->to_latex());
  }

  SECTION("asy_cost_single_node") {
    auto const p1 =
        parse_expr_antisymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(AsyCost{asy_cost_single_node(to_eval_node(p1))} ==
            AsyCost{2, 2, 3});

    auto const p2 = parse_expr_antisymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");
    auto const n2 = to_eval_node(p2);

    REQUIRE(AsyCost{asy_cost_single_node(n2)} == AsyCost{2, 3, 2});
    REQUIRE(AsyCost{asy_cost_single_node(n2.left())} == AsyCost{2, 3, 2});

    auto const p3 =
        parse_expr_antisymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const n3 = to_eval_node(p3);
    REQUIRE(AsyCost{asy_cost_single_node(n3)} == AsyCost{2, 2, 1});
    REQUIRE(AsyCost{asy_cost_single_node(n3.left())} == AsyCost{2, 3, 1});
  }

  SECTION("asy_cost") {
    auto const p1 =
        parse_expr_antisymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(asy_cost(to_eval_node(p1)) == AsyCost{2, 2, 3});

    auto const p2 = parse_expr_antisymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");

    auto const np2 = to_eval_node(p2);
    REQUIRE(asy_cost(np2) == AsyCost{2, 3, 2} + AsyCost{2, 3, 2});

    auto const p3 =
        parse_expr_antisymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const np3 = to_eval_node(p3);
    REQUIRE(asy_cost(np3) == AsyCost{2, 2, 1} + AsyCost{2, 3, 1});

    auto const t1 = parse_expr_antisymm(L"I{i1,i2,i3;a1,a2,a3}");
    auto const nt1a = to_eval_node_antisymm(t1);

    REQUIRE(asy_cost_symm_off(nt1a) == AsyCost{36, 3, 3});  // 36*O^3*V^3

    auto const nt1s = to_eval_node_symm(t1);
    REQUIRE(asy_cost_symm_off(nt1s) == AsyCost{6, 3, 3});  // 6*O^3*V^3

    auto const s1 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::symm);
    auto const ns1 = to_eval_node(s1);
    REQUIRE(asy_cost(ns1) == AsyCost{{1, 4}, 2, 2});  // 1/4 * O^2V^2

    auto const s2 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::antisymm);
    auto const ns2 = to_eval_node(s2);
    REQUIRE(asy_cost(ns2) == AsyCost{{1, 2}, 2, 2});  // 1/2 * O^2V^2

    auto const s3 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::nonsymm);
    auto const ns3 = to_eval_node(s3);
    REQUIRE(asy_cost(ns3) == AsyCost{2, 2});  // O^2V^2

    auto const p4 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::symm);
    auto const np4 = to_eval_node(p4);
    REQUIRE(asy_cost(np4) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p5 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::antisymm);
    auto const np5 = to_eval_node(p5);
    REQUIRE(asy_cost(np5) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p6 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::antisymm);
    auto const np6 = to_eval_node(p6);
    REQUIRE(asy_cost(np6) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p7 = parse_expr(L"I{i1;a1} * I{i2;a2}", Symmetry::nonsymm);
    auto const np7 = to_eval_node(p7);
    REQUIRE(asy_cost(np7) == AsyCost{2, 2});  // 1/2 * 2 * O^2V^4

    auto const p8 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::nonsymm);
    auto const np8 = to_eval_node(p8);
    REQUIRE(asy_cost(np8) == AsyCost{2, 2, 4});  // 2 * O^2V^4
  }
}
