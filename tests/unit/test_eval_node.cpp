#include "catch.hpp"

#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/eval_seq.hpp>
#include <SeQuant/core/parse_expr.hpp>

// validates if x is constructible from tspec using parse_expr
auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
  return x.to_latex() == sequant::parse_expr_asymm(tspec)->to_latex();
};

TEST_CASE("TEST EVAL_NODE", "[EvalNode]") {
  using namespace sequant;
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("product") {
    // 1/16 * (A * B) * C
    const auto p1 = parse_expr_asymm(
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
    auto const sum1 = parse_expr_asymm(
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
    const auto p1 = parse_expr_asymm(
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

    auto const p2 = parse_expr_asymm(L"1/4 * g_{i2,i1}^{a1,a2}");

    auto const n2 = to_eval_node(p2);
    REQUIRE(to_expr(n2)->to_latex() == p2->to_latex());
  }

  SECTION("linearize_eval_node") {
    const auto p1 = parse_expr_asymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");
    REQUIRE(linearize_eval_node(to_eval_node(p1))->to_latex() ==
            p1->to_latex());

    auto const p2 = parse_expr_asymm(L"1/4 * g_{i2,i1}^{a1,a2}");
    REQUIRE(linearize_eval_node(to_eval_node(p2))->to_latex() ==
            parse_expr_asymm(L"1/4 * g_{i2,i1}^{a1,a2}")->to_latex());
  }

  SECTION("asy_cost_single_node") {
    auto const p1 =
        parse_expr_asymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(asy_cost_single_node(to_eval_node(p1)) == AsyCost{2, 3});

    auto const p2 = parse_expr_asymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");
    auto const n2 = to_eval_node(p2);

    REQUIRE(asy_cost_single_node(n2) == AsyCost{3, 2});
    REQUIRE(asy_cost_single_node(n2.left()) == AsyCost{3, 2});

    auto const p3 =
        parse_expr_asymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const n3 = to_eval_node(p3);
    REQUIRE(asy_cost_single_node(n3) == AsyCost{2, 1});
    REQUIRE(asy_cost_single_node(n3.left()) == AsyCost{3, 1});
  }

  SECTION("asy_cost") {
    auto const p1 =
        parse_expr_asymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(asy_cost(to_eval_node(p1)) == AsyCost{2, 3});

    auto const p2 = parse_expr_asymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");

    auto const n2 = to_eval_node(p2);
    REQUIRE(asy_cost(n2) == AsyCost{3, 2} + AsyCost{3, 2});

    auto const p3 =
        parse_expr_asymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const n3 = to_eval_node(p3);
    REQUIRE(asy_cost(n3) == AsyCost{2, 1} + AsyCost{3, 1});
  }
}
