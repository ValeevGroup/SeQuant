#include "catch.hpp"

#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/parse_expr.hpp>

namespace {

// validates if x is constructible from tspec using parse_expr
auto validate_tensor = [](const auto& x, std::wstring_view tspec) -> bool {
  return x.to_latex() ==
         sequant::parse_expr(tspec, sequant::Symmetry::antisymm)->to_latex();
};

auto eval_node(sequant::ExprPtr const& expr) {
  return sequant::to_eval_node<sequant::EvalExpr>(expr);
}

enum struct NodePos {
  L,  // Left
  R   // Right
};

sequant::EvalExpr node(sequant::EvalNode<sequant::EvalExpr> const& n,
                       std::initializer_list<NodePos> pos) {
  if (pos.size() == 0) return *n;
  auto n_ = n;
  for (auto p : pos) {
    assert(!n_.leaf() && "Accessing child of leaf!");
    n_ = p == NodePos::L ? n_.left() : n_.right();
  }

  return *n_;
}

std::wstring tikz(sequant::EvalNode<sequant::EvalExpr> const& n) noexcept {
  return n.tikz<std::wstring>(
      [](auto&& n) { return L"$" + n->expr()->to_latex() + L"$"; },
      [](auto&&) { return L""; });
}

}  // namespace

TEST_CASE("TEST EVAL_NODE", "[EvalNode]") {
  using namespace sequant;
  auto L = NodePos::L;
  auto R = NodePos::R;

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

    auto node1 = eval_node(p1);

    REQUIRE(validate_tensor(node(node1, {}).as_tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(node(node1, {R}).as_constant() == Constant{rational{1, 16}});

    REQUIRE(
        validate_tensor(node(node1, {L}).as_tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(
        validate_tensor(node(node1, {L, L}).as_tensor(), L"I_{a1,a2}^{a3,a4}"));

    REQUIRE(
        validate_tensor(node(node1, {L, R}).as_tensor(), L"t_{a3,a4}^{i1,i2}"));

    REQUIRE(validate_tensor(node(node1, {L, L, L}).as_tensor(),
                            L"g_{i3,i4}^{a3,a4}"));

    REQUIRE(validate_tensor(node(node1, {L, L, R}).as_tensor(),
                            L"t_{a1,a2}^{i3,i4}"));

    // 1/16 * A * (B * C)
    auto node2p = Product{p1->as<Product>().scalar(), {}};
    node2p.append(1, p1->at(0), Product::Flatten::No);
    node2p.append(1, ex<Product>(Product{p1->at(1), p1->at(2)}),
                  Product::Flatten::No);

    auto const node2 = eval_node(ex<Product>(node2p));

    REQUIRE(validate_tensor(node(node2, {}).as_tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(
        validate_tensor(node(node2, {L}).as_tensor(), L"I_{a1,a2}^{i1,i2}"));

    REQUIRE(node(node2, {R}).as_constant() == Constant{rational{1, 16}});

    REQUIRE(
        validate_tensor(node(node2, {L, L}).as_tensor(), L"g{i3,i4; a3,a4}"));

    REQUIRE(validate_tensor(node(node2, {L, R}).as_tensor(),
                            L"I{a1,a2,a3,a4;i3,i4,i1,i2}"));

    REQUIRE(validate_tensor(node(node2, {L, R, L}).as_tensor(),
                            L"t{a1,a2;i3,i4}"));

    REQUIRE(validate_tensor(node(node2, {L, R, R}).as_tensor(),
                            L"t{a3,a4;i1,i2}"));
  }

  SECTION("sum") {
    auto const sum1 = parse_expr_antisymm(
        L"X^{i1,i2}_{a1,a2} "
        L"+ Y^{i1, i2}_{a1,a2}"
        L"+ g_{i3,a1}^{i1,i2} * t_{a2}^{i3}");

    auto const node1 = eval_node(sum1);
    REQUIRE(node1->op_type() == EvalOp::Sum);
    REQUIRE(node1.left()->op_type() == EvalOp::Sum);
    REQUIRE(validate_tensor(node1.left()->as_tensor(), L"I^{i1,i2}_{a1,a2}"));
    REQUIRE(validate_tensor(node1.left().left()->as_tensor(),
                            L"X^{i1,i2}_{a1,a2}"));
    REQUIRE(validate_tensor(node1.left().right()->as_tensor(),
                            L"Y^{i1,i2}_{a1,a2}"));

    REQUIRE(node1.right()->op_type() == EvalOp::Prod);
    REQUIRE(
        (validate_tensor(node1.right()->as_tensor(), L"I_{a2,a1}^{i1,i2}") ||
         validate_tensor(node1.right()->as_tensor(), L"I_{a1,a2}^{i2,i1}")));
    REQUIRE(validate_tensor(node1.right().left()->as_tensor(),
                            L"g_{i3,a1}^{i1,i2}"));
    REQUIRE(
        validate_tensor(node1.right().right()->as_tensor(), L"t_{a2}^{i3}"));
  }

  SECTION("to_expr") {
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");

    auto p1_after = ex<Product>(rational{1,16}, ExprPtrList{});
    p1_after->as<Product>().append(
        1, ex<Product>(ExprPtrList{p1->at(0)->clone(), p1->at(1)->clone()}),
        Product::Flatten::No);
    p1_after->as<Product>().append(1, p1->at(2)->clone(), Product::Flatten::No);

    auto const n1 = eval_node(p1);

    REQUIRE(to_expr(n1)->to_latex() == p1_after->to_latex());

    auto const p2 = parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}");

    auto const n2 = eval_node(p2);
    REQUIRE(to_expr(n2)->to_latex() == p2->to_latex());
  }

  SECTION("linearize_eval_node") {
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");
    REQUIRE(linearize_eval_node(eval_node(p1))->to_latex() == p1->to_latex());

    auto const p2 = parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}");
    REQUIRE(linearize_eval_node(eval_node(p2))->to_latex() ==
            parse_expr_antisymm(L"1/4 * g_{i2,i1}^{a1,a2}")->to_latex());
  }

  SECTION("asy_cost_single_node") {
    auto const p1 =
        parse_expr_antisymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(AsyCost{asy_cost_single_node(eval_node(p1))} == AsyCost{2, 2, 3});

    auto const p2 = parse_expr_antisymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");
    auto const n2 = eval_node(p2);

    REQUIRE(AsyCost{asy_cost_single_node(n2)} == AsyCost{2, 3, 2});
    REQUIRE(AsyCost{asy_cost_single_node(n2.left())} == AsyCost{2, 3, 2});

    auto const p3 =
        parse_expr_antisymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const n3 = eval_node(p3);
    REQUIRE(AsyCost{asy_cost_single_node(n3)} == AsyCost{2, 2, 1});
    REQUIRE(AsyCost{asy_cost_single_node(n3.left())} == AsyCost{2, 3, 1});
  }

  SECTION("asy_cost") {
    auto const p1 =
        parse_expr_antisymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(asy_cost(eval_node(p1)) == AsyCost{2, 2, 3});

    auto const p2 = parse_expr_antisymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");

    auto const np2 = eval_node(p2);
    REQUIRE(asy_cost(np2) == AsyCost{2, 3, 2} + AsyCost{2, 3, 2});

    auto const p3 =
        parse_expr_antisymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const np3 = eval_node(p3);
    REQUIRE(asy_cost(np3) == AsyCost{2, 2, 1} + AsyCost{2, 3, 1});

    auto const s1 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::symm);
    auto const ns1 = eval_node(s1);
    REQUIRE(asy_cost(ns1) == AsyCost{{1, 4}, 2, 2});  // 1/4 * O^2V^2

    auto const s2 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::antisymm);
    auto const ns2 = eval_node(s2);
    REQUIRE(asy_cost(ns2) == AsyCost{{1, 2}, 2, 2});  // 1/2 * O^2V^2

    auto const s3 =
        parse_expr(L"I{i1,i2;a1,a2} + I{i1,i2;a1,a2}", Symmetry::nonsymm);
    auto const ns3 = eval_node(s3);
    REQUIRE(asy_cost(ns3) == AsyCost{2, 2});  // O^2V^2

    auto const p4 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::symm);
    auto const np4 = eval_node(p4);
    REQUIRE(asy_cost(np4) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p5 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::antisymm);
    auto const np5 = eval_node(p5);
    REQUIRE(asy_cost(np5) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p6 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::antisymm);
    auto const np6 = eval_node(p6);
    REQUIRE(asy_cost(np6) == AsyCost{{1, 2}, 2, 4});  // 1/4 * 2 * O^2V^4

    auto const p7 = parse_expr(L"I{i1;a1} * I{i2;a2}", Symmetry::nonsymm);
    auto const np7 = eval_node(p7);
    REQUIRE(asy_cost(np7) == AsyCost{2, 2});  // 1/2 * 2 * O^2V^4

    auto const p8 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::nonsymm);
    auto const np8 = eval_node(p8);
    REQUIRE(asy_cost(np8) == AsyCost{2, 2, 4});  // 2 * O^2V^4
  }
}
