#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/asy_cost.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <cassert>
#include <initializer_list>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

namespace {

auto eval_node(sequant::ExprPtr const& expr) {
  return sequant::eval_node<sequant::EvalExpr>(expr);
}

enum struct Npos {
  L,  // Left
  R   // Right
};

sequant::EvalExpr node(sequant::EvalNode<sequant::EvalExpr> const& n,
                       std::initializer_list<Npos> pos) {
  if (pos.size() == 0) return *n;
  auto n_ = n;
  for (auto p : pos) {
    assert(!n_.leaf() && "Accessing child of leaf!");
    n_ = p == Npos::L ? n_.left() : n_.right();
  }

  return *n_;
}

[[maybe_unused]] std::wstring tikz(
    sequant::EvalNode<sequant::EvalExpr> const& n) noexcept {
  return n.tikz(
      [](auto&& n) { return L"$" + n->expr()->to_latex() + L"$"; });
}

}  // namespace

TEST_CASE("TEST EVAL_NODE", "[EvalNode]") {
  using namespace sequant;
  auto L = Npos::L;
  auto R = Npos::R;

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

    REQUIRE_THAT(node(node1, {}).as_tensor(), EquivalentTo("I{a1,a2;i1,i2}:A"));

    REQUIRE(node(node1, {R}).as_constant() == Constant{rational{1, 16}});

    REQUIRE_THAT(node(node1, {L}).as_tensor(),
                 EquivalentTo("I{a1,a2;i1,i2}:A"));

    REQUIRE_THAT(node(node1, {L, L}).as_tensor(),
                 EquivalentTo("I{a1,a2;a3,a4}:A"));

    REQUIRE_THAT(node(node1, {L, R}).as_tensor(),
                 EquivalentTo("t{a3,a4;i1,i2}:A"));

    REQUIRE_THAT(node(node1, {L, L, L}).as_tensor(),
                 EquivalentTo("g{i3,i4;a3,a4}:A"));

    REQUIRE_THAT(node(node1, {L, L, R}).as_tensor(),
                 EquivalentTo("t{a1,a2;i3,i4}:A"));

    // 1/16 * A * (B * C)
    auto node2p = Product{p1->as<Product>().scalar(), {}};
    node2p.append(1, p1->at(0), Product::Flatten::No);
    node2p.append(1, ex<Product>(Product{p1->at(1), p1->at(2)}),
                  Product::Flatten::No);

    auto const node2 = eval_node(ex<Product>(node2p));

    REQUIRE_THAT(node(node2, {}).as_tensor(), EquivalentTo("I{a1,a2;i1,i2}:N"));

    REQUIRE_THAT(node(node2, {L}).as_tensor(),
                 EquivalentTo("I{a1,a2;i1,i2}:N"));

    REQUIRE(node(node2, {R}).as_constant() == Constant{rational{1, 16}});

    REQUIRE_THAT(node(node2, {L, L}).as_tensor(),
                 EquivalentTo("g{i3,i4; a3,a4}:A"));

    REQUIRE_THAT(node(node2, {L, R}).as_tensor(),
                 EquivalentTo("I{a1,a2,a3,a4;i3,i4,i1,i2}:A"));

    REQUIRE_THAT(node(node2, {L, R, L}).as_tensor(),
                 EquivalentTo("t{a1,a2;i3,i4}:A"));

    REQUIRE_THAT(node(node2, {L, R, R}).as_tensor(),
                 EquivalentTo("t{a3,a4;i1,i2}:A"));
  }

  SECTION("sum") {
    auto const sum1 = parse_expr_antisymm(
        L"X^{i1,i2}_{a1,a2} "
        L"+ Y^{i1, i2}_{a1,a2}"
        L"+ g_{i3,a1}^{i1,i2} * t_{a2}^{i3}");

    auto const node1 = eval_node(sum1);
    REQUIRE(node1->op_type() == EvalOp::Sum);
    REQUIRE(node1.left()->op_type() == EvalOp::Sum);
    REQUIRE_THAT(node1.left()->as_tensor(), EquivalentTo("I{a1,a2;i1,i2}:A"));
    REQUIRE_THAT(node1.left().left()->as_tensor(),
                 EquivalentTo("X{a1,a2;i1,i2}:A"));
    REQUIRE_THAT(node1.left().right()->as_tensor(),
                 EquivalentTo("Y{a1,a2;i1,i2}:A"));

    REQUIRE(node1.right()->op_type() == EvalOp::Prod);
    REQUIRE_THAT(node1.right()->as_tensor(),
                 EquivalentTo("I{a1,a2;i1,i2}:N-C-N"));
    REQUIRE_THAT(node1.right().left()->as_tensor(),
                 EquivalentTo("g{i3,a1;i1,i2}:A"));
    REQUIRE_THAT(node1.right().right()->as_tensor(),
                 EquivalentTo("t{a2;i3}:A"));
  }

  SECTION("variable") {
    auto prod1 = parse_expr(L"a * b * c");
    auto node1 = eval_node(prod1);

    REQUIRE(node1->op_type() == EvalOp::Prod);
    REQUIRE(node1.left()->op_type() == EvalOp::Prod);

    REQUIRE(node(node1, {}).is_variable());
    REQUIRE(node(node1, {L}).is_variable());
    REQUIRE(node(node1, {R}).is_variable());
    REQUIRE(node(node1, {R}).as_variable() == Variable{L"c"});
    REQUIRE(node(node1, {L, R}).as_variable() == Variable{L"b"});
    REQUIRE(node(node1, {L, L}).as_variable() == Variable{L"a"});

    auto sum1 = parse_expr(L"a + b + c");
    auto node2 = eval_node(sum1);

    REQUIRE(node2->op_type() == EvalOp::Sum);
    REQUIRE(node2.left()->op_type() == EvalOp::Sum);

    REQUIRE(node(node2, {}).is_variable());
    REQUIRE(node(node2, {L}).is_variable());
    REQUIRE(node(node2, {R}).is_variable());
    REQUIRE(node(node2, {R}).as_variable() == Variable{L"c"});
    REQUIRE(node(node2, {L, R}).as_variable() == Variable{L"b"});
    REQUIRE(node(node2, {L, L}).as_variable() == Variable{L"a"});

    auto prod2 = parse_expr(L"a * t{i1;a1}");
    auto node3 = eval_node(prod2);
    REQUIRE_THAT(node(node3, {}).as_tensor(), EquivalentTo("I{i1;a1}"));
    REQUIRE_THAT(node(node3, {R}).as_tensor(), EquivalentTo("t{i1;a1}"));
    REQUIRE(node(node3, {L}).as_variable() == Variable{L"a"});
  }

  SECTION("to_expr") {
    const auto p1 = parse_expr_antisymm(
        L"1/16 "
        L"* g_{i3, i4}^{a3, a4}"
        L"* t_{a1, a2}^{i3, i4}"
        L"* t_{a3,a4}^{i1,i2}");

    auto p1_after = ex<Product>(rational{1, 16}, ExprPtrList{});
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
    REQUIRE(linearize_eval_node(eval_node(p2))->to_latex() == p2->to_latex());
  }

  SECTION("single node") {
    auto asy_cost_single_node = FlopsWithSymm{};
    auto const p1 =
        parse_expr_antisymm(L"g_{i2, a1}^{a2, a3} * t_{a2, a3}^{i1, i2}");
    REQUIRE(asy_cost_single_node(eval_node(p1)) == AsyCost{2, 2, 3});

    auto const p2 = parse_expr_antisymm(
        L"g_{i2,i3}^{a2,a3} * t_{a2}^{i1} * t_{a1,a3}^{i2,i3}");
    auto const n2 = eval_node(p2);

    REQUIRE(asy_cost_single_node(n2) == AsyCost{2, 3, 2});
    REQUIRE(asy_cost_single_node(n2.left()) == AsyCost{2, 3, 2});

    auto const p3 =
        parse_expr_antisymm(L"g_{i2,i3}^{i1,a2} * t_{a2}^{i2} * t_{a1}^{i3}");
    auto const n3 = eval_node(p3);
    REQUIRE(asy_cost_single_node(n3) == AsyCost{2, 2, 1});
    REQUIRE(asy_cost_single_node(n3.left()) == AsyCost{2, 3, 1});
  }

  SECTION("asy_cost") {
    auto asy_cost = [](auto const& n) {
      return sequant::asy_cost(n, FlopsWithSymm{});
    };
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
    REQUIRE(asy_cost(np7) == AsyCost{{1, 2}, 2, 2});  // 1/2 * O^2V^2

    auto const p8 =
        parse_expr(L"I{i1,i2;a3,a4} * I{a3,a4;a1,a2}", Symmetry::nonsymm);
    auto const np8 = eval_node(p8);
    REQUIRE(asy_cost(np8) == AsyCost{2, 2, 4});  // 2 * O^2V^4
  }

  SECTION("minimum storage") {
    auto p1 = parse_expr(L"2 * g{i2,a1;a2,a3} * t{a2,a3;i2,i1}");
    auto const n1 = eval_node(p1);
    // evaluation happens in two steps.
    // g and t are contracted to give an intermediate I{a1;i1}
    // I is then scaled by 2 to give the final result I{a1;i1}
    // - contraction of g and t requires the storage of g, plus that of t
    //   plus that of I. Which is OV^3 + O^2V^2 + OV.
    // - scaling of I requires OV + OV = 2 * OV.
    // The maximum of the two steps is OV^3 + O^2V^2 + OV, which is the minimum
    // storage requirement.
    REQUIRE(min_storage(n1) == AsyCost{1, 3} + AsyCost{2, 2} + AsyCost{1, 1});

    auto p2 = parse_expr(L"1/2 * (g{a1,a2; a3,a4} t{a3;i1}) t{a4;i2}");
    auto const n2 = eval_node(p2);
    REQUIRE(min_storage(n2) == AsyCost{0, 4} + AsyCost{1, 3} + AsyCost{1, 1});
  }
}
