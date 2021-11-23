#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>

#include "catch.hpp"

TEST_CASE("TEST_EVAL_EXPR", "[EvalExpr]") {
  using namespace std::string_literals;
  using sequant::EvalExpr;
  using namespace sequant;
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto parse_expr_antisymm = [](auto const& xpr){return parse_expr(xpr, Symmetry::antisymm);};

  SECTION("Constructors") {
    auto t1 = parse_expr_antisymm(L"t_{i1, i2}^{a1, a2}");

    REQUIRE_NOTHROW(EvalExpr{t1->as<sequant::Tensor>()});

    auto p1 = parse_expr_antisymm(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};

    REQUIRE_NOTHROW(EvalExpr{c2, c3, EvalOp::Prod});
  }

  SECTION("EvalExpr::EvalOp types") {
    auto t1 = parse_expr_antisymm(L"t_{i1, i2}^{a1, a2}");

    auto x1 = EvalExpr(t1->as<Tensor>());

    REQUIRE(x1.op() == EvalOp::Id);

    auto p1 = parse_expr_antisymm(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};
    const auto& c4 = EvalExpr{c2, c3, EvalOp::Prod};
    REQUIRE(c4.op() == EvalOp::Prod);

    const auto c5 =
        EvalExpr{parse_expr_antisymm(L"I^{i3,a1}_{i1,i2}")->as<Tensor>()};
    const auto& c6 = EvalExpr{c2, c5, EvalOp::Sum};

    REQUIRE(c6.op() == EvalOp::Sum);

    auto x2 = EvalExpr(parse_expr_antisymm(L"A_{a_1, a_2}^{i_1, i_2}")->as<Tensor>());
    REQUIRE(x2.op() == EvalOp::Id);

    auto p2 = parse_expr_antisymm(L"A{a1,a2;i1,i2} * I{a1,a2;i1,i2}");
    auto const x3 = EvalExpr{p2->at(0)->as<Tensor>()};
    auto const x4 = EvalExpr{p2->at(1)->as<Tensor>()};
    auto const x5 = EvalExpr{x3, x4, EvalOp::Antisymm};

    REQUIRE(x5.op() == EvalOp::Antisymm);

    auto p3 = parse_expr(L"S{a1,a2;i1,i2}:S * I{a1,a2;i1,i2}", Symmetry::nonsymm);
    auto const x6 = EvalExpr{p3->at(0)->as<Tensor>()};
    auto const x7 = EvalExpr{p3->at(1)->as<Tensor>()};
    auto const x8 = EvalExpr{x6, x7, EvalOp::Symm};
    REQUIRE(x8.op() == EvalOp::Symm);
  }

  SECTION("Sequant expression") {
    const auto& str_t1 = L"g_{a1,a2}^{a3,a4}";
    const auto& str_t2 = L"t_{a3,a4}^{i1,i2}";
    const auto& t1 = parse_expr_antisymm(str_t1);

    const auto& t2 = parse_expr_antisymm(str_t2);

    const auto& x1 = EvalExpr{t1->as<Tensor>()};
    const auto& x2 = EvalExpr{t2->as<Tensor>()};

    REQUIRE(*t1 == x1.tensor());
    REQUIRE(*t2 == x2.tensor());

    const auto& x3 = EvalExpr{x1, x2, EvalOp::Prod};

    REQUIRE_NOTHROW(x3.tensor());

    const auto& prod_indices =
        x3.tensor().const_braket() |
        ranges::views::transform([](const auto& x) { return x.label(); }) |
        ranges::to<container::set<std::wstring_view>>;

    const auto& expected_indices =
        std::initializer_list<std::wstring_view>{L"i_1", L"i_2", L"a_1",
                                                 L"a_2"} |
        ranges::to<container::set<std::wstring_view>>;

    REQUIRE(x3.op() == EvalOp::Prod);

    REQUIRE(prod_indices == expected_indices);

    const auto t4 = parse_expr_antisymm(L"g_{i3,i4}^{a3,a4}")->as<Tensor>();
    const auto t5 =
        parse_expr_antisymm(L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}")->as<Tensor>();

    const auto& x45 = EvalExpr{EvalExpr{t4}, EvalExpr{t5}, EvalOp::Prod};
    const auto& x54 = EvalExpr{EvalExpr{t5}, EvalExpr{t4}, EvalOp::Prod};

    REQUIRE(x45.tensor().to_latex() ==
            parse_expr_antisymm(L"I_{a1,a2}^{i1,i2}")->to_latex());
    REQUIRE(x45.tensor().to_latex() == x54.tensor().to_latex());
  }

  SECTION("Hash value") {
    const auto t1 = parse_expr_antisymm(L"t_{i1}^{a1}")->as<Tensor>();
    const auto t2 = parse_expr_antisymm(L"t_{i2}^{a2}")->as<Tensor>();
    const auto t3 = parse_expr_antisymm(L"t_{i1,i2}^{a1,a2}")->as<Tensor>();

    const auto& x1 = EvalExpr{t1};
    const auto& x2 = EvalExpr{t2};

    const auto& x12 = EvalExpr{x1, x2, EvalOp::Prod};
    const auto& x21 = EvalExpr{x2, x1, EvalOp::Prod};

    REQUIRE(x1.hash() == x2.hash());
    REQUIRE(x12.hash() == x21.hash());

    const auto& x3 = EvalExpr{t3};
    // const auto& x123 = EvalExpr{x12, x3};

    REQUIRE_FALSE(x1.hash() == x3.hash());
    REQUIRE_FALSE(x12.hash() == x3.hash());
  }

  SECTION("Symmetry of product") {
    // whole bra <-> ket contraction between two antisymmetric tensors
    const auto t1 = parse_expr_antisymm(L"g_{i3,i4}^{i1,i2}")->as<Tensor>();
    const auto t2 = parse_expr_antisymm(L"t_{a1,a2}^{i3,i4}")->as<Tensor>();

    const auto x12 = EvalExpr{EvalExpr{t1}, EvalExpr{t2}, EvalOp::Prod};

    REQUIRE(x12.tensor().symmetry() == Symmetry::antisymm);

    // whole bra <-> ket contraction between two symmetric tensors
    const auto t3 =
        parse_expr(L"g_{i3,i4}^{i1,i2}", Symmetry::symm)->as<Tensor>();
    const auto t4 =
        parse_expr(L"t_{a1,a2}^{i3,i4}", Symmetry::symm)->as<Tensor>();

    const auto x34 = EvalExpr{EvalExpr{t3}, EvalExpr{t4}, EvalOp::Prod};

    REQUIRE(x34.tensor().symmetry() == Symmetry::symm);

    // outer product of the same tensor
    const auto t5 = parse_expr(L"f_{i1}^{a1}", Symmetry::nonsymm)->as<Tensor>();
    const auto t6 = parse_expr(L"f_{i2}^{a2}", Symmetry::nonsymm)->as<Tensor>();

    const auto& x56 = EvalExpr{EvalExpr{t5}, EvalExpr{t6}, EvalOp::Prod};

    REQUIRE(x56.tensor().symmetry() == Symmetry::antisymm);

    // contraction of some indices from a bra to a ket
    const auto t7 = parse_expr_antisymm(L"g_{a1,a2}^{i1,a3}")->as<Tensor>();
    const auto t8 = parse_expr_antisymm(L"t_{a3}^{i2}")->as<Tensor>();

    const auto x78 = EvalExpr{EvalExpr{t7}, EvalExpr{t8}, EvalOp::Prod};

    REQUIRE(x78.tensor().symmetry() == Symmetry::nonsymm);
  }

  SECTION("Symmetry of sum") {
    auto tensor = [](Symmetry s) {
      return parse_expr(L"I_{i1,i2}^{a1,a2}", s)->as<Tensor>();
    };

    auto symmetry = [](const EvalExpr& x) { return x.tensor().symmetry(); };

    auto imed = [](const Tensor& t1, const Tensor& t2) {
      return EvalExpr{EvalExpr{t1}, EvalExpr{t2}, EvalOp::Sum};
    };

    const auto t1 = tensor(Symmetry::antisymm);
    const auto t2 = tensor(Symmetry::antisymm);

    const auto t3 = tensor(Symmetry::symm);
    const auto t4 = tensor(Symmetry::symm);

    const auto t5 = tensor(Symmetry::nonsymm);
    const auto t6 = tensor(Symmetry::nonsymm);

    // sum of two antisymm tensors.
    REQUIRE(symmetry(imed(t1, t2)) == Symmetry::antisymm);

    // sum of one antisymm and one symmetric tensors
    REQUIRE(symmetry(imed(t1, t3)) == Symmetry::symm);

    // sum of two symmetric tensors
    REQUIRE(symmetry(imed(t3, t4)) == Symmetry::symm);

    // sum of one antisymmetric and one one nonsymmetric tensors
    REQUIRE(symmetry(imed(t1, t5)) == Symmetry::nonsymm);

    // sum of one symmetric and one nonsymmetric tensors
    REQUIRE(symmetry(imed(t3, t5)) == Symmetry::nonsymm);

    // sum of two nonsymmetric tensors
    REQUIRE(symmetry(imed(t5, t6)) == Symmetry::nonsymm);
  }
}
