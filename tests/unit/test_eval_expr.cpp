#include <SeQuant/core/sequant.hpp>
#include <SeQuant/domain/utils/eval_expr.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>

#include "catch.hpp"

TEST_CASE("TEST_EVAL_EXPR", "[eval_expr]") {
  using namespace std::string_literals;
  using sequant::utils::eval_expr;
  using sequant::utils::parse_expr;
  using namespace sequant;
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto ex_asym = [](std::wstring_view spec) {
    return utils::parse_expr(spec, Symmetry::antisymm);
  };

  SECTION("Constructors") {
    auto t1 = ex_asym(L"t_{i1, i2}^{a1, a2}");

    REQUIRE_NOTHROW(eval_expr{t1->as<sequant::Tensor>()});

    auto p1 = ex_asym(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = eval_expr{p1->at(0)->as<Tensor>()};
    const auto& c3 = eval_expr{p1->at(1)->as<Tensor>()};

    REQUIRE_NOTHROW(eval_expr{c2, c3});
  }

  SECTION("Eval ops_met") {
    auto t1 = ex_asym(L"t_{i1, i2}^{a1, a2}");

    auto x1 = eval_expr(t1->as<Tensor>());

    REQUIRE(x1.op() == eval_expr::eval_op::Id);

    auto p1 = ex_asym(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = eval_expr{p1->at(0)->as<Tensor>()};
    const auto& c3 = eval_expr{p1->at(1)->as<Tensor>()};
    const auto& c4 = eval_expr{c2, c3};
    REQUIRE(c4.op() == eval_expr::eval_op::Prod);

    const auto c5 = eval_expr{ex_asym(L"I^{i3,a1}_{i1,i2}")->as<Tensor>()};
    const auto& c6 = eval_expr{c2, c5};

    REQUIRE(c6.op() == eval_expr::eval_op::Sum);
  }

  SECTION("Sequant expression") {
    const auto& str_t1 = L"g_{a1,a2}^{a3,a4}";
    const auto& str_t2 = L"t_{a3,a4}^{i1,i2}";
    const auto& t1 = ex_asym(str_t1);

    const auto& t2 = ex_asym(str_t2);

    const auto& x1 = eval_expr{t1->as<Tensor>()};
    const auto& x2 = eval_expr{t2->as<Tensor>()};

    REQUIRE(*t1 == x1.tensor());
    REQUIRE(*t2 == x2.tensor());

    const auto& x3 = eval_expr{x1, x2};

    REQUIRE_NOTHROW(x3.tensor());

    const auto& prod_indices =
        x3.tensor().const_braket() |
        ranges::views::transform([](const auto& x) { return x.label(); }) |
        ranges::to<container::set<std::wstring_view>>;

    const auto& expected_indices =
        std::initializer_list<std::wstring_view>{L"i_1", L"i_2", L"a_1",
                                                 L"a_2"} |
        ranges::to<container::set<std::wstring_view>>;

    REQUIRE(x3.op() == eval_expr::eval_op::Prod);

    REQUIRE(prod_indices == expected_indices);

    const auto t4 = ex_asym(L"g_{i3,i4}^{a3,a4}")->as<Tensor>();
    const auto t5 = ex_asym(L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}")->as<Tensor>();

    const auto& x45 = eval_expr{eval_expr{t4}, eval_expr{t5}};
    const auto& x54 = eval_expr{eval_expr{t5}, eval_expr{t4}};

    REQUIRE(x45.tensor().to_latex() ==
            ex_asym(L"I_{a1,a2}^{i1,i2}")->to_latex());
    REQUIRE(x45.tensor().to_latex() == x54.tensor().to_latex());
  }

  SECTION("Hash value") {
    const auto t1 = ex_asym(L"t_{i1}^{a1}")->as<Tensor>();
    const auto t2 = ex_asym(L"t_{i2}^{a2}")->as<Tensor>();
    const auto t3 = ex_asym(L"t_{i1,i2}^{a1,a2}")->as<Tensor>();

    const auto& x1 = eval_expr{t1};
    const auto& x2 = eval_expr{t2};

    const auto& x12 = eval_expr{x1, x2};
    const auto& x21 = eval_expr{x2, x1};

    REQUIRE(x1.hash() == x2.hash());
    REQUIRE(x12.hash() == x21.hash());

    const auto& x3 = eval_expr{t3};
    const auto& x123 = eval_expr{x12, x3};

    REQUIRE_FALSE(x1.hash() == x3.hash());
    REQUIRE_FALSE(x12.hash() == x3.hash());
  }

  SECTION("Symmetry of product") {
    // whole bra <-> ket contraction between two antisymmetric tensors
    const auto t1 = ex_asym(L"g_{i3,i4}^{i1,i2}")->as<Tensor>();
    const auto t2 = ex_asym(L"t_{a1,a2}^{i3,i4}")->as<Tensor>();

    const auto x12 = eval_expr{eval_expr{t1}, eval_expr{t2}};

    REQUIRE(x12.tensor().symmetry() == Symmetry::antisymm);

    // whole bra <-> ket contraction between two symmetric tensors
    const auto t3 =
        parse_expr(L"g_{i3,i4}^{i1,i2}", Symmetry::symm)->as<Tensor>();
    const auto t4 =
        parse_expr(L"t_{a1,a2}^{i3,i4}", Symmetry::symm)->as<Tensor>();

    const auto x34 = eval_expr{eval_expr{t3}, eval_expr{t4}};

    REQUIRE(x34.tensor().symmetry() == Symmetry::symm);

    // outer product of the same tensor
    const auto t5 = parse_expr(L"f_{i1}^{a1}", Symmetry::nonsymm)->as<Tensor>();
    const auto t6 = parse_expr(L"f_{i2}^{a2}", Symmetry::nonsymm)->as<Tensor>();

    const auto& x56 = eval_expr{eval_expr{t5}, eval_expr{t6}};

    REQUIRE(x56.tensor().symmetry() == Symmetry::symm);

    // contraction of some indices from a bra to a ket
    const auto t7 = ex_asym(L"g_{a1,a2}^{i1,a3}")->as<Tensor>();
    const auto t8 = ex_asym(L"t_{a3}^{i2}")->as<Tensor>();

    const auto x78 = eval_expr{eval_expr{t7}, eval_expr{t8}};

    REQUIRE(x78.tensor().symmetry() == Symmetry::nonsymm);
  }

  SECTION("Symmetry of sum") {
    auto tensor = [](Symmetry s) {
      return parse_expr(L"I_{i1,i2}^{a1,a2}", s)->as<Tensor>();
    };

    auto symmetry = [](const eval_expr& x) { return x.tensor().symmetry(); };

    auto imed = [](const Tensor& t1, const Tensor& t2) {
      return eval_expr{eval_expr{t1}, eval_expr{t2}};
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

  SECTION("Canonicalization") {
    auto evxpr1 = eval_expr(
        parse_expr(L"g_{i1,i2}^{a1,a2}", Symmetry::antisymm)->as<Tensor>());
    auto evxpr2 = eval_expr(
        parse_expr(L"g_{i2,i1}^{a1,a2}", Symmetry::antisymm)->as<Tensor>());

    REQUIRE(evxpr1.tensor() == evxpr2.tensor());
    REQUIRE_FALSE(evxpr1.scalar() == evxpr2.scalar());
    REQUIRE(evxpr1.hash() == evxpr2.hash());
  }
}
