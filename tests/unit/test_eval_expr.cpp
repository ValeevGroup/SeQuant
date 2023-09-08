#include "catch.hpp"

#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>

#include <string_view>
#include <initializer_list>
#include <memory>
#include <string>

#include <range/v3/all.hpp>

template <typename Iterable1, typename Iterable2>
bool same_index_labels(Iterable1 const& indices1,
                       Iterable2 const& indices2) noexcept {
  return ranges::equal(indices1, indices2,
                       [](auto const& idx1, auto const& idx2) {
                         return idx1.label() == idx2.label();
                       });
}

namespace sequant {
Tensor parse_tensor(std::wstring_view tnsr, Symmetry s = Symmetry::nonsymm) {
  return parse_expr(tnsr, s)->as<Tensor>();
}

Constant parse_constant(std::wstring_view c) {
  return parse_expr(c)->as<Constant>();
}

// Variable parse_variable(std::string_view v) {
//   //
// }
}  // namespace sequant

TEST_CASE("TEST_EVAL_EXPR", "[EvalExpr]") {
  using namespace std::string_literals;
  using sequant::EvalExpr;
  using namespace sequant;
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("Constructors") {
    auto t1 = parse_tensor(L"t_{i1, i2}^{a1, a2}");

    REQUIRE_NOTHROW(EvalExpr{t1});

    auto p1 = parse_expr(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};

    REQUIRE_NOTHROW(EvalExpr{c2, c3, EvalOp::Prod});

    REQUIRE_NOTHROW(EvalExpr{Variable{L"λ"}});
  }

  SECTION("EvalExpr::EvalOp types") {
    auto t1 = parse_tensor(L"t_{i1, i2}^{a1, a2}");

    auto x1 = EvalExpr(t1);

    REQUIRE(x1.op_type() == EvalOp::Id);

    auto p1 = parse_expr(L"g_{i3,a1}^{i1,i2} * t_{a2}^{a3}");

    const auto& c2 = EvalExpr{p1->at(0)->as<Tensor>()};
    const auto& c3 = EvalExpr{p1->at(1)->as<Tensor>()};
    const auto& c4 = EvalExpr{c2, c3, EvalOp::Prod};
    REQUIRE(c4.op_type() == EvalOp::Prod);

    const auto c5 = EvalExpr{parse_tensor(L"I^{i3,a1}_{i1,i2}")};
    const auto& c6 = EvalExpr{c2, c5, EvalOp::Sum};

    REQUIRE(c6.op_type() == EvalOp::Sum);

    auto x2 = EvalExpr(parse_expr(L"1/2")->as<Constant>());
    REQUIRE(x2.op_type() == EvalOp::Id);

    REQUIRE(EvalExpr{Variable{L"λ"}}.op_type() == EvalOp::Id);
  }

  SECTION("ResultType types") {
    auto T = [](std::wstring_view xpr) { return EvalExpr{parse_tensor(xpr)}; };

    auto C = [](std::wstring_view xpr) {
      return EvalExpr{parse_constant(xpr)};
    };

    auto result_type = [](EvalExpr const& left,   //
                          EvalExpr const& right,  //
                          EvalOp op) -> ResultType {
      return EvalExpr{left, right, op}.result_type();
    };

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                T(L"Y{i1;a1}"),  //
                EvalOp::Sum      //
                ) == ResultType::Tensor);

    auto x = EvalExpr{T(L"X{i1;a1}"), T(L"Y{a1;i1}"), EvalOp::Prod};

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                T(L"Y{a1;i1}"),  //
                EvalOp::Prod     //
                ) == ResultType::Scalar);

    REQUIRE(result_type(                //
                T(L"X{i1,i2; a3,a4}"),  //
                T(L"Y{a3,a4; a1,a2}"),  //
                EvalOp::Prod            //
                ) == ResultType::Tensor);

    REQUIRE(result_type(         //
                T(L"X{i1;a1}"),  //
                C(L"2.5"),       //
                EvalOp::Prod     //
                ) == ResultType::Tensor);

    REQUIRE(result_type(      //
                C(L"1.5"),    //
                C(L"2.5"),    //
                EvalOp::Prod  //
                ) == ResultType::Scalar);

    REQUIRE(result_type(     //
                C(L"1.5"),   //
                C(L"2.5"),   //
                EvalOp::Sum  //
                ) == ResultType::Scalar);
  }

  SECTION("Sequant expression") {
    const auto& str_t1 = L"g_{a1,a2}^{a3,a4}";
    const auto& str_t2 = L"t_{a3,a4}^{i1,i2}";
    const auto& t1 = parse_expr(str_t1);

    const auto& t2 = parse_expr(str_t2);

    const auto& x1 = EvalExpr{t1->as<Tensor>()};
    const auto& x2 = EvalExpr{t2->as<Tensor>()};

    REQUIRE(*t1 == x1.expr()->as<Tensor>());
    REQUIRE(*t2 == x2.expr()->as<Tensor>());

    const auto& x3 = EvalExpr{x1, x2, EvalOp::Prod};

    REQUIRE_NOTHROW(x3.expr()->as<Tensor>());

    const auto& prod_indices =
        x3.expr()->as<Tensor>().const_braket() |
        ranges::views::transform([](const auto& x) { return x.label(); }) |
        ranges::to<container::set<std::wstring_view>>;

    const auto& expected_indices =
        std::initializer_list<std::wstring_view>{L"i_1", L"i_2", L"a_1",
                                                 L"a_2"} |
        ranges::to<container::set<std::wstring_view>>;

    REQUIRE(x3.op_type() == EvalOp::Prod);

    REQUIRE(prod_indices == expected_indices);

    const auto t4 = parse_tensor(L"g_{i3,i4}^{a3,a4}");
    const auto t5 = parse_tensor(L"I_{a1,a2,a3,a4}^{i1,i2,i3,i4}");

    const auto& x45 = EvalExpr{EvalExpr{t4}, EvalExpr{t5}, EvalOp::Prod};
    const auto& x54 = EvalExpr{EvalExpr{t5}, EvalExpr{t4}, EvalOp::Prod};

    REQUIRE(x45.to_latex() == parse_expr(L"I_{a1,a2}^{i1,i2}")->to_latex());
    REQUIRE(x45.to_latex() == x54.to_latex());
  }

  SECTION("Hash value") {
    const auto t1 = parse_tensor(L"t_{i1}^{a1}", Symmetry::antisymm);
    const auto t2 = parse_tensor(L"t_{i2}^{a2}", Symmetry::antisymm);
    const auto t3 = parse_tensor(L"t_{i1,i2}^{a1,a2}", Symmetry::antisymm);

    const auto& x1 = EvalExpr{t1};
    const auto& x2 = EvalExpr{t2};

    const auto& x12 = EvalExpr{x1, x2, EvalOp::Prod};
    const auto& x21 = EvalExpr{x2, x1, EvalOp::Prod};

    REQUIRE(x1.hash_value() == x2.hash_value());
    REQUIRE(x12.hash_value() == x21.hash_value());

    const auto& x3 = EvalExpr{t3};

    REQUIRE_FALSE(x1.hash_value() == x3.hash_value());
    REQUIRE_FALSE(x12.hash_value() == x3.hash_value());
  }

  SECTION("Symmetry of product") {
    // whole bra <-> ket contraction between two antisymmetric tensors
    const auto t1 = parse_tensor(L"g_{i3,i4}^{i1,i2}", Symmetry::antisymm);
    const auto t2 = parse_tensor(L"t_{a1,a2}^{i3,i4}", Symmetry::antisymm);

    const auto x12 = EvalExpr{EvalExpr{t1}, EvalExpr{t2}, EvalOp::Prod};

    REQUIRE(x12.expr()->as<Tensor>().symmetry() == Symmetry::antisymm);

    // whole bra <-> ket contraction between two symmetric tensors
    const auto t3 =
        parse_expr(L"g_{i3,i4}^{i1,i2}", Symmetry::symm)->as<Tensor>();
    const auto t4 =
        parse_expr(L"t_{a1,a2}^{i3,i4}", Symmetry::symm)->as<Tensor>();

    const auto x34 = EvalExpr{EvalExpr{t3}, EvalExpr{t4}, EvalOp::Prod};

    REQUIRE(x34.expr()->as<Tensor>().symmetry() == Symmetry::symm);

    // outer product of the same tensor
    const auto t5 = parse_expr(L"f_{i1}^{a1}", Symmetry::nonsymm)->as<Tensor>();
    const auto t6 = parse_expr(L"f_{i2}^{a2}", Symmetry::nonsymm)->as<Tensor>();

    const auto& x56 = EvalExpr{EvalExpr{t5}, EvalExpr{t6}, EvalOp::Prod};

    REQUIRE(x56.expr()->as<Tensor>().symmetry() == Symmetry::antisymm);

    // contraction of some indices from a bra to a ket
    const auto t7 = parse_tensor(L"g_{a1,a2}^{i1,a3}", Symmetry::antisymm);
    const auto t8 = parse_tensor(L"t_{a3}^{i2}", Symmetry::antisymm);

    const auto x78 = EvalExpr{EvalExpr{t7}, EvalExpr{t8}, EvalOp::Prod};

    REQUIRE(x78.expr()->as<Tensor>().symmetry() == Symmetry::nonsymm);

    // whole bra <-> ket contraction between symmetric and antisymmetric tensors
    auto const t9 =
        parse_expr(L"g_{a1,a2}^{a3,a4}", Symmetry::antisymm)->as<Tensor>();
    auto const t10 =
        parse_expr(L"t_{a3,a4}^{i1,i2}", Symmetry::symm)->as<Tensor>();
    auto const x910 = EvalExpr{EvalExpr{t9}, EvalExpr{t10}, EvalOp::Prod};
    REQUIRE(x910.expr()->as<Tensor>().symmetry() == Symmetry::symm);
  }

  SECTION("Symmetry of sum") {
    auto tensor = [](Symmetry s) {
      return parse_expr(L"I_{i1,i2}^{a1,a2}", s)->as<Tensor>();
    };

    auto symmetry = [](const EvalExpr& x) {
      return x.expr()->as<Tensor>().symmetry();
    };

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

    // sum of an antisymmetric and a nonsymmetric tensors
    REQUIRE(symmetry(imed(t1, t5)) == Symmetry::nonsymm);

    // sum of one symmetric and one nonsymmetric tensors
    REQUIRE(symmetry(imed(t3, t5)) == Symmetry::nonsymm);

    // sum of two nonsymmetric tensors
    REQUIRE(symmetry(imed(t5, t6)) == Symmetry::nonsymm);
  }

  SECTION("Inner and outer indices") {
    auto same_inner_outer = [](Tensor const& t, IndexList inner,
                               IndexList outer) -> bool {
      auto [in, out] = EvalExpr{t}.inner_outer_indices();
      return same_index_labels(in, inner) && same_index_labels(out, outer);
    };

    auto c1 = parse_expr(L"3/2")->as<Constant>();

    REQUIRE(EvalExpr{c1}.inner_outer_indices().inner.empty());
    REQUIRE(EvalExpr{c1}.inner_outer_indices().outer.empty());

    auto f1 = parse_expr(L"f{i2;a2}")->as<Tensor>();
    REQUIRE(same_inner_outer(f1,                              //
                             {Index{L"i_2"}, Index{L"a_2"}},  // inner
                             {}                               // outer
                             ));

    auto f2 = parse_expr(L"f{i2;a2<i1,i2>}")->as<Tensor>();
    REQUIRE(same_inner_outer(f2,                             //
                             {Index{L"a_2"}},                // inner
                             {Index{L"i_1"}, Index{L"i_2"}}  // outer
                             ));

    auto t1 = parse_expr(L"t{a2,a3;i1,i2}")->as<Tensor>();
    REQUIRE(same_inner_outer(
        t1,                                                            //
        {Index{L"a_2"}, Index{L"a_3"}, Index{L"i_1"}, Index{L"i_2"}},  // inner
        {}                                                             // outer
        ));
    auto t2 = parse_expr(L"t{a2<i1,i2>,a3<i1,i2>;i1,i2}")->as<Tensor>();
    REQUIRE(same_inner_outer(t2,                              //
                             {Index{L"a_2"}, Index{L"a_3"}},  // inner
                             {Index{L"i_1"}, Index{L"i_2"}}   // outer
                             ));
  }

  SECTION("Debug") {
    auto t1 =
        EvalExpr{parse_expr(L"O{a_1<i_1,i_2>;a_1<i_3,i_2>}", Symmetry::nonsymm)
                     ->as<Tensor>()};
    auto t2 =
        EvalExpr{parse_expr(L"O{a_2<i_1,i_2>;a_2<i_3,i_2>}", Symmetry::nonsymm)
                     ->as<Tensor>()};

    REQUIRE_NOTHROW(EvalExpr{t1, t2, EvalOp::Prod});
  }
}
