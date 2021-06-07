#include "catch.hpp"

#include <tiledarray.h>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval.hpp>
#include <SeQuant/domain/eval/eval_ta.hpp>

#include <cstdlib>
#include <string>
#include <vector>

template <typename T>
struct tensor_yield {
  T const& t_vo;
  T const& t_vvoo;
  T const& f_ov;
  T const& g_oovv;
  T const& operator()(sequant::Tensor const& texpr) {
    if (texpr.rank() == 1)
      return texpr.label() == L"t" ? t_vo : f_ov;
    else
      return texpr.label() == L"t" ? t_vvoo : g_oovv;
  }
};

TEST_CASE("TEST_EVAL_USING_TA", "[eval]") {
  using ranges::views::transform;
  using sequant::parse_expr_asymm;
  using sequant::to_eval_node;
  using TA::TArrayD;

  // tnsr is assumed to be single tiled
  auto norm = [](TArrayD const& tnsr) { return tnsr.find(0).get().norm(); };

  std::srand(2021);

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  auto& world = TA::get_default_world();

  auto make_rand_tensor = [&world](auto const& tr) {
    auto res = TArrayD{world, tr};
    // res.fill_random();
    res.init_elements([](auto const&) {
      return static_cast<double>(std::rand()) / RAND_MAX;
    });
    return res;
  };

  const size_t nocc = 2, nvirt = 20;
  //
  auto const tr1o = TA::TiledRange1{0, nocc};
  auto const tr1v = TA::TiledRange1{0, nvirt};

  // clang-format off
  TA::TiledRange tr_ov{tr1o, tr1v},
                 tr_vo{tr1v, tr1o},
                 tr_oovv{tr1o, tr1o, tr1v, tr1v},
                 tr_vvoo{tr1v, tr1v, tr1o, tr1o};
  // clang-format on

  auto const t_vo = make_rand_tensor(tr_vo);
  auto const t_vvoo = make_rand_tensor(tr_vvoo);
  auto const f_ov = make_rand_tensor(tr_ov);
  auto const g_oovv = make_rand_tensor(tr_oovv);

  auto yielder = tensor_yield<TArrayD>{t_vo, t_vvoo, f_ov, g_oovv};

  // nominal(empty) cache manager
  auto manager = sequant::utils::cache_manager<TArrayD>{{}, {}};

  auto eval_bnode = [&yielder, &manager](auto const& node) {
    auto inst = sequant::eval::eval_instance_ta{node};
    return inst.evaluate(yielder, manager);
  };

  auto eval_bnode_sym = [&yielder, &manager](auto const& node) {
    auto inst = sequant::eval::eval_instance_ta{node};
    return inst.evaluate_symm(yielder, manager);
  };

  auto eval_bnode_asym = [&yielder, &manager](auto const& node) {
    auto inst = sequant::eval::eval_instance_ta{node};
    return inst.evaluate_asymm(yielder, manager);
  };

  SECTION("summation") {
    auto expr1 = parse_expr_asymm(L"t_{a1}^{i1} + f_{i1}^{a1}");

    auto sum1_man = TArrayD{};
    sum1_man("0,1") = t_vo("0,1") + f_ov("1,0");

    auto sum1_eval = eval_bnode(to_eval_node(expr1));

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 = parse_expr_asymm(L"2 * t_{a1}^{i1} + 1.5 * f_{i1}^{a1}");

    auto sum2_man = TArrayD{};
    sum2_man("0,1") = 2 * t_vo("0,1") + 1.5 * f_ov("1,0");

    auto sum2_eval = eval_bnode(to_eval_node(expr2));

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 =
        parse_expr_asymm(L"1/2.0 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}");
    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") =
        1 / 2.0 * g_oovv("i2,i4,a2,a4") * t_vvoo("a1,a2,i1,i2");
    auto prod1_eval = eval_bnode(to_eval_node(expr1));

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_expr_asymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}");

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 * g_oovv("i3,i4,a3,a4") *
                               t_vvoo("a2,a4,i1,i2") * t_vvoo("a1,a3,i3,i4");

    auto prod2_eval = eval_bnode(to_eval_node(expr2));

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
  }

  SECTION("sum and product") {
    auto expr1 = parse_expr_asymm(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2} ");

    auto man1 = TArrayD{};
    man1("a1,a2,i1,i2") = -1.0 / 4 * g_oovv("i3,i4,a3,a4") *
                              t_vvoo("a2,a4,i1,i2") * t_vvoo("a1,a3,i3,i4") +
                          1.0 / 16 * g_oovv("i3,i4,a3,a4") *
                              t_vvoo("a1,a2,i3,i4") * t_vvoo("a3,a4,i1,i2");

    auto eval1 = eval_bnode(to_eval_node(expr1));

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Antisymmetrization") {
    auto man1 = TArrayD{};
    man1("0,1,2,3") = g_oovv("0,1,2,3") - g_oovv("1,0,2,3") +
                      g_oovv("1,0,3,2") - g_oovv("0,1,3,2");

    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    auto expr1 = parse_expr_asymm(L"0.5 * g_{i1, i2}^{a1, a2}");

    auto eval1 = eval_bnode_asym(to_eval_node(expr1));

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Symmetrization") {
    auto man1 = TArrayD{};
    man1("0,1,2,3") = g_oovv("0,1,2,3") + g_oovv("1,0,3,2");
    man1("0,1,2,3") = 0.5 * man1("0,1,2,3");

    auto expr1 = parse_expr_asymm(L"0.5 * g_{i1, i2}^{a1, a2}");

    auto eval1 = eval_bnode_sym(to_eval_node(expr1));

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }
}
