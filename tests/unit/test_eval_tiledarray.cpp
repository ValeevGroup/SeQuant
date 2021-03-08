#include "catch.hpp"

#include <tiledarray.h>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/eval_tiledarray.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>
#include <cstdlib>
#include <string>
#include <vector>

TEST_CASE("TEST_EVALUATE_TILEDARRAY", "[evaluate]") {
  using namespace sequant;
  using eval::antisymmetrize;
  using eval::evaluate;
  using eval::symmetrize;
  using ranges::views::transform;
  using TA::TArrayD;
  using utils::parse_expr;

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

  auto norm = [](TArrayD const& tensor) {
    using ranges::views::iota;

    auto annot =
        eval::detail::ords_to_annot(iota(size_t{0}, tensor.range().rank()));

    return sqrt(tensor(annot).dot(tensor(annot)));
  };

  const size_t nocc = 2, nvirt = 20;
  //
  auto const tr1o = TA::TiledRange1{0, nocc};
  auto const tr1v = TA::TiledRange1{0, nvirt};

  TA::TiledRange tr_ov{tr1o, tr1v};
  TA::TiledRange tr_vo{tr1v, tr1o};
  TA::TiledRange tr_oovv{tr1o, tr1o, tr1v, tr1v};
  TA::TiledRange tr_vvoo{tr1v, tr1v, tr1o, tr1o};

  auto const t_vo = make_rand_tensor(tr_vo);
  auto const t_vvoo = make_rand_tensor(tr_vvoo);
  auto const f_ov = make_rand_tensor(tr_ov);
  auto const g_oovv = make_rand_tensor(tr_oovv);

  auto tensor_yield = [&t_vo, &f_ov, &t_vvoo, &g_oovv](Tensor const& texpr) {
    if (texpr.rank() == 1)
      return texpr.label() == L"t" ? t_vo : f_ov;
    else
      return texpr.label() == L"t" ? t_vvoo : g_oovv;
  };

  SECTION("summation") {
    auto expr1 = parse_expr(L"t_{a1}^{i1} + f_{i1}^{a1}", Symmetry::antisymm);

    auto sum1_man = TArrayD{};
    sum1_man("0,1") = t_vo("0,1") + f_ov("1,0");

    auto sum1_eval = evaluate<TArrayD>(expr1, false, tensor_yield);

    REQUIRE(norm(sum1_man) == Approx(norm(sum1_eval)));

    auto expr2 =
        parse_expr(L"2 * t_{a1}^{i1} + 1.5 * f_{i1}^{a1}", Symmetry::antisymm);

    auto sum2_man = TArrayD{};
    sum2_man("0,1") = 2 * t_vo("0,1") + 1.5 * f_ov("1,0");

    auto sum2_eval = evaluate<TArrayD>(expr2, false, tensor_yield);

    REQUIRE(norm(sum2_man) == Approx(norm(sum2_eval)));
  }

  SECTION("product") {
    auto expr1 = parse_expr(L"1/2.0 * g_{i2,i4}^{a2,a4} * t_{a1,a2}^{i1,i2}",
                            Symmetry::antisymm);
    TArrayD prod1_man{};
    prod1_man("i4,a1,a4,i1") =
        1 / 2.0 * g_oovv("i2,i4,a2,a4") * t_vvoo("a1,a2,i1,i2");
    auto prod1_eval = evaluate<TArrayD>(expr1, false, tensor_yield);

    REQUIRE(norm(prod1_man) == Approx(norm(prod1_eval)));

    auto expr2 = parse_expr(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}",
        Symmetry::antisymm);

    auto prod2_man = TArrayD{};
    prod2_man("a1,a2,i1,i2") = -1 / 4.0 * g_oovv("i3,i4,a3,a4") *
                               t_vvoo("a2,a4,i1,i2") * t_vvoo("a1,a3,i3,i4");
    auto prod2_eval =
        evaluate<TArrayD>(expr2, true, tensor_yield);  // optimize = true

    REQUIRE(norm(prod2_man) == Approx(norm(prod2_eval)));
  }

  SECTION("sum and product") {
    auto expr1 = parse_expr(
        L"-1/4 * g_{i3,i4}^{a3,a4} * t_{a2,a4}^{i1,i2} * t_{a1,a3}^{i3,i4}"
        " + "
        " 1/16 * g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2} ",
        Symmetry::antisymm);

    auto man1 = TArrayD{};
    man1("a1,a2,i1,i2") = -1.0 / 4 * g_oovv("i3,i4,a3,a4") *
                              t_vvoo("a2,a4,i1,i2") * t_vvoo("a1,a3,i3,i4") +
                          1.0 / 16 * g_oovv("i3,i4,a3,a4") *
                              t_vvoo("a1,a2,i3,i4") * t_vvoo("a3,a4,i1,i2");

    auto eval1 = evaluate<TArrayD>(expr1, true, tensor_yield);

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Antisymmetrization") {
    auto const& tnsr1 = g_oovv;
    auto man1 = TArrayD{};
    man1("0,1,2,3") = tnsr1("0,1,2,3") - tnsr1("1,0,2,3") + tnsr1("1,0,3,2") -
                      tnsr1("0,1,3,2");

    auto eval1 = antisymmetrize(tnsr1);

    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }

  SECTION("Symmetrization") {
    auto const& tnsr1 = g_oovv;
    auto man1 = TArrayD{};
    man1("0,1,2,3") = tnsr1("0,1,2,3") + tnsr1("1,0,3,2");

    auto eval1 = symmetrize(tnsr1);
    REQUIRE(norm(man1) == Approx(norm(eval1)));
  }
}
