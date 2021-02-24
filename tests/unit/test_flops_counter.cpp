#include "catch.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/binary_node.hpp>
#include <SeQuant/domain/utils/eval_seq.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <SeQuant/domain/utils/flops_counter.hpp>

size_t evaluate_flops(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    size_t nocc, size_t nvirt) {
  return node.evaluate(sequant::utils::flops_counter{nocc, nvirt});
}

size_t evaluate_flops(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    std::pair<size_t, size_t> nov) {
  return evaluate_flops(node, std::get<0>(nov), std::get<1>(nov));
}

TEST_CASE("TEST_OPS_COUNTER", "[flops_counter]") {
  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  using namespace sequant;
  using utils::binarize_expr;
  using utils::eval_expr;
  using utils::eval_seq;
  using utils::parse_expr;

  const std::pair<size_t, size_t> no_lt_nv = {2, 3};
  const std::pair<size_t, size_t> no_gt_nv = {3, 2};
  const std::pair<size_t, size_t> no_eq_nv = {2, 2};

  SECTION("Identity operation") {
    auto t1 = parse_expr(L"g_{i3,i4}^{a3,a4}", Symmetry::nonsymm);
    auto tree = binarize_expr(t1);

    REQUIRE(evaluate_flops(tree, no_lt_nv) == 0);
    REQUIRE(evaluate_flops(tree, no_gt_nv) == 0);
    REQUIRE(evaluate_flops(tree, no_eq_nv) == 0);
  }

  SECTION("Summation operation") {
    const auto srange1 = parse_expr(
        L"I1_{i1,i2}^{a1,a2}"
        " + I2_{i1,i2}^{a1,a2}",
        Symmetry::nonsymm);

    const auto tree1 = binarize_expr(srange1);

    auto flops1 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);
      return nocc * nocc * nvirt * nvirt;
    };

    REQUIRE(evaluate_flops(tree1, no_lt_nv) == flops1(no_lt_nv));
    REQUIRE(evaluate_flops(tree1, no_gt_nv) == flops1(no_gt_nv));
    REQUIRE(evaluate_flops(tree1, no_eq_nv) == flops1(no_eq_nv));

    const auto srange2 = parse_expr(
        L"I1_{i1,i2}^{a1,a2}"
        " + I2_{i1,i2}^{a1,a2}"
        " + I3_{i1,i2}^{a1,a2}",
        Symmetry::nonsymm);

    const auto tree2 = binarize_expr(srange2);

    auto flops2 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);
      return 2 * nocc * nocc * nvirt * nvirt;
    };

    REQUIRE(evaluate_flops(tree2, no_lt_nv) == flops2(no_lt_nv));
    REQUIRE(evaluate_flops(tree2, no_gt_nv) == flops2(no_gt_nv));
    REQUIRE(evaluate_flops(tree2, no_eq_nv) == flops2(no_eq_nv));

    // todo
    // sum with scaling
  }

  SECTION("Product operation") {
    const auto prod1 = parse_expr(
        L"g_{i3,i4}^{a3,a4}"
        " * t_{a1,a2}^{i3,i4}",
        Symmetry::nonsymm);

    auto flops1 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);
      return (size_t)(std::pow(nocc, 2) * std::pow(nvirt, 4));
    };

    const auto tree1 = binarize_expr(prod1);

    REQUIRE(evaluate_flops(tree1, no_lt_nv) == flops1(no_lt_nv));
    REQUIRE(evaluate_flops(tree1, no_gt_nv) == flops1(no_gt_nv));
    REQUIRE(evaluate_flops(tree1, no_eq_nv) == flops1(no_eq_nv));

    auto const prod2_0 = ex<Product>(Product{
        parse_expr(L"g_{i3,i4}^{a3,a4} * t_{a1,a2}^{i3,i4}",
                   Symmetry::nonsymm),  // fac
        parse_expr(L"t_{a3,a4}^{i1,i2}",
                   Symmetry::nonsymm)  // fac
    });

    auto flops2_0 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);

      // prod2_0 implies:
      //
      //   ( g_(i3, i4)^(a3, a4) * t_(a1, a2)^(i3, i4) ) * t_(a3, a4)^(i1, i2)
      //     =>
      //     I_(a1, a2)^(a3, a4)                 : flops nocc^2 * nvirt^4
      //            * t_(a3, a4)^(i1, i2)
      //        => I_(a1, a2)^(i1, i2)           : flops nocc^2 * nvirt^4
      //
      return (size_t)(2 * std::pow(nocc, 2) * std::pow(nvirt, 4));
    };

    const auto tree2_0 = binarize_expr(prod2_0);

    REQUIRE(evaluate_flops(tree2_0, no_lt_nv) == flops2_0(no_lt_nv));
    REQUIRE(evaluate_flops(tree2_0, no_gt_nv) == flops2_0(no_gt_nv));
    REQUIRE(evaluate_flops(tree2_0, no_eq_nv) == flops2_0(no_eq_nv));

    auto prod2_1 = Product{};
    prod2_1.append(parse_expr(L"g_{i3,i4}^{a3,a4}", Symmetry::nonsymm));
    prod2_1.append(parse_expr(L"t_{a1,a2}^{i3,i4} * t_{a3,a4}^{i1,i2}",
                              Symmetry::nonsymm));

    auto flops2_1 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);

      // prod2_1 implies:
      //
      //   g_(i3, i4)^(a3, a4) * ( t_(a1, a2)^(i3, i4) * t_(a3, a4)^(i1, i2) )
      //     => g_(i3, i4)^(a3, a4) * I_(a1,a2,a3,a4)^(i3,i4,i1,i2)
      //                                         : flops nocc^4 * nvirt^4
      //        => I_(a1, a2)^(i1, i2)           : flops nocc^4 * nvirt^4
      //
      return (size_t)(2 * std::pow(nocc, 4) * std::pow(nvirt, 4));
    };

    const auto tree2_1 = binarize_expr(ex<Product>(prod2_1));

    REQUIRE(evaluate_flops(tree2_1, no_lt_nv) == flops2_1(no_lt_nv));
    REQUIRE(evaluate_flops(tree2_1, no_gt_nv) == flops2_1(no_gt_nv));
    REQUIRE(evaluate_flops(tree2_1, no_eq_nv) == flops2_1(no_eq_nv));
  }
}
