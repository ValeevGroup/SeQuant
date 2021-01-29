#include "catch.hpp"

#include <SeQuant/core/expr.hpp>
#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/binary_expr.hpp>
#include <SeQuant/domain/utils/eval_sequence.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <SeQuant/domain/utils/flops_counter.hpp>

// std::wostream& operator<<(
//     std::wostream& os,
//     sequant::utils::binary_expr<sequant::utils::eval_expr>::node_ptr const&
//         node) {
//   os << "\n";
//
//   sequant::utils::digraph_binary_expr<sequant::utils::eval_expr>(
//       os, node, [](const auto& x) {
//         return L"\"$" + x->data().seq_expr()->to_latex() + L"$\"";
//       });
//
//   os << "\n";
//   return os;
// }

size_t evaluate_flops(sequant::utils::binary_expr<
                          sequant::utils::eval_expr>::node_ptr const& node,
                      size_t nocc, size_t nvirt) {
  return sequant::utils::evaluate_binary_expr<sequant::utils::eval_expr>(
      node, sequant::utils::flops_counter{nocc, nvirt});
}

size_t evaluate_flops(sequant::utils::binary_expr<
                          sequant::utils::eval_expr>::node_ptr const& node,
                      std::pair<size_t, size_t> nov) {
  return evaluate_flops(node, std::get<0>(nov), std::get<1>(nov));
}

auto tnsr_nsym = [](std::wstring_view spec) {
  return sequant::utils::parse_expr(spec, sequant::Symmetry::nonsymm)
      ->as<sequant::Tensor>();
};

TEST_CASE("TEST_OPS_COUNTER", "[flops_counter]") {
  using namespace sequant;
  using utils::binarize_evxpr_range;
  using utils::binarize_flat_prod;
  using utils::eval_expr;
  using utils::eval_sequence;
  using utils::parse_expr;

  const std::pair<size_t, size_t> no_lt_nv = {2, 3};
  const std::pair<size_t, size_t> no_gt_nv = {3, 2};
  const std::pair<size_t, size_t> no_eq_nv = {2, 2};

  SECTION("Identity operation") {
    const Tensor t1 = tnsr_nsym(L"g_(i3,i4)^(a3,a4)");
    auto tree =
        binarize_evxpr_range(ranges::views::single(eval_expr{eval_expr{t1}}));

    REQUIRE(evaluate_flops(tree, no_lt_nv) == 0);
    REQUIRE(evaluate_flops(tree, no_gt_nv) == 0);
    REQUIRE(evaluate_flops(tree, no_eq_nv) == 0);
  }

  // SECTION("Scaling operation") {
  //   const auto prod = parse_expr(L"1/2 * g_(i1,i2)^(a1,a2)",
  //   Symmetry::antisymm)
  //                         ->as<Product>();
  //
  //   const auto tree = binarize_flat_prod{prod}(eval_sequence<size_t>{0});
  //
  //   auto flops = [](const std::pair<size_t, size_t>& ov) {
  //     size_t nocc = std::get<0>(ov);
  //     size_t nvirt = std::get<1>(ov);
  //     return nocc * nocc * nvirt * nvirt;
  //   };
  //
  //   REQUIRE(evaluate_flops(tree, no_lt_nv) == flops(no_lt_nv));
  //   REQUIRE(evaluate_flops(tree, no_gt_nv) == flops(no_gt_nv));
  //   REQUIRE(evaluate_flops(tree, no_eq_nv) == flops(no_eq_nv));
  // }

  SECTION("Summation operation") {
    auto seq_node = [](std::wstring_view spec) {
      return eval_expr{tnsr_nsym(spec)};
    };

    const auto srange1 = {seq_node(L"I1_(i1,i2)^(a1,a2)"),
                          seq_node(L"I2_(i1,i2)^(a1,a2)")};

    const auto tree1 = binarize_evxpr_range(srange1);

    auto flops1 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);
      return nocc * nocc * nvirt * nvirt;
    };

    REQUIRE(evaluate_flops(tree1, no_lt_nv) == flops1(no_lt_nv));
    REQUIRE(evaluate_flops(tree1, no_gt_nv) == flops1(no_gt_nv));
    REQUIRE(evaluate_flops(tree1, no_eq_nv) == flops1(no_eq_nv));

    const auto srange2 = {seq_node(L"I1_(i1,i2)^(a1,a2)"),
                          seq_node(L"I2_(i1,i2)^(a1,a2)"),
                          seq_node(L"I3_(i1,i2)^(a1,a2)")};

    const auto tree2 = binarize_evxpr_range(srange2);

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
    using seq_t = eval_sequence<size_t>;

    const auto prod1 = parse_expr(
                           L"g_(i3,i4)^(a3,a4)"
                           " * t_(a1,a2)^(i3,i4)",
                           Symmetry::nonsymm)
                           ->as<Product>();

    auto flops1 = [](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);
      return (size_t)(std::pow(nocc, 2) * std::pow(nvirt, 4));
    };

    const auto tree1 = binarize_flat_prod{prod1}(seq_t{0, {1}});

    REQUIRE(evaluate_flops(tree1, no_lt_nv) == flops1(no_lt_nv));
    REQUIRE(evaluate_flops(tree1, no_gt_nv) == flops1(no_gt_nv));
    REQUIRE(evaluate_flops(tree1, no_eq_nv) == flops1(no_eq_nv));

    const auto prod2 = parse_expr(
                           L"g_(i3,i4)^(a3,a4)"
                           " * t_(a1,a2)^(i3,i4)"
                           " * t_(a3,a4)^(i1,i2)",
                           Symmetry::nonsymm)
                           ->as<Product>();

    const auto seq2_0 = seq_t{0, {1, 2}};

    auto flops2_0 = [&flops1](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);

      // seq2_0 implies:
      //
      //   ( g_(i3, i4)^(a3, a4) * t_(a1, a2)^(i3, i4) ) * t_(a3, a4)^(i1, i2)
      //     =>
      //     I_(a1, a2)^(a3, a4)                 : flops nocc^2 * nvirt^4
      //            * t_(a3, a4)^(i1, i2)
      //        => I_(a1, a2)^(i1, i2)           : flops nocc^2 * nvirt^4
      //
      return (size_t)(2 * std::pow(nocc, 2) * std::pow(nvirt, 4));
    };

    const auto tree2_0 = binarize_flat_prod{prod2}(seq2_0);

    REQUIRE(evaluate_flops(tree2_0, no_lt_nv) == flops2_0(no_lt_nv));
    REQUIRE(evaluate_flops(tree2_0, no_gt_nv) == flops2_0(no_gt_nv));
    REQUIRE(evaluate_flops(tree2_0, no_eq_nv) == flops2_0(no_eq_nv));

    const auto seq2_1 = seq_t{0, {seq_t{1, {2}}}};

    auto flops2_1 = [&flops1](const std::pair<size_t, size_t>& ov) {
      size_t nocc = std::get<0>(ov);
      size_t nvirt = std::get<1>(ov);

      // seq2_1 implies:
      //
      //   g_(i3, i4)^(a3, a4) * ( t_(a1, a2)^(i3, i4) * t_(a3, a4)^(i1, i2) )
      //     => g_(i3, i4)^(a3, a4) * I_(a1,a2,a3,a4)^(i3,i4,i1,i2)
      //                                         : flops nocc^4 * nvirt^4
      //        => I_(a1, a2)^(i1, i2)           : flops nocc^4 * nvirt^4
      //
      return (size_t)(2 * std::pow(nocc, 4) * std::pow(nvirt, 4));
    };

    const auto tree2_1 = binarize_flat_prod{prod2}(seq2_1);

    REQUIRE(evaluate_flops(tree2_1, no_lt_nv) == flops2_1(no_lt_nv));
    REQUIRE(evaluate_flops(tree2_1, no_gt_nv) == flops2_1(no_gt_nv));
    REQUIRE(evaluate_flops(tree2_1, no_eq_nv) == flops2_1(no_eq_nv));
  }
}
