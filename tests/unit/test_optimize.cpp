#include "catch.hpp"

#include <SeQuant/domain/factorize/optimize.hpp>
#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <iostream>

using namespace sequant;

auto get_tex = [](utils::binary_expr<utils::eval_expr>::node_ptr const& node) {
  return L"$" + node->data().seq_expr()->to_latex() + L"$";
};

auto get_hash = [](utils::binary_expr<utils::eval_expr>::node_ptr const& node) {
  return std::to_wstring(node->data().hash());
};

auto get_hash_and_tex =
    [](utils::binary_expr<utils::eval_expr>::node_ptr const& node) {
      return L"\"" + get_tex(node) + L"\\n" + get_hash(node) + L"\"";
    };

auto evxpr_range = [](const Product& prod) {
  return prod.factors() | ranges::views::transform([](const auto& expr) {
           return utils::eval_expr{expr->template as<Tensor>()};
         }) |
         ranges::to<container::svector<utils::eval_expr>>;
};

auto yield_hash = [](const auto& node, auto& container) {
  utils::visit_inorder_binary_expr<utils::eval_expr>(
      node, [&container](const auto& x) {
        if (!x->leaf()) container.insert(x->data().hash());
      });
};

TEST_CASE("TEST_OPTIMIZE", "[optimize]") {
  using factorize::single_term_opt;
  using utils::parse_expr;

  SECTION("Single term optimization") {
    const container::set<size_t> imed_hashes = {};

    const size_t nocc = 2, nvirt = 3;  // nocc < nvirt

    const auto prod1 = parse_expr(
                           L"g_(i3,i4)^(a3,a4)"     // T1
                           " * t_(a1,a2)^(i3,i4)"   // T2
                           " * t_(a3,a4)^(i1,i2)",  // T3
                           Symmetry::nonsymm)
                           ->as<Product>();
    //
    // Cost of evaluation prod1:
    //
    // ((T1 * T2) * T3)  : 2 * O^2 * V^4  best if nocc > nvirt
    //
    // ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
    //
    // ((T2 * T3) * T1)  : 2 * O^4 * V^4  worst sequence of evaluation
    //

    const auto result1 =
        single_term_opt(evxpr_range(prod1), nocc, nvirt, imed_hashes);

    REQUIRE(result1.ops < std::numeric_limits<size_t>::max());
    for (const auto& x : result1.optimal_seqs) REQUIRE(x);

    // there are no degenerate evaluation sequences
    REQUIRE(result1.optimal_seqs.size() == 1);

    const auto opt_seq1 = utils::eval_sequence<size_t>{0, {2, 1}};

    REQUIRE(*utils::binarize_prod{prod1}(opt_seq1) ==
            *result1.optimal_seqs.at(0));

    //
    const auto prod2 = parse_expr(
                           L"   g_(i3,i4)^(a3,a4)"
                           L" * t_(a3,a4)^(i1,i2)"
                           L" * t_(a1)^(i3)"
                           L" * t_(a2)^(i4)",
                           Symmetry::nonsymm)
                           ->as<Product>();

    const auto result2_naive =
        single_term_opt(evxpr_range(prod2), nocc, nvirt, {});
    // there will be two degenerate evaluation sequences for prod2 when no
    REQUIRE(result2_naive.optimal_seqs.size() == 2);

    // now let's reuse intermediates from the optimal sequence of evaluation of
    // prod1, to compute the optimal sequence for prod2

    container::set<size_t> imed_hashes_prod1{};
    yield_hash(result1.optimal_seqs.at(0), imed_hashes_prod1);

    const auto result2_discounted =
        single_term_opt(evxpr_range(prod2), nocc, nvirt, imed_hashes_prod1);
    REQUIRE(result2_discounted.ops < result2_naive.ops);
  }
}
