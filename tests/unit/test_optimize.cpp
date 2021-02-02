#include "catch.hpp"

#include <SeQuant/domain/factorize/optimize.hpp>
#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/expr_parse.hpp>
#include <iostream>

using namespace sequant;

auto get_tex = [](utils::binary_expr<utils::eval_expr>::node_ptr const& node) {
  return L"\"$" + node->data().seq_expr()->to_latex() + L"$\"";
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
                           L"g_{i3,i4}^{a3,a4}"     // T1
                           " * t_{a1,a2}^{i3,i4}"   // T2
                           " * t_{a3,a4}^{i1,i2}",  // T3
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

    REQUIRE(*utils::binarize_flat_prod{prod1}(opt_seq1) ==
            *result1.optimal_seqs.at(0));

    //
    const auto prod2 = parse_expr(
                           L"   g_{i3,i4}^{a3,a4}"
                           L" * t_{a3,a4}^{i1,i2}"
                           L" * t_{a1}^{i3}"
                           L" * t_{a2}^{i4}",
                           Symmetry::nonsymm)
                           ->as<Product>();

    const auto result2_naive =
        single_term_opt(evxpr_range(prod2), nocc, nvirt,
                        {});  // there will be two degenerate evaluation
                              // sequences for prod2 when no
    REQUIRE(result2_naive.optimal_seqs.size() == 2);

    // now let's reuse intermediates from the optimal sequence of evaluation of
    // prod1, to compute the optimal sequence for prod2
    container::set<size_t> imed_hashes_prod1{};
    yield_hash(result1.optimal_seqs.at(0), imed_hashes_prod1);

    const auto result2_discounted =
        single_term_opt(evxpr_range(prod2), nocc, nvirt, imed_hashes_prod1);
    REQUIRE(result2_discounted.ops < result2_naive.ops);

    // yet another example
    auto prod3 =
        parse_expr(L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3} * t_{a3}^{i4}",
                   Symmetry::antisymm);
    auto prod4 = parse_expr(L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3}",
                            Symmetry::antisymm);
    // we show that two the evaluation trees for prod3
    //  - one: single term optimized on prod3 alone
    //  - two: single term optimized on prod3 with the intermediate from prod4
    //  are not the same.
    auto prod3_sto =
        std::move(*(single_term_opt(prod3->as<Product>(), nocc, nvirt, {})
                        .optimal_seqs.begin()));

    // finding the intermediate from the evaluation tree of prod4
    auto prod4_sto =
        std::move(*(single_term_opt(prod4->as<Product>(), nocc, nvirt, {})
                        .optimal_seqs.begin()));

    auto imeds_prod4 = container::set<size_t>{};
    utils::visit_inorder_binary_expr<utils::eval_expr>(
        prod4_sto, [&imeds_prod4](auto const& n) {
          if (!n->leaf()) imeds_prod4.emplace(n->data().hash());
        });

    auto prod3_sto_with_imeds = std::move(
        *single_term_opt(prod3->as<Product>(), nocc, nvirt, imeds_prod4)
             .optimal_seqs.begin());

    REQUIRE_FALSE(*prod3_sto == *prod3_sto_with_imeds);
  }

  SECTION("Most expensive term") {
    using factorize::most_expensive;
    const container::set<size_t> imed_hashes = {};

    const size_t nocc = 2, nvirt = 3;  // nocc < nvirt

    const auto expr = parse_expr(
        L"f_{i3}^{i1}*t_{a1,a2}^{i2,i3}"
        L"    + "
        L"    g_{a1,a2}^{a3,a4}*t_{a3,a4}^{i1,i2}",
        Symmetry::antisymm);

    auto expensive = most_expensive(*expr, nocc, nvirt, {});

    REQUIRE(*expr->at(1) == *expensive.mets.begin()->first);
  }

  SECTION("Multiple terms Hartono") {
    // a: a1   i: i1
    // b: a2   j: i2
    // c: a3   k: i3
    // d: a4   l: i4

    auto prod1 =
        parse_expr(L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3} * t_{a3}^{i4}",
                   Symmetry::antisymm);
    auto prod2 = parse_expr(L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3}",
                            Symmetry::antisymm);

    auto terms = container::svector<ExprPtr>{prod1, prod2};

    size_t nocc = 2, nvirt = 5;
    auto opt_terms = factorize::multi_term_opt_hartono(terms, nocc, nvirt);

    auto prod2_opt = std::move(
        *(factorize::single_term_opt(prod2->as<Product>(), nocc, nvirt, {})
              .optimal_seqs.begin()));
    auto imeds = container::set<size_t>{};
    utils::visit_inorder_binary_expr<utils::eval_expr>(
        prod2_opt, [&imeds](auto const& n) {
          if (!n->leaf()) imeds.emplace(n->data().hash());
        });

    auto prod1_opt = std::move(
        *(factorize::single_term_opt(prod1->as<Product>(), nocc, nvirt, imeds)
              .optimal_seqs.begin()));

    REQUIRE(**(opt_terms.at(prod2).optimal_seqs.begin()) == *prod2_opt);
    REQUIRE(**(opt_terms.at(prod1).optimal_seqs.begin()) == *prod1_opt);
  }
}
