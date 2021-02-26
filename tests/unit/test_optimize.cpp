#include "catch.hpp"

#include <SeQuant/domain/optimize/optimize.hpp>
#include <SeQuant/domain/utils/binarize_expr.hpp>
#include <SeQuant/domain/utils/parse_expr.hpp>
#include <iostream>

using namespace sequant;

auto yield_interm_hash = [](utils::binary_node<utils::eval_expr> const& node,
                            auto& container) {
  node.visit_internal(
      [&container](const auto& x) { container.emplace(x.hash()); });
};

TEST_CASE("TEST_OPTIMIZE", "[optimize]") {
  using optimize::single_term_opt;
  using utils::binarize_expr;
  using utils::parse_expr;

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("Single term optimization") {
    const container::set<size_t> imed_hashes = {};

    const size_t nocc = 2, nvirt = 5;  // nocc < nvirt

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

    // canonicalization set false as it is unnecessary here
    const auto result1 =
        single_term_opt(prod1, nocc, nvirt, imed_hashes, false);

    REQUIRE(result1.ops < std::numeric_limits<size_t>::max());

    // there are no degenerate evaluation sequences
    REQUIRE(result1.optimal_seqs.size() == 1);

    auto prod1_opt =
        ex<Product>(Product{prod1.at(0), prod1.at(2), prod1.at(1)});

    REQUIRE(binarize_expr(prod1_opt) == result1.optimal_seqs.at(0));

    //
    const auto prod2 = parse_expr(
                           L"   g_{i3,i4}^{a3,a4}"
                           L" * t_{a3,a4}^{i1,i2}"
                           L" * t_{a1}^{i3}"
                           L" * t_{a2}^{i4}",
                           Symmetry::nonsymm)
                           ->as<Product>();

    // canon set on
    const auto result2_naive = single_term_opt(prod2, nocc, nvirt, {}, true);

    // there will be two degenerate evaluation sequences for prod2 when no
    REQUIRE(result2_naive.optimal_seqs.size() == 2);

    // now let's reuse intermediates from the optimal sequence of evaluation of
    // prod1, to compute the optimal sequence for prod2
    container::set<size_t> imed_hashes_prod1{};
    yield_interm_hash(result1.optimal_seqs.at(0), imed_hashes_prod1);

    const auto result2_discounted =
        single_term_opt(prod2, nocc, nvirt, imed_hashes_prod1, true);
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
        std::move(*(single_term_opt(prod3->as<Product>(), nocc, nvirt, {}, true)
                        .optimal_seqs.begin()));

    // finding the intermediate from the evaluation tree of prod4
    auto prod4_sto =
        std::move(*(single_term_opt(prod4->as<Product>(), nocc, nvirt, {}, true)
                        .optimal_seqs.begin()));

    auto imeds_prod4 = container::set<size_t>{};
    yield_interm_hash(prod4_sto, imeds_prod4);

    auto prod3_sto_with_imeds = std::move(
        *single_term_opt(prod3->as<Product>(), nocc, nvirt, imeds_prod4, true)
             .optimal_seqs.begin());

    REQUIRE_FALSE(*prod3_sto == *prod3_sto_with_imeds);
  }

  SECTION("Most expensive term") {
    using optimize::most_expensive;
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
    auto opt_terms = optimize::multi_term_opt_hartono(terms, nocc, nvirt);

    auto prod2_opt = std::move(*(
        optimize::single_term_opt(prod2->as<Product>(), nocc, nvirt, {}, true)
            .optimal_seqs.begin()));
    auto imeds = container::set<size_t>{};
    yield_interm_hash(prod2_opt, imeds);

    auto prod1_opt =
        std::move(*(optimize::single_term_opt(prod1->as<Product>(), nocc,
                                               nvirt, imeds, true)
                        .optimal_seqs.begin()));

    REQUIRE(**(opt_terms.at(prod2).optimal_seqs.begin()) == *prod2_opt);
    REQUIRE(**(opt_terms.at(prod1).optimal_seqs.begin()) == *prod1_opt);
  }
}
