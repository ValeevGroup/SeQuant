#include "catch.hpp"

#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <iostream>

auto yield_interm_hash = [](sequant::EvalNode const& node) {
  auto cont = sequant::container::set<sequant::EvalExpr::hash_t>{};
  node.visit_internal([&cont](const auto& n) { cont.emplace(n->hash()); });
  return cont;
};

TEST_CASE("TEST_OPTIMIZE", "[optimize]") {
  using namespace sequant;
  using optimize::single_term_opt;

  sequant::TensorCanonicalizer::register_instance(
      std::make_shared<sequant::DefaultTensorCanonicalizer>());

  SECTION("Single term optimization") {
    const auto prod1 = parse_expr_asymm(
                           L"g_{i3,i4}^{a3,a4}"     // T1
                           " * t_{a1,a2}^{i3,i4}"   // T2
                           " * t_{a3,a4}^{i1,i2}")  // T3
                           ->as<Product>();
    //
    // Cost of evaluation prod1:
    //
    // ((T1 * T2) * T3)  : 2 * O^2 * V^4  best if nocc > nvirt
    //
    // this is the one we want to find
    // ((T1 * T3) * T2)  : 2 * O^4 * V^2  best if nvirt > nocc
    //
    // ((T2 * T3) * T1)  : 2 * O^4 * V^4  worst sequence of evaluation
    //

    // canonicalization set false as it is unnecessary here
    const auto result1 = single_term_opt(prod1, false);

    REQUIRE(result1.cost != sequant::AsyCost::max());

    // there are no degenerate evaluation sequences
    REQUIRE(result1.optimal_seqs.size() == 1);

    auto prod1_opt = ex<Product>(
        sequant::ExprPtrList{prod1.at(0), prod1.at(2), prod1.at(1)});

    REQUIRE(to_eval_node(prod1_opt) == result1.optimal_seqs.at(0));

    //
    const auto prod2 = parse_expr_asymm(
                           L"   g_{i3,i4}^{a3,a4}"
                           L" * t_{a3,a4}^{i1,i2}"
                           L" * t_{a1}^{i3}"
                           L" * t_{a2}^{i4}")
                           ->as<Product>();

    // canon set on
    const auto result2_naive = single_term_opt(prod2, true);

    // there will be two degenerate evaluation sequences for prod2
    REQUIRE(result2_naive.optimal_seqs.size() == 2);

    // now let's reuse intermediates from the optimal sequence of evaluation of
    // prod1, to compute the optimal sequence for prod2
    auto imed_hashes_prod1 = yield_interm_hash(result1.optimal_seqs.at(0));

    const auto result2_discounted = single_term_opt(
        prod2,  //
        true,   // canonicalization on
        [&imed_hashes_prod1](
            auto const& n) {  // discount existing intermediate costs
          if (imed_hashes_prod1.contains(n->hash())) return false;
          imed_hashes_prod1.emplace(n->hash());
          return true;
        });
    REQUIRE(result2_discounted.cost < result2_naive.cost);

    // yet another example
    auto prod3 = parse_expr_asymm(
        L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3} * t_{a3}^{i4}");
    auto prod4 = parse_expr_asymm(L"t_{a1,a2}^{i1,i2} * g_{i2,i3}^{a2,a3}");
    // we show that two the evaluation trees for prod3
    //  - one: single term optimized on prod3 alone
    //  - two: single term optimized on prod3 with the intermediate from prod4
    //  are not the same.
    auto prod3_sto = std::move(
        *(single_term_opt(prod3->as<Product>(), true).optimal_seqs.begin()));

    // finding the intermediate from the evaluation tree of prod4
    auto prod4_sto = std::move(
        *(single_term_opt(prod4->as<Product>(), true).optimal_seqs.begin()));

    auto imeds_prod4 = yield_interm_hash(prod4_sto);

    auto prod3_sto_with_imeds = std::move(
        *single_term_opt(prod3->as<Product>(), true,
                         [&imeds_prod4](auto const& n) {
                           return !((imeds_prod4.contains(n->hash())));
                         })
             .optimal_seqs.begin());

    REQUIRE_FALSE(*prod3_sto == *prod3_sto_with_imeds);
  }
}
