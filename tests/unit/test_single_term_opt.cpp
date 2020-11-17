#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/factorize/single_term_opt.hpp>

#include "catch.hpp"

TEST_CASE("TEST_SINGLE_TERM_OPT", "[single_term_opt]") {
  using namespace sequant;
  using namespace sequant::factorize;

  auto prod1 = ex<Product>(  //
      Product{
          1. / 16.,
          {ex<Tensor>(Tensor{L"g",
                             {L"i_3", L"i_4"},
                             {L"a_3", L"a_4"},
                             Symmetry::antisymm}),  //
           ex<Tensor>(Tensor{L"t",
                             {L"a_1", L"a_2"},
                             {L"i_3", L"i_4"},
                             Symmetry::antisymm}),  //
           ex<Tensor>(Tensor{
               L"t", {L"a_3", L"a_4"}, {L"i_1", L"i_2"}, Symmetry::antisymm})}}
      //
  );

  SECTION("Optimal eval_sequence") {
    // equiv to ((0 1) 2)
    auto tree1 = eval_sequence{0, {1, 2}};

    // equiv to ((0 2) 1)
    auto tree2 = eval_sequence{0, {2, 1}};

    // equiv to ((1 2) 0)
    auto tree3 = eval_sequence{1, {2, 0}};

    // nocc (= 10) < nvirt (= 20)
    REQUIRE(OptimalRootedTree{prod1, 10, 20}.tree() == tree2);

    // nocc (= 20) > nvirt (= 10)
    REQUIRE(OptimalRootedTree{prod1, 20, 10}.tree() == tree1);

    //
    // nocc (= 10) == nvirt (= 10)
    //
    // optimal tree could be tree1 or tree2
    // but it can never be tree3 as the eval
    // sequence implied by tree3 always leads
    // to more flops than either tree1 or tree2
    // no matter nocc (>/</==) nvirt
    //
    REQUIRE_FALSE(OptimalRootedTree{prod1, 10, 10}.tree() == tree3);
  }

  SECTION("repack_product") {
    // equiv to ((0 1) 2)
    auto prod1Tree1 = eval_sequence{0, {1, 2}};

    // equiv to ((0 2) 1)
    auto prod1Tree2 = eval_sequence{0, {2, 1}};

    // equiv to ((1 2) 0)
    auto prod1Tree3 = eval_sequence{1, {2, 0}};

    const auto& prod1Repack1 = prod1;

    auto prod1Repack2 = ex<Product>(
        Product{1. / 16., {prod1->at(0), prod1->at(2), prod1->at(1)}});

    auto prod1Repack3 = ex<Product>(
        Product{1. / 16., {prod1->at(1), prod1->at(2), prod1->at(0)}});

    REQUIRE(*repack_prod(prod1, prod1Tree1) == *prod1Repack1);
    REQUIRE(*repack_prod(prod1, prod1Tree2) == *prod1Repack2);
    REQUIRE(*repack_prod(prod1, prod1Tree3) == *prod1Repack3);

    auto prod2 = ex<Product>(Product{
        -1. / 2,
        {
            ex<Tensor>(Tensor{L"g",
                              {L"i_3", L"i_4"},
                              {L"a_3", L"a_4"},
                              Symmetry::antisymm}),  //
            ex<Tensor>(
                Tensor{L"t", {L"a_1"}, {L"i_3"}, Symmetry::antisymm}),  //
            ex<Tensor>(
                Tensor{L"t", {L"a_3"}, {L"i_4"}, Symmetry::antisymm}),  //
            ex<Tensor>(Tensor{L"t",
                              {L"a_2", L"a_4"},
                              {L"i_2", L"i_3"},
                              Symmetry::antisymm})  //
        }                                           //
    });

    auto prod2Tree1 = eval_sequence{0, {1, eval_sequence{2, {3}}}};
    auto prod2Repack1 =
        ex<Product>(Product{-1. / 2, {prod2->at(0), prod2->at(1)}});

    auto& prod2Repack1_deref = prod2Repack1->as<Product>();
    // Note: appending without expansion
    prod2Repack1_deref.append(ex<Product>(Product{prod2->at(2), prod2->at(3)}));

    REQUIRE(*repack_prod(prod2, prod2Tree1) == *prod2Repack1);
  }

  SECTION("Single term optimization: exhaustive scan of eval sequence") {
    // optimization of prod1
    // equiv to ((0 1) 2)
    auto tree1 = eval_sequence{0, {1, 2}};

    // equiv to ((0 2) 1)
    auto tree2 = eval_sequence{0, {2, 1}};

    // nocc (= 2) < nvirt (= 3)
    REQUIRE(*sto_exhaustive_scan(prod1, 2, 3) == *repack_prod(prod1, tree2));

    // nocc (= 3) > nvirt (= 2)
    REQUIRE(*sto_exhaustive_scan(prod1, 3, 2) == *repack_prod(prod1, tree1));
  }
}
