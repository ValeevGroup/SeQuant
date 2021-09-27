//
// Created by Bimal Gaudel on 9/22/21.
//
#include "catch.hpp"
#include <SeQuant/core/clone.hpp>
#include <SeQuant/core/tensor.hpp>

TEST_CASE("TEST_CLONE", "[clone]") {
  using namespace sequant;
  SECTION("Tensor") {
    REQUIRE(clone(ex<Tensor>(L"t",
                     IndexList{L"i_1"},
                     IndexList{L"a_1"}))->is<Tensor>());
  }
  SECTION("Constant") {
    REQUIRE(clone(ex<Constant>(1))->is<Constant>());
  }

  SECTION("Product") {
    auto t1 = ex<Tensor>(L"g",
          IndexList{L"i_2", L"i_3"},
          IndexList{L"a_2", L"a_3"});
    auto t2 = ex<Tensor>(L"t",
                           IndexList{L"a_1", L"a_3"},
                           IndexList{L"i_2", L"i_3"});
    auto t3 = ex<Tensor>(L"t",
                           IndexList{L"a_2"},
                           IndexList{L"i_1"});
    auto prod1 = ex<Product>(ExprPtrList{t1,t2,t3});
    REQUIRE(prod1 == clone(prod1));

    auto prod2 = ex<Product>(1./2,ExprPtrList{});
    prod2->as<Product>().append(t1);
    prod2->as<Product>().append(ex<Product>(1./3, ExprPtrList{t2,t3}));

    REQUIRE(prod2->at(0)->is<Tensor>());
    REQUIRE(prod2->at(1)->is<Product>());

    REQUIRE(*clone(prod2) == *prod2);
  }

  SECTION("Sum") {
    auto g1 = ex<Tensor>(L"g",
                         IndexList{L"i_2", L"a_1"},
                         IndexList{L"i_1", L"a_2"});
    auto g2 = ex<Tensor>(L"g",
                         IndexList{L"i_2", L"i_3"},
                         IndexList{L"a_2", L"a_3"});
    auto t1 = ex<Tensor>(L"t",
                         IndexList{L"a_2"},
                         IndexList{L"i_2"});
    auto t2 = ex<Tensor>(L"t",
                         IndexList{L"a_2"},
                         IndexList{L"i_1"});
    auto t3 = ex<Tensor>(L"t",
                         IndexList{L"a_1", L"a_3"},
                         IndexList{L"i_2", L"i_3"});

    auto prod1 = ex<Product>(-1, ExprPtrList{g1, t1});
    auto prod2 = ex<Product>(-1./2, ExprPtrList{g2,t2,t3});
    auto sum = ex<Sum>(ExprPtrList{prod1, prod2});
    REQUIRE(*sum == *clone(sum));
  }
}