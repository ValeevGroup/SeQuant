//
// Created by Nakul Teke on 12/20/19.
//

#include <iostream>

#include "catch.hpp"
#include "SeQuant/domain/mbpt/spin.hpp"

TEST_CASE("Spin"){
  using namespace sequant;

  SECTION("constructors"){
    // Check {alpha, beta, null} on index, term
    Index i1(L"i_1", IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha));
    REQUIRE(i1.space()== IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha));
    REQUIRE(i1.space().qns() == IndexSpace::alpha);

    Index a1(L"a", IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta));
    REQUIRE(a1.space()== IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta));
    REQUIRE(a1.space().qns() == IndexSpace::beta);

    Index i2(L"i_2", IndexSpace::instance(IndexSpace::active_occupied));
    REQUIRE(i2.space()== IndexSpace::instance(IndexSpace::active_occupied));
    REQUIRE(i2.space().qns() == IndexSpace::nullqns);

    // remove spin


  } //  SECTION("constructors")

  SECTION("index transformations"){

    // add spin


    // remove spin

  } //  SECTION("constructors")

  SECTION("algorithms"){

    // define sequence of operators
    Tensor g = Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}},
                       {Index{L"a_1"}, Index{L"a_2"}}, Symmetry::antisymm);
    Tensor t2 = Tensor(L"t",{Index{L"a_1"}, Index{L"a_2"}}, {Index{L"i_1"}, Index{L"i_2"}}, Symmetry::antisymm);

    auto G = std::make_shared<Tensor>(g);
    auto T2 = std::make_shared<Tensor>(t2);

    Product expr;
    expr.scale(0.25);
    expr.append(1,G);
    expr.append(1,T2);
    // std::wcout << "expr: " << expr.to_latex() << std::endl;

    auto productPtr = std::make_shared<Product>(expr);
    // auto result = spintrace(productPtr);
    // std::wcout << "result:\n" <<  to_latex_align(result) << std::endl;

    Sum sum{};
    sum.append(productPtr);
    auto sumPtr = std::make_shared<Sum>(sum);
    auto result = spintrace(sumPtr);
    std::wcout << "sum result:\n" <<  to_latex_align(result) << std::endl;




    // anti-symmetrize

    // zero out terms

    // add together

  } // SECTION("algorithms")

  SECTION("expr"){
    // spin trace on constant, tensor, product and sum

  }

  SECTION("latex"){

  } // SECTION("latex")

} // TEST_CASE("Spin", "[terms]")

