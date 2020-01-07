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
    Tensor f = Tensor(L"F",{Index{L"i_1"}},{Index{L"a_1"}});
    Tensor t1 = Tensor(L"t",{Index{L"a_1"}},{Index{L"i_1"}});
    Tensor t1_2 = Tensor(L"t",{Index{L"a_2"}},{Index{L"i_2"}});

    auto F = std::make_shared<Tensor>(f);
    auto T1 = std::make_shared<Tensor>(t1);
    auto T1_2 = std::make_shared<Tensor>(t1_2);

    Tensor g = Tensor(L"g", {Index{L"i_1"}, Index{L"i_2"}},
                       {Index{L"a_1"}, Index{L"a_2"}}, Symmetry::antisymm);
    Tensor t2 = Tensor(L"t",{Index{L"a_1"}, Index{L"a_2"}}, {Index{L"i_1"}, Index{L"i_2"}}, Symmetry::antisymm);

    auto G = std::make_shared<Tensor>(g);
    auto T2 = std::make_shared<Tensor>(t2);
/*
    {
      auto exprPtr = std::make_shared<Constant>(0.25);
      auto result = spintrace(exprPtr);
      // REQUIRE(result->is<Constant>());
      std::wcout << "\nExpr:\n" << exprPtr->to_latex() << "\nSpin traced:\n" << result->to_latex() << "\n";
    }

    {
      Tensor expr = g;
      auto exprPtr = std::make_shared<Tensor>(expr);
      auto result = spintrace(exprPtr);
       // REQUIRE(result->is<Sum>());
      std::wcout << "\n" << exprPtr->to_latex() << "\nSpin traced:\n" << result->to_latex() << "\n";
    }
*/
    {
      Product expr;
      expr.append(1,F);
      expr.append(1,T1);
      auto exprPtr = std::make_shared<Product>(expr);
      {
        std::wcout << "\n" << exprPtr->to_latex() << std::endl;
         auto result = spintrace(exprPtr);
        TensorCanonicalizer::register_instance(
            std::make_shared<DefaultTensorCanonicalizer>());
        canonicalize(result);
        std::wcout << "\nSpin traced:\n" << to_latex_align(result) << std::endl;
          REQUIRE(result->is<Sum>());
      }


      Product expr2;
      expr2.scale(0.5);
      expr2.append(1,G);
      expr2.append(1,T1);
      expr2.append(1,T1_2);
      auto expr2Ptr = std::make_shared<Product>(expr2);
      {
        std::wcout << "\nExpr:\n" << expr2Ptr->to_latex() << std::endl;
        auto result = spintrace(expr2Ptr);
        REQUIRE(result->is<Sum>());
        std::wcout << "\nSpin traced:\n" << to_latex_align(result) << "\n";
      }

      Product expr3;
      expr3.scale(0.25);
      expr3.append(1,G);
      expr3.append(1,T2);
      auto expr3Ptr = std::make_shared<Product>(expr3);
      {
        std::wcout << "\nExpr:\n" << expr3Ptr->to_latex() << std::endl;
        auto result = spintrace(expr3Ptr);
        REQUIRE(result->is<Sum>());
        std::wcout << "\nSpin traced:\n" << to_latex_align(result) << "\n";
      }

      Sum sum{};
      sum.append(exprPtr);
      sum.append(expr2Ptr);
      sum.append(expr3Ptr);
      ExprPtr sumPtr (new Sum(sum));
      {
        std::wcout << "\nExpr:\n" << sumPtr->to_latex() << std::endl;
        auto result = spintrace(sumPtr);
        REQUIRE(result->is<Sum>());
        std::wcout << "\nSpin traced:\n" << to_latex_align(result) << "\n";
      }
    }

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

