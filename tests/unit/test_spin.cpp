//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.hpp"
#include "catch.hpp"

TEST_CASE("Spin Trace") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("Constant") {
    auto exprPtr = ex<Constant>(1. / 4);
    auto result = spintrace(exprPtr);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->is_atom());
    REQUIRE(to_latex(result) == L"{{{\\frac{1}{4}}}}");
  }

  SECTION("Tensor") {
    const auto expr = ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                 WstrList{L"p_3", L"p_4"}, Symmetry::antisymm);
    auto result = spintrace(expr);
//    std::wcout << "\nExpr:\n"
//               << expr->to_latex() << "\nTraced:\n"
//               << result->to_latex() << std::endl;
    REQUIRE(result->is<Sum>());
    // REQUIRE(result->size() == 2);
    REQUIRE(to_latex(result) ==
            L"{ \\left({{g^{{p_3}{p_4}}_{{p_1}{p_2}}}} - "
            L"{{g^{{p_3}{p_4}}_{{p_2}{p_1}}}}\\right) }");
  }

  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                      ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = spintrace(expr);
    canonicalize(result);
//    std::wcout << "\nExpr:\n"
//               << expr->to_latex() << "\nTraced:\n"
//               << result->to_latex();
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 1);
    REQUIRE(
        to_latex(result) ==
        L"{ \\left({{{2}}{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}}\\right) }");
  }

  SECTION("Scaled Product") {
    {
      // 1/2 * g * t1 * t1
      const auto expr =
          ex<Constant>(1. / 2) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"});
      auto result = spintrace(expr);
/*
      std::cout << "\\begin{align}\n";
      Sum sum_of_canonicalized;
      for(auto &&p : *result){
        if(p->is<Product>()){
          std::wcout << "& " << p->to_latex() <<  " \\\\ \n";
          ExprPtr pp = std::make_shared<Product>(p->as<Product>());
          canonicalize(pp);
          std::wcout << "& " << pp->to_latex() <<  " \\\\ \n";
          sum_of_canonicalized.append(pp);
        }
      }
      std::cout << "\\end{align}\n";
      ExprPtr sum_of_canonicalized_terms = std::make_shared<Sum>(sum_of_canonicalized);
      rapid_simplify(sum_of_canonicalized_terms);
      canonicalize(sum_of_canonicalized_terms);
      std::wcout << "sum_of_canonicalized_terms: " << sum_of_canonicalized_terms->to_latex() << std::endl;
*/
      canonicalize(result);
//      std::wcout << "\nExpr:\n"
//                 << expr->to_latex() << "\nTraced:\n"
//                 << result->to_latex();
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(to_latex(result) ==
              L"{ \\left( - "
              L"{{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_2}}}{t^{{i_2}}_{{"
              L"a_1}}}} + "
              L"{{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}"
              L"}_{{a_2}}}}\\right) }");
    }

    {
      // 1/4 * g * t2
      const auto expr =
          ex<Constant>(1. / 4) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm);

      auto result = spintrace(expr);
      canonicalize(result);
//      std::wcout << "\nExpr:\n"
//                 << expr->to_latex() << "\nTraced:\n"
//                 << result->to_latex() << std::endl;
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_"
              L"1}{a_2}}}} - "
              L"{{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}}"
              L"\\right) }");
    }
  }

  SECTION("Sum") {
    // f * t1 + 1/2 * g * t1 * t1 + 1/4 * g * t2
    const auto ex1 = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                     ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    const auto ex2 = ex<Constant>(1. / 2) *
                     ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                                WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
                     ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"}) *
                     ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"});
    const auto ex3 = ex<Constant>(1. / 4) *
                     ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                                WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
                     ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"},
                                WstrList{L"i_1", L"i_2"}, Symmetry::antisymm);

    auto expr = ex1 + ex2 + ex3;
    auto result = spintrace(expr);
    canonicalize(result);
//    std::wcout << "\nExpr:\n"
//               << expr->to_latex() << "\nTraced:\n"
//               << result->to_latex() << std::endl;
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 5);
    REQUIRE(
        to_latex(result) ==
        L"{ "
        L"\\left({{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_2}"
        L"}}} - {{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}} + "
        L"{{{2}}{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} - "
        L"{{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_2}}}{t^{{i_2}}_{{a_1}}}}"
        L" + "
        L"{{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}}_{{a_"
        L"2}}}}\\right) }");
  }

}  // TEST_CASE("Spin Trace")
