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
    const auto expr = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                                 WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
    auto result = spintrace(expr);
    std::wcout << "\nExpr:\n" << expr->to_latex() << "\nTraced:\n" << result->to_latex() << std::endl;
    REQUIRE(result->is<Sum>());
    // REQUIRE(result->size() == 2);
    REQUIRE(to_latex(result) ==
            L"{ \\left({{g^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
            L"{{g^{{a_1}{a_2}}_{{i_2}{i_1}}}}\\right) }");
  }

  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                      ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = spintrace(expr);
    canonicalize(result);
    std::wcout << "\nExpr:\n" << expr->to_latex() << "\nTraced:\n" << result->to_latex() << std::endl;
    REQUIRE(result->is<Sum>());
    // REQUIRE(result->size() == 2);
//    REQUIRE(to_latex(result) ==
//            L"{ \\left({{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
//            L"{{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}}\\right) }");
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
      canonicalize(result);
      std::wcout << "\nExpr:\n" << expr->to_latex() << "\nTraced:\n" << result->to_latex();
      rapid_simplify(result);
      REQUIRE(result->is<Sum>());
//      REQUIRE(result->size() == 6);

//      REQUIRE(to_latex(result) ==
//              L"{ "
//              L"\\left({{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}"
//              L"_{{a_1}}}{t^{{i_2}}_{{a_2}}}} - "
//              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
//              L"}}{t^{{i_2}}_{{a_2}}}} + "
//              L"{{{\\frac{1}{2}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
//              L"}}{t^{{i_2}}_{{a_2}}}} + "
//              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}"
//              L"}}{t^{{i_2}}_{{a_2}}}} + "
//              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}"
//              L"}}{t^{{i_2}}_{{a_2}}}} - "
//              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
//              L"}}{t^{{i_2}}_{{a_2}}}}\\right) }");
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
      std::wcout << "\nExpr:\n" << expr->to_latex() << "\nTraced:\n" << result->to_latex() << std::endl;
      REQUIRE(result->is<Sum>());
//      REQUIRE(result->size() == 12);
//      REQUIRE(to_latex(result) ==
//              L"{ "
//              L"\\left({{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{"
//              L"i_2}}_{{a_1}{a_2}}}} - "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
//              L"{a_2}{a_1}}}} - "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
//              L"{a_1}{a_2}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
//              L"{a_2}{a_1}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{"
//              L"{a_2}{a_1}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
//              L"{a_2}{a_1}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{"
//              L"{a_1}{a_2}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
//              L"{a_1}{a_2}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
//              L"{a_1}{a_2}}}} - "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
//              L"{a_2}{a_1}}}} - "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
//              L"{a_1}{a_2}}}} + "
//              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
//              L"{a_2}{a_1}}}}\\right) }");
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
    std::wcout << "\nExpr:\n" << expr->to_latex() << "\nTraced:\n" << result->to_latex() << std::endl;
    REQUIRE(result->is<Sum>());
//    REQUIRE(result->size() == 20);
//    REQUIRE(to_latex(result) ==
//            L"{ \\left({{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
//            L"{{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
//            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} - "
//            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} + "
//            L"{{{\\frac{1}{2}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} + "
//            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} + "
//            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} - "
//            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
//            L"{t^{{i_2}}_{{a_2}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_1}{a_2}}}} - "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_2}{a_1}}}} - "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
//            L"a_1}{a_2}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
//            L"a_2}{a_1}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{{"
//            L"a_2}{a_1}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_2}{a_1}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{{"
//            L"a_1}{a_2}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_1}{a_2}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_1}{a_2}}}} - "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
//            L"a_2}{a_1}}}} - "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
//            L"a_1}{a_2}}}} + "
//            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
//            L"a_2}{a_1}}}}\\right) }");
  }

}  // TEST_CASE("Spin Trace")
