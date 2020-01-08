//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.hpp"
#include "catch.hpp"

TEST_CASE("Spin Trace") {
  using namespace sequant;
#if 0
  SECTION("constructors") {
    // Check {alpha, beta, null} on index, term
    Index i1(L"i_1", IndexSpace::instance(IndexSpace::active_occupied,
                                          IndexSpace::alpha));
    REQUIRE(i1.space() == IndexSpace::instance(IndexSpace::active_occupied,
                                               IndexSpace::alpha));
    REQUIRE(i1.space().qns() == IndexSpace::alpha);

    Index a1(L"a_1", IndexSpace::instance(IndexSpace::active_occupied,
                                          IndexSpace::beta));
    REQUIRE(a1.space() == IndexSpace::instance(IndexSpace::active_occupied,
                                               IndexSpace::beta));
    REQUIRE(a1.space().qns() == IndexSpace::beta);

    Index i2(L"i_2", IndexSpace::instance(IndexSpace::active_occupied));
    REQUIRE(i2.space() == IndexSpace::instance(IndexSpace::active_occupied));
    REQUIRE(i2.space().qns() == IndexSpace::nullqns);

    Index a2(L"a_2", IndexSpace::instance(IndexSpace::active_occupied,
                                          IndexSpace::alpha));

    std::wcout << "i1: " << i1.label() << " " << "a1: " << a1.label() << "\n";
    // remove spin
    auto t1 = Tensor(L"t", {a1}, {i1});
    auto t2 = Tensor(L"t", {a2}, {i2});
    std::wcout << t1.to_latex() << " " << t2.to_latex() << "\n";

    std::cout << " t1: " << tensor_symm(t1) << " ";
    std::cout << " t2: " << tensor_symm(t2) << "\n";
    ExprPtr T1(new Tensor(t1));
    ExprPtr T2(new Tensor(t2));
    auto exprPtr = T1 * T2;
    std::wcout << "expr: " << exprPtr->to_latex() << "\n";
    remove_spin(exprPtr);
    std::wcout << "expr spin removed: " << exprPtr->to_latex() << "\n";
  }  //  SECTION("constructors")
#endif

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
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 2);
    REQUIRE(to_latex(result) ==
            L"{ \\left({{g^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
            L"{{g^{{a_1}{a_2}}_{{i_2}{i_1}}}}\\right) }");
  }

  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                      ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = spintrace(expr);
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 2);
    REQUIRE(to_latex(result) ==
            L"{ \\left({{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
            L"{{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}}\\right) }");
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
      rapid_simplify(result);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 6);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}"
              L"_{{a_1}}}{t^{{i_2}}_{{a_2}}}} - "
              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
              L"}}{t^{{i_2}}_{{a_2}}}} + "
              L"{{{\\frac{1}{2}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
              L"}}{t^{{i_2}}_{{a_2}}}} + "
              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}"
              L"}}{t^{{i_2}}_{{a_2}}}} + "
              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}"
              L"}}{t^{{i_2}}_{{a_2}}}} - "
              L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}"
              L"}}{t^{{i_2}}_{{a_2}}}}\\right) }");
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
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 12);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{"
              L"i_2}}_{{a_1}{a_2}}}} - "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
              L"{a_2}{a_1}}}} - "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
              L"{a_1}{a_2}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
              L"{a_2}{a_1}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{"
              L"{a_2}{a_1}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
              L"{a_2}{a_1}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{"
              L"{a_1}{a_2}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
              L"{a_1}{a_2}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
              L"{a_1}{a_2}}}} - "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{"
              L"{a_2}{a_1}}}} - "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
              L"{a_1}{a_2}}}} + "
              L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{"
              L"{a_2}{a_1}}}}\\right) }");
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
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 20);
    REQUIRE(to_latex(result) ==
            L"{ \\left({{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
            L"{{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} + "
            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} - "
            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} + "
            L"{{{\\frac{1}{2}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} + "
            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} + "
            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} - "
            L"{{{\\frac{1}{2}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}}_{{a_1}}}"
            L"{t^{{i_2}}_{{a_2}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_1}{a_2}}}} - "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_2}{a_1}}}} - "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
            L"a_1}{a_2}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
            L"a_2}{a_1}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{{"
            L"a_2}{a_1}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_2}{a_1}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_2}{a_1}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_2}{i_1}}_{{"
            L"a_1}{a_2}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_1}{a_2}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_1}{a_2}}}} - "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{"
            L"a_2}{a_1}}}} - "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
            L"a_1}{a_2}}}} + "
            L"{{{\\frac{1}{4}}}{g^{{a_1}{a_2}}_{{i_2}{i_1}}}{t^{{i_1}{i_2}}_{{"
            L"a_2}{a_1}}}}\\right) }");
  }

}  // TEST_CASE("Spin Trace")
