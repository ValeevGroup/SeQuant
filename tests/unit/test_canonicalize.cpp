#include <iostream>

#include "../../src/SeQuant/core/expr.hpp"
#include "../../src/SeQuant/core/expr_algorithm.hpp"
#include "../../src/SeQuant/core/tensor.hpp"
#include "catch.hpp"

TEST_CASE("Canonicalizer", "[algorithms]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("nonsymmetric tensor") {
    {
    auto op = ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                         WstrList{L"p_3", L"p_4"}, Symmetry::nonsymm);
    canonicalize(op);
    REQUIRE(to_latex(op) == L"{g^{{p_3}{p_4}}_{{p_1}{p_2}}}");
    }
    {
      auto op = ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                           WstrList{L"p_3", L"p_4"}, Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE(to_latex(op) == L"{g^{{p_4}{p_3}}_{{p_1}{p_2}}}");
    }
    {
      auto op = ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                           WstrList{L"p_4", L"p_3"}, Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE(to_latex(op) == L"{g^{{p_4}{p_3}}_{{p_1}{p_2}}}");
    }
    {
      auto op = ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                           WstrList{L"p_4", L"p_3"}, Symmetry::nonsymm);
      canonicalize(op);
      REQUIRE(to_latex(op) == L"{g^{{p_3}{p_4}}_{{p_1}{p_2}}}");
    }
  }

  SECTION("sum of products") {
    auto input = ex<Constant>(0.5) *
                     ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                WstrList{L"p_3", L"p_4"}, Symmetry::symm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm) +
                 ex<Constant>(0.5) *
                     ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                                WstrList{L"p_4", L"p_3"}, Symmetry::symm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm);
    canonicalize(input);
    REQUIRE(to_latex(input) == L"{ \\left({{g^{{p_1}{p_4}}_{{p_2}{p_3}}}{t^{{p_2}}_{{p_1}}}{t^{{p_3}}_{{p_4}}}}\\right) }");
  }
}