#include <iostream>

#include "SeQuant/core/expr.hpp"
#include "SeQuant/core/expr_algorithm.hpp"
#include "SeQuant/core/tensor.hpp"
#include "catch.hpp"

TEST_CASE("Canonicalizer", "[algorithms]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("Tensors") {
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
                                WstrList{L"p_3", L"p_4"}, Symmetry::nonsymm, BraKetSymmetry::symm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm) +
                 ex<Constant>(0.5) *
                     ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                                WstrList{L"p_4", L"p_3"}, Symmetry::nonsymm, BraKetSymmetry::symm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm);
    std::wcout <<  "\nInput: " << to_latex(input) << std::endl;
    canonicalize(input);
    std::wcout << "nonsymmetric: " << to_latex(input) << std::endl;
  }

  SECTION("asymmetric op") {
    auto input = ex<Constant>(0.5) *
        ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                   WstrList{L"p_3", L"p_4"}, Symmetry::symm, BraKetSymmetry::symm) *
        ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                   Symmetry::nonsymm) *
        ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                   Symmetry::nonsymm) +
        ex<Constant>(0.5) *
            ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                       WstrList{L"p_4", L"p_3"}, Symmetry::symm, BraKetSymmetry::symm) *
            ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                       Symmetry::nonsymm) *
            ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                       Symmetry::nonsymm);
     std::wcout << "\nInput: " << to_latex(input) << std::endl;
    canonicalize(input);
    std::wcout << "asymmetric: " << to_latex(input) << std::endl;
  }

  SECTION("antisymmetric op") {

    // CASE 1
    {
      auto input = ex<Constant>(0.5) *
          ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                     WstrList{L"p_3", L"p_4"}, Symmetry::antisymm, BraKetSymmetry::symm) *
          ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                     Symmetry::nonsymm) +
          ex<Constant>(0.5) *
              ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                         WstrList{L"p_4", L"p_3"}, Symmetry::antisymm, BraKetSymmetry::symm) *
              ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                         Symmetry::nonsymm);
      std::wcout << "\nInput: " << to_latex(input) << std::endl;
      canonicalize(input);
      // REQUIRE(to_latex(input) == L"{ \\left({{g^{{p_2}{p_4}}_{{p_1}{p_3}}}{t^{{p_1}}_{{p_2}}}{t^{{p_3}}_{{p_4}}}}\\right) }");
      std::wcout << "antisymmetric: " << to_latex(input) << std::endl;
    }

    // CASE 2
    {
      auto input = ex<Constant>(0.5) *
          ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                     WstrList{L"p_3", L"p_4"}, Symmetry::antisymm, BraKetSymmetry::symm) *
          ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}}) *
          ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}}) +
          ex<Constant>(0.5) *
              ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                         WstrList{L"p_3", L"p_4"}, Symmetry::antisymm, BraKetSymmetry::symm) *
              ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}}) *
              ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}});
      std::wcout << "\n**Input: " << to_latex(input) << std::endl;
      canonicalize(input);
      // REQUIRE(to_latex(input) == L"{ \\left({{g^{{p_2}{p_4}}_{{p_1}{p_3}}}{t^{{p_1}}_{{p_2}}}{t^{{p_3}}_{{p_4}}}}\\right) }");
      std::wcout << "antisymmetric: " << to_latex(input) << std::endl;

    }
  }

}