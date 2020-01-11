#include <iostream>

#include "../../src/SeQuant/core/expr.hpp"
#include "../../src/SeQuant/core/expr_algorithm.hpp"
#include "../../src/SeQuant/core/tensor.hpp"
#include "catch.hpp"

TEST_CASE("Canonicalizer", "[algorithms]") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("nonsymmetric op") {
    auto input = ex<Constant>(0.5) *
                     ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                WstrList{L"p_3", L"p_4"}, Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm) +
                 ex<Constant>(0.5) *
                     ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                                WstrList{L"p_4", L"p_3"}, Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                                Symmetry::nonsymm) *
                     ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                                Symmetry::nonsymm);
    std::wcout << to_latex(input) << std::endl;
    canonicalize(input);
    std::wcout << "nonsymmetric: " << to_latex(input) << std::endl;
  }

  SECTION("asymmetric op") {
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
    // std::wcout << "" << to_latex(input) << std::endl;
    canonicalize(input);
    std::wcout << "asymmetric: " << to_latex(input) << std::endl;
  }

  SECTION("antisymmetric op") {
    auto input = ex<Constant>(0.5) *
        ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                   WstrList{L"p_3", L"p_4"}, Symmetry::antisymm) *
        ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                   Symmetry::nonsymm) *
        ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                   Symmetry::nonsymm) +
        ex<Constant>(0.5) *
            ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                       WstrList{L"p_4", L"p_3"}, Symmetry::antisymm) *
            ex<Tensor>(L"t", IndexList{{L"p_3"}}, IndexList{{L"p_1"}},
                       Symmetry::nonsymm) *
            ex<Tensor>(L"t", IndexList{{L"p_4"}}, IndexList{{L"p_2"}},
                       Symmetry::nonsymm);
    // std::wcout << "" << to_latex(input) << std::endl;
    canonicalize(input);
    std::wcout << "nonsymmetric: " << to_latex(input) << std::endl;
  }

}