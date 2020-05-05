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
    {
      // CASE 1: Non-symmetric tensors
      auto input = ex<Constant>(0.5) *
                       ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                  WstrList{L"p_3", L"p_4"}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_3"}},
                                  IndexList{{L"p_1"}}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_4"}},
                                  IndexList{{L"p_2"}}, Symmetry::nonsymm) +
                   ex<Constant>(0.5) *
                       ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                                  WstrList{L"p_4", L"p_3"}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_3"}},
                                  IndexList{{L"p_1"}}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_4"}},
                                  IndexList{{L"p_2"}}, Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE(to_latex(input) ==
              L"{ \\bigl({{g^{{p_2}{p_3}}_{{p_1}{p_4}}}{t^{{p_1}}_{{p_2}}}{t^{{p_4}}_{{p_3}}}}\\bigr) }");
    }

    // CASE 2: Symmetric tensors
    {
      auto input = ex<Constant>(0.5) *
                       ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                  WstrList{L"p_3", L"p_4"}, Symmetry::symm) *
                       ex<Tensor>(L"t", IndexList{{L"p_3"}},
                                  IndexList{{L"p_1"}}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_4"}},
                                  IndexList{{L"p_2"}}, Symmetry::nonsymm) +
                   ex<Constant>(0.5) *
                       ex<Tensor>(L"g", WstrList{L"p_2", L"p_1"},
                                  WstrList{L"p_4", L"p_3"}, Symmetry::symm) *
                       ex<Tensor>(L"t", IndexList{{L"p_3"}},
                                  IndexList{{L"p_1"}}, Symmetry::nonsymm) *
                       ex<Tensor>(L"t", IndexList{{L"p_4"}},
                                  IndexList{{L"p_2"}}, Symmetry::nonsymm);
      canonicalize(input);
      REQUIRE(to_latex(input) ==
              L"{ "
              L"\\bigl({{g^{{p_1}{p_4}}_{{p_2}{p_3}}}{t^{{p_2}}_{{p_1}}}{t^{{p_"
              L"3}}_{{p_4}}}}\\bigr) }");
    }

    // Case 3: Anti-symmetric tensors
    {
      auto input =
          ex<Constant>(0.5) *
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
      canonicalize(input);
      REQUIRE(to_latex(input) ==
              L"{ "
              L"\\bigl({{\\bar{g}^{{p_1}{p_4}}_{{p_2}{p_3}}}{t^{{p_2}}_{{p_1}}}"
              L"{t^{{p_"
              L"3}}_{{p_4}}}}\\bigr) }");
    }

    // Case 4: permuted indices
    {
      auto input =
          ex<Constant>(4. / 3.) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                         WstrList{L"a_3", L"i_1"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", IndexList{{L"a_2"}}, IndexList{{L"i_3"}},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{L"a_1", L"a_3"},
                         IndexList{L"i_4", L"i_2"}, Symmetry::antisymm) -
          ex<Constant>(1. / 3.) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                         WstrList{L"i_1", L"a_3"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", IndexList{{L"a_2"}}, IndexList{{L"i_4"}},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{L"a_1", L"a_3"},
                         IndexList{L"i_3", L"i_2"}, Symmetry::antisymm);
      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE(to_latex(input) ==
              L"{ "
              "\\bigl({{\\bar{g}^{{i_1}{a_3}}_{{i_3}{i_4}}}{t^{{i_3}}_{{a_2}}}{"
              "t^{{i_2}{i_4}}_{{a_1}{a_3}}}}\\bigr) }");
    }

    // Case 4: permuted indices from CCSD R2 biorthogonal configuration
    {
      auto input =
          ex<Constant>(4. / 3.) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                         WstrList{L"a_3", L"i_1"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{{L"a_2"}}, IndexList{{L"i_3"}},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{L"a_1", L"a_3"},
                         IndexList{L"i_4", L"i_2"}, Symmetry::nonsymm) -
          ex<Constant>(1. / 3.) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                         WstrList{L"i_1", L"a_3"}, Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{{L"a_2"}}, IndexList{{L"i_4"}},
                         Symmetry::nonsymm) *
              ex<Tensor>(L"t", IndexList{L"a_1", L"a_3"},
                         IndexList{L"i_3", L"i_2"}, Symmetry::nonsymm);

      canonicalize(input);
      REQUIRE(input->size() == 1);
      REQUIRE(to_latex(input) ==
               L"{ \\bigl({{g^{{a_3}{i_1}}_{{i_3}{i_4}}}{t^{{i_3}}_{{a_2}}}{t^{{i_4}{i_2}}_{{a_1}{a_3}}}}\\bigr) }");
    }
  }
}