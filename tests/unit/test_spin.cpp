//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.cpp"
#include "catch.hpp"

TEST_CASE("Spin") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("No proto") {
    Index i1(L"i_1");
    Index i2(L"i_2");
    Index a1(L"a_1");
    Index a2(L"a_2");

    Index i3(L"i_3", IndexSpace::instance(IndexSpace::active_occupied),
             {i1, i2});
    Index a3(L"a_3", IndexSpace::instance(IndexSpace::active_occupied),
             {a1, a2});

    const auto input = ex<Tensor>(L"t", IndexList{i3}, IndexList{a3});
    // TODO: Use exception for passing
    // REQUIRE_THROWS_AS(spintrace(input),std::abort());
  }

  SECTION("Tensor: can_expand, is_tensor_spin_symm, remove_spin") {
    auto p1 = Index(L"p⁺_1",
                    IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p2 =
        Index(L"p⁻_2", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    auto p3 = Index(L"p⁺_3",
                    IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p4 =
        Index(L"p⁻_4", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));

    auto input = ex<Tensor>(L"t", IndexList{p1, p2}, IndexList{p3, p4});
    REQUIRE(can_expand(input->as<Tensor>()) == true);
    REQUIRE(is_tensor_spin_symm(input->as<Tensor>()) == true);

    auto result = remove_spin(input);
    for (auto i : result->as<Tensor>().const_braket())
      REQUIRE(i.space() ==
              IndexSpace::instance(IndexSpace::all, IndexSpace::nullqns));

    input = ex<Tensor>(L"t", IndexList{p1, p3}, IndexList{p2, p4});
    REQUIRE(can_expand(input->as<Tensor>()) == false);
    REQUIRE(is_tensor_spin_symm(input->as<Tensor>()) == false);
  }

  SECTION("Tensor: expand_antisymm") {
    // 1-body
    auto input = ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = expand_antisymm(input->as<Tensor>());
    REQUIRE(input->as<Tensor>() == result->as<Tensor>());
    REQUIRE(!result->is<Sum>());
    REQUIRE(to_latex(result) == L"{t^{{i_1}}_{{a_1}}}");

    // 1-body
    input = ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"}, Symmetry::antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(input->as<Tensor>().symmetry() == Symmetry::antisymm);
    REQUIRE(result->as<Tensor>().symmetry() == Symmetry::nonsymm);
    REQUIRE(!result->is<Sum>());
    REQUIRE(to_latex(result) == L"{t^{{i_1}}_{{a_1}}}");

    // 2-body
    input = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                       Symmetry::antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(result->is<Sum>());
    REQUIRE(to_latex(result) ==
            L"{ \\bigl({{g^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
            L"{{g^{{a_1}{a_2}}_{{i_2}{i_1}}}}\\bigr) }");

    // 3-body
    input = ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                       WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(result->is<Sum>());
    REQUIRE(to_latex(result) ==
            L"{ \\bigl({{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_3}{a_2}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_2}{a_1}{a_3}}}} + "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_2}{a_3}{a_1}}}} + "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_3}{a_1}{a_2}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_3}{a_2}{a_1}}}}\\bigr) }");
  }

  SECTION("Constant") {
    auto exprPtr = ex<Constant>(1. / 4);
    auto result = spintrace(exprPtr);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->is_atom());
    REQUIRE(to_latex(result) == L"{{{\\frac{1}{4}}}}");
  }

  SECTION("Tensor") {
    const auto expr = ex<Constant>(0.25) *
                      ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                 WstrList{L"p_3", L"p_4"}, Symmetry::antisymm);
    auto result = spintrace(expr);
    REQUIRE(result->is<Sum>());
    canonicalize(result);
    REQUIRE(result->size() == 2);
    REQUIRE(to_latex(result) ==
            L"{ \\bigl({{{2}}{g^{{p_3}{p_4}}_{{p_1}{p_2}}}} - {{{2}}{g^{{p_4}{p_3}}_{{p_1}"
            L"{p_2}}}}\\bigr) }");
  }

  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                      ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = spintrace(expr, {{L"i_1", L"a_1"}});
    canonicalize(result);
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 1);
    REQUIRE(
        to_latex(result) ==
        L"{ \\bigl({{{2}}{f^{{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}}\\bigr) }");
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
      auto result = spintrace(expr, {{L"i_1", L"a_1"}});
      canonicalize(result);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(to_latex(result) ==
              L"{ \\bigl( - {{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}}_{{a_1}}}{t^{{i_1}}_{"
          L"{a_2}}}} + {{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}}_{"
          L"{a_2}}}}\\bigr) }");
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
    auto result = ex<Constant>(0.5) * spintrace(expr);
    expand(result);
    rapid_simplify(result);
    canonicalize(result);
    REQUIRE(result->is<Sum>());
    REQUIRE(result->size() == 5);
    REQUIRE(
        to_latex(result) == L"{ \\bigl( - {{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}} + {{"
        L"{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_2}}}} + {{{2}}{f^{"
        L"{a_1}}_{{i_1}}}{t^{{i_1}}_{{a_1}}}} - {{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}"
        L"}_{{a_1}}}{t^{{i_1}}_{{a_2}}}} + {{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}"
        L"}_{{a_1}}}{t^{{i_2}}_{{a_2}}}}\\bigr) }");
  }  // Sum

  SECTION("Expand Antisymmetrizer") {
    // 0-body
    {
      auto input = ex<Constant>(1);
      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 0);
      REQUIRE(result->is_atom());

      input =
          ex<Constant>(1) * ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"},
                                       Symmetry::antisymm);
      result = expand_A_operator(input);
      REQUIRE(result->size() == 0);
      REQUIRE(result->is_atom());
    }

    // 1-body
    {
      auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"},
                              Symmetry::antisymm) *
                   ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"},
                              Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 1);
      REQUIRE(!result->is<Sum>());
    }

    // 2-body
    {
      auto input = ex<Constant>(1. / 4.) *
                   ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                              WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 1);
      REQUIRE(to_latex(result) ==
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}}");

      input = ex<Constant>(1. / 4.) *
              ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"},
                         WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
              ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                         WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
      result = expand_A_operator(input);
      REQUIRE(result->size() == 4);
      REQUIRE(result->is<Sum>());
      REQUIRE(
          to_latex(result) ==
          L"{ \\bigl({{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_2}{a_1}}_{{i_1}{i_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_2}{i_1}}}} + "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_2}{a_1}}_{{i_2}{i_1}}}}\\bigr) }");

      // 1/4 * A * g * t1 * t1
      input = ex<Constant>(1. / 4.) *
              ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"},
                         WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
              ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"},
                         WstrList{L"a_3", L"a_4"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_2"});
      result = expand_A_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_1}{a_2}}}{t^"
              L"{{i_1}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}} - "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_2}{a_1}}}{t^{{i_1}}"
              L"_{{a_3}}}{t^{{i_2}}_{{a_4}}}} - "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_1}{a_2}}}{t^{{i_2}}"
              L"_{{a_3}}}{t^{{i_1}}_{{a_4}}}} + "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_2}{a_1}}}{t^{{i_2}}"
              L"_{{a_3}}}{t^{{i_1}}_{{a_4}}}}\\bigr) }");

      // 1/4 * A * g * t1 * t1 * t1 * t1
      input = ex<Constant>(1. / 4.) *
              ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"},
                         WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                         WstrList{L"a_3", L"a_4"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_2"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"});
      result = expand_A_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(
          to_latex(result) ==
          L"{ "
          L"\\bigl({{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_"
          L"1}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_3}}_{{a_1}}}{t^{{i_4}}_{{a_2}"
          L"}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_2}}_{{"
          L"a_3}}}{t^{{i_1}}_{{a_4}}}{t^{{i_3}}_{{a_1}}}{t^{{i_4}}_{{a_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_1}}_{{"
          L"a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_3}}_{{a_2}}}{t^{{i_4}}_{{a_1}}}} + "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_2}}_{{"
          L"a_3}}}{t^{{i_1}}_{{a_4}}}{t^{{i_3}}_{{a_2}}}{t^{{i_4}}_{{a_1}}}}"
          L"\\bigr) }");
    }

    // 3-body
    {
      auto input =
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(to_latex(result) == L"{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}");

      input = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                         WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                         WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      result = expand_A_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 36);
    }

    {  // 4-body
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3", L"i_4"},
                     WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     WstrList{L"i_1", L"i_2", L"i_3", L"i_4"},
                     Symmetry::antisymm);
      auto asm_input = expand_A_operator(input);
      REQUIRE(asm_input->size() == 576);
      REQUIRE(asm_input->is<Sum>());
    }

    {  // 5-body
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                     WstrList{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                     WstrList{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                     Symmetry::antisymm);
      auto asm_input = expand_A_operator(input);
      REQUIRE(asm_input->size() == 14400);
      REQUIRE(asm_input->is<Sum>());
    }
  }

  SECTION("Expand Symmetrizer") {
    {  // 2-body
      const auto input =
          ex<Tensor>(L"S", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm);
      auto result = expand_S_operator(input);
      REQUIRE(result->size() == 2);
      REQUIRE(result->is<Sum>());
      REQUIRE(to_latex(result) ==
              L"{ \\bigl({{t^{{i_1}{i_2}}_{{a_1}{a_2}}}} + "
              L"{{t^{{i_2}{i_1}}_{{a_2}{a_1}}}}\\bigr) }");
    }

    {  // 3-body
      const auto input =
          ex<Tensor>(L"S", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result = expand_S_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 6);
      REQUIRE(to_latex(result) ==
              L"{ \\bigl({{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}} + "
              L"{{t^{{i_1}{i_3}{i_2}}_{{a_1}{a_3}{a_2}}}} + "
              L"{{t^{{i_2}{i_1}{i_3}}_{{a_2}{a_1}{a_3}}}} + "
              L"{{t^{{i_2}{i_3}{i_1}}_{{a_2}{a_3}{a_1}}}} + "
              L"{{t^{{i_3}{i_1}{i_2}}_{{a_3}{a_1}{a_2}}}} + "
              L"{{t^{{i_3}{i_2}{i_1}}_{{a_3}{a_2}{a_1}}}}\\bigr) }");
    }

    {  // 4-body
      const auto input =
          ex<Tensor>(L"S", WstrList{L"i_1", L"i_2", L"i_3", L"i_4"},
                     WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     Symmetry::nonsymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     WstrList{L"i_1", L"i_2", L"i_3", L"i_4"},
                     Symmetry::antisymm);
      auto result = expand_S_operator(input);
      REQUIRE(result->size() == 24);
      REQUIRE(result->is<Sum>());
      REQUIRE(to_latex(result) ==
              L"{ \\bigl({{t^{{i_1}{i_2}{i_3}{i_4}}_{{a_1}{a_2}{a_3}{a_4}}}} + "
              L"{{t^{{i_1}{i_2}{i_4}{i_3}}_{{a_1}{a_2}{a_4}{a_3}}}} + "
              L"{{t^{{i_1}{i_3}{i_2}{i_4}}_{{a_1}{a_3}{a_2}{a_4}}}} + "
              L"{{t^{{i_1}{i_3}{i_4}{i_2}}_{{a_1}{a_3}{a_4}{a_2}}}} + "
              L"{{t^{{i_1}{i_4}{i_2}{i_3}}_{{a_1}{a_4}{a_2}{a_3}}}} + "
              L"{{t^{{i_1}{i_4}{i_3}{i_2}}_{{a_1}{a_4}{a_3}{a_2}}}} + "
              L"{{t^{{i_2}{i_1}{i_3}{i_4}}_{{a_2}{a_1}{a_3}{a_4}}}} + "
              L"{{t^{{i_2}{i_1}{i_4}{i_3}}_{{a_2}{a_1}{a_4}{a_3}}}} + "
              L"{{t^{{i_2}{i_3}{i_1}{i_4}}_{{a_2}{a_3}{a_1}{a_4}}}} + "
              L"{{t^{{i_2}{i_3}{i_4}{i_1}}_{{a_2}{a_3}{a_4}{a_1}}}} + "
              L"{{t^{{i_2}{i_4}{i_1}{i_3}}_{{a_2}{a_4}{a_1}{a_3}}}} + "
              L"{{t^{{i_2}{i_4}{i_3}{i_1}}_{{a_2}{a_4}{a_3}{a_1}}}} + "
              L"{{t^{{i_3}{i_1}{i_2}{i_4}}_{{a_3}{a_1}{a_2}{a_4}}}} + "
              L"{{t^{{i_3}{i_1}{i_4}{i_2}}_{{a_3}{a_1}{a_4}{a_2}}}} + "
              L"{{t^{{i_3}{i_2}{i_1}{i_4}}_{{a_3}{a_2}{a_1}{a_4}}}} + "
              L"{{t^{{i_3}{i_2}{i_4}{i_1}}_{{a_3}{a_2}{a_4}{a_1}}}} + "
              L"{{t^{{i_3}{i_4}{i_1}{i_2}}_{{a_3}{a_4}{a_1}{a_2}}}} + "
              L"{{t^{{i_3}{i_4}{i_2}{i_1}}_{{a_3}{a_4}{a_2}{a_1}}}} + "
              L"{{t^{{i_4}{i_1}{i_2}{i_3}}_{{a_4}{a_1}{a_2}{a_3}}}} + "
              L"{{t^{{i_4}{i_1}{i_3}{i_2}}_{{a_4}{a_1}{a_3}{a_2}}}} + "
              L"{{t^{{i_4}{i_2}{i_1}{i_3}}_{{a_4}{a_2}{a_1}{a_3}}}} + "
              L"{{t^{{i_4}{i_2}{i_3}{i_1}}_{{a_4}{a_2}{a_3}{a_1}}}} + "
              L"{{t^{{i_4}{i_3}{i_1}{i_2}}_{{a_4}{a_3}{a_1}{a_2}}}} + "
              L"{{t^{{i_4}{i_3}{i_2}{i_1}}_{{a_4}{a_3}{a_2}{a_1}}}}\\bigr) }");
    }

    {
      const auto input = ex<Constant>(4.0) *
          ex<Tensor>(L"S", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::nonsymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"}, WstrList{L"a_4", L"a_5"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_4"}) *
          ex<Tensor>(L"t", WstrList{L"a_5"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_5", L"i_3"});
      // std::wcout << "Input: " << to_latex(input) << "\n";
      auto result = expand_S_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 6);
      result->canonicalize();
      rapid_simplify(result);
      REQUIRE(to_latex(result) == L"{ \\bigl({{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_4}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_1}}_{{a_5}}}{t^{{i_5}{i_3}}_{{a_1}{a_2}}}} + {{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_5}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_1}}_{{a_5}}}{t^{{i_3}{i_4}}_{{a_1}{a_2}}}} + {{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_4}}_{{a_1}}}{t^{{i_2}}_{{a_4}}}{t^{{i_3}}_{{a_5}}}{t^{{i_1}{i_5}}_{{a_2}{a_3}}}} + {{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_4}}_{{a_2}}}{t^{{i_3}}_{{a_4}}}{t^{{i_1}}_{{a_5}}}{t^{{i_5}{i_2}}_{{a_1}{a_3}}}} + {{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_4}}_{{a_1}}}{t^{{i_3}}_{{a_4}}}{t^{{i_2}}_{{a_5}}}{t^{{i_5}{i_1}}_{{a_2}{a_3}}}} + {{{4}}{g^{{a_4}{a_5}}_{{i_4}{i_5}}}{t^{{i_5}}_{{a_2}}}{t^{{i_3}}_{{a_4}}}{t^{{i_1}}_{{a_5}}}{t^{{i_2}{i_4}}_{{a_1}{a_3}}}}\\bigr) }");
    }
  }

  SECTION("Symmetrize expression"){
    {
      // g * t1 + g * t1
      auto input = ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"a_3"}, Symmetry::symm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_2"}) +
          ex<Tensor>(L"g", WstrList{L"a_2", L"a_1"}, WstrList{L"i_2", L"a_3"}, Symmetry::symm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});
      auto result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{S^{{a_1}{a_2}}_{{i_1}{i_2}}}{g^{{i_2}{a_3}}_{{a_1}{a_2}}}{t^{{i_1}}_{{a_3}}}}");
/*
      input = ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"a_3"}, Symmetry::symm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_2"}) +
          ex<Tensor>(L"g", WstrList{L"a_2", L"a_1"}, WstrList{L"i_2", L"a_3"}, Symmetry::symm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});
      result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, false);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{S^{{a_1}{a_2}}_{{i_1}{i_2}}}{g^{{i_2}{a_3}}_{{a_1}{a_2}}}{t^{{i_1}}_{{a_3}}}}");
*/

    }

    {
        // g * t1 * t1 * t1 + g * t1 * t1 * t1
      auto input = ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"i_1", L"a_3"}, Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_2"}) +
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"i_2", L"a_3"}, Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});

      auto result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{S^{{a_2}{a_1}}_{{i_3}{i_4}}}{g^{{i_4}{a_3}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}}_{{a_2}}}{t^{{i_3}}_{{a_3}}}}");
/*
      input = ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"i_1", L"a_3"}, Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_2"}) +
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"i_2", L"a_3"}, Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});
      result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, false);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{S^{{a_2}{a_1}}_{{i_3}{i_4}}}{g^{{i_4}{a_3}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}}_{{a_2}}}{t^{{i_3}}_{{a_3}}}}");
*/
    }

    {
      // g * t1 * t1 * t2 + g * t1 * t1 * t2
      auto input = ex<Constant>(2.0) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"}, WstrList{L"i_1", L"i_2"}) +
          ex<Constant>(2.0) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_2", L"a_4"}, WstrList{L"i_2", L"i_1"});
      // std::wcout << "inp: " << to_latex(input) << std::endl;
      auto result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{{2}}{S^{{a_1}{a_2}}_{{i_1}{i_2}}}{g^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_3}}_{{a_4}}}{t^{{i_4}}_{{a_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_3}}}}");
/*
      input = ex<Constant>(2.0) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"}, WstrList{L"i_1", L"i_2"}) +
          ex<Constant>(2.0) *
              ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},Symmetry::symm) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_4"}) *
              ex<Tensor>(L"t", WstrList{L"a_2", L"a_4"}, WstrList{L"i_2", L"i_1"});
      result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, false);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{{2}}{S^{{a_1}{a_2}}_{{i_1}{i_2}}}{g^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_3}}_{{a_4}}}{t^{{i_4}}_{{a_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_3}}}}");
*/
    }

    {
      // 3-particle operators

    }

    {
      // 4-particle operators
    }
  }

  SECTION("Transform expression") {
    // - A * g * t1
    const auto input =
        ex<Constant>(-1) *
        ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
        ex<Tensor>(L"g", WstrList{L"i_2", L"a_1"}, WstrList{L"i_1", L"a_2"},
                   Symmetry::antisymm) *
        ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"});
    auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
    expand(result);
    rapid_simplify(result);
    canonicalize(result);
    REQUIRE(
        to_latex(result) ==
        L"{ \\bigl({{{2}}{g^{{i_1}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}} - "
        L"{{g^{{a_2}{i_1}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}}\\bigr) }");

    std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                     {Index{L"i_2"}, Index{L"i_1"}}};
    auto transformed_result = transform_expression(result, idxmap);
    REQUIRE(transformed_result->is<Sum>());
    REQUIRE(transformed_result->size() == 2);
    REQUIRE(
        to_latex(transformed_result) ==
        L"{ \\bigl({{{2}}{g^{{i_2}{a_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}} - "
        L"{{g^{{a_2}{i_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}}\\bigr) }");
  }

  SECTION("CCSD R1") {
    // These terms from CCSD R1 equations
    {
      // A * f
      const auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"i_1"});
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      REQUIRE(to_latex(result) == L"{ \\bigl({{f^{{i_1}}_{{a_1}}}}\\bigr) }");
    }

    {
      // Transform indices in an expression
      // - A * g * t1
      const auto input =
          ex<Constant>(-1) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"a_1"}, WstrList{L"i_1", L"a_2"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"});
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      REQUIRE(
          to_latex(result) ==
              L"{ \\bigl({{{2}}{g^{{i_1}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}} - "
              L"{{g^{{a_2}{i_1}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}}\\bigr) }");

      std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                       {Index{L"i_2"}, Index{L"i_1"}}};
      auto transformed_result = transform_expression(result, idxmap);
      REQUIRE(transformed_result->is<Sum>());
      REQUIRE(transformed_result->size() == 2);
      REQUIRE(
          to_latex(transformed_result) ==
              L"{ \\bigl({{{2}}{g^{{i_2}{a_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}} - "
              L"{{g^{{a_2}{i_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}}\\bigr) }");
    }

    {
      // - A * f * t1
      const auto input = ex<Constant>(-1) *
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"f", WstrList{L"i_2"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_2"});
      // std::wcout << "input:3  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      // std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(
          to_latex(result) ==
              L"{ \\bigl( - {{f^{{i_1}}_{{i_2}}}{t^{{i_2}}_{{a_1}}}}\\bigr) }");
    }

    {
      // A * f * t1
      const auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"a_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_1"});
      // std::wcout << "input:4  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      // std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl({{f^{{a_2}}_{{a_1}}}{t^{{i_1}}_{{a_2}}}}\\bigr) }");
    }

    {
      // -1/2 * A * g * t2
      const auto input =
          ex<Constant>(-1. / 2.) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"i_1", L"a_2"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_2", L"i_3"},
                         Symmetry::antisymm);
      // std::wcout << "input:5  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      // std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl( - "
          L"{{{2}}{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_3}{i_2}}_{{a_1}{a_2}}"
          L"}} + "
          L"{{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_2}{i_3}}_{{a_1}{a_2}}}}"
          L"\\bigr) }");
    }

    {
      // -1/2 * A * g * t2
      const auto input =
          ex<Constant>(-0.5) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"a_1"}, WstrList{L"a_2", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2", L"a_3"}, WstrList{L"i_1", L"i_2"},
                         Symmetry::antisymm);
      //  std::wcout << "input:6  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ "
          L"\\bigl({{{2}}{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_"
          L"2}{a_3}}}} - "
          L"{{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_2}{a_3}}}}"
          L"\\bigr) }");
    }

    {
      // A * f * t2
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"f", WstrList{L"i_2"}, WstrList{L"a_2"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                         Symmetry::antisymm);
      //  std::wcout << "input:7  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(
          to_latex(result) == L"{ \\bigl({{{2}}{f^{{a_2}}_{{i_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_2}}}} - {{f^{{a_2}}_{{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}}\\bigr) }");
    }

    {
      // A * g * t1 * t1
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"a_1"}, WstrList{L"a_2", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) *
              ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});
      // std::wcout << "input:8  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ "
          L"\\bigl({{{2}}{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}{"
          L"t^{{i_1}}_{{a_3}}}} - "
          L"{{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_3}}}{t^{{i_1}}_{{"
          L"a_2}}}}\\bigr) }");
    }

    {
      // A * g * t2 * t2
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"i_1", L"a_2"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"});
      //  std::wcout << "input:9  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl({{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_3}}_{{a_2}}}} - {{{2}}{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}}_{{a_2}}}}\\bigr) }");
    }

    {
      // A * f * t1 * t1
      const auto input = ex<Constant>(-1) *
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"f", WstrList{L"i_2"}, WstrList{L"a_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_2"});
      //  std::wcout << "input:10  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl( - "
          L"{{f^{{a_2}}_{{i_2}}}{t^{{i_2}}_{{a_1}}}{t^{{i_1}}_{{a_2}}}}"
          L"\\bigr) }");
    }

    {
      // -1/2 * A * g * t1 * t2
      const auto input =
          ex<Constant>(-1. / 2.) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_2"}) *
              ex<Tensor>(L"t", WstrList{L"a_2", L"a_3"}, WstrList{L"i_1", L"i_3"},
                         Symmetry::antisymm);
      //  std::wcout << "input:11  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) == L"{ \\bigl({{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_3}{i_1}}_{{a_2}{a_3}}}} - {{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}{i_1}}_{{a_2}{a_3}}}}\\bigr) }");
    }

    {
      // -1/2 * A * g * t1 * t2
      const auto input =
          ex<Constant>(-1. / 2.) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_3"}, WstrList{L"i_2", L"i_3"},
                         Symmetry::antisymm);
      //  std::wcout << "input:12  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
//        std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl( - "
          L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_1}}_{{a_3}}}{t^{{i_3}"
          L"{i_2}}_{{a_1}{a_2}}}} + "
          L"{{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_1}}_{{a_3}}}{t^{{i_2}{i_3}"
          L"}_{{a_1}{a_2}}}}\\bigr) }");
    }

    {
      // A * g * t1 * t2
      const auto input =
          ex<Constant>(1) *
              ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
              ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"},
                         Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_3"}, WstrList{L"i_1", L"i_3"},
                         Symmetry::antisymm);
      //  std::wcout << "input:13  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ \\bigl( - {{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_2}}}{t^{{i_3}{i_1}}_{{a_1}{a_3}}}} + {{{4}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_2}}}{t^{{i_1}{i_3}}_{{a_1}{a_3}}}} + {{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_3}}}} - {{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_3}}}{t^{{i_1}{i_3}}_{{a_1}{a_2}}}}\\bigr) }");
    }

    {
      // - A * g * t1 * t1 * t1
      auto input = ex<Constant>(-1.) *
          ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"},
                     WstrList{L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"});
      //  std::wcout << "input:14  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
          L"{ "
          L"\\bigl({{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_"
          L"3}}_{{a_2}}}{t^{{i_1}}_{{a_3}}}} - "
          L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}"
          L"}_{{a_2}}}{t^{{i_1}}_{{a_3}}}}\\bigr) }");
    }
  }  // CCSD R1

  SECTION("Swap bra kets"){
    // Constant
    {
      auto input = ex<Constant>(0.5);
      auto result = swap_bra_ket(input);
      REQUIRE(result->is_atom());
      REQUIRE(result->is<Constant>());
      REQUIRE(result->to_latex() == L"{{{\\frac{1}{2}}}}");
    }

    // Tensor
    {
      auto input = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm);
      auto result = swap_bra_ket(input);
      REQUIRE(result->is_atom());
      REQUIRE(result->is<Tensor>());
      REQUIRE(result->to_latex() == L"{g^{{i_1}{i_2}}_{{a_1}{a_2}}}");
    }

    // Product
    {
      auto input = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
      auto result = swap_bra_ket(input);
      REQUIRE(result->size() == 2);
      REQUIRE(result->is<Product>());
      REQUIRE(result->to_latex() == L"{{g^{{a_5}{a_6}}_{{i_5}{i_6}}}{t^{{i_2}}_{{a_6}}}}");
    }

    // Sum
    {
      auto input = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
                   ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
      auto result = swap_bra_ket(input);
      REQUIRE(result->size() == 2);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->to_latex() == L"{ \\bigl({f^{{i_1}}_{{i_5}}} + {{g^{{a_5}{a_6}}_{{i_5}{i_6}}}{t^{{i_2}}_{{a_6}}}}\\bigr) }");
    }
  }

  SECTION("Check term"){
    { // A3 * f * t3
      auto input =  ex<Constant>(1./12) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"f", WstrList{L"i_4"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);

      std::wcout << input->size() << " " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = expand_antisymm(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      simplify(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n\n";
    }

    { // f * t3
      auto input = ex<Tensor>(L"f", WstrList{L"i_4"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                         WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);

      std::wcout << input->size() << " " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = expand_antisymm(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      simplify(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n\n";
    }

    { // A * g * t3
      auto input = ex<Constant>(-1./4) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"}, WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      std::wcout << input->size() << " " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = expand_antisymm(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      simplify(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n\n";
    }

    { // g * t3
      auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"}, WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      std::wcout << input->size() << " " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = expand_antisymm(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << result->size() << " " << to_latex(result) << "\n";
      simplify(result);
      std::wcout << result->size() << " " << to_latex(result) << "\n\n";
    }

  }


  SECTION("Open-shell spin-tracing"){
    // A * g

    // A * g * t1

    // A * f_oo * t2
    auto input = ex<Constant>(0.5) *
        ex<Tensor>(L"f", WstrList{L"i_3"}, WstrList{L"i_1"}) *
        ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"},
                         WstrList{L"i_2", L"i_3"}, Symmetry::antisymm);

    auto result = open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
    REQUIRE(result.size() == 4);
    for(auto& r_i : result){ // Check sizes for all terms
      std::wcout << to_latex(r_i) << std::endl;
    }
  }

#if 0
  SECTION("NON-ORTHOGONAL TRANSFORMATION"){
    //﻿http://dx.doi.org/10.1063/1.4907278
    // Intermediates
#if 1
    // 39
    auto W_ijam = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
        ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
    std::wcout << "39 W_ijam: " << to_latex(W_ijam) << "\n\n";

    // 36
    auto F_em = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"i_5"}) +
        (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
    std::wcout << "36 F_em: " << to_latex(F_em) << "\n\n";

    // 35
    auto F_ea = ex<Tensor>(L"f", WstrList{L"a_5"}, WstrList{L"a_1"}) -
        F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) -
        (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            ex<Tensor>(L"t", WstrList{L"i_6", L"i_5"}, WstrList{L"a_6", L"a_1"}) +
        (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm)) *
            ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_6"});
    std::wcout << "35 F_ea: " << to_latex(F_ea) << "\n\n";

    // 34
    auto F_im = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"i_5"}) +
        F_em->clone() * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) +
        (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_5"}) +
        (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_6"});
    std::wcout << "34 F_im: " << to_latex(F_im) << "\n\n";

    // 37
    auto W_ejmn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"});
    std::wcout << "37 W_ejmn: " << to_latex(W_ejmn) << "\n\n";

    // 37'
    auto W_eimn = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
    std::wcout << "W_eimn " << to_latex(W_eimn) << "\n\n";

    // 37"
//    auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) +
//        ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
    auto W_iemn = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"});
    std::wcout << "W_iemn " << to_latex(W_iemn) << "\n\n";

    // 45
//    auto W___eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) +
//        ex<Constant>(0.5) * ( ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
//                        (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}))-
//        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"});
//    std::wcout << "W___eima " << to_latex(W___eima) << "\n\n";

    // 43 // Antisymmetrized
    auto W_eima = ex<Tensor>(L"g", WstrList{L"a_5", L"i_1"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
        W_eimn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) +
        ex<Constant>(0.25) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
        ex<Constant>(0.25) *  ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);
    std::wcout << "W_eima " << to_latex(W_eima) << "\n\n";

    // 44
    auto W_iema = ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) -
        W_iemn->clone() * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"}) +
        ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_6"}) -
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"});
    std::wcout << "W_iema " << to_latex(W_iema) << "\n\n";

    // 49
    auto W_ijmn_temp = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) +
          ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) +
          ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) *
              (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}));
    auto W_ijmn_temp_c = W_ijmn_temp->clone();
    expand(W_ijmn_temp_c);
    std::wcout << "W_ijmn_temp_c: " << to_latex(W_ijmn_temp_c) << "\n\n";
    std::map<Index, Index> P_ij_mn;
    {
      Index i(L"i_1");
      Index j(L"i_2");
      Index m(L"i_5");
      Index n(L"i_6");
      P_ij_mn.emplace(std::make_pair(i, j));
      P_ij_mn.emplace(std::make_pair(j, i));
      P_ij_mn.emplace(std::make_pair(m, n));
      P_ij_mn.emplace(std::make_pair(n, m));
    }
    std::wcout << "W_ijmn_temp transformed: " << to_latex(transform_expression(W_ijmn_temp_c->clone(), P_ij_mn)) << "\n\n";

    auto W_ijmn = W_ijmn_temp_c->clone() + transform_expression(W_ijmn_temp_c->clone(), P_ij_mn);
    std::wcout << "W_ijmn " << to_latex(W_ijmn) << "\n\n";

    //
    auto temp_ab_ = W_iema->clone() * ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm);
    auto temp_ab_c = temp_ab_->clone();
    std::wcout << "temp_ab_c " << to_latex(temp_ab_c) << "\n\n";
    expand(temp_ab_c);

    std::map<Index, Index> P_ab;
    {
      Index a(L"a_1");
      Index b(L"a_2");
      P_ab.emplace(std::make_pair(a, b));
      P_ab.emplace(std::make_pair(b, a));
    }
    auto temp_ab = ex<Constant>(0.5) * temp_ab_c->clone() + transform_expression(temp_ab_c->clone(), P_ab);
    std::wcout << "temp_ab " << to_latex(temp_ab) << "\n\n";

    // CCSD Z2
    auto temp1 = ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) +
        ex<Tensor>(L"g", WstrList{L"i_1", L"a_5"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_5"}) -
        W_ijam->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}) +
        F_ea->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2"}, WstrList{L"a_5", L"a_2"}) -
        F_im->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_2"}, WstrList{L"a_1", L"a_2"}) +
        ex<Constant>(0.5) * (ex<Constant>(2.0) * W_eima->clone() - W_iema->clone()) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_5", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_2", L"a_5"}, Symmetry::nonsymm)) -
        temp_ab->clone() +
        ex<Constant>(0.5) * W_ijmn->clone() * (ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"}, Symmetry::nonsymm)) +
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) *
                     (ex<Tensor>(L"t", WstrList{L"i_1", L"i_2"}, WstrList{L"a_5", L"a_6"}, Symmetry::nonsymm) + ex<Tensor>(L"t", WstrList{L"i_1"}, WstrList{L"a_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}, Symmetry::nonsymm));

    auto temp1_c = temp1->clone();
    expand(temp1_c);
/*
    std::map<Index, Index> P_ab_ij;
    {
      Index a(L"a_1");
      Index b(L"a_2");
      Index i(L"i_1");
      Index j(L"i_2");
      P_ij_mn.emplace(std::make_pair(i, j));
      P_ij_mn.emplace(std::make_pair(j, i));
      P_ab.emplace(std::make_pair(a, b));
      P_ab.emplace(std::make_pair(b, a));
    }
*/
    auto ccsd_z2 = ex<Tensor>(L"S", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"}) * temp1_c->clone(); // + transform_expression(temp1, P_ab_ij);
    expand(ccsd_z2);
    std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";
    ccsd_z2 = swap_bra_ket(ccsd_z2);
    std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";
    canonicalize(ccsd_z2);
    rapid_simplify(ccsd_z2);
    std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";
//     ccsd_z2 = swap_bra_ket(ccsd_z2);
//     std::wcout << "Result: " << ccsd_z2->size() << "\n" << to_latex(ccsd_z2) << "\n";
#endif

    //////////////////////////////////
    ///    Required for T3 != 0    ///
    //////////////////////////////////

    auto expanded_T3 = [&] (container::svector<Index> bra, container::svector<Index> ket){
      assert(bra.size() == 3);
      assert(bra.size() == ket.size());

      auto t3_1 = Tensor(L"t", bra, ket);
      container::svector<Index> ket2 = {ket[1], ket[0], ket[2]};
      auto t3_2 = Tensor(L"t", bra, ket2);
      container::svector<Index> ket3 = {ket[2], ket[1], ket[0]};
      auto t3_3 = Tensor(L"t", bra, ket3);
      auto result = ex<Constant>(2.0) * ex<Tensor>(t3_1) - ex<Tensor>(t3_2) - ex<Tensor>(t3_3);
      return result;
    };

    // 38
    auto W_efam = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) -
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_1"});
    std::wcout << "W_efam " << to_latex(W_efam) << "\n\n";

    // 46'
    auto W___ejam = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"a_1", L"i_5"}, Symmetry::nonsymm) -
        ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"i_6", L"i_1"}, Symmetry::nonsymm);
    std::wcout << "W___ejam " << to_latex(W___ejam) << "\n\n";

    // 45'
    auto W___ejmb = ex<Tensor>(L"g", WstrList{L"a_5", L"i_2"}, WstrList{L"i_5", L"a_2"}, Symmetry::nonsymm) +
        ex<Constant>(0.5) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_6", L"i_1"}, WstrList{L"a_1", L"a_6"}, Symmetry::nonsymm)) -
        ex<Constant>(0.5) * ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_6", L"i_5"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_1", L"i_6"}, WstrList{L"a_6", L"a_1"}, Symmetry::nonsymm);
    std::wcout << "W___ejmb " << to_latex(W___ejmb) << "\n\n";

    // 41
    ExprPtr W_ejab;
    {

      auto g_tau = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) *
                   (ex<Tensor>(L"t", WstrList{L"i_2" L"i_5"}, WstrList{L"a_6", L"a_2"}) + ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"}));
      auto P_g_tau = ex<Constant>(0.5) * g_tau->clone() + transform_expression(g_tau->clone(), P_ab);

      Index a(L"a_1"),
          b(L"a_2"),
          e(L"a_5"),
          f(L"a_6"),
          i(L"i_1"),
          j(L"i_2"),
          m(L"i_5"),
          n(L"i_6");

      container::svector<Index> nmj = {n, m, j};
      container::svector<Index> fab = {f, a, b};

      W_ejab = ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) -
        W___ejam->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"}) -
        W___ejmb->clone() * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) +
        W_ejmn->clone() * (ex<Tensor>(L"t", WstrList{L"i_5", L"i_6"}, WstrList{L"a_1", L"a_2"}) + ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_1"}) * ex<Tensor>(L"t", WstrList{L"i_6"}, WstrList{L"a_2"})) +
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) * ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) +
        ex<Constant>(0.5) * (ex<Constant>(2.0) * ex<Tensor>(L"g", WstrList{L"a_6", L"a_5"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm) - ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"a_1"}, Symmetry::nonsymm)) *
            (ex<Constant>(2.0) * ex<Tensor>(L"t", WstrList{L"i_5", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) - ex<Tensor>(L"t", WstrList{L"i_2", L"i_5"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) -
                ex<Tensor>(L"t", WstrList{L"i_2"}, WstrList{L"a_6"}) * ex<Tensor>(L"t", WstrList{L"i_5"}, WstrList{L"a_2"})) -
        P_g_tau->clone() -
        ex<Tensor>(L"g", WstrList{L"a_5", L"a_6"}, WstrList{L"i_5", L"i_6"}, Symmetry::nonsymm) * expanded_T3(nmj, fab);
      std::wcout << "W_ejab " << to_latex(W_ejab) << "\n\n";
    }

    ExprPtr ccsdt_r2;
    {
      Index a(L"a_1"),
          b(L"a_2"),
          e(L"a_5"),
          f(L"a_6"),
          i(L"i_1"),
          j(L"i_2"),
          m(L"i_5"),
          n(L"i_6");

      container::svector<Index> mij = {m, i, j};
      container::svector<Index> min = {m, i, n};
      container::svector<Index> eab = {e, a, b};
      container::svector<Index> feb = {f, e, b};

//      std::wcout << "36 F_em: " << to_latex(F_em) << "\n\n";
//      std::wcout << "37 W_ejmn: " << to_latex(W_ejmn) << "\n\n";
//      std::wcout << "38 W_efam: " << to_latex(W_efam) << "\n\n";

      ccsdt_r2 = ex<Constant>(0.5) * F_em->clone() * expanded_T3(mij, eab) +
          W_efam->clone() * expanded_T3(mij, feb) -
          W_ejmn->clone() * expanded_T3(min, eab);
      std::wcout << "ccsdt_r2: " << ccsdt_r2->size() << " " << to_latex(ccsdt_r2) << "\n\n";
    }
    ccsdt_r2 = ccsd_z2 + ccsdt_r2;

    // 31
    ExprPtr ccsdt_r3;
    {
      Index a(L"a_1"),
       b(L"a_2"),
       c(L"a_3"),
       d(L"a_4"),
       e(L"a_5"),
       f(L"a_6"),

       i(L"i_1"),
       j(L"i_2"),
       k(L"i_3"),
       l(L"i_4"),
       m(L"i_5"),
       n(L"i_6");

      container::svector<Index> mjk = {m, j, k};
      container::svector<Index> ebc = {e, b, c};

      ccsdt_r3 = W_ejab->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_3"}, WstrList{L"a_5", L"a_3"}) -
          W_ijam->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_3"}, WstrList{L"a_2", L"a_3"}) +
          ex<Constant>(0.5) * F_ea->clone() * ex<Tensor>(L"t", WstrList{ L"i_1", L"i_2", L"i_3"}, WstrList{L"a_5", L"a_2", L"a_3"}) -
          ex<Constant>(0.5) * F_im->clone() * ex<Tensor>(L"t", WstrList{ L"i_5", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}) +
          ex<Constant>(0.25) * (ex<Constant>(2.0) * W_eima->clone() - W_iema->clone()) * expanded_T3(mjk, ebc);

    }
    expand(ccsdt_r3);
    canonicalize(ccsdt_r3);
    rapid_simplify(ccsdt_r3);
    std::wcout << "ccsdt_r3: " << ccsdt_r3->size() << " " << to_latex(ccsdt_r3) << "\n\n";

    ccsdt_r3 = ex<Tensor>(L"S", WstrList{L"a_1", L"a_2", L"a_3"}, WstrList{L"i_1", L"i_2", L"i_3"}) * ccsdt_r3->clone();
    expand(ccsdt_r3);
    canonicalize(ccsdt_r3);
    rapid_simplify(ccsdt_r3);
    std::wcout << "ccsdt_r3: " << ccsdt_r3->size() << " " << to_latex(ccsdt_r3) << "\n\n";

  } // NON-ORTHOGONAL TRANSFORMATION
#endif
}
