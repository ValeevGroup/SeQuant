//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.cpp"
#include "catch.hpp"

TEST_CASE("Spin") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  auto reset_idx_tags = [](ExprPtr& expr) {
    if (expr->is<Tensor>())
      ranges::for_each(expr->as<Tensor>().const_braket(),
                       [](const Index& idx) { idx.reset_tag(); });
  };

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
  SECTION("ASCII label") {
    auto p1 = Index(L"p⁺_1",
                    IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p2 =
        Index(L"p⁻_2", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    auto p3 = Index(L"p⁺_3",
                    IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p4 =
        Index(L"p⁻_4", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    auto alpha1 =
        Index(L"α⁺_1", IndexSpace::instance(IndexSpace::complete_unoccupied,
                                            IndexSpace::alpha));

    REQUIRE(p1.ascii_label() == "pa_1");
    REQUIRE(p2.ascii_label() == "pb_2");
    REQUIRE(p3.ascii_label() == "pa_3");
    REQUIRE(p4.ascii_label() == "pb_4");
    REQUIRE(alpha1.ascii_label() == "alphaa_1");
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

    auto spin_swap_tensor = swap_spin(input->as<Tensor>());
    REQUIRE(to_latex(spin_swap_tensor) == L"{t^{{p⁻_3}{p⁺_4}}_{{p⁻_1}{p⁺_2}}}");

    auto result = remove_spin(input);
    for (auto& i : result->as<Tensor>().const_braket())
      REQUIRE(i.space() ==
              IndexSpace::instance(IndexSpace::all, IndexSpace::nullqns));

    input = ex<Tensor>(L"t", IndexList{p1, p3}, IndexList{p2, p4});
    REQUIRE(to_latex(swap_spin(input)) == L"{t^{{p⁺_2}{p⁺_4}}_{{p⁻_1}{p⁻_3}}}");
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
    REQUIRE(to_latex(swap_spin(exprPtr)) == L"{{{\\frac{1}{4}}}}");
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
      REQUIRE(to_latex(result) == L"{\\bar{t}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}");

      input = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                         WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                         WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      result = expand_A_operator(input);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 36);
    }

#if 0
    // A3 * f * t3
    {
      auto input = ex<Constant>(1./12) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                           WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"f", WstrList{L"i_4"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      std::wcout << __LINE__  << " Input: " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      result->visit(reset_idx_tags);
      simplify(result);
      std::wcout << __LINE__ << " Result: " << to_latex(result) << "\n";
    }
#endif

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
              L"{ \\bigl({{\\bar{t}^{{i_1}{i_2}}_{{a_1}{a_2}}}} + "
              L"{{\\bar{t}^{{i_2}{i_1}}_{{a_2}{a_1}}}}\\bigr) }");
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
              L"{ \\bigl({{\\bar{t}^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}} + "
              L"{{\\bar{t}^{{i_1}{i_3}{i_2}}_{{a_1}{a_3}{a_2}}}} + "
              L"{{\\bar{t}^{{i_2}{i_1}{i_3}}_{{a_2}{a_1}{a_3}}}} + "
              L"{{\\bar{t}^{{i_2}{i_3}{i_1}}_{{a_2}{a_3}{a_1}}}} + "
              L"{{\\bar{t}^{{i_3}{i_1}{i_2}}_{{a_3}{a_1}{a_2}}}} + "
              L"{{\\bar{t}^{{i_3}{i_2}{i_1}}_{{a_3}{a_2}{a_1}}}}\\bigr) }");
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
              L"{ \\bigl({{\\bar{t}^{{i_1}{i_2}{i_3}{i_4}}_{{a_1}{a_2}{a_3}{a_4}}}} + "
              L"{{\\bar{t}^{{i_1}{i_2}{i_4}{i_3}}_{{a_1}{a_2}{a_4}{a_3}}}} + "
              L"{{\\bar{t}^{{i_1}{i_3}{i_2}{i_4}}_{{a_1}{a_3}{a_2}{a_4}}}} + "
              L"{{\\bar{t}^{{i_1}{i_3}{i_4}{i_2}}_{{a_1}{a_3}{a_4}{a_2}}}} + "
              L"{{\\bar{t}^{{i_1}{i_4}{i_2}{i_3}}_{{a_1}{a_4}{a_2}{a_3}}}} + "
              L"{{\\bar{t}^{{i_1}{i_4}{i_3}{i_2}}_{{a_1}{a_4}{a_3}{a_2}}}} + "
              L"{{\\bar{t}^{{i_2}{i_1}{i_3}{i_4}}_{{a_2}{a_1}{a_3}{a_4}}}} + "
              L"{{\\bar{t}^{{i_2}{i_1}{i_4}{i_3}}_{{a_2}{a_1}{a_4}{a_3}}}} + "
              L"{{\\bar{t}^{{i_2}{i_3}{i_1}{i_4}}_{{a_2}{a_3}{a_1}{a_4}}}} + "
              L"{{\\bar{t}^{{i_2}{i_3}{i_4}{i_1}}_{{a_2}{a_3}{a_4}{a_1}}}} + "
              L"{{\\bar{t}^{{i_2}{i_4}{i_1}{i_3}}_{{a_2}{a_4}{a_1}{a_3}}}} + "
              L"{{\\bar{t}^{{i_2}{i_4}{i_3}{i_1}}_{{a_2}{a_4}{a_3}{a_1}}}} + "
              L"{{\\bar{t}^{{i_3}{i_1}{i_2}{i_4}}_{{a_3}{a_1}{a_2}{a_4}}}} + "
              L"{{\\bar{t}^{{i_3}{i_1}{i_4}{i_2}}_{{a_3}{a_1}{a_4}{a_2}}}} + "
              L"{{\\bar{t}^{{i_3}{i_2}{i_1}{i_4}}_{{a_3}{a_2}{a_1}{a_4}}}} + "
              L"{{\\bar{t}^{{i_3}{i_2}{i_4}{i_1}}_{{a_3}{a_2}{a_4}{a_1}}}} + "
              L"{{\\bar{t}^{{i_3}{i_4}{i_1}{i_2}}_{{a_3}{a_4}{a_1}{a_2}}}} + "
              L"{{\\bar{t}^{{i_3}{i_4}{i_2}{i_1}}_{{a_3}{a_4}{a_2}{a_1}}}} + "
              L"{{\\bar{t}^{{i_4}{i_1}{i_2}{i_3}}_{{a_4}{a_1}{a_2}{a_3}}}} + "
              L"{{\\bar{t}^{{i_4}{i_1}{i_3}{i_2}}_{{a_4}{a_1}{a_3}{a_2}}}} + "
              L"{{\\bar{t}^{{i_4}{i_2}{i_1}{i_3}}_{{a_4}{a_2}{a_1}{a_3}}}} + "
              L"{{\\bar{t}^{{i_4}{i_2}{i_3}{i_1}}_{{a_4}{a_2}{a_3}{a_1}}}} + "
              L"{{\\bar{t}^{{i_4}{i_3}{i_1}{i_2}}_{{a_4}{a_3}{a_1}{a_2}}}} + "
              L"{{\\bar{t}^{{i_4}{i_3}{i_2}{i_1}}_{{a_4}{a_3}{a_2}{a_1}}}}\\bigr) }");
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
      auto result = factorize_S_operator(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}}, true);
      REQUIRE(result->is<Sum>() == false);
      REQUIRE(to_latex(result) == L"{{{2}}{S^{{a_1}{a_2}}_{{i_1}{i_2}}}{g^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_3}}_{{a_4}}}{t^{{i_4}}_{{a_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_3}}}}");
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

  SECTION("Closed-shell spintrace CCSD") {
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
      simplify(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      REQUIRE(
          to_latex(result) ==
              L"{ \\bigl( - {{f^{{i_1}}_{{i_2}}}{t^{{i_2}}_{{a_1}}}}\\bigr) }");
    }

    {
      // A * f * t1
      const auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"a_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_1"});
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
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
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      REQUIRE(to_latex(result) ==
          L"{ "
          L"\\bigl({{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_"
          L"3}}_{{a_2}}}{t^{{i_1}}_{{a_3}}}} - "
          L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}"
          L"}_{{a_2}}}{t^{{i_1}}_{{a_3}}}}\\bigr) }");
    }
  }  // CCSD R1

  SECTION("Closed-shell spintrace CCSDT terms"){
    { // A3 * f * t3
      auto input =  ex<Constant>(1./12) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"f", WstrList{L"i_4"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);

      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 36);
      result = expand_antisymm(result);
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      simplify(result);
      REQUIRE(result->size() == 4);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\bigl({{{2.000000}}{S^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{f^{"
              L"{i_3}}_{{i_4}}}{t^{{i_4}{i_1}{i_2}}_{{a_1}{a_2}{a_3}}}} - "
              L"{{{4.000000}}{S^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{f^{{i_3}}_"
              L"{{i_4}}}{t^{{i_1}{i_4}{i_2}}_{{a_1}{a_2}{a_3}}}} + "
              L"{{{4}}{S^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{f^{{i_3}}_{{i_4}}"
              L"}{t^{{i_1}{i_2}{i_4}}_{{a_1}{a_2}{a_3}}}} - "
              L"{{{2}}{S^{{a_1}{a_2}{a_3}}_{{i_1}{i_2}{i_3}}}{f^{{i_3}}_{{i_4}}"
              L"}{t^{{i_2}{i_1}{i_4}}_{{a_1}{a_2}{a_3}}}}\\bigr) }");
    }

    { // f * t3
      auto input = ex<Tensor>(L"f", WstrList{L"i_4"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                         WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);

      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 2);
      result = expand_antisymm(result);
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      REQUIRE(result->size() == 6);
      simplify(result);
      REQUIRE(result->size() == 6);
    }

    { // A * g * t3
      auto input = ex<Constant>(-1./4) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"}, WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 36);
      result = expand_antisymm(result);
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      REQUIRE(result->size() == 72);
      simplify(result);
      REQUIRE(result->size() == 20);
    }

    { // g * t3
      auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"}, WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"}, WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(result->size() == 2);
      result = expand_antisymm(result);
      result = closed_shell_spintrace(
          input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      REQUIRE(result->size() == 12);
      simplify(result);
      REQUIRE(result->size() == 12);
    }
  }

  SECTION("Merge P operators"){
      auto P1  = Tensor(L"P", WstrList{L"i_1", L"i_2"}, {});
      auto P2  = Tensor(L"P", {}, WstrList{L"a_1", L"a_2"});
      auto P3  = Tensor(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"});
      auto P4  = Tensor(L"P", {}, {});
      auto P12 = merge_operators(P1, P2);
      auto P34 = merge_operators(P3, P4);
      auto P11 = merge_operators(P1, P1);
      REQUIRE(to_latex(P12) == L"{P^{{a_1}{a_2}}_{{i_1}{i_2}}}");
      REQUIRE(to_latex(P34) == L"{P^{{a_1}{a_2}}_{{i_1}{i_2}}}");
      REQUIRE(to_latex(P11) == L"{P^{}_{{i_1}{i_2}{i_1}{i_2}}}");
  }

  SECTION("Permutation operators"){


      //      auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
      //                              WstrList{L"a_4", L"a_5"}, Symmetry::antisymm) *
      //                   ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"},
      //                              WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
      //                   ex<Tensor>(L"t", WstrList{L"a_2",L"a_3", L"a_5"},
      //                              WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);
      //
      //      auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
      //                              WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
      //                   ex<Tensor>(L"t", WstrList{L"a_1",L"a_2", L"a_3"},
      //                              WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);

      auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"},
                              WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
                                  ex<Tensor>(L"t", WstrList{L"a_2", L"a_3"},
                                             WstrList{L"i_3", L"i_4"}, Symmetry::antisymm);

      std::wcout << "Input: " << to_latex(input) << std::endl;
#if 0
      // Large index for dummys
      auto p2_input = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm) * input;
      std::wcout << "Input: " << to_latex(p2_input) << std::endl;
      auto p2_result = expand_P_operator(p2_input);
      std::wcout << "Result: " << to_latex(p2_result) << std::endl;
      assert(p2_result->size() == 2);
      //      assert(to_latex(p2_result) ==
      //          L"{ "
      //          L"\\bigl({{\\bar{g}^{{a_4}{a_5}}_{{i_4}{i_5}}}{\\bar{t}^{{i_1}{i_2}}_"
      //          L"{{a_1}{a_4}}}{\\bar{t}^{{i_3}{i_4}{i_5}}_{{a_2}{a_3}{a_5}}}} + "
      //          L"{{\\bar{g}^{{a_4}{a_5}}_{{i_4}{i_5}}}{\\bar{t}^{{i_2}{i_1}}_{{a_1}{"
      //          L"a_4}}}{\\bar{t}^{{i_3}{i_4}{i_5}}_{{a_2}{a_3}{a_5}}}}\\bigr) }");

      auto p3_input = (ex<Constant>(1.) - ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm)) * input;
      std::wcout << "\nInput: " << to_latex(p3_input) << std::endl;
      expand(p3_input);
      std::wcout << "\nInput: " << to_latex(p3_input) << std::endl;
      auto p3_result = expand_P_operator(p3_input, false);
      simplify(p3_input);
      std::wcout << "Result: " << to_latex(p3_result) << std::endl;

      auto p4_input = (ex<Constant>(1.) - ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm) -
          ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm)) * input;
      std::wcout << "\nInput: " << to_latex(p4_input) << std::endl;
      expand(p4_input);
      std::wcout << "\nInput: " << to_latex(p4_input) << std::endl;
      auto p4_result = expand_P_operator(p4_input, false);
      simplify(p4_input);
      std::wcout << "Result: " << to_latex(p4_result) << std::endl;

      auto p5_input = (ex<Constant>(1.) - ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm) -
          ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"}, WstrList{L"a_99", L"a_100"}, Symmetry::nonsymm)) *
              (ex<Constant>(1.) - ex<Tensor>(L"P", WstrList{L"i_98", L"i_99"}, WstrList{L"a_1", L"a_2"}, Symmetry::nonsymm) -
              ex<Tensor>(L"P", WstrList{L"i_98", L"i_99"}, WstrList{L"a_1", L"a_3"}, Symmetry::nonsymm)) * input;
      std::wcout << "\nInput: " << to_latex(p5_input) << std::endl;
      expand(p5_input);
      std::wcout << "\nInput: " << to_latex(p5_input) << std::endl;
      auto p5_result = expand_P_operator(p5_input, false);
      p5_result->visit(reset_idx_tags);
      simplify(p5_result);
      std::wcout << "Result: " << to_latex(p5_result) << std::endl;

      p5_result = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) * p5_result;
      expand(p5_result);
      p5_result = expand_A_operator(p5_result);
      p5_result->visit(reset_idx_tags);
      simplify(p5_result);
      std::wcout << "p5_result:\n" << to_latex(p5_result) << "\n";
#endif

      auto A_12 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
      auto A_23 = ex<Tensor>(L"A", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"}, Symmetry::antisymm);
      auto A2 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
      auto A3 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm);
      auto A4 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3", L"i_4"}, WstrList{L"a_1", L"a_2", L"a_3", L"a_4"}, Symmetry::antisymm);
      auto A5 = ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"},
                           WstrList{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"}, Symmetry::antisymm);

      {

        auto Avec2 = open_shell_A_op(A2->as<Tensor>());
        auto P3 = open_shell_P_op_vector(A3->as<Tensor>());
        auto Avec3 = open_shell_A_op(A3->as<Tensor>());
        assert(P3[0]->size() == 0);
        assert(P3[1]->size() == 9);
        assert(P3[2]->size() == 9);
        assert(P3[3]->size() == 0);

        auto P4 = open_shell_P_op_vector(A4->as<Tensor>());
        auto Avec4 = open_shell_A_op(A4->as<Tensor>());
        assert(P4[0]->size() == 0);
        assert(P4[1]->size() == 16);
        assert(P4[2]->size() == 36);
        assert(P4[3]->size() == 16);
        assert(P4[4]->size() == 0);
        std::wcout << "P4[2]: " << to_latex(P4[2]) << std::endl;

        auto P5 = open_shell_P_op_vector(A5->as<Tensor>());
        auto Avec5 = open_shell_A_op(A5->as<Tensor>());
        assert(P5[0]->size() == 0);
        assert(P5[1]->size() == 25);
        assert(P5[2]->size() == 100);
        assert(P5[3]->size() == 100);
        assert(P5[4]->size() == 25);
        assert(P5[5]->size() == 0);

      }

#if 0
{
        auto P13_b = Tensor(L"P", {}, WstrList{L"a_1", L"a_3"},
                            Symmetry::nonsymm);
        auto P13_k = Tensor(L"P", WstrList{L"i_1", L"i_3"}, {},
                            Symmetry::nonsymm);
        auto P12_b = Tensor(L"P", {}, WstrList{L"a_1", L"a_2"},
                            Symmetry::nonsymm);
        auto P12_k = Tensor(L"P", WstrList{L"i_1", L"i_2"}, {},
                            Symmetry::nonsymm);

        auto P23_b = Tensor(L"P", {}, WstrList{L"a_2", L"a_3"},
                            Symmetry::nonsymm);
        auto P23_k = Tensor(L"P", WstrList{L"i_2", L"i_3"}, {},
                            Symmetry::nonsymm);

        // aab expr
        auto expr = (ex<Constant>(1.) - ex<Tensor>(P13_b) - ex<Tensor>(P23_b)) * (ex<Constant>(1) - ex<Tensor>(P13_k) - ex<Tensor>(P23_k));
        // abb expr
        // auto expr = (ex<Constant>(1.) - ex<Tensor>(P12_b) - ex<Tensor>(P13_b)) * (ex<Constant>(1) - ex<Tensor>(P12_k) - ex<Tensor>(P13_k));
        std::wcout << __LINE__ << "L " << to_latex(expr) << "\n";
        simplify(expr);
        std::wcout << __LINE__ << "L " << to_latex(expr) << "\n";

        for(auto& term : *expr){
          if(term->is<Product>()){
            std::vector<Tensor> t_list;
            for(auto& t : *term){
              t_list.push_back(t->as<Tensor>());
            }
            if(t_list.size() == 2){
              term = merge_operators(t_list[0], t_list[1]);
            }
          }
        }
        std::wcout << __LINE__ << "L " << to_latex(expr) << "\n";
        auto expr1 = expr * input;
        expand(expr1);
        rapid_simplify(expr1);
        std::wcout << __LINE__ << "L " << to_latex(expr1) << "\n";
        expr1 = expand_P_operator(expr1, false);
        expr1->visit(reset_idx_tags);
        expand(expr1);
        simplify(expr1);
        std::wcout << __LINE__ << "L " << expr1->size() << " " << to_latex(expr1) << "\n";

        expr1 = A_12 * expr1;
        expand(expr1);
        canonicalize(expr1);
        std::wcout << "expr1:\n" << to_latex(expr1) << "\n";
        expr1 = expand_A_operator(expr1);
        expr1->visit(reset_idx_tags);
        simplify(expr1);
        std::wcout << "expr1:\n" << to_latex(expr1) << "\n";



        //        auto result = expand_P_operator(expr, false);
        //        simplify(result);
        //        std::wcout << __LINE__ << "L " << to_latex(result) << "\n";
      }
#endif

auto P13_b = ex<Tensor>(L"P", WstrList{},
                        WstrList{L"a_1", L"a_3"},
                        Symmetry::nonsymm);
      auto P13_k = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                              WstrList{},
                              Symmetry::nonsymm);
      auto P12_b = ex<Tensor>(L"P", WstrList{},
                              WstrList{L"a_1", L"a_2"},
                              Symmetry::nonsymm);
      auto P12_k = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                              WstrList{},
                              Symmetry::nonsymm);

      auto P23_b = ex<Tensor>(L"P", WstrList{},
                              WstrList{L"a_2", L"a_3"},
                              Symmetry::nonsymm);
      auto P23_k = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                              WstrList{},
                              Symmetry::nonsymm);

      auto P4_1313 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                                WstrList{L"a_1", L"a_3"},
                                Symmetry::nonsymm);
      auto P4_1323 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                                WstrList{L"a_2", L"a_3"},
                                Symmetry::nonsymm);
      auto P4_2313 = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                                WstrList{L"a_1", L"a_3"},
                                Symmetry::nonsymm);
      auto P4_2323 = ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                                WstrList{L"a_2", L"a_3"},
                                Symmetry::nonsymm);

      auto P4_1212 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                                WstrList{L"a_1", L"a_2"},
                                Symmetry::nonsymm);
      auto P4_1213 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_2"},
                                WstrList{L"a_1", L"a_3"},
                                Symmetry::nonsymm);
      auto P4_1312 = ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                                WstrList{L"a_1", L"a_2"},
                                Symmetry::nonsymm);

      auto p_aab = ex<Constant>(1.) - P13_b - P23_b - P13_k - P23_k +
          P4_1313 + P4_1323 + P4_2313 + P4_2323;
      std::wcout << "\np_aab: " << to_latex(p_aab) << std::endl;

      auto p_abb = ex<Constant>(1.) - P13_b - P12_b - P13_k - P12_k +
          P4_1212 + P4_1213 + P4_1312 + P4_1313;
      std::wcout << "\np_abb: " << to_latex(p_abb) << std::endl;

      auto p6_input = p_aab * input;
      expand(p6_input);
      auto p6_result = expand_P_operator(p6_input, false);
      p6_result->visit(reset_idx_tags);
      simplify(p6_result);
      std::wcout << __LINE__ << " p6_result: " << p6_result->size() << " " << to_latex(p6_result) << "\n";

      p6_result = A_12 * p6_result;
      expand(p6_result);
      canonicalize(p6_result);
      std::wcout << "p6_result:\n" << to_latex(p6_result) << "\n";
      p6_result = expand_A_operator(p6_result);
      p6_result->visit(reset_idx_tags);
      simplify(p6_result);
      std::wcout << "p6_result:\n" << to_latex(p6_result) << "\n";

      auto p7_input = p_abb * input;
      expand(p7_input);
      auto p7_result = expand_P_operator(p7_input, false);
      p7_result->visit(reset_idx_tags);
      simplify(p7_result);

      p7_result = A_23 * p7_result;
      expand(p7_result);
      p7_result = expand_A_operator(p7_result);
      p7_result->visit(reset_idx_tags);
      simplify(p7_result);
      std::wcout << "p7_result:\n" << to_latex(p7_result) << "\n";

      auto expanded_A = A3 * input;
      expanded_A = expand_A_operator(expanded_A);
      expanded_A->visit(reset_idx_tags);
      simplify(expanded_A);
      std::wcout << "expanded_A:\n" << to_latex(expanded_A) << "\n";
      assert(to_latex(p6_result) == to_latex(p7_result));
      assert(to_latex(p6_result) == to_latex(expanded_A));

    }

  SECTION("Expand P operator pair-wise") {
    auto P1 = Tensor(L"P", WstrList{L"i_1", L"i_2"}, {});
    auto P2 = Tensor(L"P", WstrList{L"i_1", L"i_2", L"i_3", L"i_4"}, {});
    auto P3 = Tensor(L"P", {}, WstrList{L"a_1", L"a_2"});
    auto P4 = Tensor(L"P", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"});
    auto P5 = Tensor(L"P", WstrList{L"i_1", L"i_2", L"i_3", L"i_4"},
                     WstrList{L"a_1", L"a_2", L"a_3", L"a_4"});

    // g* t3
    auto input = ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"},
                            WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
                 ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"},
                            WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);

    size_t n_p = 0;
    for (auto& P : {P1, P2, P3, P4, P5}) {
      auto term = ex<Tensor>(P) * input;
      expand(term);
      auto result = expand_P_operator(term);
      switch (n_p) {
        case 0 : REQUIRE(to_latex(result) ==
                  L"{ "
                  L"\\bigl({{\\bar{g}^{{i_2}{a_4}}_{{i_4}{a_1}}}{\\bar{t}^{{i_"
                  L"1}{i_3}{i_4}}_{{a_2}{a_3}{a_4}}}}\\bigr) }");
                break;
        case 1: REQUIRE(to_latex(result) ==
                  L"{ "
                  L"\\bigl({{\\bar{g}^{{i_2}{a_4}}_{{i_3}{a_1}}}{\\bar{t}^{{i_"
                  L"1}{i_4}{i_3}}_{{a_2}{a_3}{a_4}}}}\\bigr) }");
                break;
        case 2: REQUIRE(to_latex(result) ==
                  L"{ "
                  L"\\bigl({{\\bar{g}^{{i_1}{a_4}}_{{i_4}{a_2}}}{\\bar{t}^{{i_"
                  L"2}{i_3}{i_4}}_{{a_1}{a_3}{a_4}}}}\\bigr) }");
                break;
        case 3: REQUIRE(to_latex(result) ==
                  L"{ "
                  L"\\bigl({{\\bar{g}^{{i_2}{a_4}}_{{i_4}{a_2}}}{\\bar{t}^{{i_"
                  L"1}{i_3}{i_4}}_{{a_1}{a_3}{a_4}}}}\\bigr) }");
                break;
        case 4: REQUIRE(to_latex(result) ==
                  L"{ "
                  L"\\bigl({{\\bar{g}^{{i_2}{a_3}}_{{i_3}{a_2}}}{\\bar{t}^{{i_"
                  L"1}{i_4}{i_3}}_{{a_1}{a_4}{a_3}}}}\\bigr) }");
                 break;
        default: break;
      }
      ++n_p;
    }

    auto input2 = ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"},
                            WstrList{L"i_1", L"i_2"}, Symmetry::antisymm);
    auto term = ex<Tensor>(P2) * input2;
    expand(term);
    auto result = expand_P_operator(term);
    REQUIRE(to_latex(result) ==
            L"{ \\bigl({{\\bar{g}^{{i_2}{i_1}}_{{a_1}{a_2}}}}\\bigr) }");
}

#define open_shell_tests 1
#if open_shell_tests
  SECTION("Open-shell spin-tracing"){

      auto occA = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::alpha);
      auto virA = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::alpha);
      auto occB = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta);
      auto virB = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::beta);
      const auto i1A = Index(L"i⁺_1", occA);
      const auto i2A = Index(L"i⁺_2", occA);
      const auto i3A = Index(L"i⁺_3", occA);
      const auto i4A = Index(L"i⁺_4", occA);
      const auto i5A = Index(L"i⁺_5", occA);
      const auto i1B = Index(L"i⁻_1", occB);
      const auto i2B = Index(L"i⁻_2", occB);
      const auto i3B = Index(L"i⁻_3", occB);

      const auto a1A = Index(L"a⁺_1", virA);
      const auto a2A = Index(L"a⁺_2", virA);
      const auto a3A = Index(L"a⁺_3", virA);
      const auto a1B = Index(L"a⁻_1", virB);
      const auto a2B = Index(L"a⁻_2", virB);
      const auto a3B = Index(L"a⁻_3", virB);

    // Logger::get_instance().canonicalize = true;
    // Tensor canonicalize
    {
      auto t3 = ex<Tensor>(Tensor(L"t", {a3A, a2B, a2A}, {i1A, i2B, i3A}));
      auto f = ex<Tensor>(Tensor(L"f",{a1A},{a2A}));
      auto ft3 = f*t3;
      ft3->canonicalize();
      REQUIRE(to_latex(ft3) == L"{{f^{{a⁺_2}}_{{a⁺_1}}}{t^{{i⁺_3}{i⁺_1}{i⁻_2}}_{{a⁺_2}{a⁺_3}{a⁻_2}}}}");
    }

    //  g
    {
      auto input =  ex<Constant>(0.25) *
                    ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"},
                               WstrList{L"i_1", L"i_2"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
      REQUIRE(result.size() == 3);
      REQUIRE(to_latex(result[0]) == L"{{{\\frac{1}{4}}}{\\bar{g}^{{i⁺_1}{i⁺_2}}_{{a⁺_1}{a⁺_2}}}}");
      REQUIRE(to_latex(result[1]) == L"{{{\\frac{1}{4}}}{g^{{i⁺_1}{i⁻_2}}_{{a⁺_1}{a⁻_2}}}}");
      REQUIRE(to_latex(result[2]) == L"{{{\\frac{1}{4}}}{\\bar{g}^{{i⁻_1}{i⁻_2}}_{{a⁻_1}{a⁻_2}}}}");
    }

    // f_oo * t2
    {
      auto input = ex<Constant>(0.5) *
                   ex<Tensor>(L"f", WstrList{L"i_3"}, WstrList{L"i_1"}) *
                   ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"},
                              WstrList{L"i_2", L"i_3"}, Symmetry::antisymm);

      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
      REQUIRE(result.size() == 3);
      REQUIRE(to_latex(result[0]) == L"{{{\\frac{1}{2}}}{f^{{i⁺_1}}_{{i⁺_3}}}{\\bar{t}^{{i⁺_2}{i⁺_3}}_{{a⁺_1}{a⁺_2}}}}");
      REQUIRE(to_latex(result[1]) == L"{{{-\\frac{1}{2}}}{f^{{i⁺_1}}_{{i⁺_2}}}{t^{{i⁺_2}{i⁻_2}}_{{a⁺_1}{a⁻_2}}}}");
      REQUIRE(to_latex(result[2]) == L"{{{\\frac{1}{2}}}{f^{{i⁻_1}}_{{i⁻_3}}}{\\bar{t}^{{i⁻_2}{i⁻_3}}_{{a⁻_1}{a⁻_2}}}}");
    }

    // g * t1
    {
      auto input = ex<Constant>(0.5) *
          ex<Tensor>(L"g", WstrList{L"i_3", L"a_1"},
                     WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_3"}, Symmetry::nonsymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
      REQUIRE(result.size() == 3);
      REQUIRE(to_latex(result[0]) ==
              L"{{{\\frac{1}{2}}}{\\bar{g}^{{i⁺_1}{i⁺_2}}_{{i⁺_3}{a⁺_1}}}{t^{{"
              L"i⁺_3}}_{{a⁺_2}}}}");
      REQUIRE(to_latex(result[1]) ==
              L"{{{-\\frac{1}{2}}}{g^{{i⁺_1}{i⁻_2}}_{{a⁺_1}{i⁻_1}}}{t^{{i⁻_1}}_"
              L"{{a⁻_2}}}}");
      REQUIRE(to_latex(result[2]) ==
              L"{{{\\frac{1}{2}}}{\\bar{g}^{{i⁻_1}{i⁻_2}}_{{i⁻_3}{a⁻_1}}}{t^{{"
              L"i⁻_3}}_{{a⁻_2}}}}");
    }

    Logger::get_instance().canonicalize = false;
    // f * t3
    {
      auto input = ex<Constant>(1./12) *
          ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"a_4"}) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << "Input: " << to_latex(input) << "\n"
                 << "Results:\n";
      for(const auto &r : result) { std::wcout << "\t" << to_latex(r) << "\n"; }
      std::wcout << "\n";
      REQUIRE(result.size() == 4);
    }
    Logger::get_instance().canonicalize = false;

    // g*t3 (CCSDT R1 1)
    {
      auto input = ex<Constant>(0.25) *
          ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"},
                     WstrList{L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}});
      REQUIRE(result.size() == 2);
      REQUIRE(to_latex(result[0]) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{4}}}{\\bar{g}^{{a⁺_2}{a⁺_3}}_{{i⁺_2}{i⁺_3}}"
              L"}{\\bar{t}^{{i⁺_1}{i⁺_2}{i⁺_3}}_{{a⁺_1}{a⁺_2}{a⁺_3}}}} + "
              L"{{g^{{a⁺_2}{a⁻_1}}_{{i⁺_2}{i⁻_1}}}{t^{{i⁺_1}{i⁺_2}{i⁻_1}}_{{a⁺_"
              L"1}{a⁺_2}{a⁻_1}}}} - "
              L"{{g^{{a⁺_2}{a⁻_1}}_{{i⁺_2}{i⁻_1}}}{t^{{i⁺_2}{i⁺_1}{i⁻_1}}_{{a⁺_"
              L"1}{a⁺_2}{a⁻_1}}}} + "
              L"{{{\\frac{1}{2}}}{\\bar{g}^{{a⁻_1}{a⁻_2}}_{{i⁻_1}{i⁻_2}}}{t^{{"
              L"i⁺_1}{i⁻_1}{i⁻_2}}_{{a⁺_1}{a⁻_1}{a⁻_2}}}}\\bigr) }");
      REQUIRE(to_latex(result[1]) ==
              L"{ "
              L"\\bigl({{{\\frac{1}{4}}}{\\bar{g}^{{a⁻_2}{a⁻_3}}_{{i⁻_2}{i⁻_3}}"
                       L"}{\\bar{t}^{{i⁻_1}{i⁻_2}{i⁻_3}}_{{a⁻_1}{a⁻_2}{a⁻_3}}}}"
                       L" + {{{\\frac{1}{2}}}{\\bar{g}^{{a⁺_1}{a⁺_2}}_{{i⁺_1}{i⁺"
                       L"_2}}}{t^{{i⁺_1}{i⁺_2}{i⁻_1}}_{{a⁺_1}{a⁺_2}{a⁻_1}}}} + "
                       L"{{g^{{a⁺_1}{a⁻_2}}_{{i⁺_1}{i⁻_2}}}{t^{{i⁺_1}{i⁻_1}{i⁻_2"
                       L"}}_{{a⁺_1}{a⁻_1}{a⁻_2}}}} - {{g^{{a⁺_1}{a⁻_2}}_{{i⁺_1}"
                       L"{i⁻_2}}}{t^{{i⁺_1}{i⁻_2}{i⁻_1}}_{{a⁺_1}{a⁻_1}{a⁻_2}}}"
                       L"}\\bigr) }");
    }

    // g*t3 (CCSDT R2 8)
    {
      auto input = ex<Constant>(0.25) *
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                     WstrList{L"i_1", L"a_3"}, Symmetry::antisymm) *
                     ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                                WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
      std::wcout << "Input: " << to_latex(input) << "\n"
      << "Results:\n";
      for(const auto &r : result) {
        std::wcout << "\t" << to_latex(r) << "\n";
      }
      std::wcout << "\n";
      REQUIRE(result.size() == 3);
    }

    // g*t1*t3 (CCSDT R2 20)
    {
      auto input = ex<Constant>(0.25) *
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
                     WstrList{L"a_3", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}, Symmetry::nonsymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_4"},
                     WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
      std::wcout << "Input: " << to_latex(input) << "\n"
      << "Results:\n";
      for(const auto &r : result) {
        std::wcout << "\t" << to_latex(r) << "\n";
      }
      std::wcout << "\n";
      REQUIRE(result.size() == 3);
    }

    // g * t3 (CCSDT R3 1)
    {
      auto input = ex<Constant>(-0.25) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"a_1"},
                           WstrList{L"i_1", L"a_4"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_3", L"a_4"},
                     WstrList{L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << "CCSDT R3 1 Input: " << to_latex(input) << "\n"
                 << "Results:\n";
      for(const auto &r : result) {
        std::wcout << "\t" << to_latex(r) << "\n";
      }
      std::wcout << "\n";
      REQUIRE(result.size() == 4);

      auto A3_aaa = Tensor(L"A",{i1A,i2A,i3A},{a1A,a2A,a3A},Symmetry::antisymm);
      auto result2 = ex<Tensor>(A3_aaa) * result[0];
      expand(result2);
      result2 = expand_A_operator(result2);
      result2->visit(reset_idx_tags);
      canonicalize(result2);
      rapid_simplify(result2);
      std::wcout << "Result2: " << result2->size() << " " << to_latex(result2) << "\n";

    }

    // g * t3 (CCSDT R3 4)
    {
      auto input = ex<Constant>(1./24) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_3"},
                           WstrList{L"a_1", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                           WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1",L"a_2", L"a_3"},
                           WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm) +
         ex<Constant>(1./24) *
         ex<Tensor>(L"A", WstrList{L"i_2", L"i_3"},
                   WstrList{L"a_2", L"a_3"}, Symmetry::antisymm) *
        ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                   WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
       ex<Tensor>(L"t", WstrList{L"a_1",L"a_2", L"a_3"},
                   WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});

      //      auto result = expand_A_operator(input);
      //      result->visit(reset_idx_tags);
      //      canonicalize(result);
      //      rapid_simplify(result);
      std::wcout << __LINE__ << " Input: " << to_latex(input) << "\n"
      << "Result:\n";
      std::wcout << "AAB:\t" << to_latex(result[1]) << "\n";

      input = ex<Tensor>(Tensor(L"A",{i1A,i2A},{a1A,a2A},Symmetry::antisymm)) * result[1];
      expand(input);
      auto result2 = expand_A_operator(input);
      result2->visit(reset_idx_tags);
      canonicalize(result2);
      rapid_simplify(result2);
      std::wcout << "result2: " << to_latex(result2) << "\n";
    }

    // g * t3 (CCSDT R3 4)
    {
      auto input = ex<Constant>(1./24) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"},
                     WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                    WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1",L"a_2", L"a_3"},
                     WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm) +
          ex<Constant>(1./24) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_3"},
                    WstrList{L"a_1", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                   WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1",L"a_2", L"a_3"},
                    WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});

      //      auto result = expand_A_operator(input);
      //      result->visit(reset_idx_tags);
      //      canonicalize(result);
      //      rapid_simplify(result);
      std::wcout << __LINE__ << " Input: " << to_latex(input) << "\n"
                 << "Result:\n";
      std::wcout << "ABB:\t" << to_latex(result[2]) << "\n";

      // Multiply A(same_spin)
      auto occB = IndexSpace::instance(IndexSpace::active_occupied, IndexSpace::beta);
      auto virB = IndexSpace::instance(IndexSpace::active_unoccupied, IndexSpace::beta);
      const auto i2B = Index(L"i⁻_2", occB);
      const auto i3B = Index(L"i⁻_3", occB);
      const auto a3B = Index(L"a⁻_3", virB);
      const auto a2B = Index(L"a⁻_2", virB);

      input = ex<Tensor>(Tensor(L"A",{i2B,i3B},{a2B,a3B},Symmetry::antisymm)) * result[2];
      expand(input);
      auto result2 = expand_A_operator(input);
      result2->visit(reset_idx_tags);
      canonicalize(result2);
      rapid_simplify(result2);
      std::wcout << "result2: " << to_latex(result2) << "\n";


    }

    // aab: g*t3 (CCSDT R3 4)
    {
      auto A2_aab = Tensor(L"A", {i1A, i2A}, {a1A, a2A}, Symmetry::antisymm);
      auto A2_abb = Tensor(L"A", {i2B, i3B}, {a2B, a3B}, Symmetry::antisymm);

      auto g = Tensor(L"g", {i3A, i4A}, {i1A,i2A}, Symmetry::antisymm);
      auto t3 = Tensor(L"t", {a1A, a2A, a3B}, {i3A,i4A,i3B}, Symmetry::nonsymm);

      auto input = ex<Constant>(1./12) * ex<Tensor>(A2_aab) * ex<Tensor>(g) * ex<Tensor>(t3);
      auto result = expand_A_operator(input);
      result->visit(reset_idx_tags);
      canonicalize(result);
      rapid_simplify(result);
      std::wcout << "CCSDT R3 4 Input: " << to_latex(input) << "\n"
                 << "Result:\n" << to_latex(result) << "\n";

      g = Tensor(L"g", {i4A, i5A}, {i1A,i2A}, Symmetry::antisymm);
      t3 = Tensor(L"t", {a1A, a2A, a3B}, {i4A,i5A,i3B}, Symmetry::nonsymm);

      input = ex<Constant>(1./12) * ex<Tensor>(A2_aab) * ex<Tensor>(g) * ex<Tensor>(t3);
      result = expand_A_operator(input);
      result->visit(reset_idx_tags);
      canonicalize(result);
      rapid_simplify(result);
      std::wcout << "Input: " << to_latex(input) << "\n"
      << "Result:\n" << to_latex(result) << "\n";

    }

    // CCSDT R3 10 aaa, bbb
    {
      auto input = ex<Constant>(1./8) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                     WstrList{L"a_4", L"a_5"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"},
                     WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2",L"a_3", L"a_5"},
                     WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);

      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      std::wcout << __LINE__ << " Input: " << to_latex(input) << "\n"
                 << "Result_aaa:\n" << result[0]->size();
      std::wcout << "\t" << to_latex(result[0]) << "\n";

      auto A3_aaa = Tensor(L"A",{i1A,i2A,i3A},{a1A,a2A,a3A},Symmetry::antisymm);
      auto result2 = ex<Tensor>(A3_aaa) * result[0];
      expand(result2);
      result2 = expand_A_operator(result2);
      result2->visit(reset_idx_tags);
      canonicalize(result2);
      rapid_simplify(result2);
      std::wcout << "Result2: " << result2->size() << " " << to_latex(result2) << "\n";

      auto A3_bbb = Tensor(L"A",{i1B,i2B,i3B},{a1B,a2B,a3B},Symmetry::antisymm);
      auto result3 = ex<Tensor>(A3_bbb) * result[3];
      expand(result3);
      result3 = expand_A_operator(result3);
      result3->visit(reset_idx_tags);
      canonicalize(result3);
      rapid_simplify(result3);
      std::wcout << "Result3: " << result3->size() << " " << to_latex(result3) << "\n";

    }

    // CCSDT R3 10 aab
    {
      auto input = ex<Constant>(1./8) *
          ex<Tensor>(L"P", WstrList{L"i_1", L"i_3"},
                     WstrList{L"a_1", L"a_3"}, Symmetry::nonsymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                    WstrList{L"a_4", L"a_5"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"},
                     WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2",L"a_3", L"a_5"},
                     WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm) +
          ex<Constant>(1./8) *
          ex<Tensor>(L"P", WstrList{L"i_2", L"i_3"},
                     WstrList{L"a_2", L"a_3"}, Symmetry::nonsymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                     WstrList{L"a_4", L"a_5"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"},
                      WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2",L"a_3", L"a_5"},
                     WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);

      input = expand_P_operator(input);
      input->visit(reset_idx_tags);
      auto result =
          open_shell_spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      for(auto &r : result) std::cout << r->size() << " ";
      std::cout << "\n";

      auto result_aab = ex<Tensor>(Tensor(L"A", {i1A, i2A},
                                   {a1A, a2A}, Symmetry::antisymm)) * result[1];
      expand(result_aab);
      result_aab = expand_A_operator(result_aab);
      result_aab->visit(reset_idx_tags);
      canonicalize(result_aab);
      rapid_simplify(result_aab);
      REQUIRE(result_aab->size() == 18);

      auto input2 = ex<Constant>(1./8) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_4", L"i_5"},
                     WstrList{L"a_4", L"a_5"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_4"},
                     WstrList{L"i_1", L"i_2"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2",L"a_3", L"a_5"},
                     WstrList{L"i_3", L"i_4",L"i_5"}, Symmetry::antisymm);

      auto result2 = open_shell_spintrace(input2, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}, {L"i_3", L"a_3"}});
      REQUIRE(result2[1]->size() == 24);
    }

  }
#endif



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
