//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.hpp"
#include "catch.hpp"

#define MANDATORY_TESTS 1
#define CCSD_r1 1

TEST_CASE("Spin") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  SECTION("Tensor spintrace") {

    {
      const auto input = ex<Constant>(1./36.) *
          ex<Tensor>(L"A",WstrList{L"i_1", L"i_2", L"i_3"}, WstrList{L"a_1", L"a_2", L"a_3"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      // std::wcout << "input:    " << to_latex(input) << "\n";
      auto asm_input = expand_A_operator(input);
      // std::wcout << "expanded: " << to_latex(asm_input) << "\n";
      REQUIRE(asm_input->size() == 36);
    }

    {
      const auto input = ex<Constant>(1./576.) *
          ex<Tensor>(L"A",WstrList{L"i_1", L"i_2", L"i_3", L"i_4"}, WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3", L"a_4"},
                     WstrList{L"i_1", L"i_2", L"i_3", L"i_4"}, Symmetry::antisymm);
      // std::wcout << "input:    " << to_latex(input) << "\n";
      auto asm_input = expand_A_operator(input);
      // std::wcout << "expanded: " << to_latex(asm_input) << "\n";
      REQUIRE(asm_input->size() == 576);
    }

    {
      const auto input = ex<Constant>(1./14400.) *
          ex<Tensor>(L"A",WstrList{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"}, WstrList{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3", L"a_4", L"a_5"},
                     WstrList{L"i_1", L"i_2", L"i_3", L"i_4", L"i_5"}, Symmetry::antisymm);
      // std::wcout << "input:    " << to_latex(input) << "\n";
      auto asm_input = expand_A_operator(input);
      // std::wcout << "expanded: " << to_latex(asm_input) << "\n";
      REQUIRE(asm_input->size() == 14400);
    }


    {
      const auto input =
          ex<Constant>(1. / 2.) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_3", L"a_1"}, WstrList{L"a_3", L"a_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_3"}) *
          ex<Tensor>(L"t", WstrList{L"a_2", L"a_4"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm);
      // std::wcout << "input:  " << to_latex(input) << "\n";
      auto asm_input = expand_A_operator(input);
      // std::wcout << "antisy: " << to_latex(asm_input) << "\n";
      for(auto& term : *asm_input){
        std::wcout << "input:  " << to_latex(term) << "\n";
        auto result = spintrace(term, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
        expand(result);
        rapid_simplify(result);
        canonicalize(result);
        std::wcout << "result: " << to_latex(result) << "\n\n";
        break;
      }
//      auto result = spintrace(input, {{L"i_1", L"a_1"}, {L"i_2", L"a_2"}});
//      expand(result);
//      rapid_simplify(result);
//      canonicalize(result);
//      std::wcout << "result: " << to_latex(result) << "\n";
    }
#if 0
    {
      const auto input = ex<Constant>(1.0/2.0) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"}, Symmetry::antisymm);
      std::wcout << "input:  " << to_latex(input) << "\n";
      auto result = ex<Constant>(1) * spintrace(input);
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      std::wcout << "result: " << to_latex(result) << "\n";
      std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                       {Index{L"i_2"}, Index{L"i_1"}}};
      auto transformed_result = transform_expression(result, idxmap, 1.0/3.0);
      std::wcout << "transformed_result: " << to_latex(transformed_result) << "\n";
    }
#endif
    /*
        {
          const auto input = ex<Constant>(1/8) *
              ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"}, WstrList{L"a_3",
    L"a_4"}, Symmetry::antisymm) * ex<Tensor>(L"t", WstrList{L"a_3",
    L"a_4"},WstrList{L"i_1", L"i_2"}, Symmetry::antisymm); // -
    //          ex<Constant>(1/8) *
    //          ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"}, WstrList{L"a_3",
    L"a_4"}, Symmetry::antisymm) *
    //          ex<Tensor>(L"t", WstrList{L"a_3", L"a_4"},WstrList{L"i_2",
    L"i_1"}, Symmetry::antisymm); auto result = ex<Constant>(0.25) *
    spintrace(input); expand(result); rapid_simplify(result);
          canonicalize(result);
          std::wcout << "1/8*g*t2: " << to_latex(result) << "\n\n";
        }

        {
          const auto input = ex<Constant>(0.25) *
              ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1",
    L"a_2"}, Symmetry::antisymm) * ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"},
    WstrList{L"a_3", L"a_4"}, Symmetry::antisymm) * ex<Tensor>(L"t",
    WstrList{L"a_3"}, WstrList{L"i_1"}) * ex<Tensor>(L"t", WstrList{L"a_4"},
    WstrList{L"i_2"}) * ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"});
          std::wcout << "input:  " << to_latex(input) << "\n";
          auto result = expand_A_operator(input);
          std::wcout << "result: " << to_latex(result) << "\n\n";
        }
    */

#if CCSD_r1
    // These terms from CCSD R1 equations
    {
      // A * F
      const auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                         ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"i_1"});
      //  std::wcout << "input:1  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) == L"{ \\left({{f^{{i_1}}_{{a_1}}}}\\right) }");
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
      // std::wcout << "input:2  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      // std::wcout << "result: " << to_latex(result) << "\n";
      REQUIRE(
          to_latex(result) ==
          L"{ \\left({{{2}}{g^{{i_1}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}} - "
          L"{{g^{{a_2}{i_1}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}}\\right) }");

      std::map<Index, Index> idxmap = {{Index{L"i_1"}, Index{L"i_2"}},
                                       {Index{L"i_2"}, Index{L"i_1"}}};
      auto transformed_result = transform_expression(result, idxmap);
      REQUIRE(transformed_result->is<Sum>());
      REQUIRE(transformed_result->size() == 2);
      REQUIRE(
          to_latex(transformed_result) ==
          L"{ \\left({{{2}}{g^{{i_2}{a_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}} - "
          L"{{g^{{a_2}{i_2}}_{{a_1}{i_1}}}{t^{{i_1}}_{{a_2}}}}\\right) }");
      // std::wcout << "transformed_result: " << to_latex(transformed_result) <<
      // "\n\n";
    }

    {
      // - A * f * t1
      const auto input = ex<Constant>(-1) *
                         ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                         ex<Tensor>(L"f", WstrList{L"i_2"}, WstrList{L"i_1"}) *
                         ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_2"});
      //  std::wcout << "input:3  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(
          to_latex(result) ==
          L"{ \\left( - {{f^{{i_1}}_{{i_2}}}{t^{{i_2}}_{{a_1}}}}\\right) }");
    }

    {
      // A * f * t1
      const auto input = ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                         ex<Tensor>(L"f", WstrList{L"a_1"}, WstrList{L"a_2"}) *
                         ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_1"});
      //  std::wcout << "input:4  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
              L"{ \\left({{f^{{a_2}}_{{a_1}}}{t^{{i_1}}_{{a_2}}}}\\right) }");
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
      //  std::wcout << "input:5  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
              L"{ \\left( - "
              L"{{{2}}{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_3}{i_2}}_{{a_1}{a_2}}"
              L"}} + "
              L"{{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_2}{i_3}}_{{a_1}{a_2}}}}"
              L"\\right) }");
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
              L"\\left({{{2}}{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_"
              L"2}{a_3}}}} - "
              L"{{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_2}{a_3}}}}"
              L"\\right) }");
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
          to_latex(result) ==
          L"{ \\left( - {{f^{{a_2}}_{{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}} + "
          L"{{{2}}{f^{{a_2}}_{{i_2}}}{t^{{i_1}{i_2}}_{{a_1}{a_2}}}}\\right) }");
    }

    {
      // A * g * t1 * t1
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1"}, WstrList{L"a_1"}) *
          ex<Tensor>(L"g", WstrList{L"i_2", L"a_1"}, WstrList{L"a_2", L"a_3"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"});
      //  std::wcout << "input:8  " << to_latex(input) << "\n";
      auto result = ex<Constant>(0.5) * spintrace(input, {{L"i_1", L"a_1"}});
      expand(result);
      rapid_simplify(result);
      canonicalize(result);
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{2}}{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_2}}}{"
              L"t^{{i_1}}_{{a_3}}}} - "
              L"{{g^{{a_3}{a_2}}_{{a_1}{i_2}}}{t^{{i_2}}_{{a_3}}}{t^{{i_1}}_{{"
              L"a_2}}}}\\right) }");
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
              L"{ \\left( - "
              L"{{{2}}{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}"
              L"}_{{a_2}}}} + "
              L"{{g^{{a_2}{i_1}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_3}}_{{"
              L"a_2}}}}\\right) }");
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
              L"{ \\left( - "
              L"{{f^{{a_2}}_{{i_2}}}{t^{{i_2}}_{{a_1}}}{t^{{i_1}}_{{a_2}}}}"
              L"\\right) }");
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
      REQUIRE(to_latex(result) ==
              L"{ \\left( - "
              L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}"
              L"{i_1}}_{{a_2}{a_3}}}} + "
              L"{{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_3}{i_1}"
              L"}_{{a_2}{a_3}}}}\\right) }");
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
      //  std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(to_latex(result) ==
              L"{ \\left( - "
              L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_1}}_{{a_3}}}{t^{{i_3}"
              L"{i_2}}_{{a_1}{a_2}}}} + "
              L"{{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_1}}_{{a_3}}}{t^{{i_2}{i_3}"
              L"}_{{a_1}{a_2}}}}\\right) }");
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
              L"{ "
              L"\\left({{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_2}}}{t^{{i_"
              L"2}{i_1}}_{{a_1}{a_3}}}} + "
              L"{{{4}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_2}}}{t^{{i_1}"
              L"{i_3}}_{{a_1}{a_3}}}} - "
              L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_2}}}{t^{{i_3}"
              L"{i_1}}_{{a_1}{a_3}}}} - "
              L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_3}}}{t^{{i_1}"
              L"{i_3}}_{{a_1}{a_2}}}}\\right) }");
    }

    {
      // - A * g * t1 * t1 * t1
      auto input =
          ex<Constant>(-1.) *
          ex<Tensor>(L"g", WstrList{L"i_2", L"i_3"}, WstrList{L"a_2", L"a_3"},
                     Symmetry::antisymm) *
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
              L"\\left({{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_2}}_{{a_1}}}{t^{{i_"
              L"3}}_{{a_2}}}{t^{{i_1}}_{{a_3}}}} - "
              L"{{{2}}{g^{{a_2}{a_3}}_{{i_2}{i_3}}}{t^{{i_3}}_{{a_1}}}{t^{{i_2}"
              L"}_{{a_2}}}{t^{{i_1}}_{{a_3}}}}\\right) }");
    }
#endif
  }

#if MANDATORY_TESTS
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

    // 2-body
    input = ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                       Symmetry::antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(result->is<Sum>());
    REQUIRE(to_latex(result) ==
            L"{ \\left({{g^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
            L"{{g^{{a_1}{a_2}}_{{i_2}{i_1}}}}\\right) }");

    // 3-body
    input = ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                       WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
    result = expand_antisymm(input->as<Tensor>());
    REQUIRE(result->is<Sum>());
    REQUIRE(to_latex(result) ==
            L"{ \\left({{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_3}{a_2}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_2}{a_1}{a_3}}}} + "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_2}{a_3}{a_1}}}} + "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_3}{a_1}{a_2}}}} - "
            L"{{t^{{i_1}{i_2}{i_3}}_{{a_3}{a_2}{a_1}}}}\\right) }");
  }

  SECTION("Anti-symmetrizer Operator") {
    // 2-body
    {
      const auto input =
          ex<Constant>(1. / 4) * ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                                            WstrList{L"a_1", L"a_2"},
                                            Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(to_latex(result) ==
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}}");
    }

    // 2-body
    {
      const auto input =
          ex<Constant>(1. / 4) *
          ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(
          to_latex(result) ==
          L"{ \\left({{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_1}{i_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_2}{a_1}}_{{i_1}{i_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_1}{a_2}}_{{i_2}{i_1}}}} + "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_2}{a_1}}_{{i_2}{i_1}}}}\\right) }");
    }

    {
      // 1/4 * A * g * t1 * t1
      const auto input =
          ex<Constant>(1. / 4.) *
          ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"a_1", L"a_2"}, WstrList{L"a_3", L"a_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_2"});
      // std::wcout << "input:  " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      // std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_1}{a_2}}}{t^"
              L"{{i_1}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}} - "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_2}{a_1}}}{t^{{i_1}}"
              L"_{{a_3}}}{t^{{i_2}}_{{a_4}}}} - "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_1}{a_2}}}{t^{{i_2}}"
              L"_{{a_3}}}{t^{{i_1}}_{{a_4}}}} + "
              L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{a_2}{a_1}}}{t^{{i_2}}"
              L"_{{a_3}}}{t^{{i_1}}_{{a_4}}}}\\right) }");
    }

    {
      // 1/4 * A * g * t1 * t1 * t1 * t1
      const auto input =
          ex<Constant>(1. / 4.) *
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_3", L"i_4"}, WstrList{L"a_3", L"a_4"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_3"}, WstrList{L"i_1"}) *
          ex<Tensor>(L"t", WstrList{L"a_4"}, WstrList{L"i_2"}) *
          ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_3"}) *
          ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_4"});
      // std::wcout << "input:  " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      // std::wcout << "result: " << to_latex(result) << "\n\n";
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 4);
      REQUIRE(
          to_latex(result) ==
          L"{ "
          L"\\left({{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_"
          L"1}}_{{a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_3}}_{{a_1}}}{t^{{i_4}}_{{a_2}"
          L"}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_2}}_{{"
          L"a_3}}}{t^{{i_1}}_{{a_4}}}{t^{{i_3}}_{{a_1}}}{t^{{i_4}}_{{a_2}}}} - "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_1}}_{{"
          L"a_3}}}{t^{{i_2}}_{{a_4}}}{t^{{i_3}}_{{a_2}}}{t^{{i_4}}_{{a_1}}}} + "
          L"{{{\\frac{1}{4}}}{\\bar{g}^{{a_3}{a_4}}_{{i_3}{i_4}}}{t^{{i_2}}_{{"
          L"a_3}}}{t^{{i_1}}_{{a_4}}}{t^{{i_3}}_{{a_2}}}{t^{{i_4}}_{{a_1}}}}"
          L"\\right) }");
    }

    // 3-body
    {
      const auto input =
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(to_latex(result) == L"{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}");

      auto antisymm_t = expand_antisymm(input->as<Tensor>());
      //  std::wcout << "3-body result: " << to_latex(antisymm_t) << "\n";
    }

    // 3-body
    {
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      //  std::wcout << "3-body input: " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      //  std::wcout << "3-body result: " << to_latex(result) << "\n";
    }

    {
      const auto input =
          ex<Constant>(1. / 4) +
          ex<Constant>(1. / 2) *
              ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                         WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"}) *
              ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) +
          ex<Constant>(1. / 4) *
              ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"},
                         WstrList{L"a_1", L"a_2"}, Symmetry::antisymm) *
              ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"},
                         WstrList{L"i_1", L"i_2"}, Symmetry::antisymm);
      //  std::wcout << "\ninput: " << to_latex(input) << "\n";
      auto A_found = check_A_operator(input);
      if (A_found) {
        auto result = expand_A_operator(input);
        //  std::wcout << "result: " << to_latex(result) << "\n";
        auto iter = 0;
        for (auto &&summand : *result) {
          ++iter;
          //  std::wcout << "term: " << to_latex(summand) << "\n";
          auto spin_traced = spintrace(summand);
          //  std::wcout << "sptr: " << to_latex(spin_traced) << "\n";
          canonicalize(spin_traced);
          //  std::wcout << "can: " << to_latex(spin_traced) << "\n\n";
          if (iter > 0) break;
        }
        // auto spin_traced = spintrace(result);
        // //  std::wcout << "spin_traced: " << to_latex(spin_traced) << "\n";
      } else {
        // std::cout << "No A found.\n";
      }
    }

    {
      const auto input =
          ex<Constant>(1. / 4) *
          ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm);

      //  std::wcout << "\ninput: " << to_latex(input) << "\n";
      auto A_found = check_A_operator(input);
      if (A_found) {
        auto result = expand_A_operator(input);
        //  std::wcout << "result: " << to_latex(result) << "\n";
        auto iter = 0;
        for (auto &&summand : *result) {
          ++iter;
          //  std::wcout << "term: " << to_latex(summand) << "\n";
          auto spin_traced = spintrace(summand);
          //  std::wcout << "sptr: " << to_latex(spin_traced) << "\n";
          canonicalize(spin_traced);
          //  std::wcout << "can: " << to_latex(spin_traced) << "\n\n";
          if (iter > 0) break;
        }
        // auto spin_traced = spintrace(result);
        // //  std::wcout << "spin_traced: " << to_latex(spin_traced) << "\n";
      }
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
                L"{ \\left({{g^{{p_3}{p_4}}_{{p_1}{p_2}}}} - "
                L"{{g^{{p_3}{p_4}}_{{p_2}{p_1}}}}\\right) }");
  }

  SECTION("Antisymmetrizer check") {
    const auto input = ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                  WstrList{L"p_3", L"p_4"}, Symmetry::antisymm);

    auto result = expand_antisymm(input->as<Tensor>());
    //  std::wcout << "input:  " << to_latex(input) << "\n";
    //  std::wcout << "result: " << to_latex(result) << "\n";
    canonicalize(result);
    //  std::wcout << "canoni: " << to_latex(result) << "\n\n";
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
      auto result = spintrace(expr, {{L"i_1", L"a_1"}});
      canonicalize(result);
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(to_latex(result) ==
              L"{ \\left( - "
              L"{{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_2}}}{t^{{i_2}}_{{"
              L"a_1}}}} + "
              L"{{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}}_{{a_1}}}{t^{{i_2}"
              L"}_{{a_2}}}}\\right) }");
    }
#if 0
    {
      // 1/4 * g * t2
      const auto expr =
          ex<Constant>(1. / 4.) *
          ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"i_3", L"a_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm);
      std::wcout << __LINE__ << "L " << to_latex(expr) << "\n";
      auto result = spintrace(expr);  //{L"i_1", L"a_1"}
      canonicalize(result);
      std::wcout << __LINE__ << "L " << to_latex(result) << "\n";
      REQUIRE(result->is<Sum>());
      REQUIRE(result->size() == 2);
      REQUIRE(to_latex(result) ==
              L"{ "
              L"\\left({{{2}}{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_1}{i_2}}_{{a_"
              L"1}{a_2}}}} - "
              L"{{g^{{a_1}{a_2}}_{{i_1}{i_2}}}{t^{{i_2}{i_1}}_{{a_1}{a_2}}}}"
              L"\\right) }");
    }
#endif
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
#endif
}  // TEST_CASE("Spin Trace")
