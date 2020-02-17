//
// Created by Nakul Teke on 12/20/19.
//

#include "SeQuant/domain/mbpt/spin.hpp"
#include "catch.hpp"

TEST_CASE("Spin") {
  using namespace sequant;

  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

SECTION("Tensor spintrace"){

}

# if 0
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
    auto p1 =
        Index(L"p⁺_1", IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
    auto p2 =
        Index(L"p⁻_2", IndexSpace::instance(IndexSpace::all, IndexSpace::beta));
    auto p3 =
        Index(L"p⁺_3", IndexSpace::instance(IndexSpace::all, IndexSpace::alpha));
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

    // 3-body
    {
      const auto input =
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      auto result = expand_A_operator(input);
      REQUIRE(to_latex(result) == L"{t^{{i_1}{i_2}{i_3}}_{{a_1}{a_2}{a_3}}}");

      auto antisymm_t = expand_antisymm(input->as<Tensor>());
      std::wcout << "3-body result: " << to_latex(antisymm_t) << "\n";
    }

    // 3-body
    {
      const auto input =
          ex<Tensor>(L"A", WstrList{L"i_1", L"i_2", L"i_3"},
                     WstrList{L"a_1", L"a_2", L"a_3"}, Symmetry::antisymm) *
          ex<Tensor>(L"t", WstrList{L"a_1", L"a_2", L"a_3"},
                     WstrList{L"i_1", L"i_2", L"i_3"}, Symmetry::antisymm);
      std::wcout << "3-body input: " << to_latex(input) << "\n";
      auto result = expand_A_operator(input);
      std::wcout << "3-body result: " << to_latex(result) << "\n";
    }
#endif
#if 0
    {
      const auto input =
          ex<Constant>(1. / 4) +
              ex<Constant>(1. / 2) *
                  ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                             Symmetry::antisymm) *
                  ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"}) *
                  ex<Tensor>(L"t", WstrList{L"a_2"}, WstrList{L"i_2"}) +
              ex<Constant>(1. / 4) *
                  ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                             Symmetry::antisymm) *
                  ex<Tensor>(L"t", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                             Symmetry::antisymm);
      std::wcout << "\ninput: " << to_latex(input) << "\n";
      auto A_found = check_A_operator(input);
      if (A_found) {
        auto result = expand_A_operator(input);
        std::wcout << "result: " << to_latex(result) << "\n";
        auto iter = 0;
        for (auto &&summand: *result) {
          ++iter;
          std::wcout << "term: " << to_latex(summand) << "\n";
          auto spin_traced = spintrace(summand);
          std::wcout << "sptr: " << to_latex(spin_traced) << "\n";
          canonicalize(spin_traced);
          std::wcout << "can: " << to_latex(spin_traced) << "\n\n";
          if(iter>0) break;
        }
        // auto spin_traced = spintrace(result);
        // std::wcout << "spin_traced: " << to_latex(spin_traced) << "\n";
      } else {
        std::cout << "No A found.\n";
      }
    }

    {
    const auto input = ex<Constant>(1. / 4) *
        ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                   Symmetry::antisymm) *
        ex<Tensor>(L"g", WstrList{L"i_1", L"i_2"}, WstrList{L"a_1", L"a_2"},
                   Symmetry::antisymm);

    std::wcout << "\ninput: " << to_latex(input) << "\n";
    auto A_found = check_A_operator(input);
    if (A_found) {
      auto result = expand_A_operator(input);
      std::wcout << "result: " << to_latex(result) << "\n";
      auto iter = 0;
      for (auto &&summand: *result) {
        ++iter;
        std::wcout << "term: " << to_latex(summand) << "\n";
        auto spin_traced = spintrace(summand);
        std::wcout << "sptr: " << to_latex(spin_traced) << "\n";
        canonicalize(spin_traced);
        std::wcout << "can: " << to_latex(spin_traced) << "\n\n";
        if(iter>0) break;
      }
      // auto spin_traced = spintrace(result);
      // std::wcout << "spin_traced: " << to_latex(spin_traced) << "\n";
    }

  }

    {
      const auto input = ex<Constant>(-0.5) *
          ex<Tensor>(L"A", WstrList{L"a_1", L"a_2"}, WstrList{L"i_1", L"i_2"},
                     Symmetry::antisymm) *
          ex<Tensor>(L"g", WstrList{L"a_3", L"a_4"}, WstrList{L"a_1", L"a_2"},
                     Symmetry::antisymm) *
      ex<Tensor>(L"t", WstrList{L"i_2", L"i_1"}, WstrList{L"a_3", L"a_4"},
                 Symmetry::antisymm);

      std::wcout << "\ninput: " << to_latex(input) << "\n";
      auto A_found = check_A_operator(input);
      if(A_found){
        auto result = expand_A_operator(input);
        std::wcout << "result: " << to_latex(result) << "\n";
        auto iter = 0;
        for(auto&& summand: *result){
          ++iter;
          std::wcout<< "term: " << to_latex(summand) << "\n";
          auto spin_traced = spintrace(summand);
          std::wcout<< "sptr: " << to_latex(spin_traced) << "\n";
          canonicalize(spin_traced);
          std::wcout << "can: " << to_latex(spin_traced) << "\n\n";
          if(iter>0) break;
        }
        // auto spin_traced = spintrace(result);
        // std::wcout << "spin_traced: " << to_latex(spin_traced) << "\n";
      }

    }

#endif
  }

#if 0

  SECTION("Constant") {
    auto exprPtr = ex<Constant>(1. / 4);
    auto result = spintrace(exprPtr);
    REQUIRE(result->is<Constant>());
    REQUIRE(result->is_atom());
    REQUIRE(to_latex(result) == L"{{{\\frac{1}{4}}}}");
  }

  SECTION("Tensor") {
    const auto expr = ex<Constant>(0.25) * ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                                 WstrList{L"p_3", L"p_4"}, Symmetry::antisymm);
    auto result = spintrace(expr);
    REQUIRE(result->is<Sum>());
    std::wcout << "input:  " << to_latex(expr) << "\n";
    std::wcout << "result: " << to_latex(result) << "\n";
    canonicalize(result);
    std::wcout << "canoni: " << to_latex(result) << "\n\n";
//  REQUIRE(result->size() == 2);
//    REQUIRE(to_latex(result) ==
//            L"{ \\left({{g^{{p_3}{p_4}}_{{p_1}{p_2}}}} - "
//            L"{{g^{{p_3}{p_4}}_{{p_2}{p_1}}}}\\right) }");
  }

  SECTION("Antisymmetrizer check"){
    const auto input = ex<Tensor>(L"g", WstrList{L"p_1", L"p_2"},
                         WstrList{L"p_3", L"p_4"}, Symmetry::antisymm);

    auto result = expand_antisymm(input->as<Tensor>());
    std::wcout << "input:  " << to_latex(input) << "\n";
    std::wcout << "result: " << to_latex(result) << "\n";
    canonicalize(result);
    std::wcout << "canoni: " << to_latex(result) << "\n\n";
  }



  SECTION("Product") {
    const auto expr = ex<Tensor>(L"f", WstrList{L"i_1"}, WstrList{L"a_1"}) *
                      ex<Tensor>(L"t", WstrList{L"a_1"}, WstrList{L"i_1"});
    auto result = spintrace(expr);
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
      auto result = spintrace(expr);
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
