//
// Created by Eduard Valeyev on 7/19/23.
//

#include "catch.hpp"

#include "SeQuant/core/latex.hpp"

TEST_CASE("latex", "[util]") {
  using namespace sequant;

  SECTION("greek->latex") {
    using namespace std::string_literals;
    REQUIRE(greek_characters_to_latex("alpha"s) == "alpha");
    REQUIRE_THROWS_AS(greek_characters_to_latex("α"s) == "\\alpha",
                      std::invalid_argument);

    REQUIRE(greek_characters_to_latex(std::wstring(L"alpha")) == L"alpha");
    REQUIRE(greek_characters_to_latex(std::wstring(L"α")) == L"\\alpha");

    REQUIRE(greek_characters_to_latex(std::wstring(L"Alpha")) == L"Alpha");
    REQUIRE(greek_characters_to_latex(std::wstring(L"Α")) == L"A");
    REQUIRE(greek_characters_to_latex(std::wstring(L"Γ")) == L"\\Gamma");

#if __cplusplus >= 202002L
    REQUIRE(greek_characters_to_latex(std::u8string(u8"alpha")) == u8"alpha");
    REQUIRE_THROWS_AS(
        greek_characters_to_latex(std::u8string(u8"α")) == u8"\\alpha",
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        greek_characters_to_latex(std::u8string(u8"Γ")) == u8"\\Gamma",
        std::invalid_argument);

    REQUIRE(greek_characters_to_latex(std::u16string(u"alpha")) == u"alpha");
    REQUIRE(greek_characters_to_latex(std::u16string(u"α")) == u"\\alpha");
    REQUIRE(greek_characters_to_latex(std::u16string(u"Γ")) == u"\\Gamma");

    REQUIRE(greek_characters_to_latex(std::u32string(U"alpha")) == U"alpha");
    REQUIRE(greek_characters_to_latex(std::u32string(U"α")) == U"\\alpha");
    REQUIRE(greek_characters_to_latex(std::u32string(U"Γ")) == U"\\Gamma");

#endif
  }

  SECTION("\\tilde") {
    std::wstring tilde__a = L"ã";
    REQUIRE(tilde__a.size() == 2);
    REQUIRE(diactrics_to_latex(tilde__a) == L"\\tilde{a}");
    std::wstring tilde_a = L"ã";
    REQUIRE(tilde_a.size() == 1);
    REQUIRE(diactrics_to_latex(tilde_a) == L"\\tilde{a}");
    std::wstring tilde_A = L"Ã";
    REQUIRE(tilde_A.size() == 1);
    REQUIRE(diactrics_to_latex(tilde_A) == L"\\tilde{A}");

    std::wstring tilde__f = L"f̃";
    REQUIRE(tilde__f.size() == 2);
    REQUIRE(diactrics_to_latex(tilde__f) == L"\\tilde{f}");
  }
}
