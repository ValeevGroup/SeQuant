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
}
