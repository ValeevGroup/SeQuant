//
// Created by Eduard Valeyev on 7/19/23.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/meta.hpp>

#include <stdexcept>
#include <string>

TEST_CASE("latex", "[util]") {
  using namespace sequant;

  SECTION("greek->latex") {
    using namespace std::string_literals;
    REQUIRE(io::latex::greek_characters_to_string("alpha"s) == "alpha");
    REQUIRE_THROWS_AS(io::latex::greek_characters_to_string("α"s) == "\\alpha",
                      Exception);

    REQUIRE(io::latex::greek_characters_to_string(std::wstring(L"alpha")) ==
            L"alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::wstring(L"α")) ==
            L"\\alpha");

    REQUIRE(io::latex::greek_characters_to_string(std::wstring(L"Alpha")) ==
            L"Alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::wstring(L"Α")) == L"A");
    REQUIRE(io::latex::greek_characters_to_string(std::wstring(L"Γ")) ==
            L"\\Gamma");

    REQUIRE(io::latex::greek_characters_to_string(std::u8string(u8"alpha")) ==
            u8"alpha");
    REQUIRE_THROWS_AS(io::latex::greek_characters_to_string(
                          std::u8string(u8"α")) == u8"\\alpha",
                      Exception);
    REQUIRE_THROWS_AS(io::latex::greek_characters_to_string(
                          std::u8string(u8"Γ")) == u8"\\Gamma",
                      Exception);

    REQUIRE(io::latex::greek_characters_to_string(std::u16string(u"alpha")) ==
            u"alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::u16string(u"α")) ==
            u"\\alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::u16string(u"Γ")) ==
            u"\\Gamma");

    REQUIRE(io::latex::greek_characters_to_string(std::u32string(U"alpha")) ==
            U"alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::u32string(U"α")) ==
            U"\\alpha");
    REQUIRE(io::latex::greek_characters_to_string(std::u32string(U"Γ")) ==
            U"\\Gamma");
  }

  SECTION("diactrids->latex") {
    std::wstring tilde__a = L"ã";
    REQUIRE(tilde__a.size() == 2);
    REQUIRE(io::latex::diactrics_to_string(tilde__a) == L"\\tilde{a}");
    std::wstring tilde_a = L"ã";
    REQUIRE(tilde_a.size() == 1);
    REQUIRE(io::latex::diactrics_to_string(tilde_a) == L"\\tilde{a}");
    std::wstring tilde_A = L"Ã";
    REQUIRE(tilde_A.size() == 1);
    REQUIRE(io::latex::diactrics_to_string(tilde_A) == L"\\tilde{A}");

    std::wstring tilde__f = L"f̃";
    REQUIRE(tilde__f.size() == 2);
    REQUIRE(io::latex::diactrics_to_string(tilde__f) == L"\\tilde{f}");
    std::wstring grave__f = L"f̀";
    REQUIRE(grave__f.size() == 2);
    REQUIRE(io::latex::diactrics_to_string(grave__f) == L"\\grave{f}");
    std::wstring acute__f = L"f́";
    REQUIRE(acute__f.size() == 2);
    REQUIRE(io::latex::diactrics_to_string(acute__f) == L"\\acute{f}");
    std::wstring caron__f = L"f̌";
    REQUIRE(caron__f.size() == 2);
    REQUIRE(io::latex::diactrics_to_string(caron__f) == L"\\check{f}");

    // notice the size is 1 since these are precomposed characters
    std::wstring hat_A = L"Â";
    REQUIRE(hat_A.size() == 1);
    REQUIRE(io::latex::diactrics_to_string(hat_A) == L"\\hat{A}");
    std::wstring hat_S = L"Ŝ";
    REQUIRE(hat_S.size() == 1);
    REQUIRE(io::latex::diactrics_to_string(hat_S) == L"\\hat{S}");
  }

  SECTION("UTF->latex") {
    std::wstring tilde__alpha = L"α̃";
    REQUIRE(tilde__alpha.size() == 2);
    REQUIRE(io::latex::utf_to_string(tilde__alpha) == L"\\tilde{\\alpha}");
    std::wstring hat__alpha = L"α̂";
    REQUIRE(hat__alpha.size() == 2);
    REQUIRE(io::latex::utf_to_string(hat__alpha) == L"\\hat{\\alpha}");
  }
}
