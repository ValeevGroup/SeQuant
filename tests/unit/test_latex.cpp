//
// Created by Eduard Valeyev on 7/19/23.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/io/latex/latex.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <stdexcept>
#include <string>

TEST_CASE("latex", "[util]") {
  using namespace sequant;

  SECTION("greek->latex") {
    REQUIRE(io::latex::greek_characters_to_string("") == "");
    REQUIRE(io::latex::greek_characters_to_string("alpha") == "alpha");
    REQUIRE(io::latex::greek_characters_to_string("α") == "\\alpha");
    REQUIRE(io::latex::greek_characters_to_string("Alpha") == "Alpha");
    REQUIRE(io::latex::greek_characters_to_string("Α") == "A");
    REQUIRE(io::latex::greek_characters_to_string("Γ") == "\\Gamma");
    REQUIRE(io::latex::greek_characters_to_string("aγb") == "a\\gammab");
    // other non-ASCII characters (incl. ones outside of the BMP) are kept
    REQUIRE(io::latex::greek_characters_to_string("α↑𝛼ω") ==
            "\\alpha↑𝛼\\omega");
  }

  SECTION("diactrids->latex") {
    std::string tilde__a = "ã";
    REQUIRE(to_utf32(tilde__a).size() == 2);
    REQUIRE(io::latex::diactrics_to_string(tilde__a) == "\\tilde{a}");
    std::string tilde_a = "ã";
    REQUIRE(to_utf32(tilde_a).size() == 1);
    REQUIRE(io::latex::diactrics_to_string(tilde_a) == "\\tilde{a}");
    std::string tilde_A = "Ã";
    REQUIRE(to_utf32(tilde_A).size() == 1);
    REQUIRE(io::latex::diactrics_to_string(tilde_A) == "\\tilde{A}");

    std::string tilde__f = "f̃";
    REQUIRE(to_utf32(tilde__f).size() == 2);
    REQUIRE(io::latex::diactrics_to_string(tilde__f) == "\\tilde{f}");
    std::string grave__f = "f̀";
    REQUIRE(to_utf32(grave__f).size() == 2);
    REQUIRE(io::latex::diactrics_to_string(grave__f) == "\\grave{f}");
    std::string acute__f = "f́";
    REQUIRE(to_utf32(acute__f).size() == 2);
    REQUIRE(io::latex::diactrics_to_string(acute__f) == "\\acute{f}");
    std::string caron__f = "f̌";
    REQUIRE(to_utf32(caron__f).size() == 2);
    REQUIRE(io::latex::diactrics_to_string(caron__f) == "\\check{f}");
    // characters after a replacement are kept, whether replaced or not
    std::string tilde__f_g = "f̃g";
    REQUIRE(io::latex::diactrics_to_string(tilde__f_g) == "\\tilde{f}g");
    std::string tilde_a_arrow = "ã→";
    REQUIRE(io::latex::diactrics_to_string(tilde_a_arrow) == "\\tilde{a}→");
    // a non-ASCII character that follows a combined pair and ends the string
    // (a spin-cased label) is not taken to carry the pair's diacritic
    std::string tilde__mu_up = "μ̃↑";
    REQUIRE(to_utf32(tilde__mu_up).size() == 3);
    REQUIRE(io::latex::diactrics_to_string(tilde__mu_up) == "\\tilde{μ}↑");
    REQUIRE(io::latex::utf_to_string(tilde__mu_up) == "\\tilde{\\mu}↑");

    // notice the size is 1 since these are precomposed characters (i.e., single
    // code points)
    std::string hat_A = "Â";
    REQUIRE(to_utf32(hat_A).size() == 1);
    REQUIRE(io::latex::diactrics_to_string(hat_A) == "\\hat{A}");
    std::string hat_S = "Ŝ";
    REQUIRE(to_utf32(hat_S).size() == 1);
    REQUIRE(io::latex::diactrics_to_string(hat_S) == "\\hat{S}");
  }

  SECTION("UTF->latex") {
    std::string tilde__alpha = "α̃";
    REQUIRE(to_utf32(tilde__alpha).size() == 2);
    REQUIRE(io::latex::utf_to_string(tilde__alpha) == "\\tilde{\\alpha}");
    std::string hat__alpha = "α̂";
    REQUIRE(to_utf32(hat__alpha).size() == 2);
    REQUIRE(io::latex::utf_to_string(hat__alpha) == "\\hat{\\alpha}");
  }
}
