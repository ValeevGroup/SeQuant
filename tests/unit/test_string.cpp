//
// Created by Eduard Valeyev on 7/19/23.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/wstring.hpp>

#include <memory>
#include <type_traits>

TEST_CASE("string", "[util]") {
  using namespace sequant;

  SECTION("SQ_STRLIT") {
    STATIC_REQUIRE(
        std::is_same_v<decltype(SQ_STRLIT(char, "alpha")), const char(&)[6]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(wchar_t, "alpha")),
                                  const wchar_t(&)[6]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(char8_t, "alpha")),
                                  const char8_t(&)[6]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(char16_t, "alpha")),
                                  const char16_t(&)[6]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(char32_t, "alpha")),
                                  const char32_t(&)[6]>);

    STATIC_REQUIRE(
        std::is_same_v<decltype(SQ_STRLIT(char, "ɑ")), const char(&)[3]>);
    STATIC_REQUIRE(
        std::is_same_v<decltype(SQ_STRLIT(wchar_t, "ɑ")), const wchar_t(&)[2]>);
    STATIC_REQUIRE(
        std::is_same_v<decltype(SQ_STRLIT(char8_t, "ɑ")), const char8_t(&)[3]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(char16_t, "ɑ")),
                                  const char16_t(&)[2]>);
    STATIC_REQUIRE(std::is_same_v<decltype(SQ_STRLIT(char32_t, "ɑ")),
                                  const char32_t(&)[2]>);
  }
}
