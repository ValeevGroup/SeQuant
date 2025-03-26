//
// Created by Eduard Valeyev on 7/19/23.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/meta.hpp>

#define STATIC_ASSERT_AND_CHECK(X) \
  static_assert(X);                \
  CHECK(X);

TEST_CASE("meta", "[meta]") {
  using namespace sequant;
  using namespace sequant::meta;

  SECTION("castable_to_any") {
    STATIC_ASSERT_AND_CHECK((is_statically_castable_v<castable_to_any, Index>));
    STATIC_ASSERT_AND_CHECK(
        (is_statically_castable_v<range_value_t<std::array<castable_to_any, 0>>,
                                  Index>));
  }

  SECTION("std_array_size") {
    STATIC_ASSERT_AND_CHECK((std_array_size_v<std::array<int, 3>> == 3));
    STATIC_ASSERT_AND_CHECK(
        (std_array_size_v<decltype(std::array{L"a_1", L"a_12", L"a_123"})> ==
         3));
  }

  SECTION("literal_to_string") {
    STATIC_ASSERT_AND_CHECK(
        (std::is_same_v<literal_to_string_t<wchar_t[4]>, std::wstring>));
    STATIC_ASSERT_AND_CHECK(
        (std::is_same_v<literal_to_string_t<wchar_t[]>, std::wstring>));
    STATIC_ASSERT_AND_CHECK(
        (std::is_same_v<literal_to_string_t<const wchar_t[]>, std::wstring>));
    STATIC_ASSERT_AND_CHECK(
        (std::is_same_v<literal_to_string_t<const wchar_t*>, std::wstring>));
  }
}
