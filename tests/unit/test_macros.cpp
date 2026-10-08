//
// Created by Eduard Valeyev on 10/16/25.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/utility/macros.hpp>

#include <filesystem>

TEST_CASE("macros", "[elements]") {
  SECTION("SEQUANT_ASSERT") {
    std::filesystem::path this_file("tests/unit/test_macros.cpp");

    const std::string file_path = this_file.string();

    this_file.make_preferred();

    const std::string file_path_native = this_file.string();

    CAPTURE(file_path);
    CAPTURE(file_path_native);

    if (sequant::assert_behavior() == sequant::AssertBehavior::Throw) {
      try {
        // clang-format off
#line 1000  // to make sure the line number of the next line is fixed
        sequant::assert_failed("test");
        // clang-format on
        FAIL("Assert should have thrown");
      } catch (sequant::Exception& ex) {
        CAPTURE(ex.what());
        // see #line up there
        // N.B. clang <16 has std::source_location produce wrong line numbers
        // when initialized as default argument see
        // https://github.com/llvm/llvm-project/issues/56379
#if !defined(SEQUANT_CXX_COMPILER_IS_CLANG) || __clang_major__ >= 16
        bool found = std::string_view(ex.what()).find(
                         file_path + ":1000 in function ") != std::string::npos;
        found |=
            std::string_view(ex.what()).find(
                file_path_native + ":1000 in function ") != std::string::npos;
        REQUIRE(found);
#endif
      }
      try {
        // clang-format off
#line 2000  // to make sure the line number of the next line is fixed
        SEQUANT_ASSERT(1 == 0 && "1 != 0", "testing SEQUANT_ASSERT");
        // clang-format on
        FAIL("Assert should have thrown");
      } catch (sequant::Exception& ex) {
        CAPTURE(ex.what());
        // see #line up there
        // N.B. clang <16 has std::source_location produce wrong line numbers
        // when initialized as default argument see
        // https://github.com/llvm/llvm-project/issues/56379
#if defined(SEQUANT_CXX_COMPILER_IS_CLANG) && __clang_major__ >= 16
        bool found = std::string_view(ex.what()).find(
                         file_path + ":2000 in function ") != std::string::npos;
        found |=
            std::string_view(ex.what()).find(
                file_path_native + ":2000 in function ") != std::string::npos;
        REQUIRE(found);
#endif
      }
    }
  }

  SECTION("SEQUANT_ENFORCE") {
    int evaluations = 0;
    SEQUANT_ENFORCE(++evaluations == 1);
    REQUIRE(evaluations == 1);

    if (sequant::assert_behavior() != sequant::AssertBehavior::Abort) {
      REQUIRE_THROWS_AS([] { SEQUANT_ENFORCE(false); }(), sequant::Exception);
      const auto fail_with_message = [] {
      // clang-format off
#line 3000
        SEQUANT_ENFORCE(1 == 0, "invalid input");
        // clang-format on
      };
      REQUIRE_THROWS_WITH(
          fail_with_message(),
          Catch::Matchers::ContainsSubstring("invalid input") &&
              Catch::Matchers::ContainsSubstring(
                  "SEQUANT_ENFORCE(1 == 0) failed") &&
              Catch::Matchers::ContainsSubstring("test_macros.cpp:3000"));
    }
  }
}
