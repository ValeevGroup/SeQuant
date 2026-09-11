//
// Created by Eduard Valeyev on 10/16/25.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/utility/macros.hpp>

#include <filesystem>

TEST_CASE("macros", "[elements]") {
  SECTION("SEQUANT_ASSERT") {
    std::filesystem::path this_file("tests");
    this_file /= "unit";
    this_file /= "test_macros.cpp";

    this_file.make_preferred();

    const std::string this_file_path = this_file.string();

    CAPTURE(this_file_path);

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
        REQUIRE(std::string_view(ex.what()).find(this_file_path +
                                                 ":1000 in function ") !=
                std::string::npos);
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
        REQUIRE(std::string_view(ex.what()).find(this_file_path +
                                                 ":2000 in function ") !=
                std::string::npos);
#endif
      }
    }
  }
}
