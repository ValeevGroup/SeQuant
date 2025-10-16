//
// Created by Eduard Valeyev on 10/16/25.
//

#include <catch2/catch_test_macros.hpp>

#include "catch2_sequant.hpp"

#include <SeQuant/core/utility/macros.hpp>

TEST_CASE("macros", "[elements]") {
  SECTION("SEQUANT_ASSERT") {
#if SEQUANT_ASSERT_BEHAVIOR == SEQUANT_ASSERT_THROW
    try {
      // clang-format off
#line 1000  // to make sure the line number of the next line is fixed
sequant::assert_failed("test");
      // clang-format on
    } catch (sequant::Exception& ex) {
      // see #line up there
      std::cout << ex.what() << std::endl;
      REQUIRE(ex.what().find("tests/unit/test_macros.cpp:1000 in function ") !=
              std::string::npos);
    }
    try {
      // clang-format off
#line 2000  // to make sure the line number of the next line is fixed
SEQUANT_ASSERT(1 == 0 && "1 != 0", "testing SEQUANT_ASSERT");
      // clang-format on
    } catch (sequant::Exception& ex) {
      // see #line up there
      std::cout << ex.what() << std::endl;
      REQUIRE(ex.what().find("tests/unit/test_macros.cpp:2000 in function ") !=
              std::string::npos);
    }
#endif
  }
}
