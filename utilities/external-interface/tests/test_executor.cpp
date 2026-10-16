#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "executor.hpp"
#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

namespace {
json read_step(std::string id = "in") {
  return json{
      {"id", std::move(id)},
      {"kind", "read_input"},
      {"options", {{"file_path", test_support::input_path("simple.inp")}}}};
}
}  // namespace

TEST_CASE("Executor input validation", "[extint]") {
  test_support::ScopedStreamSilencer silence;
  Executor exec;

  SECTION("steps must be an array") {
    REQUIRE_THROWS_AS(exec.execute(json::object()), Exception);
  }

  SECTION("unknown step kind throws") {
    REQUIRE_THROWS_AS(exec.execute(json::array({{{"kind", "no_such_kind"}}})),
                      Exception);
  }

  SECTION("invalid inputs field type throws") {
    json steps = json::array({{{"kind", "canonicalize"}, {"inputs", 42}}});
    REQUIRE_THROWS_AS(exec.execute(steps), Exception);
  }

  SECTION("a data-consuming step with no inputs throws") {
    json steps = json::array({{{"kind", "canonicalize"}}});
    REQUIRE_THROWS_WITH(exec.execute(steps), Catch::Matchers::ContainsSubstring(
                                                 "at least one input"));
  }

  SECTION("read_input is allowed to run without inputs") {
    REQUIRE_NOTHROW(exec.execute(json::array({read_step()})));
  }
}

TEST_CASE("Executor option handling", "[extint]") {
  test_support::ScopedStreamSilencer silence;
  SECTION("required options missing throws") {
    Executor exec;
    REQUIRE_THROWS_WITH(exec.execute(json::array({{{"kind", "read_input"}}})),
                        Catch::Matchers::ContainsSubstring("requires options"));
  }

  SECTION("options given to a step that takes none throws") {
    Executor exec;
    json steps = json::array({read_step(),
                              {{"kind", "canonicalize"},
                               {"inputs", "in"},
                               {"options", {{"whatever", true}}}}});
    REQUIRE_THROWS_WITH(exec.execute(steps), Catch::Matchers::ContainsSubstring(
                                                 "does not take options"));
  }

  SECTION("inheriting from an unknown step throws") {
    Executor exec;
    json steps = json::array({{{"id", "a"},
                               {"kind", "cse"},
                               {"inputs", "x"},
                               {"inherit_options_from", "ghost"}}});
    REQUIRE_THROWS_WITH(exec.execute(steps),
                        Catch::Matchers::ContainsSubstring("unknown step"));
  }

  SECTION("inheriting across different kinds throws") {
    // The inheriting step comes first so its option resolution (and the kind
    // mismatch) is reached before the referenced step would run.
    Executor exec;
    json steps =
        json::array({{{"id", "b"},
                      {"kind", "cse"},
                      {"inputs", "x"},
                      {"inherit_options_from", "a"}},
                     {{"id", "a"}, {"kind", "optimize"}, {"inputs", "x"}}});
    REQUIRE_THROWS_WITH(exec.execute(steps), Catch::Matchers::ContainsSubstring(
                                                 "can't inherit options"));
  }

  SECTION("cyclic option inheritance throws") {
    Executor exec;
    json steps = json::array({{{"id", "a"},
                               {"kind", "cse"},
                               {"inputs", "x"},
                               {"inherit_options_from", "b"}},
                              {{"id", "b"},
                               {"kind", "cse"},
                               {"inputs", "x"},
                               {"inherit_options_from", "a"}}});
    REQUIRE_THROWS_WITH(exec.execute(steps),
                        Catch::Matchers::ContainsSubstring("Cyclic"));
  }
}

TEST_CASE("Executor step identity", "[extint]") {
  test_support::ScopedStreamSilencer silence;
  SECTION("duplicate step id throws") {
    Executor exec;
    json steps = json::array({read_step("dup"), read_step("dup")});
    REQUIRE_THROWS_WITH(exec.execute(steps), Catch::Matchers::ContainsSubstring(
                                                 "Duplicate step ID"));
  }

  SECTION("a minimal valid pipeline runs") {
    Executor exec;
    json steps = json::array({read_step(),
                              {{"kind", "canonicalize"}, {"inputs", "in"}},
                              {{"kind", "validate"}, {"inputs", "in"}}});
    REQUIRE_NOTHROW(exec.execute(steps));
  }
}
