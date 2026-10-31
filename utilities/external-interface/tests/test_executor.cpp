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

TEST_CASE("Executor releases data only after its last use", "[extint]") {
  test_support::ScopedStreamSilencer silence;
  Executor exec;

  SECTION("non-adjacent use") {
    json steps = json::array({read_step(),
                              {{"kind", "sort"}, {"inputs", "in"}},
                              {{"kind", "canonicalize"}, {"inputs", "in"}},
                              {{"kind", "validate"}, {"inputs", "in[0]"}}});
    REQUIRE_NOTHROW(exec.execute(steps));
  }

  SECTION("auto-generated step IDs") {
    json steps =
        json::array({read_step(),
                     {{"kind", "sort"}, {"inputs", "in"}},
                     {{"kind", "canonicalize"}, {"inputs", "in"}},
                     {{"kind", "validate"}, {"inputs", "step1.sort"}}});
    REQUIRE_NOTHROW(exec.execute(steps));
  }

  SECTION("inputs aliased by a later step") {
    // c2 leaves its (already canonical) input in place and merely aliases it
    json steps =
        json::array({read_step(),
                     {{"id", "c1"}, {"kind", "canonicalize"}, {"inputs", "in"}},
                     {{"id", "c2"}, {"kind", "canonicalize"}, {"inputs", "c1"}},
                     {{"kind", "validate"}, {"inputs", "c2"}}});
    REQUIRE_NOTHROW(exec.execute(steps));
  }

  SECTION("IDs referenced in options") {
    json steps = json::array(
        {read_step(),
         {{"id", "t1"}, {"kind", "to_export_tree"}, {"inputs", "in"}},
         {{"id", "t2"}, {"kind", "to_export_tree"}, {"inputs", "in"}},
         {{"kind", "dependencies"},
          {"inputs", "t2"},
          {"options", {{"deps", {{"t2", "t1"}}}}}}});
    REQUIRE_NOTHROW(exec.execute(steps));
  }

  SECTION("all data is released at the end") {
    json steps = json::array({read_step(),
                              {{"id", "s"}, {"kind", "sort"}, {"inputs", "in"}},
                              {{"kind", "validate"}, {"inputs", "s"}}});
    REQUIRE_NOTHROW(exec.execute(steps));

    // IDs remain known but their data is gone
    REQUIRE(exec.context().has_data("in.0"));
    REQUIRE(exec.context().has_data("s"));
    REQUIRE_THROWS_WITH(
        exec.context().get_data("s"),
        Catch::Matchers::ContainsSubstring("has already been released"));
    REQUIRE_THROWS_AS(exec.context().get_data("in"), Exception);
  }

  SECTION("each execution starts from scratch") {
    REQUIRE_NOTHROW(exec.execute(json::array({read_step()})));
    REQUIRE_NOTHROW(exec.execute(json::array({read_step()})));
    REQUIRE_THROWS_WITH(
        exec.execute(json::array({{{"kind", "validate"}, {"inputs", "in"}}})),
        Catch::Matchers::ContainsSubstring("No data available for ID 'in'"));
  }

  SECTION("invalid references report the step") {
    json steps =
        json::array({read_step(),
                     {{"id", "v"}, {"kind", "validate"}, {"inputs", "in[0-"}}});
    REQUIRE_THROWS_WITH(exec.execute(steps),
                        Catch::Matchers::ContainsSubstring(
                            "Invalid reference in step 'v' (validate)"));
  }
}
