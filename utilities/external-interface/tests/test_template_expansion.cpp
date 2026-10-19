#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "template_expansion.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

namespace {
json expand(const json &templates, const json &steps) {
  return expand_templates(steps, templates);
}
}  // namespace

TEST_CASE("template expansion substitutes parameters", "[extint]") {
  SECTION("a bound value keeps its JSON type") {
    const json templates = {
        {"t",
         {{"parameters",
           {{"flag", {{"default", true}}},
            {"count", {{"default", 1}}},
            {"src", {{"required", true}}}}},
          {"steps", json::array({{{"id", "s"},
                                  {"kind", "export"},
                                  {"inputs", "${src}"},
                                  {"options",
                                   {{"optimize", "${flag}"},
                                    {"min_index_id", "${count}"}}}}})}}}};

    const json steps = json::array({{{"template", "t"},
                                     {"id", "x"},
                                     {"parameters",
                                      {{"flag", false},
                                       {"count", 2},
                                       {"src", json::array({"a", "b"})}}}}});

    const json out = expand(templates, steps);
    REQUIRE(out.size() == 1);
    const json &opts = out.at(0).at("options");
    // A placeholder that is the entire value adopts the parameter's type
    // rather than being stringified.
    REQUIRE(opts.at("optimize").is_boolean());
    REQUIRE(opts.at("optimize") == false);
    REQUIRE(opts.at("min_index_id").is_number());
    REQUIRE(opts.at("min_index_id") == 2);
    REQUIRE(out.at(0).at("inputs").is_array());
    REQUIRE(out.at(0).at("inputs") == json::array({"a", "b"}));
    REQUIRE(out.at(0).at("id") == "x_s");
  }

  SECTION("a placeholder inside a string interpolates as text") {
    const json templates = {
        {"t",
         {{"parameters",
           {{"name", {{"required", true}}}, {"on", {{"default", true}}}}},
          {"steps", json::array({{{"id", "s"},
                                  {"kind", "export"},
                                  {"inputs", "in"},
                                  {"options",
                                   {{"output", "file_${name}.itfaa"},
                                    {"note", "flag=${on}"}}}}})}}}};

    const json steps = json::array(
        {{{"template", "t"}, {"id", "x"}, {"parameters", {{"name", "res"}}}}});

    const json out = expand(templates, steps);
    const json &opts = out.at(0).at("options");
    REQUIRE(opts.at("output") == "file_res.itfaa");
    REQUIRE(opts.at("note") == "flag=true");
  }

  SECTION("a null parameter omits its option") {
    const json templates = {
        {"t",
         {{"parameters", {{"replay", {{"default", nullptr}}}}},
          {"steps", json::array({{{"id", "s"},
                                  {"kind", "optimize"},
                                  {"inputs", "in"},
                                  {"options", {{"replay", "${replay}"}}}}})}}}};

    // Default (null) -> option omitted and the now-empty options dropped.
    json out =
        expand(templates, json::array({{{"template", "t"}, {"id", "x"}}}));
    REQUIRE_FALSE(out.at(0).contains("options"));

    // A value given -> option present and preserved.
    out =
        expand(templates,
               json::array({{{"template", "t"},
                             {"id", "y"},
                             {"parameters", {{"replay", {{"weight", 10}}}}}}}));
    REQUIRE(out.at(0).at("options").at("replay").at("weight") == 10);
  }
}

TEST_CASE("template expansion namespaces internal step IDs", "[extint]") {
  SECTION("internal references are namespaced, external ones are not") {
    const json templates = {
        {"t",
         {{"parameters", {{"src", {{"required", true}}}}},
          {"steps",
           json::array(
               {{{"id", "opt"}, {"kind", "optimize"}, {"inputs", "${src}"}},
                {{"id", "tree"},
                 {"kind", "to_export_tree"},
                 {"inputs", "opt"}}})}}}};

    const json steps = json::array(
        {{{"template", "t"},
          {"id", "base"},
          {"parameters", {{"src", json::array({"traced.en", "biorth"})}}}}});

    const json out = expand(templates, steps);
    REQUIRE(out.size() == 2);
    REQUIRE(out.at(0).at("id") == "base_opt");
    REQUIRE(out.at(1).at("id") == "base_tree");
    REQUIRE(out.at(1).at("inputs") == "base_opt");
    // External references (bound through parameters) are left untouched.
    REQUIRE(out.at(0).at("inputs") == json::array({"traced.en", "biorth"}));
  }

  SECTION("nested templates stack both prefixes") {
    const json templates = {
        {"inner",
         {{"parameters", {{"src", {{"required", true}}}}},
          {"steps", json::array({{{"id", "opt"},
                                  {"kind", "optimize"},
                                  {"inputs", "${src}"}}})}}},
        {"outer",
         {{"parameters", {{"src", {{"required", true}}}}},
          {"steps", json::array({{{"template", "inner"},
                                  {"id", "i"},
                                  {"parameters", {{"src", "${src}"}}}},
                                 {{"id", "tree"},
                                  {"kind", "to_export_tree"},
                                  {"inputs", "i_opt"}}})}}}};

    const json steps = json::array(
        {{{"template", "outer"}, {"id", "o"}, {"parameters", {{"src", "x"}}}}});

    const json out = expand(templates, steps);
    REQUIRE(out.size() == 2);
    REQUIRE(out.at(0).at("id") == "o_i_opt");
    REQUIRE(out.at(1).at("id") == "o_tree");
    REQUIRE(out.at(1).at("inputs") == "o_i_opt");
  }
}

TEST_CASE("template expansion rejects malformed uses", "[extint]") {
  const json templates = {
      {"t",
       {{"parameters", {{"needed", {{"required", true}}}}},
        {"steps",
         json::array(
             {{{"id", "s"}, {"kind", "validate"}, {"inputs", "in"}}})}}},
      {"interp",
       {{"parameters", {{"obj", {{"required", true}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "export"},
                                {"inputs", "in"},
                                {"options", {{"output", "a${obj}b"}}}}})}}},
      {"a", {{"steps", json::array({{{"template", "b"}, {"id", "ab"}}})}}},
      {"b", {{"steps", json::array({{{"template", "a"}, {"id", "ba"}}})}}}};

  SECTION("unknown template name") {
    REQUIRE_THROWS_AS(
        expand(templates,
               json::array({{{"template", "missing"}, {"id", "x"}}})),
        Exception);
  }

  SECTION("missing required parameter") {
    REQUIRE_THROWS_AS(
        expand(templates, json::array({{{"template", "t"}, {"id", "x"}}})),
        Exception);
  }

  SECTION("unknown parameter") {
    REQUIRE_THROWS_AS(
        expand(templates,
               json::array({{{"template", "t"},
                             {"id", "x"},
                             {"parameters", {{"needed", 1}, {"bogus", 2}}}}})),
        Exception);
  }

  SECTION("template use without an id") {
    REQUIRE_THROWS_AS(expand(templates, json::array({{{"template", "t"}}})),
                      Exception);
  }

  SECTION("a non-scalar cannot be interpolated into a string") {
    REQUIRE_THROWS_AS(
        expand(templates,
               json::array({{{"template", "interp"},
                             {"id", "x"},
                             {"parameters", {{"obj", json::object()}}}}})),
        Exception);
  }

  SECTION("cyclic template use") {
    REQUIRE_THROWS_AS(
        expand(templates, json::array({{{"template", "a"}, {"id", "x"}}})),
        Exception);
  }

  SECTION("a template use id may not contain a dot") {
    REQUIRE_THROWS_AS(
        expand(templates, json::array({{{"template", "t"},
                                        {"id", "a.b"},
                                        {"parameters", {{"needed", 1}}}}})),
        Exception);
  }
}
