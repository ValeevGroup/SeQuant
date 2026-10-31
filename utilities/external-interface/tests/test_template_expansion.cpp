#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "template_expansion.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

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

TEST_CASE("template expansion projects keys and values of an object",
          "[extint]") {
  const json templates = {
      {"t",
       {{"parameters", {{"map", {{"default", nullptr}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "export"},
                                {"inputs", "in"},
                                {"options",
                                 {{"names", "${keys(map)}"},
                                  {"vals", "${values(map)}"}}}}})}}},
      {"keys_named",
       {{"parameters", {{"keys", {{"required", true}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "export"},
                                {"inputs", "in"},
                                {"options", {{"v", "${keys}"}}}}})}}},
      {"bad",
       {{"parameters", {{"p", {{"default", nullptr}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "export"},
                                {"inputs", "in"},
                                {"options", {{"v", "${keys(p)}"}}}}})}}},
      {"embedded",
       {{"parameters", {{"p", {{"required", true}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "export"},
                                {"inputs", "in"},
                                {"options", {{"v", "x_${keys(p)}"}}}}})}}}};

  auto use = [&](const char *tmpl, const json &parameters) {
    return expand(
        templates,
        json::array(
            {{{"template", tmpl}, {"id", "x"}, {"parameters", parameters}}}));
  };

  SECTION("keys and values correspond and keep their JSON types") {
    const json map = {{"a", 1}, {"b", "str"}, {"c", json::array({true})}};
    const json opts = use("t", {{"map", map}}).at(0).at("options");

    const json &names = opts.at("names");
    const json &vals = opts.at("vals");
    REQUIRE(names.is_array());
    REQUIRE(vals.is_array());
    REQUIRE(names.size() == map.size());
    REQUIRE(vals.size() == map.size());
    for (std::size_t i = 0; i < names.size(); ++i) {
      REQUIRE(vals.at(i) == map.at(names.at(i).get<std::string>()));
    }
  }

  SECTION("an empty object gives empty arrays") {
    const json opts = use("t", {{"map", json::object()}}).at(0).at("options");
    REQUIRE(opts.at("names") == json::array());
    REQUIRE(opts.at("vals") == json::array());
  }

  SECTION("a null parameter omits the option") {
    REQUIRE_FALSE(use("t", json::object()).at(0).contains("options"));
  }

  SECTION("a parameter may itself be named keys") {
    const json opts =
        use("keys_named", {{"keys", json::array({1})}}).at(0).at("options");
    REQUIRE(opts.at("v") == json::array({1}));
  }

  SECTION("a non-object parameter is rejected") {
    REQUIRE_THROWS_AS(use("bad", {{"p", json::array({1})}}), Exception);
  }

  SECTION("an array cannot be interpolated into a string") {
    REQUIRE_THROWS_AS(use("embedded", {{"p", {{"a", 1}}}}), Exception);
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

TEST_CASE("template expansion handles conditional steps", "[extint]") {
  const json templates = json::parse(R"({
    "variant": {
      "parameters": {
        "src": { "required": true },
        "replay_src": { "default": [] },
        "replay": { "default": null }
      },
      "steps": [
        { "id": "opt", "kind": "optimize", "inputs": "${src}" },
        { "id": "opt_replay", "kind": "optimize", "inputs": "${replay_src}",
          "when": "${replay_src}", "options": { "replay": "${replay}" } },
        { "id": "tree", "kind": "to_export_tree", "inputs": "opt" },
        { "id": "tree_replay", "kind": "to_export_tree",
          "inputs": "opt_replay", "when": "${replay_src}" }
      ]
    },
    "export": {
      "parameters": { "trees": { "required": true },
                      "grouping": { "required": true } },
      "steps": [
        { "id": "deps", "kind": "dependencies", "inputs": "${trees}" },
        { "id": "cse", "kind": "cse", "inputs": "deps",
          "options": { "merge_inputs": "${values(grouping)}" } },
        { "id": "out", "kind": "export", "inputs": "cse",
          "options": { "grouping": "${grouping}" } }
      ]
    }
  })");

  auto ids = [](const json &steps) {
    std::vector<std::string> result;
    for (const json &step : steps) {
      result.push_back(step.at("id"));
    }
    return result;
  };
  using Ids = std::vector<std::string>;

  SECTION("a step is kept iff its condition is truthy") {
    json out = expand(templates, json::parse(R"([
      { "template": "variant", "id": "a", "parameters": { "src": "x" } }
    ])"));
    REQUIRE(ids(out) == Ids{"a_opt", "a_tree"});
    REQUIRE_FALSE(out.at(0).contains("when"));

    out = expand(templates, json::parse(R"([
      { "template": "variant", "id": "a", "parameters":
        { "src": "x", "replay_src": ["y"], "replay": { "weight": 3 } } }
    ])"));
    REQUIRE(ids(out) ==
            Ids{"a_opt", "a_opt_replay", "a_tree", "a_tree_replay"});
    REQUIRE(out.at(1).at("inputs") == json::array({"y"}));
    REQUIRE(out.at(1).at("options").at("replay").at("weight") == 3);
    REQUIRE_FALSE(out.at(3).contains("when"));
  }

  SECTION("a literal condition is taken as is") {
    const json lit = json::parse(R"({
      "t": { "steps": [
        { "id": "a", "kind": "validate", "inputs": "x", "when": false },
        { "id": "b", "kind": "validate", "inputs": "x", "when": true }
      ] }
    })");
    REQUIRE(ids(expand(lit, json::parse(R"([{"template":"t","id":"p"}])"))) ==
            Ids{"p_b"});
  }

  SECTION("a reference to a skipped step is namespaced, not captured") {
    const json t = json::parse(R"({
      "t": { "steps": [
        { "id": "a", "kind": "validate", "inputs": "x", "when": false },
        { "id": "b", "kind": "validate", "inputs": "a" }
      ] }
    })");
    // The top-level 'a' must not be mistaken for the skipped one
    const json out = expand(t, json::parse(R"([
      { "id": "a", "kind": "validate", "inputs": "x" },
      { "template": "t", "id": "p" }
    ])"));
    REQUIRE(ids(out) == Ids{"a", "p_b"});
    REQUIRE(out.at(1).at("inputs") == "p_a");
  }

  SECTION("optional references resolve against the expanded steps") {
    const json use_export = json::parse(R"({
      "template": "export", "id": "e", "parameters": {
        "trees": ["a_tree", "a_tree_replay?"],
        "grouping": { "Main": "a_opt", "Init": "a_opt_replay.persistent?" } }
    })");

    json steps =
        json::array({json::parse(
                         R"({ "template": "variant", "id": "a", "parameters":
                  { "src": "x" } })"),
                     use_export});
    json out = expand(templates, steps);
    REQUIRE(ids(out) == Ids{"a_opt", "a_tree", "e_deps", "e_cse", "e_out"});
    REQUIRE(out.at(2).at("inputs") == json::array({"a_tree"}));
    REQUIRE(out.at(3).at("options").at("merge_inputs") ==
            json::array({"a_opt"}));
    REQUIRE(out.at(4).at("options").at("grouping") ==
            json::parse(R"({ "Main": "a_opt" })"));

    steps[0]["parameters"]["replay_src"] = json::array({"y"});
    out = expand(templates, steps);
    REQUIRE(out.size() == 7);
    REQUIRE(out.at(4).at("inputs") == json::array({"a_tree", "a_tree_replay"}));
    REQUIRE(out.at(6).at("options").at("grouping") ==
            json::parse(R"({ "Main": "a_opt",
                             "Init": "a_opt_replay.persistent" })"));
  }

  SECTION("a step without any remaining input is removed in turn") {
    const json out = expand(json::object(), json::parse(R"([
      { "id": "a", "kind": "validate", "inputs": "x" },
      { "id": "b", "kind": "cse", "inputs": ["gone?", "nope.sub?"] },
      { "id": "c", "kind": "export", "inputs": "b?",
        "inherit_options_from": "b?" },
      { "id": "d", "kind": "export", "inputs": "a",
        "inherit_options_from": "c?" }
    ])"));
    REQUIRE(ids(out) == Ids{"a", "d"});
    REQUIRE_FALSE(out.at(1).contains("inherit_options_from"));
  }

  SECTION("a reference without a marker is left alone") {
    const json out = expand(json::object(), json::parse(R"([
      { "id": "a", "kind": "cse", "inputs": ["missing"] }
    ])"));
    REQUIRE(out.at(0).at("inputs") == json::array({"missing"}));
  }
}
