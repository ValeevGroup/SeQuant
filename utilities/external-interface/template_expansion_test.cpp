#include "template_expansion.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <iostream>
#include <string>

using nlohmann::json;
using sequant::Exception;
using sequant::util::extint::expand_templates;

namespace {

int failures = 0;

void check(bool condition, const std::string &what) {
  if (!condition) {
    std::cerr << "FAILED: " << what << std::endl;
    ++failures;
  }
}

template <typename Callable>
void check_throws(Callable &&callable, const std::string &what) {
  try {
    std::forward<Callable>(callable)();
    std::cerr << "FAILED: expected an exception: " << what << std::endl;
    ++failures;
  } catch (const Exception &) {
    // expected
  }
}

json expand(const json &templates, const json &steps) {
  return expand_templates(steps, templates);
}

void test_value_substitution_preserves_types() {
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

  const json steps = json::array(
      {{{"template", "t"},
        {"id", "x"},
        {"parameters",
         {{"flag", false}, {"count", 2}, {"src", json::array({"a", "b"})}}}}});

  const json out = expand(templates, steps);
  check(out.size() == 1, "single step expanded");
  const json &opts = out.at(0).at("options");
  check(opts.at("optimize").is_boolean() && opts.at("optimize") == false,
        "boolean parameter stays a boolean");
  check(opts.at("min_index_id").is_number() && opts.at("min_index_id") == 2,
        "number parameter stays a number");
  check(out.at(0).at("inputs").is_array() &&
            out.at(0).at("inputs") == json::array({"a", "b"}),
        "array parameter is spliced in as an array");
  check(out.at(0).at("id") == "x_s", "internal id is namespaced");
}

void test_interpolation() {
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
  check(opts.at("output") == "file_res.itfaa", "string interpolation");
  check(opts.at("note") == "flag=true", "boolean interpolates as text");
}

void test_null_omission() {
  const json templates = {
      {"t",
       {{"parameters", {{"replay", {{"default", nullptr}}}}},
        {"steps", json::array({{{"id", "s"},
                                {"kind", "optimize"},
                                {"inputs", "in"},
                                {"options", {{"replay", "${replay}"}}}}})}}}};

  // Default (null) -> option omitted and the now-empty options dropped.
  json out = expand(templates, json::array({{{"template", "t"}, {"id", "x"}}}));
  check(!out.at(0).contains("options"),
        "an options object emptied by a null parameter is dropped");

  // A value given -> option present and preserved.
  out = expand(templates,
               json::array({{{"template", "t"},
                             {"id", "y"},
                             {"parameters", {{"replay", {{"weight", 10}}}}}}}));
  check(out.at(0).at("options").at("replay").at("weight") == 10,
        "a given object parameter is preserved");
}

void test_internal_reference_namespacing() {
  const json templates = {{"t",
                           {{"parameters", {{"src", {{"required", true}}}}},
                            {"steps", json::array({{{"id", "opt"},
                                                    {"kind", "optimize"},
                                                    {"inputs", "${src}"}},
                                                   {{"id", "tree"},
                                                    {"kind", "to_export_tree"},
                                                    {"inputs", "opt"}}})}}}};

  const json steps = json::array(
      {{{"template", "t"},
        {"id", "base"},
        {"parameters", {{"src", json::array({"traced.en", "biorth"})}}}}});

  const json out = expand(templates, steps);
  check(out.size() == 2, "two steps expanded");
  check(out.at(0).at("id") == "base_opt", "first internal id namespaced");
  check(out.at(1).at("id") == "base_tree", "second internal id namespaced");
  check(out.at(1).at("inputs") == "base_opt",
        "internal reference is namespaced to the producing step");
  // External references (parameters) are left untouched.
  check(out.at(0).at("inputs") == json::array({"traced.en", "biorth"}),
        "external references are left untouched");
}

void test_nested_templates() {
  const json templates = {
      {"inner",
       {{"parameters", {{"src", {{"required", true}}}}},
        {"steps",
         json::array(
             {{{"id", "opt"}, {"kind", "optimize"}, {"inputs", "${src}"}}})}}},
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
  check(out.size() == 2, "nested expansion yields two steps");
  check(out.at(0).at("id") == "o_i_opt", "nested id gets both prefixes");
  check(out.at(1).at("id") == "o_tree", "outer id namespaced");
  check(out.at(1).at("inputs") == "o_i_opt",
        "reference to a nested step is fully namespaced");
}

void test_errors() {
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

  check_throws(
      [&] {
        expand(templates,
               json::array({{{"template", "missing"}, {"id", "x"}}}));
      },
      "unknown template");
  check_throws(
      [&] {
        expand(templates, json::array({{{"template", "t"}, {"id", "x"}}}));
      },
      "missing required parameter");
  check_throws(
      [&] {
        expand(templates,
               json::array({{{"template", "t"},
                             {"id", "x"},
                             {"parameters", {{"needed", 1}, {"bogus", 2}}}}}));
      },
      "unknown parameter");
  check_throws(
      [&] {
        expand(templates, json::array({{{"template", "t"}}}));
      },
      "template use without id");
  check_throws(
      [&] {
        expand(templates,
               json::array({{{"template", "interp"},
                             {"id", "x"},
                             {"parameters", {{"obj", json::object()}}}}}));
      },
      "non-scalar interpolation");
  check_throws(
      [&] {
        expand(templates, json::array({{{"template", "a"}, {"id", "x"}}}));
      },
      "cyclic template use");
  check_throws(
      [&] {
        expand(templates, json::array({{{"template", "t"},
                                        {"id", "a.b"},
                                        {"parameters", {{"needed", 1}}}}}));
      },
      "template use id containing a dot");
}

}  // namespace

int main() {
  test_value_substitution_preserves_types();
  test_interpolation();
  test_null_omission();
  test_internal_reference_namespacing();
  test_nested_templates();
  test_errors();

  if (failures > 0) {
    std::cerr << failures << " check(s) failed" << std::endl;
    return 1;
  }

  std::cout << "All template expansion checks passed" << std::endl;
  return 0;
}
