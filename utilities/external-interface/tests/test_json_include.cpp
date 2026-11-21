#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "json_include.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;
using Catch::Matchers::ContainsSubstring;

namespace {

namespace fs = std::filesystem;

struct TempDir {
  fs::path path;

  TempDir()
      : path(fs::temp_directory_path() /
             ("extint_include_" + std::to_string(counter()++))) {
    fs::remove_all(path);
    fs::create_directories(path);
  }
  ~TempDir() { fs::remove_all(path); }

  void write(const std::string &name, const std::string &content) const {
    fs::create_directories((path / name).parent_path());
    std::ofstream(path / name) << content;
  }

  static int &counter() {
    static int c = 0;
    return c;
  }
};

}  // namespace

TEST_CASE("JSON includes", "[extint]") {
  TempDir dir;

  SECTION("no includes") {
    dir.write("main.json", R"({"a": [1, 2], /* comment */ "b": {"c": 3}})");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"a": [1, 2], "b": {"c": 3}})"));
  }

  SECTION("value replacement") {
    dir.write("main.json", R"({"x": {"$include": "x.json"}, "y": 1})");
    dir.write("x.json", R"([1, {"k": "v"}])");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"x": [1, {"k": "v"}], "y": 1})"));
  }

  SECTION("object merge") {
    dir.write("main.json", R"({"t": {"$include": "t.json", "local": 2}})");
    dir.write("t.json", R"({"shared": 1})");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"t": {"shared": 1, "local": 2}})"));

    dir.write("t.json", R"({"local": 1})");
    REQUIRE_THROWS_WITH(load_json_with_includes(dir.path / "main.json"),
                        ContainsSubstring("both locally"));

    dir.write("t.json", R"([1])");
    REQUIRE_THROWS_AS(load_json_with_includes(dir.path / "main.json"),
                      Exception);
  }

  SECTION("array splicing") {
    dir.write("main.json", R"({"a": [0, {"$include": "a.json"}, 3]})");
    dir.write("a.json", R"([1, 2])");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"a": [0, 1, 2, 3]})"));

    dir.write("a.json", R"({"o": 1})");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"a": [0, {"o": 1}, 3]})"));
  }

  SECTION("nested includes resolve relative to the including file") {
    dir.write("main.json", R"({"$include": "sub/a.json"})");
    dir.write("sub/a.json", R"({"b": {"$include": "deeper/b.json"}})");
    dir.write("sub/deeper/b.json", R"("leaf")");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"({"b": "leaf"})"));
  }

  SECTION("errors") {
    dir.write("main.json", R"({"x": {"$include": "missing.json"}})");
    REQUIRE_THROWS_AS(load_json_with_includes(dir.path / "main.json"),
                      Exception);

    dir.write("main.json", R"({"x": {"$include": 1}})");
    REQUIRE_THROWS_AS(load_json_with_includes(dir.path / "main.json"),
                      Exception);

    dir.write("main.json", R"({"$include": "a.json"})");
    dir.write("a.json", R"({"$include": "main.json"})");
    REQUIRE_THROWS_WITH(load_json_with_includes(dir.path / "main.json"),
                        ContainsSubstring("Cyclic"));
  }

  SECTION("the same file may be included repeatedly") {
    dir.write("main.json",
              R"([{"$include": "s.json"}, {"$include": "s.json"}])");
    dir.write("s.json", R"("s")");
    REQUIRE(load_json_with_includes(dir.path / "main.json") ==
            json::parse(R"(["s", "s"])"));
  }
}
