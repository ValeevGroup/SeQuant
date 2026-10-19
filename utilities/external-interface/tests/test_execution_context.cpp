#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "execution_context.hpp"
#include "processing_data.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <string>
#include <vector>

using namespace sequant;
using namespace sequant::util::extint;

namespace {
ProcessingData expr_data() { return ProcessingData{ExpressionData{}}; }
ProcessingData tree_data() { return ProcessingData{ExportTreeData{}}; }
}  // namespace

TEST_CASE("ExecutionContext::is_valid_id", "[extint]") {
  SECTION("plain ids") {
    REQUIRE(ExecutionContext::is_valid_id("foo", false));
    REQUIRE(ExecutionContext::is_valid_id("foo.0", false));
    REQUIRE(ExecutionContext::is_valid_id("foo_bar.baz", false));
    REQUIRE_FALSE(ExecutionContext::is_valid_id("", false));
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo!", false));
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo bar", false));
  }

  SECTION("selectors require allow_selectors") {
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo[0-2]", false));
    REQUIRE(ExecutionContext::is_valid_id("foo[0-2]", true));
    REQUIRE(ExecutionContext::is_valid_id("foo[0,2,name]", true));
  }

  SECTION("malformed selectors") {
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo[0-2", true));
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo[]", true));
    REQUIRE_FALSE(ExecutionContext::is_valid_id("foo[0-1-2]", true));
  }
}

TEST_CASE("ExecutionContext::expand_id", "[extint]") {
  SECTION("no selector") {
    REQUIRE(ExecutionContext::expand_id("foo") ==
            std::vector<std::string>{"foo"});
    REQUIRE(ExecutionContext::expand_id("foo.0") ==
            std::vector<std::string>{"foo.0"});
  }

  SECTION("ranges and lists") {
    REQUIRE(ExecutionContext::expand_id("foo[0-2]") ==
            std::vector<std::string>{"foo.0", "foo.1", "foo.2"});
    REQUIRE(ExecutionContext::expand_id("foo[0,2]") ==
            std::vector<std::string>{"foo.0", "foo.2"});
    REQUIRE(ExecutionContext::expand_id("foo[0-1,name,5]") ==
            std::vector<std::string>{"foo.0", "foo.1", "foo.name", "foo.5"});
  }

  SECTION("reversed range is normalised") {
    REQUIRE(ExecutionContext::expand_id("foo[2-0]") ==
            std::vector<std::string>{"foo.0", "foo.1", "foo.2"});
  }

  SECTION("invalid id throws") {
    REQUIRE_THROWS_AS(ExecutionContext::expand_id("foo!"), Exception);
  }
}

TEST_CASE("ExecutionContext data storage", "[extint]") {
  ExecutionContext ctx;

  SECTION("set_data then retrieve") {
    ctx.set_data("input", 0, expr_data());

    REQUIRE(ctx.has_data("input.0"));
    REQUIRE_FALSE(ctx.has_data("input.1"));
    REQUIRE(ctx.dataset_size("input.0") == 1);

    auto data = ctx.get_data("input.0");
    REQUIRE(data.size() == 1);
  }

  SECTION("duplicate data id throws") {
    ctx.set_data("input", 0, expr_data());
    REQUIRE_THROWS_AS(ctx.set_data("input", 0, expr_data()), Exception);
  }

  SECTION("invalid prefix throws") {
    REQUIRE_THROWS_AS(ctx.set_data("in put", 0, expr_data()), Exception);
  }

  SECTION("empty counters / dataset throw") {
    REQUIRE_THROWS_AS(ctx.set_data("input", std::vector<std::size_t>{},
                                   std::vector<ProcessingData>{expr_data()}),
                      Exception);
    REQUIRE_THROWS_AS(ctx.set_data("input", std::vector<std::size_t>{0},
                                   std::vector<ProcessingData>{}),
                      Exception);
  }

  SECTION("dataset_size on invalid id throws") {
    REQUIRE_THROWS_AS(ctx.dataset_size("foo!"), Exception);
  }
}

TEST_CASE("ExecutionContext aliasing", "[extint]") {
  ExecutionContext ctx;
  ctx.set_data("input", 0, expr_data());
  ctx.set_data("input", 1, expr_data());

  SECTION("alias to an existing id") {
    ctx.add_data_alias("input.0", "first");
    REQUIRE(ctx.has_data("first"));
    REQUIRE(ctx.ids_are_equivalent("input.0", "first"));
    REQUIRE_FALSE(ctx.ids_are_equivalent("input.1", "first"));
  }

  SECTION("group alias referencing several data objects") {
    ctx.add_data_alias(std::vector<std::string>{"input.0", "input.1"}, "both");
    REQUIRE(ctx.dataset_size("both") == 2);
  }

  SECTION("aliasing a non-existent id throws") {
    REQUIRE_THROWS_AS(ctx.add_data_alias("input.9", "x"), Exception);
  }

  SECTION("alias colliding with a data id throws") {
    REQUIRE_THROWS_AS(ctx.add_data_alias("input.0", "input.1"), Exception);
  }
}

TEST_CASE("convert_data enforces the data kind", "[extint]") {
  ProcessingData expr = expr_data();
  ProcessingData tree = tree_data();

  REQUIRE_NOTHROW(convert_data<ExpressionData>(expr));
  REQUIRE_NOTHROW(convert_data<ExportTreeData>(tree));

  REQUIRE_THROWS_AS(convert_data<ExportTreeData>(expr), Exception);
  REQUIRE_THROWS_WITH(
      convert_data<ExportTreeData>(expr),
      Catch::Matchers::ContainsSubstring("Can't convert ExpressionData"));
}
