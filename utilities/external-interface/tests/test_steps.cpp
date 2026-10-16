#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "canonicalize_step.hpp"
#include "dependency_step.hpp"
#include "read_input_step.hpp"
#include "simplify_step.hpp"
#include "to_export_tree_step.hpp"
#include "validate_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

TEST_CASE("canonicalize step", "[extint]") {
  ExecutionContext ctx;
  ctx.set_data(
      "in", 0,
      ProcessingData{tests::expr_data({L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});

  SECTION("produces an ExpressionData output") {
    CanonicalizeStep step;
    REQUIRE(step.run("canon", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("canon.0"));
    auto out = ctx.get_data("canon.0");
    REQUIRE(out.size() == 1);
    REQUIRE_NOTHROW(convert_data<ExpressionData>(out.front().data.get()));
  }

  SECTION("propagates named input aliases onto the output") {
    ctx.add_data_alias("in.0", "in.energy");
    CanonicalizeStep step;
    REQUIRE(step.run("canon", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("canon.0"));
    REQUIRE(ctx.has_data("canon.energy"));
  }

  SECTION("rejects options") {
    CanonicalizeStep step;
    REQUIRE_THROWS_AS(step.set_options(json::object()), Exception);
  }

  SECTION("rejects ExportTree input") {
    ctx.set_data(
        "tree", 0,
        ProcessingData{tests::tree_data({L"X{a1;i1} = + f{a1;i1}:A-H-S"})});
    CanonicalizeStep step;
    REQUIRE_THROWS_AS(step.run("canon", ctx, {"tree.0"}), Exception);
  }
}

TEST_CASE("validate step", "[extint]") {
  SECTION("valid expression passes") {
    ExecutionContext ctx;
    ctx.set_data(
        "in", 0,
        ProcessingData{tests::expr_data({L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});
    ValidateStep step;
    REQUIRE(step.run("v", ctx, {"in.0"}) == 0);
  }

  SECTION("inconsistent external indices are rejected") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{tests::expr_data(
                     {L"R1{a1;i1} = + f{a1;i1}:A-H-S - f{a1;a2}:A-H-S"})});
    ValidateStep step;
    REQUIRE_THROWS_AS(step.run("v", ctx, {"in.0"}), Exception);
  }
}

TEST_CASE("to_export_tree step", "[extint]") {
  ExecutionContext ctx;
  ctx.set_data("in", 0,
               ProcessingData{tests::expr_data(
                   {L"R1{a1;i1} = + f{a1;i1}:A-H-S",
                    L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});

  ToExportTreeStep step;
  REQUIRE(step.run("tree", ctx, {"in.0"}) == 1);

  auto out = ctx.get_data("tree.0");
  REQUIRE(out.size() == 1);
  const ExportTreeData &data =
      convert_data<ExportTreeData>(out.front().data.get());
  REQUIRE(data.entries.size() == 2);
}

TEST_CASE("dependencies step", "[extint]") {
  SECTION("option validation") {
    DependencyStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(step.set_options(json::object()),
                        Catch::Matchers::ContainsSubstring("at least one"));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("autodetect records producer->consumer dependencies") {
    // Entry 0 produces X; entry 1 (R) consumes X as a leaf -> 1 depends on 0.
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{tests::tree_data(
                     {L"X{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
                      L"R{a1;i1} = + X{a2;i1}:A-H-S f{a1;a2}:A-H-S"})});

    DependencyStep step;
    step.set_options(json{{"autodetect", true}});
    step.run("deps", ctx, {"in.0"});

    const ExportTreeData &data =
        convert_data<ExportTreeData>(ctx.get_data("in.0").front().data.get());
    REQUIRE(data.entries.size() == 2);

    const std::size_t producer_id = data.entries[0].tree->id();
    REQUIRE(data.entries[0].dependencies.empty());
    REQUIRE(std::ranges::find(data.entries[1].dependencies, producer_id) !=
            data.entries[1].dependencies.end());
  }
}
