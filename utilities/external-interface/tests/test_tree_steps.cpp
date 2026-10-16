#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "cse_step.hpp"
#include "export_step.hpp"
#include "index_batching_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <utility>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

TEST_CASE("cse step", "[extint]") {
  SECTION("option validation") {
    CSEStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"min_usage", "two"}}), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"merge_inputs", 1}}), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("a tree with no common subexpression is passed through") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{tests::tree_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    CSEStep step;
    REQUIRE(step.run("cse", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("cse.0"));
    REQUIRE_NOTHROW(
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get()));
  }
}

TEST_CASE("batch_indices step", "[extint]") {
  SECTION("option validation") {
    IndexBatchingStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"min_unbatched", "two"}}),
                      Exception);
    REQUIRE_THROWS_WITH(
        step.set_options(json{{"selection_strategy", "middle"}}),
        Catch::Matchers::ContainsSubstring("strategy"));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("a result with too few indices is left unbatched") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{tests::tree_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    IndexBatchingStep
        step;  // min_unbatched defaults to 2, result has 2 indices
    REQUIRE(step.run("b", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("b.0"));
  }
}

TEST_CASE("export step", "[extint]") {
  SECTION("option validation") {
    // A fresh step per case: set_options accumulates into members, so reusing
    // one instance would let an earlier call satisfy a later "missing X" check.
    REQUIRE_THROWS_AS(ExportStep{}.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(ExportStep{}.set_options(json{{"output", "out.itfaa"}}),
                        Catch::Matchers::ContainsSubstring("language"));
    REQUIRE_THROWS_WITH(ExportStep{}.set_options(json{{"language", "itf"}}),
                        Catch::Matchers::ContainsSubstring("output"));
    REQUIRE_THROWS_AS(ExportStep{}.set_options(
                          json{{"language", 5}, {"output", "out.itfaa"}}),
                      Exception);
    REQUIRE_THROWS_AS(
        ExportStep{}.set_options(json{
            {"language", "itf"}, {"output", "out.itfaa"}, {"unknown", true}}),
        Exception);
  }

  SECTION("a dependency cycle is rejected when ordering the results") {
    // Build two independent trees, then wire a mutual dependency between them.
    ExportTreeData two = tests::tree_data(
        {L"X{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"Y{a1;i1} = + g{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    REQUIRE(two.entries.size() == 2);

    ExportTreeData data_a;
    data_a.entries.push_back(two.entries[0]);
    ExportTreeData data_b;
    data_b.entries.push_back(two.entries[1]);

    const std::size_t id_a = data_a.entries[0].tree->id();
    const std::size_t id_b = data_b.entries[0].tree->id();
    data_a.entries[0].dependencies = {id_b};
    data_b.entries[0].dependencies = {id_a};

    ExecutionContext ctx;
    ctx.set_data("A", 0, ProcessingData{std::move(data_a)});
    ctx.set_data("B", 0, ProcessingData{std::move(data_b)});

    ExportStep step;
    step.set_options(json{{"language", "itf"},
                          {"output", "cycle_test_output.itfaa"},
                          {"meta", json::object()}});
    REQUIRE_THROWS_WITH(
        step.run("exp", ctx, {"A.0", "B.0"}),
        Catch::Matchers::ContainsSubstring("Impossible dependencies"));
  }
}
