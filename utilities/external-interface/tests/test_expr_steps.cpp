#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "density_fitting_step.hpp"
#include "filter_step.hpp"
#include "optimization_step.hpp"
#include "output_step.hpp"
#include "projection_step.hpp"
#include "simplify_step.hpp"
#include "spintracing_step.hpp"
#include "substitute_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

namespace {
void seed_expr(ExecutionContext &ctx,
               std::initializer_list<std::string_view> defs) {
  ctx.set_data("in", 0, ProcessingData{test_support::expr_data(defs)});
}
}  // namespace

TEST_CASE("simplify step", "[extint]") {
  ExecutionContext ctx;
  seed_expr(ctx, {"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});

  SECTION("produces an ExpressionData output") {
    SimplifyStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("s.0"));
    REQUIRE_NOTHROW(
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get()));
  }

  SECTION("rejects options") {
    SimplifyStep step;
    REQUIRE_THROWS_AS(step.set_options(json::object()), Exception);
  }
}

TEST_CASE("substitute step", "[extint]") {
  SECTION("option validation") {
    SubstituteStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    // substitutions is mandatory
    REQUIRE_THROWS_WITH(
        step.set_options(json{{"tensor_equality_mode", "block"}}),
        Catch::Matchers::ContainsSubstring("mandatory"));
    REQUIRE_THROWS_AS(
        step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}},
                              {"tensor_equality_mode", "nonsense"}}),
        Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("a matching substitution rewrites the expression") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("sub.0"));
  }

  SECTION("an input left unchanged is dropped") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"q{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 0);
    REQUIRE_FALSE(ctx.has_data("sub.0"));
  }

  SECTION("result_relabeling renames the result") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}},
                          {"result_relabeling", {{"R1", "S1"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    const ExpressionData &out =
        convert_data<ExpressionData>(ctx.get_data("sub.0").front().data.get());
    REQUIRE(out.expressions.size() == 1);
    REQUIRE(out.expressions.front().label() == std::string("S1"));
  }
}

TEST_CASE("spintracing step", "[extint]") {
  SECTION("option validation") {
    SpintracingStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(step.set_options(json{{"algorithm", "bogus"}}),
                        Catch::Matchers::ContainsSubstring("Unknown value"));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
    REQUIRE_NOTHROW(step.set_options(json{{"algorithm", "rigorous"}}));
  }

  SECTION("rigorous spintrace yields ExpressionData") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SpintracingStep step;
    REQUIRE(step.run("st", ctx, {"in.0"}) == 1);
    REQUIRE_NOTHROW(
        convert_data<ExpressionData>(ctx.get_data("st.0").front().data.get()));
  }
}

TEST_CASE("project step", "[extint]") {
  SECTION("option validation") {
    ProjectionStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(step.set_options(json{{"method", "primitive"}}),
                        Catch::Matchers::ContainsSubstring("Invalid value"));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
    REQUIRE_NOTHROW(step.set_options(json{{"method", "biorthogonal"}}));
  }
}

TEST_CASE("density_fitting step", "[extint]") {
  SECTION("option validation") {
    DensityFittingStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    // auxiliary_space is required
    REQUIRE_THROWS_WITH(step.set_options(json{{"integral_label", "g"}}),
                        Catch::Matchers::ContainsSubstring("auxiliary_space"));
    REQUIRE_THROWS_AS(step.set_options(json{{"auxiliary_space", 5}}),
                      Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }
}

TEST_CASE("output step", "[extint]") {
  SECTION("option validation") {
    OutputStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(
        step.set_options(json{{"format", "yaml"}}),
        Catch::Matchers::ContainsSubstring("Invalid output format"));
    REQUIRE_THROWS_AS(step.set_options(json{{"annotate_symmetry", "yes"}}),
                      Exception);
    REQUIRE_NOTHROW(step.set_options(json{{"format", "serialize"}}));
  }

  SECTION("produces no output data") {
    test_support::ScopedStreamSilencer silence;
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    OutputStep step;
    REQUIRE(step.run("o", ctx, {"in.0"}) == 0);
  }
}

TEST_CASE("filter step", "[extint]") {
  SECTION("option validation") {
    FilterStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    // groups is mandatory
    REQUIRE_THROWS_WITH(step.set_options(json{{"keep_empty", true}}),
                        Catch::Matchers::ContainsSubstring("mandatory"));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("splits terms into labelled groups") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    FilterStep step;
    step.set_options(json{
        {"groups",
         {{"f_terms", {{"rules", {{{"type", "contains"}, {"label", "f"}}}}}},
          {"g_terms",
           {{"rules", {{{"type", "contains"}, {"label", "g"}}}}}}}}});
    // keep_empty defaults to true -> both groups are emitted
    REQUIRE(step.run("flt", ctx, {"in.0"}) == 2);
  }

  SECTION("keep_empty=false drops empty groups") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    FilterStep step;
    step.set_options(json{
        {"keep_empty", false},
        {"groups",
         {{"f_terms", {{"rules", {{{"type", "contains"}, {"label", "f"}}}}}},
          {"g_terms",
           {{"rules", {{{"type", "contains"}, {"label", "g"}}}}}}}}});
    REQUIRE(step.run("flt", ctx, {"in.0"}) == 1);
  }

  SECTION("label patterns operate on (Unicode) characters") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + λ⁺{a1;i1}:A-H-S"});
    FilterStep step;
    step.set_options(json{
        {"keep_empty", false},
        {"groups",
         {{"lambda", {{"rules", {{{"type", "contains"}, {"label", "λ."}}}}}},
          {"greek",
           {{"rules", {{{"type", "contains"}, {"label", "[αβγ]."}}}}}}}}});
    REQUIRE(step.run("flt", ctx, {"in.0"}) == 1);

    REQUIRE_THROWS_WITH(
        step.set_options(json{
            {"groups",
             {{"greek",
               {{"rules", {{{"type", "contains"}, {"label", "[α-ω]"}}}}}}}}}),
        Catch::Matchers::ContainsSubstring("Unsupported label regex"));
  }
}

TEST_CASE("optimize step", "[extint]") {
  SECTION("option validation") {
    OptimizationStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_WITH(step.set_options(json{{"objective", "bogus"}}),
                        Catch::Matchers::ContainsSubstring("Invalid value"));
    REQUIRE_THROWS_AS(step.set_options(json{{"cse", "bogus"}}), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
    REQUIRE_THROWS_WITH(step.set_options(json{{"replay", {{"weight", 2}}}}),
                        Catch::Matchers::ContainsSubstring("mandatory"));
  }

  SECTION("factorizes a contraction into ExpressionData") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    OptimizationStep step;
    REQUIRE(step.run("opt", ctx, {"in.0"}) == 1);
    REQUIRE_NOTHROW(
        convert_data<ExpressionData>(ctx.get_data("opt.0").front().data.get()));
  }

  SECTION("replay splits each input into persistent and volatile parts") {
    ExecutionContext ctx;
    seed_expr(ctx, {"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    OptimizationStep step;
    step.set_options(json{{"replay", {{"volatile", "T.*"}}}});
    REQUIRE(step.run("opt", ctx, {"in.0"}) == 2);
    REQUIRE(ctx.has_data("opt.persistent"));
    REQUIRE(ctx.has_data("opt.volatile"));
  }
}
