#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "density_fitting_step.hpp"
#include "filter_step.hpp"
#include "optimization_step.hpp"
#include "output_step.hpp"
#include "prettify_step.hpp"
#include "projection_step.hpp"
#include "simplify_step.hpp"
#include "sort_step.hpp"
#include "spintracing_step.hpp"
#include "substitute_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/io/serialization/serialization.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <nlohmann/json.hpp>

#include <string>
#include <string_view>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

namespace {
void seed_expr(ExecutionContext &ctx,
               std::initializer_list<std::wstring_view> defs) {
  ctx.set_data("in", 0, ProcessingData{test_support::expr_data(defs)});
}
}  // namespace

TEST_CASE("simplify step", "[extint]") {
  ExecutionContext ctx;
  seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});

  SECTION("produces an ExpressionData output") {
    SimplifyStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("s.0"));
    REQUIRE_NOTHROW(
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get()));
  }

  SECTION("input is consumed if it is not needed later") {
    ctx.set_needed_later_predicate([](std::string_view) { return false; });
    SimplifyStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE_THROWS_AS(ctx.get_data("in.0"), Exception);
    REQUIRE_NOTHROW(ctx.get_data("s.0"));
  }

  SECTION("input is left intact if it is needed later") {
    ctx.set_needed_later_predicate([](std::string_view) { return true; });
    SimplifyStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE_NOTHROW(ctx.get_data("in.0"));
    REQUIRE_NOTHROW(ctx.get_data("s.0"));
  }

  SECTION("rejects options") {
    SimplifyStep step;
    REQUIRE_THROWS_AS(step.set_options(json::object()), Exception);
  }
}

TEST_CASE("prettify step", "[extint]") {
  ExecutionContext ctx;

  SECTION("uses the conventions of mbpt") {
    // t precedes R among the cardinal tensor labels of mbpt
    seed_expr(ctx, {L"E = R{i1;a1}:N-N-N t{a1;i1}:N-N-N"});
    PrettifyStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    const auto &out =
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get());
    REQUIRE(out.expressions.size() == 1);
    REQUIRE(out.expressions.front() ==
            deserialize<ResultExpr>(L"E = t{a1;i1}:N-N-N R{i1;a1}:N-N-N"));
  }

  SECTION("recanonicalizes unless told not to") {
    seed_expr(ctx, {L"E = f{i3;a3}:N-N-N t{a3;i3}:N-N-N"});
    {
      PrettifyStep step;
      REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
      const auto &out =
          convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get());
      const auto &in =
          convert_data<ExpressionData>(ctx.get_data("in.0").front().data.get());
      // the dummy indices are renamed
      REQUIRE(out.expressions.front() != in.expressions.front());
    }
    {
      PrettifyStep step;
      step.set_options(json{{"recanonicalize", false}});
      REQUIRE(step.run("t", ctx, {"in.0"}) == 1);
      const auto &out =
          convert_data<ExpressionData>(ctx.get_data("t.0").front().data.get());
      const auto &in =
          convert_data<ExpressionData>(ctx.get_data("in.0").front().data.get());
      REQUIRE(out.expressions == in.expressions);
    }
  }

  SECTION("option validation") {
    PrettifyStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"recanonicalize", "yes"}}),
                      Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }
}

TEST_CASE("sort step", "[extint]") {
  ExecutionContext ctx;

  SECTION("sorts the expressions lexicographically") {
    seed_expr(ctx, {L"R1{a1;i1} = + g{a1;i1}:A-H-S + f{a1;i1}:A-H-S"});
    SortStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("s.0"));
    const auto &out =
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get());
    REQUIRE(out.expressions.size() == 1);
    REQUIRE(out.expressions.front() ==
            deserialize<ResultExpr>(
                L"R1{a1;i1} = + f{a1;i1}:A-H-S + g{a1;i1}:A-H-S"));
  }

  SECTION("already sorted input is passed through") {
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S + g{a1;i1}:A-H-S"});
    SortStep step;
    // the unchanged input is aliased as the output
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("s.0"));
    const auto &out =
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get());
    const auto &in =
        convert_data<ExpressionData>(ctx.get_data("in.0").front().data.get());
    REQUIRE(out.expressions == in.expressions);
  }

  SECTION("input is consumed if it is not needed later") {
    seed_expr(ctx, {L"R1{a1;i1} = + g{a1;i1}:A-H-S + f{a1;i1}:A-H-S"});
    ctx.set_needed_later_predicate([](std::string_view) { return false; });
    SortStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    REQUIRE_THROWS_AS(ctx.get_data("in.0"), Exception);
    const auto &out =
        convert_data<ExpressionData>(ctx.get_data("s.0").front().data.get());
    REQUIRE(out.expressions.front() ==
            deserialize<ResultExpr>(
                L"R1{a1;i1} = + f{a1;i1}:A-H-S + g{a1;i1}:A-H-S"));
  }

  SECTION("input is left intact if it is needed later") {
    seed_expr(ctx, {L"R1{a1;i1} = + g{a1;i1}:A-H-S + f{a1;i1}:A-H-S"});
    ctx.set_needed_later_predicate([](std::string_view) { return true; });
    SortStep step;
    REQUIRE(step.run("s", ctx, {"in.0"}) == 1);
    const auto &in =
        convert_data<ExpressionData>(ctx.get_data("in.0").front().data.get());
    REQUIRE(in.expressions.front() ==
            deserialize<ResultExpr>(
                L"R1{a1;i1} = + g{a1;i1}:A-H-S + f{a1;i1}:A-H-S"));
  }

  SECTION("rejects options") {
    SortStep step;
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
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("sub.0"));
  }

  SECTION("block-mode substitution renames only the targeted blocks") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a2;i2} = + f{p1;a1} T1s{a1;p1} f{i1;p2} T1s{p2;i1}"
                    L" T1s{a2;i2}"});
    SubstituteStep step;
    step.set_options(
        json{{"substitutions",
              {{"T1s{a1;p1}", "T1{a1;p1}"}, {"T1s{p1;i1}", "T1{p1;i1}"}}},
             {"tensor_equality_mode", "block"}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    const ExpressionData &out =
        convert_data<ExpressionData>(ctx.get_data("sub.0").front().data.get());
    const std::wstring str =
        io::serialization::to_string(out.expressions.front());
    INFO(toUtf8(str));
    REQUIRE(str.find(L"T1{a_1;p_1}") != std::wstring::npos);
    REQUIRE(str.find(L"T1{p_2;i_1}") != std::wstring::npos);
    REQUIRE(str.find(L"T1s{a_2;i_2}") != std::wstring::npos);
    REQUIRE(str.find(L"T1s{a_1;p_1}") == std::wstring::npos);
    REQUIRE(str.find(L"T1s{p_2;i_1}") == std::wstring::npos);
  }

  SECTION("an input left unchanged is dropped") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"q{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 0);
    REQUIRE_FALSE(ctx.has_data("sub.0"));
  }

  SECTION("keep_unchanged emits unchanged inputs") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"q{a1;i1}", "X{a1;i1}"}}},
                          {"keep_unchanged", true}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("sub.0"));
  }

  SECTION("an unchanged input is dropped even if it is not needed later") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    ctx.set_needed_later_predicate([](std::string_view) { return false; });
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"q{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 0);
    REQUIRE_FALSE(ctx.has_data("sub.0"));
  }

  SECTION("a changed input is consumed if it is not needed later") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    ctx.set_needed_later_predicate([](std::string_view) { return false; });
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    REQUIRE_NOTHROW(ctx.get_data("sub.0"));
  }

  SECTION("result_relabeling renames the result") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    SubstituteStep step;
    step.set_options(json{{"substitutions", {{"f{a1;i1}", "X{a1;i1}"}}},
                          {"result_relabeling", {{"R1", "S1"}}}});
    REQUIRE(step.run("sub", ctx, {"in.0"}) == 1);
    const ExpressionData &out =
        convert_data<ExpressionData>(ctx.get_data("sub.0").front().data.get());
    REQUIRE(out.expressions.size() == 1);
    REQUIRE(out.expressions.front().label() == std::wstring(L"S1"));
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
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
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
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
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
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    FilterStep step;
    step.set_options(json{
        {"groups",
         {{"f_terms", {{"rules", {{{"type", "contains"}, {"label", "f"}}}}}},
          {"g_terms",
           {{"rules", {{{"type", "contains"}, {"label", "g"}}}}}}}}});
    // keep_empty defaults to true -> both groups are emitted
    REQUIRE(step.run("flt", ctx, {"in.0"}) == 2);
  }

  SECTION("index count rules") {
    // f{a1;i1} has 2 indices, T2 has 4
    auto run = [](const json &rule) {
      ExecutionContext ctx;
      seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
      FilterStep step;
      step.set_options(
          json{{"keep_empty", false},
               {"groups", {{"g", {{"rules", json::array({rule})}}}}}});
      return step.run("flt", ctx, {"in.0"});
    };

    REQUIRE(run({{"type", "index_count"}, {"count", 4}}) == 1);
    REQUIRE(run({{"type", "index_count"}, {"count", 3}}) == 0);
    REQUIRE(run({{"type", "index_count"}, {"count", 3}, {"negate", true}}) ==
            1);
    REQUIRE(run({{"type", "index_spaces"}, {"spaces", {{"i", 2}, {"a", 2}}}}) ==
            1);
    REQUIRE(run({{"type", "index_spaces"}, {"spaces", {{"i", 3}}}}) == 0);
    // bounds and label restriction apply to a single tensor
    REQUIRE(run({{"type", "index_spaces"},
                 {"spaces", {{"a", {{"min", 2}}}}}}) == 1);
    REQUIRE(run({{"type", "index_spaces"},
                 {"spaces", {{"a", {{"min", 3}}}}}}) == 0);
    REQUIRE(run({{"type", "index_spaces"},
                 {"spaces", {{"a", {{"max", 1}}}}},
                 {"label", "f"}}) == 1);
    REQUIRE(run({{"type", "index_spaces"},
                 {"spaces", {{"a", {{"min", 2}}}}},
                 {"label", "f"}}) == 0);
    REQUIRE(run({{"type", "index_spaces"},
                 {"spaces", {{"a", {{"min", 2}}}}},
                 {"label", "T2"}}) == 1);
    REQUIRE_THROWS_AS(
        run({{"type", "index_spaces"}, {"spaces", {{"a", {{"foo", 1}}}}}}),
        Exception);
    REQUIRE_THROWS_AS(
        run({{"type", "index_count"}, {"count", 2}, {"label", "f"}}),
        Exception);
    REQUIRE_THROWS_AS(run({{"type", "index_spaces"},
                           {"spaces", {{"a", {{"min", 3}, {"max", 1}}}}}}),
                      Exception);
    REQUIRE_THROWS_AS(run({{"type", "index_count"}, {"count", -1}}), Exception);
    REQUIRE_THROWS_AS(
        run({{"type", "index_spaces"}, {"spaces", json::object()}}), Exception);
  }

  SECTION("keep_empty=false drops empty groups") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{a1;i1}:A-H-S"});
    FilterStep step;
    step.set_options(json{
        {"keep_empty", false},
        {"groups",
         {{"f_terms", {{"rules", {{{"type", "contains"}, {"label", "f"}}}}}},
          {"g_terms",
           {{"rules", {{{"type", "contains"}, {"label", "g"}}}}}}}}});
    REQUIRE(step.run("flt", ctx, {"in.0"}) == 1);
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
    seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    OptimizationStep step;
    REQUIRE(step.run("opt", ctx, {"in.0"}) == 1);
    REQUIRE_NOTHROW(
        convert_data<ExpressionData>(ctx.get_data("opt.0").front().data.get()));
  }

  SECTION("input handling") {
    const std::wstring_view definition =
        L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S";

    ExecutionContext copied;
    seed_expr(copied, {definition});
    copied.set_needed_later_predicate([](std::string_view) { return true; });
    ExecutionContext moved;
    seed_expr(moved, {definition});
    moved.set_needed_later_predicate([](std::string_view) { return false; });

    OptimizationStep step;
    REQUIRE(step.run("opt", copied, {"in.0"}) == 1);
    OptimizationStep other_step;
    REQUIRE(other_step.run("opt", moved, {"in.0"}) == 1);

    REQUIRE_NOTHROW(copied.get_data("in.0"));
    REQUIRE_THROWS_AS(moved.get_data("in.0"), Exception);
    REQUIRE(convert_data<ExpressionData>(
                copied.get_data("opt.0").front().data.get()) ==
            convert_data<ExpressionData>(
                moved.get_data("opt.0").front().data.get()));
  }

  SECTION("replay splits each input into persistent and volatile parts") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    OptimizationStep step;
    step.set_options(json{{"replay", {{"volatile", "T.*"}}}});
    REQUIRE(step.run("opt", ctx, {"in.0"}) == 2);
    REQUIRE(ctx.has_data("opt.persistent"));
    REQUIRE(ctx.has_data("opt.volatile"));
  }

  SECTION("replay errors name the failing input") {
    ExecutionContext ctx;
    ctx.set_data("tree", 0,
                 ProcessingData{test_support::tree_data(
                     {L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});
    OptimizationStep step;
    step.set_options(json{{"replay", {{"volatile", "T.*"}}}});
    REQUIRE_THROWS_WITH(step.run("opt", ctx, {"tree.0"}),
                        Catch::Matchers::ContainsSubstring("tree.0"));
  }

  SECTION("replay consumes inputs that are not needed later") {
    ExecutionContext ctx;
    seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    ctx.add_data_alias("in.0", "in.energy");
    ctx.set_needed_later_predicate([](std::string_view) { return false; });
    OptimizationStep step;
    step.set_options(json{{"replay", {{"volatile", "T.*"}}}});
    REQUIRE(step.run("opt", ctx, {"in.0"}) == 2);
    REQUIRE_THROWS_AS(ctx.get_data("in.0"), Exception);
    REQUIRE(ctx.has_data("opt.persistent"));
    REQUIRE(ctx.has_data("opt.volatile"));
    // The aliases of the consumed input are carried over
    REQUIRE(ctx.has_data("opt.energy"));
    REQUIRE(ctx.has_data("opt.energy.persistent"));
    REQUIRE(ctx.has_data("opt.energy.volatile"));
  }
}
