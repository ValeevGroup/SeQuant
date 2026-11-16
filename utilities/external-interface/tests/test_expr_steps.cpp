#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "density_fitting_step.hpp"
#include "filter_step.hpp"
#include "optimization_step.hpp"
#include "output_step.hpp"
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

  SECTION("deprecated objectives are rejected") {
    for (const char *objective :
         {"DenseSize", "DensePeakSize", "DensePeakSizeBatched"}) {
      OptimizationStep step;
      REQUIRE_THROWS_WITH(step.set_options(json{{"objective", objective}}),
                          Catch::Matchers::ContainsSubstring("Invalid value"));
    }
    for (const char *objective : {"DenseFLOPs", "DenseTimeSpace"}) {
      OptimizationStep step;
      REQUIRE_NOTHROW(step.set_options(json{{"objective", objective}}));
    }
  }

  SECTION("options that depend on the objective") {
    // DenseTimeSpace is the default objective
    {
      OptimizationStep step;
      REQUIRE_THROWS_WITH(
          step.set_options(json{{"intermediate_size_penalty", 1.0}}),
          Catch::Matchers::ContainsSubstring("DenseFLOPs"));
    }
    {
      OptimizationStep step;
      REQUIRE_THROWS_WITH(step.set_options(json{{"cse", "subnet"}}),
                          Catch::Matchers::ContainsSubstring("DenseFLOPs"));
    }
    {
      OptimizationStep step;
      REQUIRE_NOTHROW(step.set_options(json{{"cse", "none"}}));
    }
    {
      OptimizationStep step;
      REQUIRE_NOTHROW(
          step.set_options(json{{"objective", "DenseFLOPs"},
                                {"cse", "subnet"},
                                {"intermediate_size_penalty", 1.0}}));
    }
  }

  SECTION("batch options") {
    const json batch = {
        {"peak_threshold", 1e9}, {"contracted", {"a"}}, {"target_size", 100}};
    {
      OptimizationStep step;
      REQUIRE_NOTHROW(step.set_options(
          json{{"objective", "DenseTimeSpaceBatched"}, {"batch", batch}}));
    }
    {
      OptimizationStep step;
      REQUIRE_THROWS_WITH(
          step.set_options(json{{"objective", "DenseTimeSpaceBatched"}}),
          Catch::Matchers::ContainsSubstring("requires the 'batch'"));
    }
    {
      OptimizationStep step;
      REQUIRE_THROWS_WITH(
          step.set_options(json{{"batch", batch}}),
          Catch::Matchers::ContainsSubstring("DenseTimeSpaceBatched"));
    }

    auto with_batch = [](json b) {
      return json{{"objective", "DenseTimeSpaceBatched"}, {"batch", b}};
    };
    for (const char *mandatory : {"peak_threshold", "target_size"}) {
      json incomplete = batch;
      incomplete.erase(mandatory);
      OptimizationStep step;
      REQUIRE_THROWS_WITH(step.set_options(with_batch(incomplete)),
                          Catch::Matchers::ContainsSubstring("mandatory"));
    }
    {
      json no_space = batch;
      no_space.erase("contracted");
      OptimizationStep step;
      REQUIRE_THROWS_AS(step.set_options(with_batch(no_space)), Exception);
    }
    {
      json missing_size = batch;
      missing_size["target_size"] = {{"i", 10}};
      OptimizationStep step;
      REQUIRE_THROWS_WITH(step.set_options(with_batch(missing_size)),
                          Catch::Matchers::ContainsSubstring("'a'"));
    }
    {
      json external_only = batch;
      external_only.erase("contracted");
      external_only["external"] = {"a"};
      OptimizationStep step;
      REQUIRE_THROWS_WITH(
          step.set_options(with_batch(external_only)),
          Catch::Matchers::ContainsSubstring("spectator_indices"));
      external_only["spectator_indices"] = true;
      OptimizationStep other;
      REQUIRE_NOTHROW(other.set_options(with_batch(external_only)));
    }
    {
      json spectators = batch;
      spectators["spectator_indices"] = true;
      OptimizationStep step;
      REQUIRE_THROWS_AS(step.set_options(with_batch(spectators)), Exception);
    }
    {
      json persistent = batch;
      persistent["persistent_only"] = true;
      OptimizationStep step;
      REQUIRE_THROWS_WITH(step.set_options(with_batch(persistent)),
                          Catch::Matchers::ContainsSubstring("replay"));
      json options = with_batch(persistent);
      options["replay"] = {{"volatile", "T2"}};
      OptimizationStep other;
      REQUIRE_NOTHROW(other.set_options(options));
    }
    {
      json unknown_space = batch;
      unknown_space["contracted"] = {"no_such_space"};
      ExecutionContext ctx;
      seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
      OptimizationStep step;
      step.set_options(with_batch(unknown_space));
      REQUIRE_THROWS_WITH(step.run("opt", ctx, {"in.0"}),
                          Catch::Matchers::ContainsSubstring("no_such_space"));
    }
    {
      json bad = batch;
      bad["unknown"] = true;
      OptimizationStep step;
      REQUIRE_THROWS_AS(step.set_options(with_batch(bad)), Exception);
    }
    {
      json bad = batch;
      bad["peak_threshold"] = -1;
      OptimizationStep step;
      REQUIRE_THROWS_AS(step.set_options(with_batch(bad)), Exception);
    }
  }

  SECTION("factorizes with each objective") {
    for (const json &options : {json{{"objective", "DenseFLOPs"}},
                                json{{"objective", "DenseTimeSpace"}},
                                json{{"objective", "DenseTimeSpaceBatched"},
                                     {"batch",
                                      {{"peak_threshold", 1e9},
                                       {"contracted", {"i"}},
                                       {"target_size", 10}}}}}) {
      ExecutionContext ctx;
      seed_expr(ctx, {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S "
                      L"T1{a2;i2}:A-H-S"});
      OptimizationStep step;
      step.set_options(options);
      REQUIRE(step.run("opt", ctx, {"in.0"}) == 1);
      REQUIRE_NOTHROW(convert_data<ExpressionData>(
          ctx.get_data("opt.0").front().data.get()));
    }
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
