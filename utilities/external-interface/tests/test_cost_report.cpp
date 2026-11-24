#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "cost_report_step.hpp"

#include "temp_dir.hpp"
#include "test_extint.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

namespace {

/// Reads the whole file; the stream is closed on return, so that the file can
/// be removed afterwards (Windows refuses to delete an open file)
std::string read_file(const std::filesystem::path &path) {
  std::ifstream in(path);
  return std::string(std::istreambuf_iterator<char>(in), {});
}

}  // namespace

using test_support::TempDir;

TEST_CASE("cost_report step", "[extint]") {
  SECTION("option validation") {
    auto rejects = [](const json &options, const char *message) {
      CostReportStep step;
      REQUIRE_THROWS_WITH(step.set_options(options),
                          Catch::Matchers::ContainsSubstring(message));
    };
    rejects(json::array(), "JSON object");
    rejects(json{{"bogus", 1}}, "Unknown option");
    rejects(json{{"output", ""}}, "non-empty");
    rejects(json{{"top_n", -1}}, "non-negative");
    rejects(json{{"space_symbols", {{"i", "OO"}}}}, "single character");
    rejects(json{{"space_symbols", {{"i", ""}}}}, "single character");
    // Two code points, each a single non-ASCII character ("ΚΚ")
    rejects(json{{"space_symbols", {{"i", "\xce\x9a\xce\x9a"}}}},
            "single character");
    // One code point beyond U+FFFF ("𝒪", four bytes in UTF-8), which a
    // 16-bit wchar_t holds as two units
    CHECK_NOTHROW(CostReportStep{}.set_options(
        json{{"space_symbols", {{"i", "\xf0\x9d\x92\xaa"}}}}));
    rejects(json{{"cache", {{"min_repeats", 2}}}}, "mandatory");
    rejects(json{{"cache", {{"volatile", "T2"}, {"bogus", 1}}}},
            "Unknown option");
  }

  SECTION("names each section after the expression's result") {
    const TempDir tmp;
    const std::filesystem::path report = tmp.path / "report.md";
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    ctx.set_data("in", 1,
                 ProcessingData{test_support::expr_data(
                     {L"R2{a1;i1} = + f{a1;i1}:A-H-S"})});
    ctx.add_data_alias("in.0", "in.singles");

    CostReportStep step;
    step.set_options(json{{"output", report.string()}});
    // in.0 is reachable under two names but analyzed once; otherwise R1 would
    // be a duplicate
    REQUIRE(step.run("cost", ctx, {"in.0", "in.singles", "in.1"}) == 0);

    const std::string text = read_file(report);
    REQUIRE_THAT(text, Catch::Matchers::ContainsSubstring("\n## R1\n"));
    REQUIRE_THAT(text, Catch::Matchers::ContainsSubstring("\n## R2\n"));
  }

  SECTION("spin-annotated indices take their space's symbol") {
    const TempDir tmp;
    const std::filesystem::path report = tmp.path / "report.md";
    ExecutionContext ctx;
    ctx.set_data(
        "in", 0,
        ProcessingData{test_support::expr_data(
            {L"R1{a↑1;i↑1} = + g{i↑2,i↑3;a↑2,a↑3}:A-H-S T1{a↑2;i↑2}:A-H-S "
             L"T2{a↑1,a↑3;i↑1,i↑3}:A-H-S"})});
    CostReportStep step;
    step.set_options(json{{"output", report.string()},
                          {"space_symbols", {{"i", "O"}, {"a", "V"}}}});
    REQUIRE(step.run("cost", ctx, {"in.0"}) == 0);

    const std::string text = read_file(report);
    REQUIRE_THAT(text, Catch::Matchers::ContainsSubstring("| OV"));
  }

  SECTION("a space symbol may be any single Unicode character") {
    const TempDir tmp;
    const std::filesystem::path report = tmp.path / "report.md";
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + g{i2,i3;a2,a3}:A-H-S T1{a2;i2}:A-H-S "
                      L"T2{a1,a3;i1,i3}:A-H-S"})});
    CostReportStep step;
    // Greek capital omicron and kappa, two bytes each in UTF-8
    step.set_options(
        json{{"output", report.string()},
             {"space_symbols", {{"i", "\xce\x9f"}, {"a", "\xce\x9a"}}}});
    REQUIRE(step.run("cost", ctx, {"in.0"}) == 0);

    const std::string text = read_file(report);
    REQUIRE_THAT(text, Catch::Matchers::ContainsSubstring(
                           "| \xce\x9f\xce\x9a |"));  // "| ΟΚ |"
  }

  SECTION("space_symbols keys must be registered index spaces") {
    // Checked when the step runs, against the spaces registered then
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});
    const TempDir tmp;
    const std::filesystem::path &dir = tmp.path;
    auto rejects_key = [&](const std::string &key, const char *reason) {
      CostReportStep step;
      CHECK_NOTHROW(
          step.set_options(json{{"output", (dir / "report.md").string()},
                                {"space_symbols", {{"i", "O"}, {key, "Z"}}}}));
      CHECK_THROWS_WITH(step.run("cost", ctx, {"in.0"}),
                        Catch::Matchers::ContainsSubstring(key) &&
                            Catch::Matchers::ContainsSubstring(reason));
    };
    rejects_key("no_such_space", "Unknown index space");
    // A registered spin-annotated label is rejected too: symbols are looked
    // up by the label without its spin annotation
    rejects_key("i\xe2\x86\x91", "spin annotation");  // "i↑" in UTF-8
    CHECK_FALSE(std::filesystem::exists(dir / "report.md"));
  }

  SECTION("creates the output's directory and dumps trees next to it") {
    // The output directory does not exist yet, so the trees can only land
    // there if it is created before they are written
    const TempDir tmp;
    const std::filesystem::path out_dir = tmp.path / "nested";
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    CostReportStep step;
    step.set_options(json{{"output", (out_dir / "report.md").string()},
                          {"dump_tree", true}});
    REQUIRE(step.run("cost", ctx, {"in.0"}) == 0);

    const std::string text = read_file(out_dir / "report.md");
    const std::string tree = read_file(out_dir / "report.R1.tree.txt");
    REQUIRE_THAT(text, Catch::Matchers::ContainsSubstring("\n## R1\n"));
    REQUIRE_THAT(tree, Catch::Matchers::ContainsSubstring("T2"));
  }

  SECTION("writes no tree dumps if the report cannot be written") {
    const TempDir tmp;
    const std::filesystem::path &dir = tmp.path;
    // A directory where the report should go makes opening it fail
    std::filesystem::create_directories(dir / "report.md");
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    CostReportStep step;
    step.set_options(
        json{{"output", (dir / "report.md").string()}, {"dump_tree", true}});
    CHECK_THROWS_WITH(step.run("cost", ctx, {"in.0"}),
                      Catch::Matchers::ContainsSubstring("Cannot open"));
    CHECK_FALSE(std::filesystem::exists(dir / "report.R1.tree.txt"));
  }

  SECTION("costs the factorization the input carries") {
    const TempDir tmp;
    const std::filesystem::path &dir = tmp.path;
    ExecutionContext ctx;
    // Not left-deep: the two inner products are independent intermediates
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R{a1;a2} = (A{a1;i1} B{i1;i2}) (C{i2;i3} D{i3;a2})"})});
    CostReportStep step;
    step.set_options(
        json{{"output", (dir / "report.md").string()}, {"dump_tree", true}});
    REQUIRE(step.run("cost", ctx, {"in.0"}) == 0);

    const std::string tree = read_file(dir / "report.R.tree.txt");
    REQUIRE(tree ==
            "((A(a_1,i_1) * B(i_1,i_2)) * (C(i_2,i_3) * D(i_3,a_2)))\n");
  }

  SECTION("ranks and truncates the tables") {
    const TempDir tmp;
    const std::filesystem::path report = tmp.path / "report.md";
    ExecutionContext ctx;
    // Intermediates a×i and the larger a×x (x is i ∪ a), both of shape XX as no
    // space_symbols are given
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R{a1;i1} = (A{a1;a2} B{a2;i2}) C{i2;i1} + "
                      L"(D{a1;i2} E{i2;x1}) F{x1;i1}"})});
    auto run = [&](int top_n) {
      CostReportStep step;
      step.set_options(json{{"output", report.string()}, {"top_n", top_n}});
      REQUIRE(step.run("cost", ctx, {"in.0"}) == 0);
      return read_file(report);
    };
    using Catch::Matchers::ContainsSubstring;

    const std::string all = run(2);
    INFO(all);
    const auto ax = all.find("| 1 | I(a_1,x_1) |");
    const auto ai = all.find("| 2 | I(a_1,i_2) |");
    REQUIRE(ax != std::string::npos);
    REQUIRE(ai != std::string::npos);
    // The census reports a shape at the size of its largest intermediate
    const auto row_mb = [&](std::size_t row) {
      // "| rank | result | spaces | MB |"
      std::size_t pos = row;
      for (int bar = 0; bar < 4; ++bar) pos = all.find('|', pos) + 1;
      return all.substr(pos, all.find('|', pos) - pos);
    };
    REQUIRE_THAT(all, ContainsSubstring("| XX |" + row_mb(ax) + "| 2 | 2 |"));

    const std::string top = run(1);
    REQUIRE_THAT(top, ContainsSubstring("| 1 | I(a_1,x_1) |"));
    REQUIRE_THAT(top, !ContainsSubstring("| 2 | I("));
  }

  SECTION("rejects several expressions under one name") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::expr_data(
                     {L"R1{a1;i1} = + f{a1;i1}:A-H-S",
                      L"R1{a1;i1} = + g{a1;i1}:A-H-S"})});
    CostReportStep step;
    REQUIRE_THROWS_WITH(step.run("cost", ctx, {"in.0"}),
                        Catch::Matchers::ContainsSubstring("'R1'"));
  }

  SECTION("names must also be distinct, safe file names for tree dumps") {
    const TempDir tmp;
    auto run = [&tmp](std::vector<ResultExpr> exprs, bool dump_tree) {
      ExecutionContext ctx;
      ExpressionData data;
      data.expressions = std::move(exprs);
      ctx.set_data("in", 0, ProcessingData{std::move(data)});
      CostReportStep step;
      step.set_options(json{{"output", (tmp.path / "report.md").string()},
                            {"dump_tree", dump_tree}});
      step.run("cost", ctx, {"in.0"});
    };
    const ResultExpr r = deserialize<ResultExpr>(L"R1{a1;i1} = f{a1;i1}");
    ResultExpr lower = r, separator = r;
    lower.set_label(L"r1");
    separator.set_label(L"../R1");
    const ResultExpr unlabeled(
        bra(r.bra()), ket(r.ket()), aux(r.aux()), r.symmetry(),
        r.braket_symmetry(), r.column_symmetry(), std::nullopt, r.expression());

    // Checked before anything is written
    REQUIRE_THROWS_WITH(run({r, lower}, true),
                        Catch::Matchers::ContainsSubstring("ignoring case"));
    REQUIRE_THROWS_WITH(run({separator}, true),
                        Catch::Matchers::ContainsSubstring("path separator"));
    REQUIRE_THROWS_WITH(
        run({unlabeled, unlabeled}, false),
        Catch::Matchers::ContainsSubstring("without a result label"));
    ResultExpr placeholder = r;
    placeholder.set_label(L"Unnamed");
    REQUIRE_THROWS_WITH(run({placeholder, unlabeled}, false),
                        Catch::Matchers::ContainsSubstring("'Unnamed'"));
    CHECK_FALSE(std::filesystem::exists(tmp.path / "report.md"));
  }

  SECTION("errors name the failing input") {
    ExecutionContext ctx;
    ctx.set_data("tree", 0,
                 ProcessingData{test_support::tree_data(
                     {L"R1{a1;i1} = + f{a1;i1}:A-H-S"})});
    CostReportStep step;
    REQUIRE_THROWS_WITH(step.run("cost", ctx, {"tree.0"}),
                        Catch::Matchers::ContainsSubstring("tree.0"));
  }
}
