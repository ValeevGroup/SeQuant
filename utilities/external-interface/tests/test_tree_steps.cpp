#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "cse_step.hpp"
#include "export_step.hpp"
#include "index_batching_step.hpp"

#include "test_extint.hpp"

#include <SeQuant/core/utility/exception.hpp>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace sequant;
using namespace sequant::util::extint;
using json = nlohmann::json;

TEST_CASE("cse step", "[extint]") {
  SECTION("option validation") {
    CSEStep step;
    REQUIRE_THROWS_AS(step.set_options(json::array()), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"min_usage", "two"}}), Exception);
    REQUIRE_THROWS_AS(step.set_options(json{{"merge_inputs", 1}}), Exception);
    REQUIRE_THROWS_AS(
        step.set_options(json{{"merge_inputs", json::array({1})}}), Exception);
    REQUIRE_THROWS_AS(
        step.set_options(json{{"merge_inputs", json::array({""})}}), Exception);
    REQUIRE_NOTHROW(step.set_options(json{{"merge_inputs", true}}));
    REQUIRE_NOTHROW(
        step.set_options(json{{"merge_inputs", json::array({"ecc", "res"})}}));
    REQUIRE_THROWS_AS(step.set_options(json{{"unknown", true}}), Exception);
  }

  SECTION("a tree with no common subexpression is passed through") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::tree_data(
                     {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"})});
    CSEStep step;
    REQUIRE(step.run("cse", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("cse.0"));
    REQUIRE_NOTHROW(
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get()));
  }

  SECTION("an input without entries is passed through") {
    ExecutionContext ctx;
    ctx.set_data("in", 0, ProcessingData{ExportTreeData{}});

    CSEStep step;
    REQUIRE_NOTHROW(step.run("cse", ctx, {"in.0"}));
    REQUIRE(ctx.has_data("cse.0"));
    REQUIRE(
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get())
            .entries.empty());
  }

  SECTION("CSE retains the IDs of the original results") {
    // Both results are the same contraction, so CSE rewrites their roots
    ExportTreeData data = test_support::tree_data(
        {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    REQUIRE(data.entries.size() == 2);
    const std::vector<std::size_t> original_ids = {data.entries[0].tree->id(),
                                                   data.entries[1].tree->id()};

    ExecutionContext ctx;
    ctx.set_data("in", 0, ProcessingData{std::move(data)});

    CSEStep step;
    REQUIRE(step.run("cse", ctx, {"in.0"}) == 1);

    const ExportTreeData &result =
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get());
    REQUIRE(result.entries.size() > 2);

    for (std::size_t id : original_ids) {
      REQUIRE(std::ranges::any_of(result.entries, [id](const auto &entry) {
        return entry.tree->id() == id;
      }));
    }
  }

  SECTION("merge_inputs with group IDs merges only within each group") {
    const std::wstring expr =
        L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S";
    ExecutionContext ctx;
    for (std::size_t i = 0; i < 3; ++i) {
      ctx.set_data("in", i, ProcessingData{test_support::tree_data({expr})});
    }
    ctx.add_data_alias(std::vector<std::string_view>{"in.0", "in.1"}, "in.grp");

    CSEStep step;
    step.set_options(json{{"merge_inputs", json::array({"grp"})}});
    // One merged output for the group plus the individually processed input
    REQUIRE(step.run("cse", ctx, {"in.0", "in.1", "in.2"}) == 2);
    REQUIRE(ctx.has_data("cse.grp"));
    REQUIRE(ctx.dataset_size("cse.grp") == 1);
    REQUIRE(ctx.dataset_size("cse[0-1]") == 2);

    SECTION("overlapping groups are rejected") {
      ctx.add_data_alias("in.1", "in.other");
      CSEStep overlap;
      overlap.set_options(
          json{{"merge_inputs", json::array({"grp", "other"})}});
      REQUIRE_THROWS_AS(overlap.run("cse2", ctx, {"in.0", "in.1", "in.2"}),
                        Exception);
    }
  }

  SECTION("dependencies on the original results survive CSE") {
    // Two identical contractions share a common subexpression, so CSE
    // introduces an intermediate (the result grows by one entry).
    ExportTreeData data = test_support::tree_data(
        {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    REQUIRE(data.entries.size() == 2);

    // Wire a dependency onto each original result. The CSE step copies these
    // IDs through verbatim, so a sentinel that is not a real root id is fine.
    constexpr std::size_t dep_id = 123456;
    for (ExportTreeData::Entry &entry : data.entries) {
      entry.dependencies = {dep_id};
    }

    ExecutionContext ctx;
    ctx.set_data("in", 0, ProcessingData{std::move(data)});

    CSEStep step;
    REQUIRE(step.run("cse", ctx, {"in.0"}) == 1);
    REQUIRE(ctx.has_data("cse.0"));

    const ExportTreeData &result =
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get());

    // CSE must actually have fired, otherwise this section would not exercise
    // the dependency-carrying code path.
    REQUIRE(result.entries.size() > 2);

    std::size_t original_count = 0;
    for (const ExportTreeData::Entry &entry : result.entries) {
      if (entry.overwrite_previous.value_or(false)) {
        // Freshly introduced intermediate: no prior dependencies to carry.
        continue;
      }

      ++original_count;
      REQUIRE(std::ranges::find(entry.dependencies, dep_id) !=
              entry.dependencies.end());
    }

    REQUIRE(original_count == 2);
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
                 ProcessingData{test_support::tree_data(
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

  SECTION("irreps option validation") {
    auto with_irreps = [](json irreps) {
      return json{{"language", "itf"},
                  {"output", "out.itfaa"},
                  {"irreps", std::move(irreps)}};
    };
    auto rejects = [&](json irreps, std::string_view message) {
      REQUIRE_THROWS_WITH(
          ExportStep{}.set_options(with_irreps(std::move(irreps))),
          Catch::Matchers::ContainsSubstring(std::string(message)));
    };

    REQUIRE_NOTHROW(ExportStep{}.set_options(with_irreps(
        {{"group", "D2h"},
         {"symbols", {"Gp", "Gq"}},
         {"relations", {"Gp*Gq = B1u"}},
         {"default", "error"},
         {"tensors",
          {{"f", "Ag"}, {"T2{a1,a2;i1,i2}", {"Gp", "B1u*Gq"}}, {"X", "*"}}}})));
    REQUIRE_NOTHROW(ExportStep{}.set_options(
        with_irreps({{"group", "symbolic"}, {"tensors", {{"t", "Gx*Gy"}}}})));
    REQUIRE_NOTHROW(ExportStep{}.set_options(with_irreps(
        {{"group",
          {{"totally_symmetric", "A1"},
           {"products",
            {{"A2*A2", "A1"}, {"A2*E", "E"}, {"E*E", {"A1", "A2", "E"}}}}}},
         {"tensors", {{"t", "E"}}}})));

    rejects(json::array(), "Expected an object");
    rejects({{"tensors", json::object()}}, "'group' entry is mandatory");
    rejects({{"group", "D3h"}}, "Unknown point group");
    rejects({{"group", "D2h"}, {"unknown", 1}}, "Unknown key");
    // Symbols must be declared for concrete point groups (catches typos)
    rejects({{"group", "D2h"}, {"tensors", {{"t", "B1"}}}}, "Unknown IRREP");
    rejects({{"group", "symbolic"},
             {"symbols", {"Gp"}},
             {"tensors", {{"t", "Gq"}}}},
            "Unknown IRREP");
    rejects({{"group", "D2h"}, {"relations", {"B1u"}}}, "lhs = rhs");
    rejects({{"group", "D2h"}, {"relations", {"B1u = B2u"}}},
            "two different concrete IRREPs");
    rejects({{"group", "D2h"}, {"tensors", {{"t", json::array()}}}},
            "non-empty array");
    rejects({{"group", "D2h"}, {"tensors", {{"t", 5}}}}, "string");
    rejects(
        {{"group", {{"totally_symmetric", "A1"}, {"products", json::object()}}},
         {"tensors", {{"t", "E"}}}},
        "Unknown IRREP 'E'");
    rejects(
        {{"group", {{"totally_symmetric", "A1"}}}, {"relations", {"A1 = A1"}}},
        "only supported");
    rejects({{"group", {{"products", json::object()}}}}, "totally_symmetric");
  }

  SECTION("IRREPs are emitted in the generated code") {
    ExportTreeData data = test_support::tree_data(
        {L"R{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i3}:A-H-S "
         L"Y{i3;i2}:N-N-N"});

    ExecutionContext ctx;
    ctx.set_data("in", 0, ProcessingData{std::move(data)});

    const std::filesystem::path out =
        std::filesystem::temp_directory_path() / "extint_irreps_test.itfaa";

    json options{{"language", "itf"},
                 {"output", out.string()},
                 {"irreps",
                  {{"group", "D2h"},
                   {"symbols", {"Gp"}},
                   {"default", "error"},
                   {"tensors",
                    {{"f", "Ag"}, {"T2", "Gp"}, {"Y{i1;i2}", {"Ag", "B1u"}}}}}},
                 {"meta",
                  {{"index_spaces",
                    {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                     {"a", {{"name", "External"}, {"tag", "e"}}}}},
                   {"irrep_names", {{"Gp", "pert"}}}}}};

    ExportStep step;
    step.set_options(options);
    {
      test_support::ScopedStreamSilencer silencer;
      step.run("exp", ctx, {"in.0"});
    }

    std::ifstream stream(out);
    REQUIRE(stream.good());
    std::stringstream buffer;
    buffer << stream.rdbuf();
    const std::string code = buffer.str();
    stream.close();
    std::filesystem::remove(out);

    // Gp x Ag x {Ag, B1u} = {Gp, B1u*Gp} where Gp is renamed to 'pert'
    REQUIRE_THAT(
        code, Catch::Matchers::ContainsSubstring("R:ec // irrep:pert,B1u*Gp"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("f:ce // irrep:Ag"));
    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring("T2:eecc // irrep:pert"));
    REQUIRE_THAT(
        code, Catch::Matchers::ContainsSubstring("!Create{type:plain; irrep:"));

    // Tensors without IRREP specification are rejected if requested
    options["irreps"]["tensors"].erase("Y{i1;i2}");
    ExecutionContext ctx2;
    ctx2.set_data("in", 0,
                  ProcessingData{test_support::tree_data(
                      {L"R{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i3}:A-H-S "
                       L"Y{i3;i2}:N-N-N"})});
    ExportStep strict;
    strict.set_options(options);
    test_support::ScopedStreamSilencer silencer;
    REQUIRE_THROWS_WITH(strict.run("exp", ctx2, {"in.0"}),
                        Catch::Matchers::ContainsSubstring("No IRREPs"));
    std::filesystem::remove(out);
  }

  SECTION("IRREP names are looked up by canonical IRREP") {
    ExecutionContext ctx;
    ctx.set_data("in", 0,
                 ProcessingData{test_support::tree_data(
                     {L"R{a1;i1} = + T2{a1,a2;i1,i2}:A-H-S "
                      L"Y{i2;a2}:N-N-N"})});

    const std::filesystem::path out =
        std::filesystem::temp_directory_path() / "extint_irrep_names.itfaa";

    // Gp is the totally symmetric IRREP by virtue of the relation
    json options{{"language", "itf"},
                 {"output", out.string()},
                 {"irreps",
                  {{"group", "D2h"},
                   {"symbols", {"Gp"}},
                   {"relations", {"Gp = Ag"}},
                   {"default", "error"},
                   {"tensors", {{"T2", "Gp"}, {"Y", "Ag"}}}}},
                 {"meta",
                  {{"index_spaces",
                    {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                     {"a", {{"name", "External"}, {"tag", "e"}}}}},
                   {"irrep_names", {{"Gp", "sym"}}}}}};

    ExportStep step;
    step.set_options(options);
    {
      test_support::ScopedStreamSilencer silencer;
      step.run("exp", ctx, {"in.0"});
    }

    std::ifstream stream(out);
    REQUIRE(stream.good());
    std::stringstream buffer;
    buffer << stream.rdbuf();
    const std::string code = buffer.str();
    stream.close();
    std::filesystem::remove(out);

    REQUIRE_THAT(code,
                 Catch::Matchers::ContainsSubstring("T2:eecc // irrep:sym"));
    REQUIRE_THAT(code, Catch::Matchers::ContainsSubstring("R:ec // irrep:sym"));
  }

  SECTION("tensors are matched by label patterns and blocks") {
    // @returns The code generated for a fixed expression with the given
    // IRREP specification of the tensors
    auto export_code = [](json tensors) {
      ExecutionContext ctx;
      ctx.set_data("in", 0,
                   ProcessingData{test_support::tree_data(
                       {L"R{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i3}:A-H-S "
                        L"Y{i3;i2}:N-N-N"})});

      const std::filesystem::path out = std::filesystem::temp_directory_path() /
                                        "extint_irreps_matching_test.itfaa";

      ExportStep step;
      step.set_options(
          json{{"language", "itf"},
               {"output", out.string()},
               {"irreps",
                {{"group", "D2h"},
                 {"symbols", {"Gp"}},
                 {"default", "error"},
                 {"tensors", std::move(tensors)}}},
               {"meta",
                {{"index_spaces",
                  {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                   {"a", {{"name", "External"}, {"tag", "e"}}}}}}}});

      test_support::ScopedStreamSilencer silencer;
      step.run("exp", ctx, {"in.0"});

      std::ifstream stream(out);
      std::stringstream buffer;
      buffer << stream.rdbuf();
      stream.close();
      std::filesystem::remove(out);

      return buffer.str();
    };

    using Catch::Matchers::ContainsSubstring;

    // Labels are full-match regular expressions
    std::string code = export_code({{"[fY]", "Ag"}, {"T.", "Gp"}});
    REQUIRE_THAT(code, ContainsSubstring("Y:cc // irrep:Ag"));
    REQUIRE_THAT(code, ContainsSubstring("T2:eecc // irrep:Gp"));
    REQUIRE_THAT(code, ContainsSubstring("R:ec // irrep:Gp"));
    REQUIRE_THROWS_WITH(export_code({{"[fY]", "Ag"}, {"T", "Gp"}}),
                        ContainsSubstring("No IRREPs"));

    // A literal label takes precedence over patterns
    code = export_code({{"[fY]", "Ag"}, {"T.", "Gp"}, {"T2", "B1u"}});
    REQUIRE_THAT(code, ContainsSubstring("T2:eecc // irrep:B1u"));

    // Tensor blocks take precedence over labels
    code = export_code(
        {{"[fY]", "Ag"}, {"T2", "Gp"}, {"Y{i1;i2}", {"Ag", "B1u"}}});
    REQUIRE_THAT(code, ContainsSubstring("Y:cc // irrep:Ag,B1u"));

    // Several patterns may only match if they agree
    REQUIRE_NOTHROW(export_code({{"[fY]", "Ag"}, {"T.", "Gp"}, {".2", "Gp"}}));
    REQUIRE_THROWS_WITH(
        export_code({{"[fY]", "Ag"}, {"T.", "Gp"}, {".2", "B1u"}}),
        ContainsSubstring("Ambiguous"));

    REQUIRE_THROWS_WITH(export_code({{"T[", "Gp"}}),
                        ContainsSubstring("Invalid label regex"));
  }

  SECTION("a dependency cycle is rejected when ordering the results") {
    // Build two independent trees, then wire a mutual dependency between them.
    ExportTreeData two = test_support::tree_data(
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

  SECTION(
      "CSE intermediates reused across inputs are allocated, not accumulated") {
    // Common-subexpression elimination runs once per input (per equation) and
    // restarts its intermediate-name counter each time, so two unrelated
    // intermediates from different inputs can end up sharing a name once the
    // inputs are exported into a single section. Such a collision must not turn
    // the second definition into a load-and-accumulate onto the first (which
    // would silently add the first intermediate's value into the second).
    //
    // Here two independent inputs each contain the *same* common subexpression
    // used by two results, so CSE factors it out in both and -- because the
    // counter restarts -- gives both intermediates the same name. We then
    // export both inputs into one section and require that the shared
    // intermediate is allocated afresh for each definition.
    ExportTreeData data_a = test_support::tree_data(
        {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});
    ExportTreeData data_b = test_support::tree_data(
        {L"R3{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"R4{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});

    ExecutionContext ctx;
    ctx.set_data("A", 0, ProcessingData{std::move(data_a)});
    ctx.set_data("B", 0, ProcessingData{std::move(data_b)});

    CSEStep cse;
    REQUIRE(cse.run("cse", ctx, {"A.0", "B.0"}) == 2);
    REQUIRE(ctx.has_data("cse.0"));
    REQUIRE(ctx.has_data("cse.1"));

    // Each input gained an intermediate definition (CSE fired on both).
    const ExportTreeData &res_a =
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get());
    const ExportTreeData &res_b =
        convert_data<ExportTreeData>(ctx.get_data("cse.1").front().data.get());
    REQUIRE(res_a.entries.size() == 3);
    REQUIRE(res_b.entries.size() == 3);

    const std::filesystem::path out = std::filesystem::temp_directory_path() /
                                      "extint_cse_reuse_alloc_test.itfaa";

    ExportStep step;
    step.set_options(
        json{{"language", "itf"},
             // optimize at default (true), matching the real pipeline
             {"output", out.string()},
             {"meta",
              {{"index_spaces",
                {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                 {"a", {{"name", "External"}, {"tag", "e"}}}}}}}});

    {
      test_support::ScopedStreamSilencer silencer;
      step.run("exp", ctx, {"cse.0", "cse.1"});
    }

    std::ifstream stream(out);
    REQUIRE(stream.good());
    std::stringstream buffer;
    buffer << stream.rdbuf();
    const std::string code = buffer.str();
    stream.close();
    std::filesystem::remove(out);

    // The shared intermediate is defined once per input. Each definition must
    // be a fresh allocation. The bug signature is a definition that first
    // *loads* the stored intermediate and then accumulates onto it, i.e. a
    //     load  CSEx:blk[...]
    //     ...
    //     .CSEx:blk[...] += ...
    // sequence (reusing the previous, unrelated intermediate's value). Count
    // such accumulate-definitions; there must be none.
    std::vector<std::string> lines;
    for (std::size_t pos = 0; pos <= code.size();) {
      std::size_t eol = code.find('\n', pos);
      if (eol == std::string::npos) eol = code.size();
      lines.emplace_back(code.substr(pos, eol - pos));
      pos = eol + 1;
    }

    const auto nameblock = [](std::string_view line,
                              std::string_view prefix) -> std::string {
      if (!line.starts_with(prefix)) return {};
      line.remove_prefix(prefix.size());
      std::size_t bracket = line.find('[');
      if (bracket == std::string_view::npos) return {};
      return std::string(line.substr(0, bracket));
    };

    std::size_t accumulate_definitions = 0;
    std::size_t alloc_definitions = 0;
    for (std::size_t i = 0; i < lines.size(); ++i) {
      // A CSE result definition line: ".CSEx:blk[...] += ..."
      const std::string def = nameblock(lines[i], ".CSE");
      if (def.empty()) continue;
      alloc_definitions += 1;  // provisionally; corrected below if it's a load

      // Look back a few lines for how this intermediate was brought into scope.
      bool loaded = false;
      bool allocated = false;
      for (std::size_t j = i; j-- > 0 && i - j <= 6;) {
        if (nameblock(lines[j], "load CSE") == def) {
          loaded = true;
          break;
        }
        if (nameblock(lines[j], "alloc CSE") == def) {
          allocated = true;
          break;
        }
      }
      if (loaded && !allocated) ++accumulate_definitions;
    }

    REQUIRE(alloc_definitions >= 2);  // CSE fired on both inputs
    CHECK(accumulate_definitions == 0);
  }

  SECTION("a CSE intermediate that is a sum is not zeroed between summands") {
    // Both results contain the same two-term sum, so CSE factors it out as an
    // intermediate whose summands all accumulate into that one tensor.
    ExportTreeData data = test_support::tree_data(
        {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"
         L" + g{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"
         L" + h{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
         L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"
         L" + g{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"
         L" + k{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});

    ExecutionContext ctx;
    ctx.set_data("in", 0, ProcessingData{std::move(data)});

    CSEStep cse;
    REQUIRE(cse.run("cse", ctx, {"in.0"}) == 1);

    const ExportTreeData &res =
        convert_data<ExportTreeData>(ctx.get_data("cse.0").front().data.get());
    REQUIRE(res.entries.size() > 2);  // CSE fired

    const std::filesystem::path out = std::filesystem::temp_directory_path() /
                                      "extint_cse_sum_zeroing_test.itfaa";

    ExportStep step;
    step.set_options(json{{"language", "itf"},
                          {"output", out.string()},
                          {"meta",
                           {{"index_spaces",
                             {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                              {"a", {{"name", "External"}, {"tag", "e"}}}}}}}});

    {
      test_support::ScopedStreamSilencer silencer;
      step.run("exp", ctx, {"cse.0"});
    }

    std::ifstream stream(out);
    REQUIRE(stream.good());
    std::stringstream buffer;
    buffer << stream.rdbuf();
    const std::string code = buffer.str();
    stream.close();
    std::filesystem::remove(out);

    // The zeroing of a freshly created intermediate is the allocation's; an
    // explicit reset in the middle of the accumulation discards the summands
    // added so far.
    CAPTURE(code);
    CHECK(code.find("+=") != std::string::npos);
    CHECK(code.find(":= 0") == std::string::npos);
  }

  SECTION("names of CSE intermediates are recycled on request") {
    // Two common subexpressions with compatible shapes whose lifetimes don't
    // overlap: the second can reuse the name (and storage) of the first
    auto generate = [](const json &recycle) {
      ExportTreeData data = test_support::tree_data(
          {L"R1{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
           L"R2{a1;i1} = + f{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
           L"R3{a1;i1} = + g{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S",
           L"R4{a1;i1} = + g{i2;a2}:A-H-S T2{a1,a2;i1,i2}:A-H-S"});

      ExecutionContext ctx;
      ctx.set_data("in", 0, ProcessingData{std::move(data)});

      CSEStep cse;
      REQUIRE(cse.run("cse", ctx, {"in.0"}) == 1);

      const std::filesystem::path out = std::filesystem::temp_directory_path() /
                                        "extint_cse_name_recycling_test.itfaa";

      json options{{"language", "itf"},
                   {"output", out.string()},
                   {"meta",
                    {{"index_spaces",
                      {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                       {"a", {{"name", "External"}, {"tag", "e"}}}}}}}};
      if (!recycle.is_null()) {
        options["recycle_intermediates"] = recycle;
      }

      ExportStep step;
      step.set_options(options);

      {
        test_support::ScopedStreamSilencer silencer;
        step.run("exp", ctx, {"cse.0"});
      }

      std::ifstream stream(out);
      REQUIRE(stream.good());
      std::stringstream buffer;
      buffer << stream.rdbuf();
      stream.close();
      std::filesystem::remove(out);

      return buffer.str();
    };

    const std::string plain = generate(nullptr);
    CAPTURE(plain);
    CHECK(plain.find("CSE1") != std::string::npos);
    CHECK(plain.find("CSE2") != std::string::npos);
    CHECK(generate(false) == plain);

    for (const json &recycle :
         {json(true), json{{"shape", "block"}}, json{{"shape", "multiset"}}}) {
      const std::string recycled = generate(recycle);
      CAPTURE(recycle.dump(), recycled);
      CHECK(recycled.find("CSE1") != std::string::npos);
      CHECK(recycled.find("CSE2") == std::string::npos);
      // Every definition of the intermediate is a fresh allocation
      std::size_t num_allocs = 0;
      for (std::size_t pos = recycled.find("alloc CSE1");
           pos != std::string::npos;
           pos = recycled.find("alloc CSE1", pos + 1)) {
        ++num_allocs;
      }
      CHECK(num_allocs == 2);
    }

    ExportStep step;
    CHECK_THROWS(step.set_options(json{{"language", "itf"},
                                       {"output", "out.itfaa"},
                                       {"recycle_intermediates", 1}}));
    CHECK_THROWS(step.set_options(
        json{{"language", "itf"},
             {"output", "out.itfaa"},
             {"recycle_intermediates", {{"shape", "triangle"}}}}));
  }

  SECTION("recycling with IRREPs uses results of earlier groups") {
    // P is computed in its own group and used by the CSE intermediates of a
    // later group. Its IRREPs are not specified but have to be determined from
    // its definition.
    auto generate = [](bool recycle) {
      ExecutionContext ctx;
      ctx.set_data("P", 1,
                   ProcessingData{test_support::tree_data(
                       {L"P{a1;i1} = + f{a1;a2}:N-N-N t{a2;i1}:N-N-N"})});
      ctx.set_data("in", 0,
                   ProcessingData{test_support::tree_data(
                       {L"R1{a1;i1} = + P{a1;i2}:N-N-N g{i2;i1}:N-N-N",
                        L"R2{a1;i1} = + P{a1;i2}:N-N-N g{i2;i1}:N-N-N",
                        L"R3{a1;i1} = + P{a1;i2}:N-N-N h{i2;i1}:N-N-N",
                        L"R4{a1;i1} = + P{a1;i2}:N-N-N h{i2;i1}:N-N-N"})});

      CSEStep cse;
      REQUIRE(cse.run("cse", ctx, {"in.0"}) == 1);

      const std::filesystem::path out = std::filesystem::temp_directory_path() /
                                        "extint_recycling_irreps_test.itfaa";

      json options{{"language", "itf"},
                   {"output", out.string()},
                   // Note: groups are matched by the IDs without step prefix
                   {"grouping", {{"Init", "1"}, {"Residual", "0"}}},
                   {"irreps",
                    {{"group", "D2h"},
                     {"default", "error"},
                     {"tensors",
                      {{"f", "Ag"}, {"t", "B1u"}, {"g", "Ag"}, {"h", "Ag"}}}}},
                   {"recycle_intermediates", recycle},
                   {"meta",
                    {{"index_spaces",
                      {{"i", {{"name", "Closed"}, {"tag", "c"}}},
                       {"a", {{"name", "External"}, {"tag", "e"}}}}}}}};

      ExportStep step;
      step.set_options(options);
      {
        test_support::ScopedStreamSilencer silencer;
        step.run("exp", ctx, {"P.1", "cse.0"});
      }

      std::ifstream stream(out);
      REQUIRE(stream.good());
      std::stringstream buffer;
      buffer << stream.rdbuf();
      stream.close();
      std::filesystem::remove(out);

      return buffer.str();
    };

    const std::string plain = generate(false);
    CAPTURE(plain);
    REQUIRE(plain.find("CSE2") != std::string::npos);

    const std::string recycled = generate(true);
    CAPTURE(recycled);
    CHECK(recycled.find("CSE1") != std::string::npos);
    CHECK(recycled.find("CSE2") == std::string::npos);
  }
}
