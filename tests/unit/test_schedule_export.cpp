#include <catch2/catch_test_macros.hpp>

#include "test_export.hpp"

#include <SeQuant/core/eval/cell_table.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/export/memory_model.hpp>
#include <SeQuant/core/export/schedule_export.hpp>
#include <SeQuant/core/export/text_generator.hpp>
#include <SeQuant/core/export/value_resolver.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <string>
#include <unordered_map>
#include <vector>

using namespace sequant;

namespace {

/// A trivial ValueResolver backed by a fixed map, for exercising
/// ScheduleWalkingGenerationVisitor without a real RichSchedule/
/// ValueNodeMap pipeline (that bridge is assembled in migration step 5).
class MapValueResolver final : public ValueResolver {
 public:
  void add(std::size_t value_id, EvalExpr node) {
    m_nodes.emplace(value_id, std::move(node));
  }

  const EvalExpr &node_of(std::size_t value_id) const override {
    return m_nodes.at(value_id);
  }

 private:
  std::unordered_map<std::size_t, EvalExpr> m_nodes;
};

/// First word of each non-empty line of TextGenerator's output, i.e. the
/// sequence of semantic verbs (Create/Load/Compute/Unload/Delete/Persist/...)
/// -- checking this rather than the full text keeps the test independent of
/// exact index-label formatting.
std::vector<std::string> verbs(const std::string &generated) {
  std::vector<std::string> out;
  std::size_t pos = 0;
  while (pos < generated.size()) {
    std::size_t nl = generated.find('\n', pos);
    std::string line = generated.substr(pos, nl - pos);
    if (!line.empty()) {
      std::size_t sp = line.find(' ');
      out.push_back(line.substr(0, sp));
    }
    if (nl == std::string::npos) break;
    pos = nl + 1;
  }
  return out;
}

}  // namespace

TEST_CASE("ScheduleWalkingGenerationVisitor", "[export]") {
  auto resetter = to_export_context();

  SECTION("binary product: create, load operands, compute, unload, persist") {
    // R{i1,i2} = A{i1,a1} B{a1,i2}: a genuine binary contraction, matching
    // this first cut's supported scope.
    auto tree =
        binarize(deserialize<ResultExpr>(L"R{i1;i2} = A{i1;a1} B{a1;i2}"));
    REQUIRE(!tree.leaf());

    MapValueResolver resolver;
    resolver.add(0, *tree.left());
    resolver.add(1, *tree.right());
    resolver.add(2, *tree);

    eval::CellTable table;
    eval::TableCell a, b, r;
    a.value_id = 0;
    a.production.kind = eval::ProductionKind::Leaf;
    b.value_id = 1;
    b.production.kind = eval::ProductionKind::Leaf;
    r.value_id = 2;
    r.production.kind = eval::ProductionKind::Build;
    r.life = 0;  // forest root: no consumer.
    table.cells = {a, b, r};
    eval::Read ra, rb;
    ra.consumer = 2;
    ra.source = 0;
    rb.consumer = 2;
    rb.source = 1;
    table.reads = {ra, rb};

    eval::OrderedSchedule schedule;
    schedule.root.steps.push_back(eval::BuildStep{2});
    schedule.num_values = 3;

    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;

    ScheduleWalkingGenerationVisitor<TextGeneratorContext> visitor(
        table, schedule, resolver, gen, ctx);
    visitor.run();

    const std::string code = gen.get_generated_code();
    CHECK(code.find("Create R[") != std::string::npos);
    CHECK(code.find("Load A[") != std::string::npos);
    CHECK(code.find("Load B[") != std::string::npos);
    CHECK(code.find("Compute R[") != std::string::npos);
    CHECK(code.find("+= A[") != std::string::npos);
    CHECK(code.find("Unload A[") != std::string::npos);
    CHECK(code.find("Unload B[") != std::string::npos);
    CHECK(code.find("Persist R[") != std::string::npos);

    const auto seq = verbs(code);
    REQUIRE(seq.size() == 7);
    CHECK(seq.at(0) == "Create");
    CHECK(seq.at(1) == "Load");
    CHECK(seq.at(2) == "Load");
    CHECK(seq.at(3) == "Compute");
    CHECK(seq.at(4) == "Unload");
    CHECK(seq.at(5) == "Unload");
    CHECK(seq.at(6) == "Persist");
  }

  SECTION("shared intermediate: destroy fires only on the draining read") {
    // I is read by both R1 and R2; only the second (draining) read should
    // destroy it, the first should merely leave it resident (no generator
    // call at all between its production and that first read).
    auto i_tree =
        binarize(deserialize<ResultExpr>(L"I{i1;a1} = A{i1;a1} + B{i1;a1}"));
    auto r1_tree =
        binarize(deserialize<ResultExpr>(L"R1{i1;a1} = I{i1;a1} + C{i1;a1}"));
    auto r2_tree =
        binarize(deserialize<ResultExpr>(L"R2{i1;a1} = I{i1;a1} + D{i1;a1}"));

    MapValueResolver resolver;
    resolver.add(0, *i_tree.left());    // A
    resolver.add(1, *i_tree.right());   // B
    resolver.add(2, *i_tree);           // I
    resolver.add(3, *r1_tree.right());  // C
    resolver.add(4, *r1_tree);          // R1
    resolver.add(5, *r2_tree.right());  // D
    resolver.add(6, *r2_tree);          // R2

    eval::CellTable table;
    table.cells.resize(7);
    auto leaf = [](std::size_t vid) {
      eval::TableCell c;
      c.value_id = vid;
      c.production.kind = eval::ProductionKind::Leaf;
      return c;
    };
    table.cells.at(0) = leaf(0);
    table.cells.at(1) = leaf(1);
    table.cells.at(2).value_id = 2;
    table.cells.at(2).production.kind = eval::ProductionKind::Build;
    table.cells.at(2).life = 2;  // read by both R1 and R2
    table.cells.at(3) = leaf(3);
    table.cells.at(4).value_id = 4;
    table.cells.at(4).production.kind = eval::ProductionKind::Build;
    table.cells.at(4).life = 0;
    table.cells.at(5) = leaf(5);
    table.cells.at(6).value_id = 6;
    table.cells.at(6).production.kind = eval::ProductionKind::Build;
    table.cells.at(6).life = 0;

    auto read = [](eval::CellId consumer, eval::CellId source) {
      eval::Read r;
      r.consumer = consumer;
      r.source = source;
      return r;
    };
    table.reads = {read(2, 0), read(2, 1), read(4, 2),
                   read(4, 3), read(6, 2), read(6, 5)};

    eval::OrderedSchedule schedule;
    schedule.root.steps = {eval::BuildStep{2}, eval::BuildStep{4},
                           eval::BuildStep{6}};
    schedule.num_values = 7;

    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    ScheduleWalkingGenerationVisitor<TextGeneratorContext> visitor(
        table, schedule, resolver, gen, ctx);
    visitor.run();

    const std::string code = gen.get_generated_code();
    // I is never unloaded or deleted between its production and R1's use.
    const auto iPos = code.find("Create I[");
    const auto r1ComputePos = code.find("Compute R1[");
    REQUIRE(iPos != std::string::npos);
    REQUIRE(r1ComputePos != std::string::npos);
    const std::string between = code.substr(iPos, r1ComputePos - iPos);
    CHECK(between.find("Unload I[") == std::string::npos);
    CHECK(between.find("Delete I[") == std::string::npos);

    // I is destroyed (not merely unloaded) exactly once, after R2's compute
    // (the draining read), and never unloaded at all.
    CHECK(code.find("Unload I[") == std::string::npos);
    const auto deletePos = code.find("Delete I[");
    const auto r2ComputePos = code.find("Compute R2[");
    REQUIRE(deletePos != std::string::npos);
    REQUIRE(r2ComputePos != std::string::npos);
    CHECK(deletePos > r2ComputePos);
  }

  SECTION("Stack-model generator: a naturally Stack cell is unaffected") {
    // Same fixture as the first section; under Stack, R's only children are
    // leaves and R is a life==0 root, so natural_memory_model() never
    // touches it -- the generated sequence must be identical to a
    // RandomAccess-model generator's.
    auto tree =
        binarize(deserialize<ResultExpr>(L"R{i1;i2} = A{i1;a1} B{a1;i2}"));
    MapValueResolver resolver;
    resolver.add(0, *tree.left());
    resolver.add(1, *tree.right());
    resolver.add(2, *tree);

    eval::CellTable table;
    eval::TableCell a, b, r;
    a.value_id = 0;
    a.production.kind = eval::ProductionKind::Leaf;
    b.value_id = 1;
    b.production.kind = eval::ProductionKind::Leaf;
    r.value_id = 2;
    r.production.kind = eval::ProductionKind::Build;
    r.life = 0;
    table.cells = {a, b, r};
    eval::Read ra, rb;
    ra.consumer = 2;
    ra.source = 0;
    rb.consumer = 2;
    rb.source = 1;
    table.reads = {ra, rb};

    eval::OrderedSchedule schedule;
    schedule.root.steps.push_back(eval::BuildStep{2});

    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen(MemoryModel::Stack);
    ScheduleWalkingGenerationVisitor<TextGeneratorContext> visitor(
        table, schedule, resolver, gen, ctx);
    visitor.run();

    const auto seq = verbs(gen.get_generated_code());
    const std::vector<std::string> expected = {
        "Create", "Load", "Load", "Compute", "Unload", "Unload", "Persist"};
    CHECK(seq == expected);
  }

  SECTION(
      "Stack-model generator: a shared (natural RandomAccess) cell is "
      "persisted once and reloaded independently at each read") {
    auto i_tree =
        binarize(deserialize<ResultExpr>(L"I{i1;a1} = A{i1;a1} + B{i1;a1}"));
    auto r1_tree =
        binarize(deserialize<ResultExpr>(L"R1{i1;a1} = I{i1;a1} + C{i1;a1}"));
    auto r2_tree =
        binarize(deserialize<ResultExpr>(L"R2{i1;a1} = I{i1;a1} + D{i1;a1}"));

    MapValueResolver resolver;
    resolver.add(0, *i_tree.left());
    resolver.add(1, *i_tree.right());
    resolver.add(2, *i_tree);
    resolver.add(3, *r1_tree.right());
    resolver.add(4, *r1_tree);
    resolver.add(5, *r2_tree.right());
    resolver.add(6, *r2_tree);

    eval::CellTable table;
    table.cells.resize(7);
    auto leaf = [](std::size_t vid) {
      eval::TableCell c;
      c.value_id = vid;
      c.production.kind = eval::ProductionKind::Leaf;
      return c;
    };
    table.cells.at(0) = leaf(0);
    table.cells.at(1) = leaf(1);
    table.cells.at(2).value_id = 2;
    table.cells.at(2).production.kind = eval::ProductionKind::Build;
    table.cells.at(2).life = 2;
    table.cells.at(3) = leaf(3);
    table.cells.at(4).value_id = 4;
    table.cells.at(4).production.kind = eval::ProductionKind::Build;
    table.cells.at(4).life = 0;
    table.cells.at(5) = leaf(5);
    table.cells.at(6).value_id = 6;
    table.cells.at(6).production.kind = eval::ProductionKind::Build;
    table.cells.at(6).life = 0;

    auto read = [](eval::CellId consumer, eval::CellId source) {
      eval::Read r;
      r.consumer = consumer;
      r.source = source;
      return r;
    };
    table.reads = {read(2, 0), read(2, 1), read(4, 2),
                   read(4, 3), read(6, 2), read(6, 5)};

    eval::OrderedSchedule schedule;
    schedule.root.steps = {eval::BuildStep{2}, eval::BuildStep{4},
                           eval::BuildStep{6}};

    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen(MemoryModel::Stack);
    ScheduleWalkingGenerationVisitor<TextGeneratorContext> visitor(
        table, schedule, resolver, gen, ctx);
    visitor.run();

    const std::string code = gen.get_generated_code();
    auto count = [&](const std::string &needle) {
      std::size_t n = 0, pos = 0;
      while ((pos = code.find(needle, pos)) != std::string::npos) {
        ++n;
        pos += needle.size();
      }
      return n;
    };

    // Persisted exactly once, immediately after its own production, well
    // before either reader.
    CHECK(count("Persist I[") == 1);
    const auto persistPos = code.find("Persist I[");
    const auto r1CreatePos = code.find("Create R1[");
    REQUIRE(persistPos != std::string::npos);
    REQUIRE(r1CreatePos != std::string::npos);
    CHECK(persistPos < r1CreatePos);

    // Reloaded independently at each of its two reads.
    CHECK(count("Load I[") == 2);
    // Unloaded after its own production and after the first (non-draining)
    // read; destroyed exactly once, on the second (draining) read.
    CHECK(count("Unload I[") == 2);
    CHECK(count("Delete I[") == 1);
    const auto deletePos = code.find("Delete I[");
    const auto r2ComputePos = code.find("Compute R2[");
    REQUIRE(deletePos != std::string::npos);
    REQUIRE(r2ComputePos != std::string::npos);
    CHECK(deletePos > r2ComputePos);
  }

  SECTION("a batched schedule (nested ScopeBlock) is rejected") {
    eval::CellTable table;
    eval::OrderedSchedule schedule;
    schedule.root.steps.push_back(eval::ScopeBlock{});

    MapValueResolver resolver;
    TextGeneratorContext ctx;
    TextGenerator<TextGeneratorContext> gen;
    ScheduleWalkingGenerationVisitor<TextGeneratorContext> visitor(
        table, schedule, resolver, gen, ctx);

    CHECK_THROWS_AS(visitor.run(), Exception);
  }
}
