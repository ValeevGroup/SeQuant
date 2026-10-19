#include <SeQuant/core/eval/schedule_document.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <memory>
#include <string>

namespace {

using namespace sequant::eval::schedule;

Document make_doc() {
  Document d;
  d.fingerprint.roots = {"R1"};
  d.fingerprint.mode_order = {"aux", "occ"};
  TermRecord t;
  t.key = "K";
  t.factors = {"g{i_1;a_1}", "t{a_1;i_1}"};
  t.scalar = "1/2";
  TreeShape tree;
  tree.left = std::make_unique<TreeShape>(TreeShape{0});
  tree.right = std::make_unique<TreeShape>(TreeShape{1});
  t.tree = std::move(tree);
  t.nodes = {
      {{{"i_1", ModeKind::External}}, {{"i_1", ModeKind::External}}, false}};
  d.roots = {{"R1", {std::move(t)}}};
  return d;
}

// replaces the first occurrence of `from` in the serialized document
std::string mutate(Document const& d, std::string const& from,
                   std::string const& to) {
  auto text = to_json(d);
  auto const pos = text.find(from);
  REQUIRE(pos != std::string::npos);
  text.replace(pos, from.size(), to);
  return text;
}

}  // namespace

TEST_CASE("schedule document round-trips through JSON", "[schedule-document]") {
  auto const d = make_doc();
  auto const text = to_json(d);
  REQUIRE(text.find("\"format\": \"sequant.schedule/1\"") != std::string::npos);
  auto const back = from_json(text);
  REQUIRE(compare(d, back).empty());
  REQUIRE(back.roots[0].terms[0].nodes[0].axes[0].label == "i_1");
  REQUIRE(summary_line(d) == "roots=1 terms=1 values=- blocks=- depth=-");
}

TEST_CASE("schedule document rejects bad shapes", "[schedule-document]") {
  using Catch::Matchers::ContainsSubstring;
  auto const d = make_doc();

  // (a) ordinal out of range
  {
    Document e = make_doc();
    e.roots[0].terms[0].tree->right->leaf = 7;
    REQUIRE_THROWS_WITH(from_json(to_json(e)),
                        ContainsSubstring("K") && ContainsSubstring("ordinal"));
  }

  // (b) two inner nodes but one `nodes` entry
  {
    Document e = make_doc();
    auto& t = e.roots[0].terms[0];
    t.factors.push_back("f{a_1;i_1}");
    TreeShape inner;
    inner.left = std::move(t.tree->left);
    inner.right = std::move(t.tree->right);
    TreeShape top;
    top.left = std::make_unique<TreeShape>(std::move(inner));
    top.right = std::make_unique<TreeShape>(TreeShape{2});
    t.tree = std::move(top);
    REQUIRE_THROWS_WITH(from_json(to_json(e)), ContainsSubstring("nodes"));
  }

  // (c) axes label absent from factors
  {
    Document e = make_doc();
    e.roots[0].terms[0].nodes[0].axes[0].label = "x_9";
    REQUIRE_THROWS_WITH(from_json(to_json(e)), ContainsSubstring("x_9"));
  }

  // (d) wrong format
  REQUIRE_THROWS_WITH(
      from_json(mutate(d, "sequant.schedule/1", "sequant.schedule/0")),
      ContainsSubstring("format"));
}

TEST_CASE("schedule document compare reports one difference per changed field",
          "[schedule-document]") {
  auto const d = make_doc();
  auto d2 = d;
  d2.roots[0].terms[0].nodes[0].order_aware = true;
  auto diffs = compare(d, d2);
  REQUIRE(diffs.size() == 1);
  REQUIRE(diffs[0].path == "roots[R1].terms[K].nodes[0].order_aware");

  // a run section on one side only is not a difference
  auto d3 = d;
  d3.run = RunRecord{};
  d3.producer_version = "other";
  REQUIRE(compare(d, d3).empty());
}

namespace {

Document make_dag_doc() {
  Document d = make_doc();
  DagRecord dag;
  for (std::string id : {"v0", "v1"}) {
    ValueRecord v;
    v.id = id;
    v.expr = "g*t";
    v.hash = "h" + id;
    v.key = "k" + id;
    v.leaf = id == "v0";
    v.result_modes = {"i_1"};
    v.home_modes = {"i_1"};
    v.enclosing_modes = {"a_1"};
    v.divergent_modes = {"a_1"};
    v.home_block = 1;
    v.home_depth = 1;
    v.persistent = true;
    v.first_use = 3;
    v.last_use = 5;
    v.build_site = {"occ"};
    v.home_floor = {"occ"};
    v.forced_split_axes = {"aux"};
    v.per_axis = {{"i_1", "LoopLocal"}, {"a_1", "Reduction"}};
    OccurrenceRecord o;
    o.term = "K";
    o.node = 0;
    o.consumer = 1;
    o.home = {"i_1"};
    o.enclosing = {"a_1"};
    o.opens = {{"i_1", ModeKind::External}};
    v.occurrences = {o};
    v.siblings = {id == "v0" ? "v1" : "v0"};
    dag.values.push_back(std::move(v));
  }
  BlockRecord root;
  root.id = 0;
  root.kind = "root";
  root.steps = {{std::string("v0"), std::nullopt}, {std::nullopt, 1}};
  BlockRecord inner;
  inner.id = 1;
  inner.parent = 0;
  inner.axis = "i_1";
  inner.kind = "external";
  inner.depth = 1;
  inner.slot = 2;
  inner.altitude = 1;
  inner.pass = 1;
  inner.steps = {{std::string("v1"), std::nullopt}};
  inner.outputs = {{"v1", "AccumulateSum"}};
  dag.blocks = {root, inner};
  LoopRecord l0{"occ", 0, "external"}, l1{"aux", 1, "contracted"};
  dag.loops = {l0, l1};
  dag.loop_order = {{l0, l1, 2}};
  dag.edges = {{"v0", "v1", true}};
  d.dag = std::move(dag);
  RunRecord r;
  r.extents = {{"occ", 10}};
  r.batch_targets = {{"occ", 5}};
  r.n_batches = {{1, 2}};
  r.footprint = {{"v0", 100}};
  r.builds = {{"v1", 4}};
  r.planner = {{"peak", "1e9"}};
  d.run = std::move(r);
  return d;
}

}  // namespace

TEST_CASE("schedule document round-trips the dag and run sections",
          "[schedule-document]") {
  auto const d = make_dag_doc();
  auto const back = from_json(to_json(d));
  REQUIRE(compare(d, back).empty());
  REQUIRE(back.run.has_value());
  REQUIRE(back.run->n_batches.at(1) == 2);
  REQUIRE(back.run->planner.at("peak") == "1e9");
  REQUIRE(back.dag->values[0].occurrences[0].consumer == 1);
  REQUIRE(back.dag->blocks[1].outputs[0].second == "AccumulateSum");
  REQUIRE(std::get<2>(back.dag->loop_order[0]) == 2);
  REQUIRE(summary_line(d) == "roots=1 terms=1 values=2 blocks=2 depth=1");
}

TEST_CASE("schedule document compare covers the dag", "[schedule-document]") {
  auto const d = make_dag_doc();
  auto d2 = d;
  d2.dag->values[1].home_block = 0;
  auto diffs = compare(d, d2);
  REQUIRE(diffs.size() == 1);
  REQUIRE(diffs[0].path == "dag.values[v1].home_block");

  auto d3 = d;
  d3.dag->edges.push_back({"v1", "v0", false});
  diffs = compare(d, d3);
  REQUIRE(diffs.size() == 1);
  REQUIRE(diffs[0].path.rfind("dag.edges", 0) == 0);

  // ignored value fields
  auto d4 = d;
  d4.dag->values[0].hash = "other";
  d4.dag->values[0].last_use = 99;
  REQUIRE(compare(d, d4).empty());
}
