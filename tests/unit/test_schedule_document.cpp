#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/schedule_pin.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

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

// ---- term key, term record, pinned rebuild ---------------------------------

namespace {

sequant::ExprPtr pin_parse(std::wstring_view s) {
  return sequant::deserialize(s, {.def_perm_symm = sequant::Symmetry::Nonsymm});
}

/// extents: occupied 10, everything else 100
std::size_t pin_extent(sequant::Index const& ix) {
  return ix.space().base_key() == L"i" ? 10 : 100;
}

/// batched optimize of a canonical term; returns the optimized expression and
/// its annotations (captured via term_batch_axes)
std::pair<sequant::ExprPtr,
          sequant::container::vector<sequant::NodeBatchAnnotation>>
pin_optimize(sequant::ExprPtr const& canon) {
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"a";
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  policy.batch_target_size = [](sequant::Index const&) -> std::size_t {
    return 10;
  };
  policy.peak_threshold = 1e5;
  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = pin_extent;
  opts.batch_policy = policy;
  opts.term_batch_axes = axes_map;
  auto optimized = sequant::optimize(canon, opts);
  sequant::container::vector<sequant::NodeBatchAnnotation> annotations;
  if (auto it = axes_map->find(optimized.get()); it != axes_map->end())
    annotations = it->second;
  return {optimized, annotations};
}

}  // namespace

TEST_CASE("term key is stable under relabeling and factor order",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto a = pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_3,i_4}");
  auto b = pin_parse(L"t{a_3,a_4;i_3,i_4} * g{i_1,i_2;a_3,a_4}");
  auto const a_text = sequant::serialize(a);
  REQUIRE(term_key(canonical_term(a)) == term_key(canonical_term(b)));
  // canonical_term works on a clone
  REQUIRE(sequant::serialize(a) == a_text);
  // different externals order
  auto c = pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_4,i_3}");
  REQUIRE(term_key(canonical_term(a)) != term_key(canonical_term(c)));
  // factors of a canonical Product in order
  REQUIRE(canonical_factors(canonical_term(a)).size() == 2);
}

TEST_CASE("record_term then rebuild_pinned_term reproduces the optimized tree",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto canon = canonical_term(
      pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1;i_3} * t{a_2;i_4}"));
  auto [optimized, annotations] = pin_optimize(canon);
  REQUIRE(annotations.size() == 2);
  // the policy forces slicing somewhere: the record carries real annotations
  bool any_axes = false;
  for (auto const& n : annotations) any_axes = any_axes || !n.axes.empty();
  REQUIRE(any_axes);

  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.key == term_key(canon));
  REQUIRE(rec.factors.size() == 3);
  REQUIRE(rec.scalar == "1");
  REQUIRE(rec.tree.has_value());
  REQUIRE(!rec.tree->leaf);
  REQUIRE(rec.nodes.size() == 2);

  auto pinned = rebuild_pinned_term(rec, canon);
  REQUIRE(*pinned.expr == *optimized);  // same nested Product, same labels
  REQUIRE(pinned.nodes.size() == annotations.size());
  for (std::size_t k = 0; k < annotations.size(); ++k) {
    REQUIRE(pinned.nodes[k].axes == annotations[k].axes);
    REQUIRE(pinned.nodes[k].opened_here == annotations[k].opened_here);
    REQUIRE(pinned.nodes[k].order_aware == annotations[k].order_aware);
  }

  // a record survives re-canonicalizing a relabeled copy of the term
  auto canon2 = canonical_term(
      pin_parse(L"t{a_4;i_4} * g{i_1,i_2;a_3,a_4} * t{a_3;i_3}"));
  REQUIRE(term_key(canon2) == rec.key);
  auto pinned2 = rebuild_pinned_term(rec, canon2);
  REQUIRE(*pinned2.expr == *optimized);
}

TEST_CASE(
    "record_term handles a lone tensor, a scalar, identical factors, and "
    "rejects mixed products",
    "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;

  SECTION("lone tensor") {
    auto canon = canonical_term(pin_parse(L"t{a_1;i_1}"));
    REQUIRE(canonical_factors(canon).size() == 1);
    auto [optimized, annotations] = pin_optimize(canon);
    auto rec = record_term(canon, optimized, annotations);
    REQUIRE(rec.tree.has_value());
    REQUIRE(rec.tree->leaf == std::optional<std::size_t>{0});
    REQUIRE(rec.nodes.empty());
    REQUIRE(rec.scalar == "1");
    auto pinned = rebuild_pinned_term(rec, canon);
    REQUIRE(pinned.expr->is<sequant::Tensor>());
    REQUIRE(*pinned.expr == *canon);
    REQUIRE(pinned.nodes.empty());
  }

  SECTION("scalar") {
    auto canon =
        canonical_term(sequant::ex<sequant::Constant>(sequant::rational{1, 2}));
    REQUIRE(canonical_factors(canon).empty());
    auto rec = record_term(canon, canon, {});
    REQUIRE(!rec.tree.has_value());
    REQUIRE(rec.scalar == "1/2");
    REQUIRE(rec.factors.empty());
    auto pinned = rebuild_pinned_term(rec, canon);
    REQUIRE(*pinned.expr == *canon);
    // a negative scalar round-trips too
    auto neg = canonical_term(
        sequant::ex<sequant::Constant>(sequant::rational{-3, 4}));
    REQUIRE(*rebuild_pinned_term(record_term(neg, neg, {}), neg).expr == *neg);
  }

  SECTION("identical factors and a fractional scalar") {
    auto canon = canonical_term(
        pin_parse(L"1/2 * t{a_1;i_1} * t{a_2;i_2} * f{i_1,i_2;a_1,a_2}"));
    REQUIRE(canonical_factors(canon).size() == 3);
    auto [optimized, annotations] = pin_optimize(canon);
    auto rec = record_term(canon, optimized, annotations);
    REQUIRE(rec.scalar == "1/2");
    REQUIRE(rec.nodes.size() == 2);
    auto pinned = rebuild_pinned_term(rec, canon);
    REQUIRE(*pinned.expr == *optimized);
    REQUIRE(pinned.nodes.size() == annotations.size());
    for (std::size_t k = 0; k < annotations.size(); ++k)
      REQUIRE(pinned.nodes[k].axes == annotations[k].axes);
  }

  SECTION("mixed product") {
    auto canon =
        canonical_term(pin_parse(L"t{a_1;i_1} * (f{i_1;a_1} + g{i_1;a_1})"));
    REQUIRE(canon->is<sequant::Product>());
    REQUIRE_THROWS_WITH(record_term(canon, canon, {}),
                        ContainsSubstring("only tensor products"));
    sequant::eval::schedule::TermRecord rec;
    rec.key = term_key(canon);
    rec.scalar = "1";
    REQUIRE_THROWS_WITH(rebuild_pinned_term(rec, canon),
                        ContainsSubstring("only tensor products"));
  }
}

TEST_CASE("rebuild_pinned_term rejects foreign factors and unknown axis labels",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;
  auto canon_a =
      canonical_term(pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_3,i_4}"));
  auto canon_b = canonical_term(pin_parse(L"f{i_1;a_1} * t{a_1;i_2}"));
  auto [optimized, annotations] = pin_optimize(canon_a);
  auto rec = record_term(canon_a, optimized, annotations);
  REQUIRE(rec.nodes.size() == 1);

  REQUIRE_THROWS_WITH(rebuild_pinned_term(rec, canon_b),
                      ContainsSubstring(rec.key));

  auto bad = rec;
  bad.nodes[0].axes.push_back(
      {"z_9", sequant::eval::schedule::ModeKind::Contracted});
  REQUIRE_THROWS_WITH(rebuild_pinned_term(bad, canon_a),
                      ContainsSubstring("z_9"));
}

TEST_CASE(
    "record_term gives every inner node of a pure tensor tree one annotation",
    "[schedule-document]") {
  using namespace sequant::opt::schedule;
  // i_2 is contracted between X and Y only: the (X Y) subtree is fully
  // contracted unless binarize is told i_2 is external; the optimizer still
  // emits one annotation per inner node
  auto canon = canonical_term(pin_parse(L"X{;;i_2} * Y{;;i_2} * Z{a_1;i_1}"));
  REQUIRE(canonical_factors(canon).size() == 3);
  auto [optimized, annotations] = pin_optimize(canon);
  REQUIRE(annotations.size() == 2);
  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.nodes.size() == 2);
  auto pinned = rebuild_pinned_term(rec, canon);
  REQUIRE(*pinned.expr == *optimized);
  REQUIRE(pinned.nodes.size() == 2);
  for (std::size_t k = 0; k < annotations.size(); ++k)
    REQUIRE(pinned.nodes[k].axes == annotations[k].axes);
}

TEST_CASE("record_term accepts an empty annotation vector",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto canon = canonical_term(
      pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1;i_3} * t{a_2;i_4}"));
  auto axes_map = std::make_shared<std::unordered_map<
      sequant::Expr const*,
      sequant::container::vector<sequant::NodeBatchAnnotation>>>();
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseFLOPs;
  opts.idx_to_extent = pin_extent;
  opts.term_batch_axes = axes_map;
  auto optimized = sequant::optimize(canon, opts);
  sequant::container::vector<sequant::NodeBatchAnnotation> annotations;
  if (auto it = axes_map->find(optimized.get()); it != axes_map->end())
    annotations = it->second;
  REQUIRE(annotations.empty());
  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.tree.has_value());
  REQUIRE(!rec.tree->leaf);
  REQUIRE(rec.nodes.empty());
  auto pinned = rebuild_pinned_term(rec, canon);
  REQUIRE(*pinned.expr == *optimized);
  REQUIRE(pinned.nodes.empty());
}

TEST_CASE(
    "schedule pin rejects variables, bad trees and miscounted annotations",
    "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;

  SECTION("a Variable factor") {
    auto canon =
        canonical_term(sequant::ex<sequant::Product>(sequant::ExprPtrList{
            sequant::ex<sequant::Variable>(L"x"), pin_parse(L"f{i_1;a_1}"),
            pin_parse(L"t{a_1;i_2}")}));
    REQUIRE(canon->is<sequant::Product>());
    REQUIRE_THROWS_WITH(record_term(canon, canon, {}),
                        ContainsSubstring("only tensor products"));
    sequant::eval::schedule::TermRecord rec;
    rec.key = term_key(canon);
    rec.scalar = "1";
    REQUIRE_THROWS_WITH(rebuild_pinned_term(rec, canon),
                        ContainsSubstring("only tensor products"));
  }

  auto canon = canonical_term(
      pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1;i_3} * t{a_2;i_4}"));
  auto [optimized, annotations] = pin_optimize(canon);
  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.nodes.size() == 2);

  SECTION("a tree that uses an ordinal twice") {
    auto bad = rec;
    // replace the first leaf found by the ordinal of another leaf
    std::vector<sequant::eval::schedule::TreeShape*> leaves;
    auto collect = [&](auto&& self,
                       sequant::eval::schedule::TreeShape& t) -> void {
      if (t.leaf) {
        leaves.push_back(&t);
        return;
      }
      self(self, *t.left);
      self(self, *t.right);
    };
    collect(collect, *bad.tree);
    REQUIRE(leaves.size() == 3);
    leaves[0]->leaf = leaves[1]->leaf;
    REQUIRE_THROWS_WITH(rebuild_pinned_term(bad, canon),
                        ContainsSubstring(rec.key) &&
                            ContainsSubstring("out of range or repeated"));
  }

  SECTION("an annotation count that does not match the tree") {
    auto extra = annotations;
    extra.emplace_back();
    REQUIRE_THROWS_WITH(record_term(canon, optimized, extra),
                        ContainsSubstring("3 annotations supplied") &&
                            ContainsSubstring("2 tensor-times-tensor nodes"));
    auto bad = rec;
    bad.nodes.pop_back();
    REQUIRE_THROWS_WITH(rebuild_pinned_term(bad, canon),
                        ContainsSubstring(rec.key) &&
                            ContainsSubstring("1 annotations supplied"));
  }
}
