#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/schedule_pin.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

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

// ---- pin and record sessions through optimize_result ------------------------

namespace {

using AxesMap = std::unordered_map<
    sequant::Expr const*,
    sequant::container::vector<sequant::NodeBatchAnnotation>>;

/// batched options as pin_optimize uses them, with a fresh term_batch_axes map
sequant::OptimizeOptions session_opts() {
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
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = pin_extent;
  opts.batch_policy = policy;
  opts.term_batch_axes = std::make_shared<AxesMap>();
  return opts;
}

/// a Sum of separately parsed summands (no simplification across them)
sequant::ExprPtr session_sum(std::vector<std::wstring_view> const& terms) {
  std::vector<sequant::ExprPtr> summands;
  for (auto t : terms) summands.push_back(pin_parse(t));
  return sequant::ex<sequant::Sum>(summands);
}

/// the 4-summand root used by the round-trip tests (the last one a lone tensor)
sequant::ExprPtr session_root() {
  return session_sum({L"g{i_1,i_2;a_3,a_4} * t{a_3;i_3} * t{a_4;i_4}",
                      L"1/2 * g{i_1,i_2;a_3,a_4} * t{a_3,a_4;i_3,i_4}",
                      L"f{i_1;a_3} * t{a_3,a_4;i_3,i_4} * f{i_2;a_4}",
                      L"X{i_1,i_2;i_3,i_4}"});
}

void require_same_annotations(
    sequant::container::vector<sequant::NodeBatchAnnotation> const& a,
    sequant::container::vector<sequant::NodeBatchAnnotation> const& b) {
  REQUIRE(a.size() == b.size());
  for (std::size_t k = 0; k < a.size(); ++k) {
    REQUIRE(a[k].axes == b[k].axes);
    REQUIRE(a[k].opened_here == b[k].opened_here);
    REQUIRE(a[k].order_aware == b[k].order_aware);
  }
}

std::vector<std::string> record_keys(RootRecord const& r) {
  std::vector<std::string> keys;
  for (auto const& t : r.terms) keys.push_back(t.key);
  return keys;
}

}  // namespace

TEST_CASE("record_root then pin reproduces the optimized Sum and annotations",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto const sum = session_root();
  auto opts = session_opts();
  auto rr = record_root("R", sum, opts);
  REQUIRE(rr.record.label == "R");
  REQUIRE(rr.record.terms.size() == 4);
  REQUIRE(rr.optimized->is<sequant::Sum>());
  REQUIRE(rr.optimized->as<sequant::Sum>().size() == 4);
  auto const& recorded_axes = (*opts.term_batch_axes)[rr.optimized.get()];
  REQUIRE(!recorded_axes.empty());

  PinSession session(rr.record);
  auto opts2 = opts;
  opts2.pinned = &session;
  opts2.term_batch_axes = std::make_shared<AxesMap>();
  auto pinned = sequant::optimize_result(sum, opts2).expr;
  REQUIRE(*pinned == *rr.optimized);  // same trees, same order, same labels
  require_same_annotations((*opts2.term_batch_axes)[pinned.get()],
                           recorded_axes);
  REQUIRE_NOTHROW(session.finish());
}

TEST_CASE("pinning is strict", "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;
  auto const sum = session_root();
  auto opts = session_opts();
  auto rr = record_root("R", sum, opts);
  REQUIRE(rr.record.terms.size() == 4);

  SECTION("a missing record") {
    auto doc = rr.record;
    auto const dropped = doc.terms[1].key;
    doc.terms.erase(doc.terms.begin() + 1);
    PinSession session(doc);
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(
        sequant::optimize_result(sum, opts2),
        ContainsSubstring("missing") && ContainsSubstring(dropped));
  }

  SECTION("an unused record") {
    auto other = record_root(
        "O", session_sum({L"h{i_1;a_3} * t{a_3;i_2}", L"Y{i_1;i_2}"}), opts);
    auto doc = rr.record;
    doc.terms.push_back(other.record.terms[0]);
    PinSession session(doc);
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(sequant::optimize_result(sum, opts2),
                        ContainsSubstring("unused") &&
                            ContainsSubstring(other.record.terms[0].key));
  }

  SECTION("a duplicate term") {
    auto twice =
        session_sum({L"g{i_1;a_3} * t{a_3;i_2}", L"g{i_1;a_4} * t{a_4;i_2}"});
    REQUIRE(twice->as<sequant::Sum>().size() == 2);
    auto once = record_root("D", session_sum({L"g{i_1;a_3} * t{a_3;i_2}"}),
                            session_opts());
    PinSession session(once.record);
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(sequant::optimize_result(twice, opts2),
                        ContainsSubstring("duplicate") &&
                            ContainsSubstring(once.record.terms[0].key));
  }

  SECTION("a mixed product") {
    auto mixed = sequant::ex<sequant::Sum>(sequant::ExprPtrList{
        pin_parse(L"t{a_1;i_1} * (f{i_1;a_1} + g{i_1;a_1})"),
        pin_parse(L"h{i_1;a_1}")});
    PinSession session(RootRecord{"M", {}});
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(sequant::optimize_result(mixed, opts2),
                        ContainsSubstring("only tensor products"));
    auto opts3 = opts;
    RecordSession recording("M");
    opts3.recording = &recording;
    REQUIRE_THROWS_WITH(sequant::optimize_result(mixed, opts3),
                        ContainsSubstring("only tensor products"));
  }

  SECTION("several failures are listed in one exception") {
    auto doc = rr.record;
    auto const dropped = doc.terms[0].key;
    doc.terms.erase(doc.terms.begin());
    auto other = record_root("O", session_sum({L"h{i_1;a_3} * t{a_3;i_2}"}),
                             session_opts());
    doc.terms.push_back(other.record.terms[0]);
    PinSession session(doc);
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(sequant::optimize_result(sum, opts2),
                        ContainsSubstring("missing") &&
                            ContainsSubstring(dropped) &&
                            ContainsSubstring("unused") &&
                            ContainsSubstring(other.record.terms[0].key));
  }
}

TEST_CASE("per-summand pinning through one session", "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;
  auto const sum = session_root();
  auto opts = session_opts();
  auto rr = record_root("R", sum, opts);
  auto const& whole = rr.optimized->as<sequant::Sum>();

  // the whole-Sum pin, for comparison
  PinSession whole_session(rr.record);
  auto opts_whole = opts;
  opts_whole.pinned = &whole_session;
  auto pinned_whole = sequant::optimize_result(sum, opts_whole).expr;

  PinSession session(rr.record);
  auto opts2 = opts;
  opts2.pinned = &session;
  opts2.term_batch_axes = std::make_shared<AxesMap>();
  std::vector<sequant::ExprPtr> by_position(rr.record.terms.size());
  for (auto const& s : sum->as<sequant::Sum>().summands()) {
    auto optimized = sequant::optimize_result(s, opts2).expr;
    auto const pos = session.document_position(term_key(canonical_term(s)));
    REQUIRE(pos.has_value());
    REQUIRE(!by_position[*pos]);
    by_position[*pos] = optimized;
  }
  REQUIRE_NOTHROW(session.finish());
  for (std::size_t k = 0; k < by_position.size(); ++k) {
    REQUIRE(*by_position[k] == *whole.summand(k));
    REQUIRE(*by_position[k] == *pinned_whole->as<sequant::Sum>().summand(k));
  }
  REQUIRE(!session.document_position("not a key").has_value());

  // an incomplete per-summand pass is caught by finish()
  PinSession partial(rr.record);
  auto opts3 = opts;
  opts3.pinned = &partial;
  sequant::optimize_result(sum->as<sequant::Sum>().summand(0), opts3);
  REQUIRE_THROWS_WITH(partial.finish(), ContainsSubstring("unused"));

  // per-summand recording: records in call order
  RecordSession recording("R");
  auto opts4 = opts;
  opts4.recording = &recording;
  // each result is dropped right away (a realistic caller): its address may be
  // reused by the next result, and no record may be lost
  for (auto const& s : sum->as<sequant::Sum>().summands())
    sequant::optimize_result(s, opts4);
  auto per_summand = recording.finish();
  REQUIRE(per_summand.terms.size() == sum->as<sequant::Sum>().size());
  REQUIRE(per_summand.label == "R");
  REQUIRE(per_summand.terms.size() == 4);
  for (std::size_t k = 0; k < 4; ++k)
    REQUIRE(per_summand.terms[k].key ==
            term_key(canonical_term(sum->as<sequant::Sum>().summand(k))));
}

TEST_CASE("pinned order is the document order, not the clustered order",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  // summands 0 and 2 share the g*t contraction; reorder clusters them
  auto const sum =
      session_sum({L"g{i_1,i_2;a_3,a_4} * t{a_3,a_4;i_3,i_4} * f{i_5;i_6}",
                   L"h{i_1,i_2;i_3,i_4} * k{i_5;i_6}",
                   L"g{i_1,i_2;a_3,a_4} * t{a_3,a_4;i_3,i_4} * s{i_5;i_6}",
                   L"Z{i_1,i_2,i_5;i_3,i_4,i_6}"});
  auto opts = session_opts();
  auto rr = record_root("R", sum, opts);
  REQUIRE(rr.record.terms.size() == 4);
  std::vector<std::string> input_keys;
  for (auto const& s : sum->as<sequant::Sum>().summands())
    input_keys.push_back(term_key(canonical_term(s)));
  // reorder did change the order
  REQUIRE(record_keys(rr.record) != input_keys);

  auto doc = rr.record;
  std::reverse(doc.terms.begin(), doc.terms.end());
  PinSession session(doc);
  auto opts2 = opts;
  opts2.pinned = &session;
  opts2.term_batch_axes = std::make_shared<AxesMap>();
  auto pinned = sequant::optimize_result(sum, opts2).expr;
  REQUIRE(pinned->is<sequant::Sum>());
  auto const& ps = pinned->as<sequant::Sum>();
  REQUIRE(ps.size() == 4);
  auto const& rs = rr.optimized->as<sequant::Sum>();
  for (std::size_t k = 0; k < 4; ++k)
    REQUIRE(*ps.summand(k) == *rs.summand(3 - k));

  // the combined annotations follow the pinned order too
  RecordSession recording("R");
  auto opts3 = opts2;
  opts3.recording = &recording;
  opts3.term_batch_axes = std::make_shared<AxesMap>();
  PinSession session2(doc);
  opts3.pinned = &session2;
  auto pinned2 = sequant::optimize_result(sum, opts3).expr;
  std::vector<sequant::Expr const*> order;
  for (auto const& s : pinned2->as<sequant::Sum>().summands())
    order.push_back(s.get());
  REQUIRE(record_keys(recording.finish(order)) == record_keys(doc));
}

// ---- fingerprint ------------------------------------------------------------

namespace {

/// a two-term forest of dry-run trees over session_opts() (batchable: space a,
/// contracted); volatile leaves: label t
std::vector<sequant::eval::dryrun::EvalNodeDryRun> fingerprint_forest(
    sequant::OptimizeOptions const& opts) {
  using sequant::eval::dryrun::EvalExprDryRun;
  std::vector<sequant::eval::dryrun::EvalNodeDryRun> forest;
  for (auto s : {L"g{i_1,i_2;a_3,a_4} * t{a_3;i_3} * t{a_4;i_4}",
                 L"f{i_1;a_3} * t{a_3,a_4;i_3,i_4} * f{i_2;a_4}"}) {
    auto canon = sequant::opt::schedule::canonical_term(pin_parse(s));
    auto optimized = sequant::optimize(canon, opts);
    sequant::BinarizationOptions bopts;
    if (auto it = opts.term_batch_axes->find(optimized.get());
        it != opts.term_batch_axes->end())
      bopts.node_batch_axes = it->second;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    forest.push_back(sequant::binarize<EvalExprDryRun>(optimized, {}, bopts));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }
  return forest;
}

}  // namespace

TEST_CASE("fingerprint captures DAG-shaping inputs and nothing size-dependent",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto opts = session_opts();
  opts.batch_policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  auto const forest = fingerprint_forest(opts);
  auto fp = fingerprint_of(forest, {"R"}, opts, {L"a", L"i"});
  REQUIRE(fp.roots == std::vector<std::string>{"R"});
  REQUIRE(fp.batchable.at("a") == std::pair{true, false});
  REQUIRE(fp.batchable.at("i") == std::pair{false, false});
  REQUIRE(fp.volatile_leaves == std::vector<std::string>{"t"});
  REQUIRE(fp.mode_order == std::vector<std::string>{"a", "i"});
  REQUIRE(fp.objective_function == "DenseTimeSpaceBatched");
  REQUIRE(fp.reorder == "Reorder");
  REQUIRE(fp.scheduler == "forest_descent");
  REQUIRE(fp.cse.at("subnet") == "false");
  REQUIRE(fp.context.spaces == std::vector<std::string>{"a", "i"});
  REQUIRE(fp.context.strict_braket_symmetry ==
          sequant::get_default_context().assert_strict_braket_symmetry());

  auto opts2 = opts;
  opts2.idx_to_extent = [](sequant::Index const&) -> std::size_t { return 3; };
  opts2.batch_policy.batch_target_size = [](sequant::Index const&) {
    return std::size_t{7};
  };
  opts2.batch_policy.peak_threshold = 1e9;
  REQUIRE_NOTHROW(check_fingerprint(
      fp, fingerprint_of(forest, {"R"}, opts2, {L"a", L"i"})));

  auto planner = planner_record(opts2);
  REQUIRE(planner.at("batch_target_size") == "set");
  REQUIRE(planner.at("idx_to_extent") == "set");
  REQUIRE(planner.at("inner_pow") == "unset");
  REQUIRE(planner.at("peak_threshold") == std::to_string(1e9));
  REQUIRE(planner.count("volatile_weight") == 1);
  REQUIRE(planner.count("footprint_weight") == 1);
  REQUIRE(planner.count("peak_flops_tolerance") == 1);
  REQUIRE(planner.count("prune_outer_products") == 1);
  REQUIRE(planner.count("roofline.machine_balance") == 1);
  REQUIRE(planner.count("roofline.fast_mem_elems") == 1);
  REQUIRE(planner.count("roofline.block_tiles") == 1);
  REQUIRE(planner.count("roofline.block_prefactor") == 1);
}

TEST_CASE("check_fingerprint names every differing field",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  using Catch::Matchers::ContainsSubstring;
  auto opts = session_opts();
  auto const forest = fingerprint_forest(opts);
  auto a = fingerprint_of(forest, {"R"}, opts, {L"a", L"i"});
  REQUIRE_NOTHROW(check_fingerprint(a, a));

  auto b = a;
  b.persistent_only = !a.persistent_only;
  b.mode_order = {"i"};
  b.batchable["a"] = {false, true};
  REQUIRE_THROWS_WITH(check_fingerprint(a, b),
                      ContainsSubstring("fingerprint.persistent_only") &&
                          ContainsSubstring("fingerprint.mode_order") &&
                          ContainsSubstring("fingerprint.batchable.a"));

  // a different Context makes the live fingerprint differ
  auto ctx = sequant::get_default_context_snapshot();
  ctx.set(a.context.strict_braket_symmetry
              ? sequant::AssertStrictBraKetSymmetry::No
              : sequant::AssertStrictBraKetSymmetry::Yes);
  auto resetter = sequant::set_scoped_default_context(std::move(ctx));
  auto live = fingerprint_of(forest, {"R"}, opts, {L"a", L"i"});
  REQUIRE_THROWS_WITH(
      check_fingerprint(a, live),
      ContainsSubstring("fingerprint.context.strict_braket_symmetry"));
}
