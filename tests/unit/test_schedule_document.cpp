#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/batch_policy.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>

#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/legality.hpp>
#include <SeQuant/core/eval/node_batch_annotation.hpp>
#include <SeQuant/core/eval/op_counts.hpp>
#include <SeQuant/core/eval/ordered_executor.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/schedule_dag_record.hpp>
#include <SeQuant/core/eval/schedule_document.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/optimize/optimize.hpp>
#include <SeQuant/core/optimize/options.hpp>
#include <SeQuant/core/optimize/schedule_pin.hpp>
#include <SeQuant/core/utility/exception.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/space_qns.hpp>
#include "csv_size_fixtures.hpp"
#include "tensor_power_fixtures.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using namespace sequant::eval::schedule;
using sequant::tests::tensor_power::lax_braket_context;

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
    REQUIRE_THROWS_WITH(
        from_json(to_json(e)),
        ContainsSubstring("K") && ContainsSubstring("out of range"));
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

  // (b') leaves that are not a permutation of the factors: an ordinal used
  // twice, a factor unused -- caught by the reader, not only at pin time
  {
    Document e = make_doc();
    e.roots[0].terms[0].tree->right->leaf = 0;
    REQUIRE_THROWS_WITH(
        from_json(to_json(e)),
        ContainsSubstring("K") && ContainsSubstring("repeated"));
    Document f = make_doc();
    f.roots[0].terms[0].factors.push_back("f{a_1;i_1}");
    f.roots[0].terms[0].nodes.clear();
    REQUIRE_THROWS_WITH(from_json(to_json(f)),
                        ContainsSubstring("K") &&
                            ContainsSubstring("does not use every factor"));
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

TEST_CASE("schedule document reader caps the nesting depth",
          "[schedule-document]") {
  using Catch::Matchers::ContainsSubstring;
  // deep enough to overflow the stack if it were parsed to the end
  std::string deep(100000, '[');
  deep += std::string(100000, ']');
  REQUIRE_THROWS_WITH(from_json(deep), ContainsSubstring("nesting"));
  REQUIRE_THROWS_WITH(run_from_json(deep), ContainsSubstring("nesting"));
  // the cap leaves room for any term: this tree is 200 levels deep
  Document d = make_doc();
  auto& t = d.roots[0].terms[0];
  t.nodes.clear();
  t.factors.resize(201, "x{i_1;a_1}");
  TreeShape tree{0};
  for (std::size_t k = 1; k <= 200; ++k) {
    TreeShape up;
    up.left = std::make_unique<TreeShape>(std::move(tree));
    up.right = std::make_unique<TreeShape>(TreeShape{k});
    tree = std::move(up);
  }
  t.tree = std::move(tree);
  REQUIRE_NOTHROW(from_json(to_json(d)));
}

TEST_CASE("schedule document reader accepts a tree-only term",
          "[schedule-document]") {
  // a non-batched objective records the tree without node annotations
  Document d = make_doc();
  d.roots[0].terms[0].nodes.clear();
  Document back;
  REQUIRE_NOTHROW(back = from_json(to_json(d)));
  REQUIRE(back.roots[0].terms[0].nodes.empty());
  REQUIRE(compare(d, back).empty());
}

TEST_CASE("schedule document reader matches axis labels as whole indices",
          "[schedule-document]") {
  using Catch::Matchers::ContainsSubstring;
  // an axis with two proto indices, spelled as the V1 serializer spells it
  {
    Document d = make_doc();
    auto& t = d.roots[0].terms[0];
    t.factors = {"g{i_1,i_2;a_1<i_1,i_2>}:N-N-S", "t{a_1<i_1,i_2>;i_1}:N-N-S"};
    t.nodes = {{{{"a_1<i_1,i_2>", ModeKind::Contracted}}, {}, false}};
    REQUIRE_NOTHROW(from_json(to_json(d)));
    // a plain index of the node is still found
    t.nodes = {{{{"i_2", ModeKind::External}}, {}, false}};
    REQUIRE_NOTHROW(from_json(to_json(d)));
  }
  // a basis instance suffix is not an index: "2" in a_1<;2> or a_1<i_1;2>
  {
    Document d = make_doc();
    auto& t = d.roots[0].terms[0];
    t.factors = {"g{i_1;a_1<;2>}:N-N-S", "t{a_1<i_1;2>;i_1}:N-N-S"};
    t.nodes = {{{{"a_1<;2>", ModeKind::Contracted}}, {}, false}};
    REQUIRE_NOTHROW(from_json(to_json(d)));
    t.nodes = {{{{"a_1<i_1;2>", ModeKind::Contracted}}, {}, false}};
    REQUIRE_NOTHROW(from_json(to_json(d)));
    t.nodes = {{{{"2", ModeKind::Contracted}}, {}, false}};
    REQUIRE_THROWS_WITH(
        from_json(to_json(d)),
        ContainsSubstring("'2'") && ContainsSubstring("does not occur"));
  }
  // i_1 is not an index of a subtree whose only occupied index is i_10
  {
    Document d = make_doc();
    auto& t = d.roots[0].terms[0];
    t.factors = {"g{i_10;a_1}:N-N-S", "t{a_1;i_10}:N-N-S"};
    t.nodes = {{{{"i_1", ModeKind::External}}, {}, false}};
    REQUIRE_THROWS_WITH(
        from_json(to_json(d)),
        ContainsSubstring("'i_1'") && ContainsSubstring("does not occur"));
  }
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
    o.roles = {{"i_1", "Lockstep"}};
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
  REQUIRE(back.dag->values[0].occurrences[0].roles.size() == 1);
  REQUIRE(back.dag->values[0].occurrences[0].roles[0].role == "Lockstep");
  // a role is one of the names the writer produces
  REQUIRE_THROWS_WITH(from_json(mutate(d, "\"Lockstep\"", "\"Sideways\"")),
                      Catch::Matchers::ContainsSubstring("unknown role"));
  REQUIRE_THROWS_WITH(from_json(mutate(d, "\"LoopLocal\"", "\"Local\"")),
                      Catch::Matchers::ContainsSubstring("unknown role"));
  REQUIRE(back.dag->blocks[1].outputs[0].second == "AccumulateSum");
  REQUIRE(std::get<2>(back.dag->loop_order[0]) == 2);
  REQUIRE(summary_line(d) == "roots=1 terms=1 values=2 blocks=2 depth=1");
  // depth is the nesting of the block tree, not the blocks' loop-group depth
  auto d2 = d;
  d2.dag->blocks[1].depth = 5;
  REQUIRE(summary_line(d2) == "roots=1 terms=1 values=2 blocks=2 depth=1");
  d2.dag->blocks[1].parent.reset();
  REQUIRE(summary_line(d2) == "roots=1 terms=1 values=2 blocks=2 depth=0");
}

TEST_CASE("run section round-trips alone and JSON is made script-safe",
          "[schedule-document]") {
  using Catch::Matchers::ContainsSubstring;
  auto const d = make_dag_doc();
  REQUIRE(d.run);
  auto const back = run_from_json(to_json(*d.run));
  REQUIRE(back.n_batches == d.run->n_batches);
  REQUIRE(back.footprint == d.run->footprint);
  REQUIRE(back.builds == d.run->builds);
  REQUIRE(back.planner == d.run->planner);
  REQUIRE_THROWS_AS(run_from_json("[]"), sequant::Exception);
  REQUIRE_THROWS_AS(run_from_json("not json"), sequant::Exception);
  REQUIRE_THROWS_AS(run_from_json(R"({"extents": {}})"), sequant::Exception);
  // a block id is a whole decimal number: no sign, space or suffix
  auto const run_with_key = [](std::string const& key) {
    return R"({"extents": {}, "batch_targets": {}, "n_batches": {")" + key +
           R"(": 3}, "footprint": {}, "builds": {}, "planner": {}})";
  };
  REQUIRE_NOTHROW(run_from_json(run_with_key("12")));
  for (auto const* key : {"-1", "+1", "12abc", " 7", "", "0x1"})
    REQUIRE_THROWS_WITH(run_from_json(run_with_key(key)),
                        ContainsSubstring("not a block id"));
  // an int field must fit an int
  auto const d2 = make_dag_doc();
  REQUIRE_THROWS_WITH(
      from_json(mutate(d2, "\"slot\": 2", "\"slot\": 1099511627776")),
      ContainsSubstring("slot") && ContainsSubstring("fit"));
  REQUIRE_THROWS_WITH(
      from_json(mutate(d2, "\"home_depth\": 1", "\"home_depth\": -3000000000")),
      ContainsSubstring("home_depth") && ContainsSubstring("fit"));

  // a label containing "</script>" must not be able to end the element the
  // viewer inlines the document into
  auto doc = make_doc();
  doc.roots[0].label = "R</script><b>x\xE2\x80\xA8y\xE2\x80\xA9z";
  auto const safe = escape_for_html_script(to_json(doc));
  REQUIRE(safe.find('<') == std::string::npos);
  REQUIRE(safe.find("\xE2\x80\xA8") == std::string::npos);
  REQUIRE(safe.find("\xE2\x80\xA9") == std::string::npos);
  REQUIRE(safe.find("\\u003c/script>") != std::string::npos);
  REQUIRE(safe.find("\\u2028") != std::string::npos);
  REQUIRE(safe.find("\\u2029") != std::string::npos);
  // the escaped text denotes the same document
  auto const parsed = from_json(safe);
  REQUIRE(compare(parsed, doc).empty());
  REQUIRE(parsed.roots[0].label == doc.roots[0].label);
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

namespace {

/// L·M·pow(D,-1), D totally symmetric, L and M antisymmetric; i_1 sits in two
/// bra slots, so callers scope AssertStrictBraKetSymmetry::No
sequant::ExprPtr power_term() {
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  return antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
         antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
         pow(power_base(), rational{-1});
}

/// a Sum of power_term() and a plain tensor product
sequant::ExprPtr power_root() {
  return sequant::ex<sequant::Sum>(sequant::ExprPtrList{
      power_term(), pin_parse(L"g{i_1,i_2;a_1,a_2} * t{a_1,a_2;i_1,i_2}")});
}

}  // namespace

TEST_CASE("record_term then rebuild_pinned_term reproduce a TensorPower term",
          "[schedule-document][tensor-power]") {
  using namespace sequant::opt::schedule;
  auto _ = sequant::set_scoped_default_context(lax_braket_context());
  auto canon = canonical_term(power_term());
  REQUIRE(canon->is<sequant::Product>());
  REQUIRE(canonical_factors(canon).size() == 3);
  auto [optimized, annotations] = pin_optimize(canon);
  REQUIRE(annotations.size() == 2);

  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.factors.size() == 3);
  REQUIRE(rec.nodes.size() == 2);
  auto pinned = rebuild_pinned_term(rec, canon);
  REQUIRE(*pinned.expr == *optimized);
  REQUIRE(pinned.nodes.size() == annotations.size());
  for (std::size_t k = 0; k < annotations.size(); ++k) {
    REQUIRE(pinned.nodes[k].axes == annotations[k].axes);
    REQUIRE(pinned.nodes[k].opened_here == annotations[k].opened_here);
  }

  // a lone power is a term too
  auto lone = canonical_term(
      sequant::pow(sequant::ex<sequant::Tensor>(L"D", sequant::bra{L"i_1"},
                                                sequant::ket{L"a_1"}),
                   sequant::rational{-1}));
  REQUIRE(lone->is<sequant::TensorPower>());
  auto lone_rec = record_term(lone, lone, {});
  REQUIRE(lone_rec.tree.has_value());
  REQUIRE(lone_rec.tree->leaf == std::optional<std::size_t>{0});
  REQUIRE(*rebuild_pinned_term(lone_rec, lone).expr == *lone);
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

TEST_CASE("record_root then pin reproduce a Sum with a TensorPower term",
          "[schedule-document][tensor-power]") {
  using namespace sequant::opt::schedule;
  auto _ = sequant::set_scoped_default_context(lax_braket_context());
  auto const sum = power_root();
  auto opts = session_opts();
  auto rr = record_root("E", sum, opts);
  REQUIRE(rr.record.terms.size() == 2);
  for (auto const& t : rr.record.terms) REQUIRE(t.pinnable);
  auto const& recorded_axes = (*opts.term_batch_axes)[rr.optimized.get()];

  PinSession session(rr.record);
  auto opts2 = opts;
  opts2.pinned = &session;
  opts2.term_batch_axes = std::make_shared<AxesMap>();
  auto pinned = sequant::optimize_result(sum, opts2).expr;
  REQUIRE(*pinned == *rr.optimized);
  require_same_annotations((*opts2.term_batch_axes)[pinned.get()],
                           recorded_axes);
  REQUIRE_NOTHROW(session.finish());

  sequant::BinarizationOptions bopts;
  bopts.node_batch_axes = recorded_axes;
  REQUIRE_NOTHROW(sequant::binarize(
      sequant::ResultExpr{sequant::Variable{L"E"}, pinned}, bopts));
}

namespace {

/// records \p sum under batched and unbatched objectives, writes the record
/// to JSON, reads it back, pins from it and checks the pinned run reproduces
/// the recorded document
void require_json_pin_round_trip(sequant::ExprPtr const& sum) {
  using namespace sequant::opt::schedule;
  for (auto const objective :
       {sequant::ObjectiveFunction::DenseTimeSpaceBatched,
        sequant::ObjectiveFunction::DenseFLOPs}) {
    CAPTURE(static_cast<int>(objective));
    auto opts = session_opts();
    opts.objective_function = objective;
    auto rr = record_root("R", sum, opts);
    bool const batched =
        objective == sequant::ObjectiveFunction::DenseTimeSpaceBatched;
    bool any_nodes = false;
    for (auto const& t : rr.record.terms)
      any_nodes = any_nodes || !t.nodes.empty();
    REQUIRE(any_nodes == batched);

    Document doc;
    doc.fingerprint.roots = {"R"};
    doc.roots = {rr.record};
    auto const loaded = from_json(to_json(doc));
    REQUIRE(compare(doc, loaded).empty());

    // pin from the loaded record and record what the pinned run produced
    PinSession session(loaded.roots[0]);
    auto opts2 = opts;
    opts2.pinned = &session;
    opts2.term_batch_axes = std::make_shared<AxesMap>();
    auto rr2 = record_root("R", sum, opts2);
    REQUIRE_NOTHROW(session.finish());
    REQUIRE(*rr2.optimized == *rr.optimized);

    Document doc2;
    doc2.fingerprint.roots = {"R"};
    doc2.roots = {rr2.record};
    REQUIRE(compare(loaded, from_json(to_json(doc2))).empty());
  }
}

}  // namespace

TEST_CASE("a pin loaded from JSON reproduces the recorded document",
          "[schedule-document]") {
  require_json_pin_round_trip(session_root());
}

TEST_CASE("a pin loaded from JSON reproduces a TensorPower term",
          "[schedule-document][tensor-power]") {
  auto _ = sequant::set_scoped_default_context(lax_braket_context());
  require_json_pin_round_trip(power_root());
}

TEST_CASE("an axis with several proto indices survives JSON and pins",
          "[schedule-document]") {
  using namespace sequant::opt::schedule;
  auto canon = canonical_term(
      pin_parse(L"g{i_1,i_2;a_1<i_1,i_2>} * t{a_1<i_1,i_2>;i_1}"));
  REQUIRE(canonical_factors(canon).size() == 2);
  auto [optimized, annotations] = pin_optimize(canon);
  // annotate the single contraction node with the PNO index
  std::optional<sequant::Index> pno;
  for (auto const& f : canonical_factors(canon))
    for (auto const& ix : f->as<sequant::Tensor>().const_indices())
      if (ix.proto_indices().size() == 2) pno = ix;
  REQUIRE(pno.has_value());
  annotations.assign(1, sequant::NodeBatchAnnotation{});
  annotations[0].axes.emplace_back(*pno, sequant::BatchModeType::Contracted);
  auto rec = record_term(canon, optimized, annotations);
  REQUIRE(rec.nodes.size() == 1);
  auto const& label = rec.nodes[0].axes[0].label;
  CAPTURE(label, rec.factors);
  REQUIRE(label.find(',') != std::string::npos);
  REQUIRE(label.find(", ") == std::string::npos);  // the V1 spelling

  Document doc;
  doc.fingerprint.roots = {"R"};
  doc.roots = {{"R", {rec}}};
  Document loaded;
  REQUIRE_NOTHROW(loaded = from_json(to_json(doc)));
  auto pinned = rebuild_pinned_term(loaded.roots[0].terms[0], canon);
  REQUIRE(*pinned.expr == *optimized);
  REQUIRE(pinned.nodes.size() == 1);
  REQUIRE(pinned.nodes[0].axes == annotations[0].axes);
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

  SECTION("externals labeled differently than in the record") {
    // same terms, externals i_1,i_2 -> i_5,i_6: dummies are canonicalized
    // away, externals are not, so no key matches and the message says so
    auto relabeled =
        session_sum({L"g{i_5,i_6;a_3,a_4} * t{a_3;i_3} * t{a_4;i_4}",
                     L"1/2 * g{i_5,i_6;a_3,a_4} * t{a_3,a_4;i_3,i_4}",
                     L"f{i_5;a_3} * t{a_3,a_4;i_3,i_4} * f{i_6;a_4}",
                     L"X{i_5,i_6;i_3,i_4}"});
    PinSession session(rr.record);
    auto opts2 = opts;
    opts2.pinned = &session;
    REQUIRE_THROWS_WITH(sequant::optimize_result(relabeled, opts2),
                        ContainsSubstring("missing") &&
                            ContainsSubstring("no term matched") &&
                            ContainsSubstring("external indices"));
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
    // a recording optimizes the mixed product as usual and lists it as not
    // pinnable (no tree, no nodes), so the root's term set stays complete
    auto opts3 = opts;
    RecordSession recording("M");
    opts3.recording = &recording;
    auto const optimized = sequant::optimize_result(mixed, opts3).expr;
    REQUIRE(optimized->is<sequant::Sum>());
    REQUIRE(optimized->as<sequant::Sum>().size() == 2);
    auto const record = recording.finish();
    REQUIRE(record.terms.size() == 2);
    TermRecord const* unpinnable = nullptr;
    for (auto const& t : record.terms)
      if (!t.pinnable) {
        REQUIRE(!unpinnable);
        unpinnable = &t;
        REQUIRE(!t.tree);
        REQUIRE(t.nodes.empty());
        REQUIRE(t.factors.size() == 2);
      }
    REQUIRE(unpinnable);
    // the flag survives JSON, and a pin refuses the term
    Document doc;
    doc.fingerprint.roots = {"M"};
    doc.roots = {record};
    auto const loaded = from_json(to_json(doc));
    REQUIRE(compare(doc, loaded).empty());
    REQUIRE_THROWS_WITH(
        rebuild_pinned_term(
            *unpinnable, canonical_term(mixed->as<sequant::Sum>().summand(0))),
        ContainsSubstring("not pinnable"));
    PinSession session2(loaded.roots[0]);
    auto opts4 = opts;
    opts4.pinned = &session2;
    REQUIRE_THROWS_WITH(sequant::optimize_result(mixed, opts4),
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
  // a finished session is single-use: a further run through it is refused
  REQUIRE_THROWS_WITH(sequant::optimize_result(sum, opts2),
                      ContainsSubstring("single-use"));
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
  REQUIRE(fp.batchable.size() == 2);  // exactly the bases present
  REQUIRE(fp.context.strict_braket_symmetry ==
          sequant::get_default_context().assert_strict_braket_symmetry());
  REQUIRE(!fp.context.canonicalization_method.empty());

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

TEST_CASE("fingerprint sees the base of a TensorPower leaf",
          "[schedule-document][tensor-power]") {
  using namespace sequant::opt::schedule;
  using sequant::eval::dryrun::EvalExprDryRun;
  auto _ = sequant::set_scoped_default_context(lax_braket_context());
  auto opts = session_opts();
  opts.batch_policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"D";
  };
  auto optimized = sequant::optimize(canonical_term(power_term()), opts);
  sequant::BinarizationOptions bopts;
  if (auto it = opts.term_batch_axes->find(optimized.get());
      it != opts.term_batch_axes->end())
    bopts.node_batch_axes = it->second;
  std::vector<sequant::eval::dryrun::EvalNodeDryRun> forest;
  forest.push_back(sequant::binarize<EvalExprDryRun>(
      sequant::ResultExpr{sequant::Variable{L"E"}, optimized}, bopts));
  sequant::eval::schedule::Fingerprint fp;
  REQUIRE_NOTHROW(fp = fingerprint_of(forest, {"E"}, opts, {L"a", L"i"}));
  REQUIRE(fp.volatile_leaves == std::vector<std::string>{"D"});
  REQUIRE(fp.batchable.size() == 2);
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
  b.context.canonicalization_method = "other";
  b.context.ignore_named_index_labels = !a.context.ignore_named_index_labels;
  REQUIRE_THROWS_WITH(
      check_fingerprint(a, b),
      ContainsSubstring("fingerprint.persistent_only") &&
          ContainsSubstring("fingerprint.mode_order") &&
          ContainsSubstring("fingerprint.batchable.a") &&
          ContainsSubstring("fingerprint.context.canonicalization_method") &&
          ContainsSubstring("fingerprint.context.ignore_named_index_labels"));

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

namespace {

/// the pipeline of record_root -> binarize -> compute_dag_boulevard ->
/// analyze_legality -> build_ordered_schedule (MPQC's driver shape), over one
/// tree per root
struct DagFixture {
  std::vector<sequant::eval::dryrun::EvalNodeDryRun> forest;
  std::vector<RootRecord> roots;
  sequant::eval::RichSchedule rich;
  sequant::eval::LegalitySchedule legality;
  sequant::eval::OrderedSchedule ordered;
};

DagFixture dag_fixture(
    std::vector<std::pair<std::string, sequant::ExprPtr>> const& roots,
    sequant::BatchPolicy const& policy) {
  using sequant::eval::dryrun::EvalExprDryRun;
  DagFixture fx;
  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = pin_extent;
  opts.batch_policy = policy;
  opts.term_batch_axes = std::make_shared<AxesMap>();
  for (auto const& [label, expr] : roots) {
    auto rr = sequant::opt::schedule::record_root(label, expr, opts);
    sequant::BinarizationOptions bopts;
    if (auto it = opts.term_batch_axes->find(rr.optimized.get());
        it != opts.term_batch_axes->end())
      bopts.node_batch_axes = it->second;
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    fx.forest.push_back(
        sequant::binarize<EvalExprDryRun>(rr.optimized, {}, bopts));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
    fx.roots.push_back(std::move(rr.record));
  }
  sequant::eval::dryrun::SizeRegime const regime;
  sequant::eval::dryrun::CostModel const cm{regime};
  fx.rich = sequant::eval::compute_dag_boulevard(fx.forest, cm,
                                                 policy.batch_target_size);
  fx.legality = sequant::eval::analyze_legality(fx.rich, fx.forest, policy);
  fx.ordered =
      sequant::eval::build_ordered_schedule(fx.rich, fx.legality, policy, {});
  return fx;
}

/// occupied external modes batchable, nothing contracted; slices at any size
sequant::BatchPolicy external_occ_policy() {
  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const&) {
    return false;
  };
  policy.is_batchable_external_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"i";
  };
  policy.batch_target_size = [](sequant::Index const&) -> std::size_t {
    return 2;
  };
  policy.batch_spectator_indices = true;
  policy.peak_threshold = 1;
  return policy;
}

/// two roots sharing X = g*t; R1 contracts X's ket occ, R2 its bra occ, so
/// each root can batch X only on the other one
DagFixture two_root_fixture() {
  return dag_fixture(
      {{"R1",
        session_sum({L"g{i_1;a_1} * t{a_1;i_2} * f{i_2;i_3}", L"h{i_1;i_3}"})},
       {"R2",
        session_sum({L"g{i_1;a_1} * t{a_1;i_2} * f{i_3;i_1}", L"h{i_3;i_2}"})}},
      external_occ_policy());
}

/// what a value's expr spells it from: a leaf's tensor label, else the
/// operands
std::string rhs_of(ValueRecord const& v) {
  auto const eq = v.expr.find(" = ");
  REQUIRE(eq != std::string::npos);
  return v.expr.substr(eq + 3);
}

/// the products of g and t (the X of the fixture), by id
std::vector<ValueRecord> x_values(DagRecord const& dag) {
  std::vector<ValueRecord> xs;
  for (auto const& v : dag.values) {
    if (v.leaf) continue;
    auto const rhs = rhs_of(v);
    if ((rhs.starts_with("g") && rhs.find(" * t") != std::string::npos) ||
        (rhs.starts_with("t") && rhs.find(" * g") != std::string::npos))
      xs.push_back(v);
  }
  std::ranges::sort(xs, {}, &ValueRecord::id);
  return xs;
}

Document doc_with(std::vector<RootRecord> const& roots, DagRecord dag) {
  Document d;
  for (auto const& r : roots) d.fingerprint.roots.push_back(r.label);
  d.roots = roots;
  d.dag = std::move(dag);
  return d;
}

}  // namespace

TEST_CASE(
    "record_dag gives a node sliced by two loops two values with each other "
    "as siblings",
    "[schedule-document]") {
  auto const fx = two_root_fixture();
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  auto const xs = x_values(dag);
  REQUIRE(xs.size() == 2);
  REQUIRE(xs[0].siblings == std::vector<std::string>{xs[1].id});
  REQUIRE(xs[1].siblings == std::vector<std::string>{xs[0].id});
  REQUIRE(xs[0].home_block != xs[1].home_block);
  for (auto const& v : xs) REQUIRE(v.home_modes.size() == 1);  // neither whole
  REQUIRE(dag.blocks[0].kind == "root");
  REQUIRE(std::ranges::count_if(dag.blocks, [](auto const& b) {
            return b.kind == "external";
          }) == 2);
  // both roots are two-summand Sums: one partial sum each, the root itself
  std::set<std::string> ids;
  for (auto const& v : dag.values) ids.insert(v.id);
  REQUIRE(ids.count("R1+1") == 1);
  REQUIRE(ids.count("R2+1") == 1);
  REQUIRE(!dag.loops.empty());

  // a leaf homed at the root block feeding X[i] inside its loop block is read
  // in slices; X[i] feeding its consumer in the same block is not
  auto const edge = [&](std::string const& from, std::string const& to) {
    auto const it = std::ranges::find_if(
        dag.edges, [&](auto const& e) { return e.from == from && e.to == to; });
    REQUIRE(it != dag.edges.end());
    return *it;
  };
  std::map<std::string, ValueRecord const*> by_id;
  for (auto const& v : dag.values) by_id[v.id] = &v;
  auto const& x = xs[0];
  REQUIRE(dag.blocks[x.home_block].parent == std::size_t{0});
  std::optional<std::string> leaf_in, consumer;
  for (auto const& e : dag.edges) {
    if (e.to == x.id && by_id.at(e.from)->leaf) leaf_in = e.from;
    if (e.from == x.id) consumer = e.to;
  }
  REQUIRE(leaf_in);
  REQUIRE(consumer);
  REQUIRE(by_id.at(*leaf_in)->home_block == 0);
  REQUIRE(by_id.at(*consumer)->home_block == x.home_block);
  REQUIRE(edge(*leaf_in, x.id).sliced);
  REQUIRE(!edge(x.id, *consumer).sliced);
}

TEST_CASE("record_dag spells every value as its name and what makes it",
          "[schedule-document]") {
  auto const fx = two_root_fixture();
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  std::set<std::string> names;
  for (auto const& v : dag.values) {
    CAPTURE(v.id, v.name, v.expr);
    REQUIRE(names.insert(v.name).second);
    // the name opens expr, bare or in the form the value is produced in
    REQUIRE(v.expr.starts_with(v.name));
    auto const after = v.expr.substr(v.name.size());
    REQUIRE((after.starts_with(" = ") || after.starts_with("[")));
    auto const rhs = rhs_of(v);
    if (v.leaf) {
      REQUIRE(std::set<std::string>{"f", "g", "h", "t"}.contains(rhs));
    } else {
      REQUIRE((rhs.find(" * ") != std::string::npos ||
               rhs.find(" + ") != std::string::npos));
    }
  }
  // the whole X of R1's loop, assembled from per-batch partials, by name
  auto const xs = x_values(dag);
  REQUIRE(!xs.empty());
  for (auto const& x : xs) REQUIRE(x.name.starts_with("R"));
}

TEST_CASE("record_dag takes persistence from the cell table",
          "[schedule-document]") {
  using sequant::eval::dryrun::EvalNodeDryRun;
  auto const fx = two_root_fixture();
  // nothing volatile: a value survives an evaluation only when nothing reads
  // it again, i.e. the roots; every leaf and every X is read again or bound
  auto const dag0 =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  for (auto const& v : dag0.values)
    REQUIRE(v.persistent == (v.id == "R1+1" || v.id == "R2+1"));
  // t volatile: whatever holds t (t, X = g*t, the roots) cannot persist; the
  // leaves read by a value that holds t (g by X, f and h above it) persist
  auto const dag1 =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots,
                 [](EvalNodeDryRun const& n) {
                   return n->is_tensor() && n->as_tensor().label() == L"t";
                 });
  std::size_t leaves = 0;
  for (auto const& v : dag1.values) {
    if (v.leaf) {
      ++leaves;
      REQUIRE(v.persistent == (rhs_of(v) != "t"));
    } else {
      REQUIRE(!v.persistent);
    }
  }
  REQUIRE(leaves == 4);  // g, t, f, h
}

TEST_CASE("record_dag persistence is the flag of the assembled form",
          "[schedule-document]") {
  using sequant::eval::dryrun::EvalNodeDryRun;
  // a contracted (a-batched) root: X = g * h is a partial sum per a-batch,
  // assembled whole at the root scope, then multiplied by the volatile t
  auto const fx = dag_fixture(
      {{"R1", pin_parse(L"g{i_1,i_2;a_3,a_4} * h{a_3,a_4;i_3,i_4} * "
                        L"t{i_3,i_4;i_5,i_6}")}},
      session_opts().batch_policy);
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots,
                 [](EvalNodeDryRun const& n) {
                   return n->is_tensor() && n->as_tensor().label() == L"t";
                 });
  std::map<std::string, ValueRecord const*> by_id;
  for (auto const& v : dag.values) by_id[v.id] = &v;
  // the partial escapes its loop block as an AccumulateSum: the value is
  // homed in the block (built per batch) and its whole, read by the
  // t-holding product, is what the executor keeps across evaluations
  // (both a-loops are batched, so the partial escapes two nested blocks)
  std::set<std::string> assembled;
  for (auto const& b : dag.blocks)
    for (auto const& [id, kind] : b.outputs)
      if (kind == "AccumulateSum") assembled.insert(id);
  REQUIRE(assembled.size() == 1);
  auto const& x = *by_id.at(*assembled.begin());
  CAPTURE(x.expr, x.home_block);
  REQUIRE(x.home_block != 0);
  REQUIRE(x.persistent);
  // g and h feed only the partial, whose whole holds no t: not persistent
  for (auto const& v : dag.values)
    if (v.leaf && (rhs_of(v) == "g" || rhs_of(v) == "h"))
      REQUIRE(!v.persistent);
}

TEST_CASE("record_dag ids resolve for single-summand roots and scaled terms",
          "[schedule-document]") {
  auto const fx = dag_fixture(
      {{"R1", pin_parse(L"1/2 * g{i_1;a_1} * t{a_1;i_2} * f{i_2;i_3}")},
       {"R2", pin_parse(L"g{i_1;a_1} * t{a_1;i_2} * f{i_2;i_3}")}},
      session_opts().batch_policy);
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  REQUIRE(fx.roots.size() == 2);
  auto const& scaled = fx.roots[0].terms.at(0).key;
  auto const& plain = fx.roots[1].terms.at(0).key;
  REQUIRE(scaled != plain);

  std::size_t nscale = 0;
  for (auto const& v : dag.values) {
    REQUIRE(!v.id.starts_with("R1+"));
    REQUIRE(!v.id.starts_with("R2+"));
    if (v.id.ends_with("#scale")) {
      ++nscale;
      REQUIRE(v.id == scaled + "#scale");
    }
  }
  REQUIRE(nscale == 1);

  // the scaled term's tree ordinals are the plain term's: drop the occurrences
  // of the scale node and of the Constant leaf it multiplies by
  std::map<std::string, std::set<std::size_t>> ordinals;
  for (auto const& v : dag.values)
    for (auto const& o : v.occurrences) ordinals[o.term].insert(o.node);
  std::set<std::size_t> extra;
  for (auto const& v : dag.values)
    if (v.id == scaled + "#scale" || v.id == scaled + "#scalar")
      for (auto const& o : v.occurrences)
        if (o.term == scaled) extra.insert(o.node);
  REQUIRE(extra.size() == 2);
  for (auto n : extra) REQUIRE(ordinals[scaled].erase(n) == 1);
  REQUIRE(ordinals[scaled] == ordinals[plain]);
  REQUIRE(ordinals[plain] == std::set<std::size_t>{0, 1, 2, 3, 4});
}

TEST_CASE("record_dag is deterministic and compare sees no difference",
          "[schedule-document]") {
  auto const fx = two_root_fixture();
  auto const dag1 =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  auto const dag2 =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  REQUIRE(compare(doc_with(fx.roots, dag1), doc_with(fx.roots, dag2)).empty());
  // two independent optimize -> binarize -> DAG pipelines on the same input
  auto const fx2 = two_root_fixture();
  auto const dag3 =
      record_dag(fx2.rich, fx2.legality, fx2.ordered, fx2.forest, fx2.roots);
  REQUIRE(compare(doc_with(fx.roots, dag1), doc_with(fx2.roots, dag3)).empty());
  REQUIRE(to_json(doc_with(fx.roots, dag1)) ==
          to_json(doc_with(fx.roots, dag2)));

  std::set<std::string> ids;
  for (auto const& v : dag1.values) REQUIRE(ids.insert(v.id).second);
  REQUIRE(ids.size() == fx.rich.cells.size());
  for (auto const& e : dag1.edges) {
    REQUIRE(ids.count(e.from) == 1);
    REQUIRE(ids.count(e.to) == 1);
  }
  for (auto const& b : dag1.blocks) {
    for (auto const& s : b.steps) {
      REQUIRE(s.build.has_value() != s.block.has_value());
      if (s.build) REQUIRE(ids.count(*s.build) == 1);
      if (s.block) REQUIRE(*s.block < dag1.blocks.size());
    }
    for (auto const& [id, kind] : b.outputs) REQUIRE(ids.count(id) == 1);
  }
}

TEST_CASE("record_dag names an adjoint factor like a leaf, its operands below",
          "[schedule-document]") {
  auto const fx =
      dag_fixture({{"R1", pin_parse(L"g{i_1;a_1} * t⁺{a_1;i_2} * f{i_2;i_3}")},
                   {"R2", pin_parse(L"g{i_1;a_1} * s{a_1;i_2} * f{i_2;i_3}")}},
                  session_opts().batch_policy);
  DagRecord dag;
  REQUIRE_NOTHROW(
      dag = record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots));
  auto const& adj = fx.roots[0].terms.at(0).key;
  auto const& twin = fx.roots[1].terms.at(0).key;

  std::map<std::string, std::set<std::size_t>> ordinals;
  std::set<std::string> ids;
  for (auto const& v : dag.values) {
    ids.insert(v.id);
    for (auto const& o : v.occurrences) ordinals[o.term].insert(o.node);
  }
  REQUIRE(ordinals[adj] == ordinals[twin]);
  REQUIRE(ordinals[adj] == std::set<std::size_t>{0, 1, 2, 3, 4});
  std::size_t nsub = 0;
  for (auto const& id : ids)
    if (id.starts_with(adj + "#") && (id.ends_with("/0") || id.ends_with("/1")))
      ++nsub;
  REQUIRE(nsub == 2);
  // the adjoint value reads its one operand, marked "⁺"
  std::size_t nadjoint = 0;
  for (auto const& v : dag.values)
    if (v.expr.ends_with("⁺") && v.expr.find(" * ") == std::string::npos)
      ++nadjoint;
  REQUIRE(nadjoint > 0);
}

TEST_CASE("record_dag names a power factor like a leaf, its operands below",
          "[schedule-document][tensor-power]") {
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  auto _ = set_scoped_default_context(lax_braket_context());
  // the twin replaces pow(D,-1) by a plain tensor of D's shape
  auto const twin_expr =
      antisymm(L"L", bra{L"a_1", L"a_2"}, ket{L"i_1", L"i_2"}) *
      antisymm(L"M", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}) *
      ex<Tensor>(L"N", bra{L"i_1", L"i_2"}, ket{L"a_1", L"a_2"}, symmetric());
  auto const fx = dag_fixture({{"R1", power_term()}, {"R2", twin_expr}},
                              occ_contracted_policy());
  DagRecord dag;
  REQUIRE_NOTHROW(
      dag = record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots));
  auto const& pw = fx.roots[0].terms.at(0).key;
  auto const& twin = fx.roots[1].terms.at(0).key;

  std::map<std::string, std::set<std::size_t>> ordinals;
  std::set<std::string> ids;
  for (auto const& v : dag.values) {
    ids.insert(v.id);
    for (auto const& o : v.occurrences) ordinals[o.term].insert(o.node);
  }
  REQUIRE(ordinals[pw] == ordinals[twin]);
  REQUIRE(ordinals[pw] == std::set<std::size_t>{0, 1, 2, 3, 4});
  std::size_t nsub = 0;
  for (auto const& id : ids)
    if (id.starts_with(pw + "#") && (id.ends_with("/0") || id.ends_with("/1")))
      ++nsub;
  REQUIRE(nsub == 2);
  // the power is spelled as its operand raised to its exponent, and its id
  // is the ordinal its sub-positions hang off
  auto const power = std::ranges::find_if(dag.values, [](auto const& v) {
    return v.expr.find(" ^ ") != std::string::npos;
  });
  REQUIRE(power != dag.values.end());
  CHECK(rhs_of(*power).starts_with("D"));
  CHECK(rhs_of(*power).ends_with(" ^ -1"));
  CHECK(ids.count(power->id + "/0") == 1);
  CHECK(ids.count(power->id + "/1") == 1);

  // the dag section round-trips through JSON
  auto const doc = doc_with(fx.roots, dag);
  REQUIRE(compare(doc, from_json(to_json(doc))).empty());
}

TEST_CASE("record_dag names the nodes below a unary node's inlined operand",
          "[schedule-document][inline]") {
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  using sequant::eval::dryrun::EvalExprDryRun;
  // pow(D,-1) inside ℓ raises the inlined broadcast D: a Sum subtree below
  // the Power node
  container::svector<ResultExpr> const sys{D_broadcast(), ell_definition(),
                                           consumer_E()};
  auto inlined = binarize<EvalExprDryRun>(sys);
  std::vector<eval::dryrun::EvalNodeDryRun> forest{inlined.trees[2]};
  stamp_occ_loop(forest[0]);
  auto const policy = occ_contracted_policy();
  eval::dryrun::SizeRegime const regime;
  eval::dryrun::CostModel const cm{regime};
  auto const rich =
      eval::compute_dag_boulevard(forest, cm, policy.batch_target_size);
  auto const legality = eval::analyze_legality(rich, forest, policy);
  auto const ordered = eval::build_ordered_schedule(rich, legality, policy, {});
  std::vector<RootRecord> const roots{
      RootRecord{.label = "E", .terms = {TermRecord{.key = "E"}}}};

  DagRecord dag;
  REQUIRE_NOTHROW(dag = record_dag(rich, legality, ordered, forest, roots));
  std::set<std::string> ids;
  for (auto const& v : dag.values) REQUIRE(ids.insert(v.id).second);
  REQUIRE(ids.size() == rich.cells.size());
  // the power and its operand, the inlined D, keep their ids; D's own nodes
  // hang below the operand
  auto const power = std::ranges::find_if(dag.values, [](auto const& v) {
    return v.expr.find(" ^ ") != std::string::npos;
  });
  REQUIRE(power != dag.values.end());
  // (the exponent -1 is first met as a prefactor inside D, which names it)
  REQUIRE(ids.count(power->id + "/0") == 1);
  std::size_t nbelow = 0;
  for (auto const& id : ids)
    if (id.starts_with(power->id + "/0.")) ++nbelow;
  REQUIRE(nbelow > 0);

  auto const doc = doc_with(roots, dag);
  REQUIRE(compare(doc, from_json(to_json(doc))).empty());
}

TEST_CASE("record_dag records an inlined head like the expression it replaces",
          "[schedule-document][tensor-power][inline]") {
  using namespace sequant;
  using namespace sequant::tests::tensor_power;
  using sequant::eval::dryrun::EvalExprDryRun;
  using sequant::eval::dryrun::EvalNodeDryRun;
  std::vector<RootRecord> const roots{
      RootRecord{.label = "E", .terms = {TermRecord{.key = "E"}}}};
  auto const policy = occ_contracted_policy();
  auto const record = [&](EvalNodeDryRun node) {
    stamp_occ_loop(node);
    std::vector<EvalNodeDryRun> const forest{node};
    eval::dryrun::SizeRegime const regime;
    eval::dryrun::CostModel const cm{regime};
    auto const rich =
        eval::compute_dag_boulevard(forest, cm, policy.batch_target_size);
    auto const legality = eval::analyze_legality(rich, forest, policy);
    auto const ordered =
        eval::build_ordered_schedule(rich, legality, policy, {});
    DagRecord dag;
    REQUIRE_NOTHROW(dag = record_dag(rich, legality, ordered, forest, roots));
    REQUIRE(dag.values.size() == rich.cells.size());
    return dag;
  };

  // E = L · ℓ with ℓ = M ⊙ pow(D,-1) inlined, and E = L · (M ⊙ pow(D,-1))
  auto const sys = binarize<EvalExprDryRun>(
      std::vector<ResultExpr>{ell_definition(), consumer_E()});
  REQUIRE(sys.roots().size() == 1);
  auto const dag = record(sys.roots()[0].get());
  auto const direct = [&] {
    auto _ = set_scoped_default_context(lax_braket_context());
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
    return record(binarize<EvalExprDryRun>(hadamard_power_term()));
    SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  }();

  // the same values, named and spelled alike
  REQUIRE(dag.values.size() == 7);
  REQUIRE(dag.values.size() == direct.values.size());
  for (std::size_t v = 0; v != dag.values.size(); ++v) {
    CAPTURE(v);
    CHECK(dag.values[v].id == direct.values[v].id);
    CHECK(dag.values[v].expr == direct.values[v].expr);
  }

  auto const doc = doc_with(roots, dag);
  REQUIRE(compare(doc, from_json(to_json(doc))).empty());
}

namespace {

/// an index of each space the two-root fixture uses
sequant::Index const& occ_index() {
  static sequant::Index const ix{L"i_1"};
  return ix;
}
sequant::Index const& virt_index() {
  static sequant::Index const ix{L"a_1"};
  return ix;
}

/// dry-run evaluation of a DagFixture through the ordered multi-root entry
/// point, as the ordered executor tests drive it (fresh cache each call)
sequant::container::svector<sequant::ResultPtr> dry_run(
    DagFixture const& fx, RunRecord* run = nullptr,
    std::function<std::string(std::size_t)> const& id_of = {}) {
  using sequant::eval::dryrun::EvalNodeDryRun;
  sequant::eval::dryrun::SizeRegime regime;
  // the sizes the fixture was optimized with (pin_extent)
  regime.space_extent = {
      {std::wstring{occ_index().space().base_key()}, pin_extent(occ_index())},
      {std::wstring{virt_index().space().base_key()},
       pin_extent(virt_index())}};
  auto const cm =
      std::make_shared<sequant::eval::dryrun::CostModel const>(regime);
  sequant::eval::dryrun::DryRunLeafEvaluator const leaf_eval{cm};
  auto cache = sequant::CacheManager<EvalNodeDryRun>::empty();
  auto const aops = sequant::eval::dryrun::make_dryrun_array_ops(cm);
  cache.set_array_ops(&aops);
  std::function<std::size_t(sequant::Index const&)> const target =
      external_occ_policy().batch_target_size;
  return sequant::eval::evaluate_ordered_multiroot(
      fx.forest, fx.ordered, fx.rich,
      sequant::container::svector<std::string>(fx.forest.size()), leaf_eval,
      cache, target, {}, {}, run, id_of);
}

}  // namespace

TEST_CASE("ordered dry run fills the run record", "[schedule-document]") {
  auto const fx = two_root_fixture();
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  // record_dag keeps one value record per engine value, in value order
  REQUIRE(dag.values.size() == fx.rich.cells.size());
  auto const id_of = [&](std::size_t vid) { return dag.values.at(vid).id; };

  RunRecord run;
  auto const results = dry_run(fx, &run, id_of);
  REQUIRE(results.size() == 2);

  // the fixture's external occ extent and its batch target
  std::size_t const occ_extent = pin_extent(occ_index());
  std::size_t const occ_target =
      external_occ_policy().batch_target_size(occ_index());
  REQUIRE(occ_target < occ_extent);  // a genuine batching

  REQUIRE(run.n_batches.size() == 2);  // two external blocks
  // each external block loops over occ: ceil(extent(occ) / target(occ))
  for (auto const& [blk, n] : run.n_batches) {
    REQUIRE(blk < dag.blocks.size());
    REQUIRE(dag.blocks[blk].kind == "external");
    REQUIRE(n == (occ_extent + occ_target - 1) / occ_target);
  }
  auto const xs = x_values(dag);
  REQUIRE(xs.size() == 2);
  for (auto const& x : xs) {
    REQUIRE(run.footprint.count(x.id) == 1);
    REQUIRE(run.footprint.at(x.id) > 0);
    REQUIRE(run.builds.at(x.id) == run.n_batches.at(x.home_block));
    // X is a transient of its loop block (never assembled): its complete
    // value is the largest slice it was ever built as
    REQUIRE(run.size.at(x.id) == run.footprint.at(x.id));
  }
  // the value each block scatters into its whole is larger than the largest
  // per-batch slice that was written into it
  std::size_t scattered = 0;
  for (auto const& b : dag.blocks)
    for (auto const& [id, kind] : b.outputs)
      if (kind == "AccumulateScatter") {
        ++scattered;
        REQUIRE(run.size.at(id) > run.footprint.at(id));
      }
  REQUIRE(scattered == 2);
  // the caller's sections are left alone
  REQUIRE(run.extents.empty());
  REQUIRE(run.batch_targets.empty());
  REQUIRE(run.planner.empty());
  // every recorded id names a non-leaf value of the DAG record
  std::map<std::string, ValueRecord const*> by_id;
  for (auto const& v : dag.values) by_id[v.id] = &v;
  for (auto const& [id, n] : run.builds) {
    REQUIRE(by_id.count(id) == 1);
    REQUIRE(!by_id.at(id)->leaf);
    REQUIRE(n > 0);
    REQUIRE(run.footprint.count(id) == 1);
  }
}

/// Not run by default. Regenerates utilities/schedule-view/example.json, the
/// two-root fixture (R1, R2 sharing X = g*t, each batched on the other's
/// occupied index) as a complete document; the output path is taken from the
/// environment variable SEQUANT_SCHEDULE_EXAMPLE_OUT. The `run` section is the
/// one the ordered dry-run evaluation records (n_batches, footprint, builds);
/// extents, batch targets and planner settings are filled in from the
/// fixture's own policy and pinned extents.
TEST_CASE("write the viewer example document", "[.write-example]") {
  char const* const out_path = std::getenv("SEQUANT_SCHEDULE_EXAMPLE_OUT");
  REQUIRE(out_path != nullptr);

  auto const fx = two_root_fixture();
  Document doc;
  doc.producer_version = "example";
  doc.producer_date = "2026-10-19";

  sequant::OptimizeOptions opts;  // as dag_fixture() set it up
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = pin_extent;
  opts.batch_policy = external_occ_policy();
  doc.fingerprint = sequant::opt::schedule::fingerprint_of(
      fx.forest, {"R1", "R2"}, opts, {L"a", L"i"});
  doc.roots = fx.roots;
  doc.dag = record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);

  auto const id_of = [&](std::size_t vid) {
    return doc.dag->values.at(vid).id;
  };
  RunRecord run;
  dry_run(fx, &run, id_of);
  run.extents = {{"i", pin_extent(occ_index())},
                 {"a", pin_extent(virt_index())}};
  run.batch_targets = {
      {"i", external_occ_policy().batch_target_size(occ_index())}};
  run.planner = {{"objective_function", "DenseTimeSpaceBatched"}};
  doc.run = std::move(run);

  std::ofstream out(out_path, std::ios::binary);
  REQUIRE(out.good());
  out << to_json(doc);
}

TEST_CASE("ordered executor is unchanged without a run record",
          "[schedule-document]") {
  auto const fx = two_root_fixture();
  auto const dag =
      record_dag(fx.rich, fx.legality, fx.ordered, fx.forest, fx.roots);
  auto const id_of = [&](std::size_t vid) { return dag.values.at(vid).id; };

  auto const plain = dry_run(fx);
  auto const plain_counts = sequant::eval::detail::ordered_op_counts();
  RunRecord run;
  auto const recorded = dry_run(fx, &run, id_of);
  auto const recorded_counts = sequant::eval::detail::ordered_op_counts();

  REQUIRE(plain.size() == recorded.size());
  for (std::size_t i = 0; i != plain.size(); ++i) {
    REQUIRE(plain[i]);
    REQUIRE(recorded[i]);
    REQUIRE(plain[i]->size_in_bytes() == recorded[i]->size_in_bytes());
  }
  REQUIRE(plain_counts.builds == recorded_counts.builds);
  REQUIRE(plain_counts.assembles == recorded_counts.assembles);
  // the record's build tally is the executor's own build count
  std::size_t nbuilds = 0;
  for (auto const& [id, n] : run.builds) nbuilds += n;
  REQUIRE(nbuilds == recorded_counts.builds);
}

// ---------------------------------------------------------------------------
// size independence on a realistic input. The helpers below are copies of the
// witness fixture in test_ordered_executor.cpp (sizeind_ prefix; the sibling
// test files duplicate their fixtures rather than share a header, see
// tests/unit/CMakeLists.txt).
// ---------------------------------------------------------------------------

namespace {

// Problem size, provenance and occupancy adjustment: csv_size_fixtures.hpp.
using SizeIndProblemSize = sequant::tests::csv_sizes::CsvProblemSize;
inline constexpr auto const& kSizeIndWater20 =
    sequant::tests::csv_sizes::kWater20_pVDZF12;

/// every extent doubled
SizeIndProblemSize sizeind_doubled(SizeIndProblemSize p) {
  p.mu_tilde *= 2;
  p.aux *= 2;
  p.i_occ *= 2;
  return p;
}

std::string sizeind_slurp(std::string const& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

sequant::eval::dryrun::SizeRegime sizeind_df_regime(
    SizeIndProblemSize const& p) {
  return sequant::tests::csv_sizes::df_size_regime(p);
}

sequant::ExprPtr sizeind_flatten_product(sequant::ExprPtr const& e) {
  if (!e->is<sequant::Product>()) return e;
  auto const& p = e->as<sequant::Product>();
  return sequant::ex<sequant::Product>(p.scalar(), p.factors(),
                                       sequant::Product::Flatten::Yes);
}

/// the first \p nterms summands of the DF CCSD doubles residual, as one Sum
sequant::ExprPtr sizeind_residual_sum(std::size_t nterms) {
  auto const body = sizeind_slurp(std::string(SEQUANT_UNIT_TESTS_SOURCE_DIR) +
                                  "/data/csv_ccsd_doubles_residual_df.txt");
  REQUIRE(!body.empty());
  std::string line = body.substr(0, body.find('\n'));
  auto expr = sequant::deserialize<sequant::ExprPtr>(line);
  REQUIRE(static_cast<bool>(expr));
  REQUIRE(expr->is<sequant::Sum>());
  auto const& summands = expr->as<sequant::Sum>().summands();
  REQUIRE(summands.size() >= nterms);
  std::vector<sequant::ExprPtr> flat;
  for (std::size_t s = 0; s < nterms; ++s)
    flat.push_back(sizeind_flatten_product(summands[s]));
  return sequant::ex<sequant::Sum>(flat);
}

/// optimize (record, or replay \p pinned) -> binarize -> DAG stages ->
/// record_dag -> fingerprint_of, over the single root "R" = \p sum
Document build_document(sequant::ExprPtr const& sum,
                        SizeIndProblemSize const& problem,
                        std::size_t batch_target, RootRecord const* pinned) {
  using sequant::eval::dryrun::EvalExprDryRun;
  auto const regime = sizeind_df_regime(problem);
  auto const cm =
      std::make_shared<sequant::eval::dryrun::CostModel const>(regime);

  sequant::BatchPolicy policy;
  policy.is_batchable_contracted_index = [](sequant::Index const& ix) {
    return ix.space().base_key() == L"Κ";
  };
  policy.is_batchable_external_index = [](sequant::Index const&) {
    return false;
  };
  policy.batch_spectator_indices = false;
  policy.batch_target_size = [batch_target](sequant::Index const&) {
    return batch_target;
  };
  policy.is_volatile_leaf = [](sequant::Tensor const& t) {
    return t.label() == L"t";
  };
  policy.accumulation_factor = 1.0;
  policy.persistent_only = false;
  policy.peak_threshold = 1e11;

  sequant::OptimizeOptions opts;
  opts.objective_function = sequant::ObjectiveFunction::DenseTimeSpaceBatched;
  opts.idx_to_extent = regime.idx_to_extent();
  opts.inner_pow = regime.inner_pow_fn();
  opts.batch_policy = policy;
  opts.volatile_weight = 20.0;
  opts.roofline.machine_balance = 200.0;
  opts.roofline.fast_mem_elems = 1000000.0;
  opts.term_batch_axes = std::make_shared<AxesMap>();

  // a pinned run is recorded too, so its roots are what the pin produced
  std::optional<sequant::opt::schedule::PinSession> session;
  if (pinned) {
    session.emplace(*pinned);
    opts.pinned = &*session;
  }
  auto rr = sequant::opt::schedule::record_root("R", sum, opts);
  if (session) session->finish();
  sequant::ExprPtr const optimized = rr.optimized;
  RootRecord const record = std::move(rr.record);

  sequant::BinarizationOptions bopts;
  if (auto it = opts.term_batch_axes->find(optimized.get());
      it != opts.term_batch_axes->end())
    bopts.node_batch_axes = it->second;
  std::vector<sequant::eval::dryrun::EvalNodeDryRun> forest;
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  forest.push_back(sequant::binarize<EvalExprDryRun>(optimized, {}, bopts));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END

  auto rich = sequant::eval::compute_dag_boulevard(forest, *cm,
                                                   policy.batch_target_size);
  auto legality = sequant::eval::analyze_legality(rich, forest, policy);
  auto ordered =
      sequant::eval::build_ordered_schedule(rich, legality, policy, {});

  Document doc;
  doc.fingerprint =
      sequant::opt::schedule::fingerprint_of(forest, {"R"}, opts, {L"Κ"});
  doc.roots = {record};
  doc.dag = record_dag(rich, legality, ordered, forest, doc.roots);
  return doc;
}

}  // namespace

TEST_CASE("a pinned schedule is identical across two problem sizes",
          "[schedule-document]") {
  auto ctx = sequant::get_default_context_snapshot();
  ctx.set_first_dummy_index_ordinal(1000000);
  REQUIRE(ctx.index_basis_registry() != nullptr);
  auto isr = std::make_shared<sequant::IndexBasisRegistry>(
      *ctx.index_basis_registry());
  sequant::mbpt::add_pao_basis(isr, sequant::mbpt::Spin::any);
  sequant::mbpt::add_df_spaces(isr);
  ctx.set(isr);
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx));

  // number of leading residual summands: N = 8 already makes the planner's
  // choice differ between the regimes (see the final REQUIRE)
  constexpr std::size_t kNTerms = 8;
  auto const sum = sizeind_residual_sum(kNTerms);

  // regime A = water-20 pVDZ-F12; regime B = every extent doubled, batch
  // target doubled
  auto const A = kSizeIndWater20;
  auto const B = sizeind_doubled(kSizeIndWater20);
  auto docA = build_document(sum, A, 256, nullptr);  // records
  auto docB_free = build_document(sum, B, 512, nullptr);
  // the pin goes through the document's JSON form, as a user's file would
  auto const docA_loaded = from_json(to_json(docA));
  auto docB_pinned = build_document(sum, B, 512, &docA_loaded.roots[0]);

  REQUIRE(compare(docA, docB_pinned).empty());  // roots + dag identical
  // the fingerprint is keyed the way the DAG spells its modes: by
  // IndexBasis::base_key(), which for the named PAO basis is its name
  REQUIRE(docA.fingerprint.batchable.count("μ̃") == 1);
  for (auto const& v : docA.dag->values)
    for (auto const& m : v.result_modes)
      REQUIRE(docA.fingerprint.batchable.count(m.substr(m.rfind(':') + 1)) ==
              1);
  REQUIRE(!compare(docA, docB_free).empty());  // the planner changed its mind
  REQUIRE_NOTHROW(sequant::opt::schedule::check_fingerprint(
      docA.fingerprint, docB_pinned.fingerprint));
}
