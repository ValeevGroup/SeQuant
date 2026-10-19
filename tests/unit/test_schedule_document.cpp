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
