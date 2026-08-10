// Task 1 of the whole-scope batched DAG execution design (see
// doc/dev/specs/2026-08-10-whole-scope-batched-dag-execution-design.md): pins
// build_scope_schedule (SeQuant/core/eval/scope_schedule.hpp) -- the narrow
// scope-tree data structure + builder consumed by the later Task-2+ executor.
// No execution here, only the placement projection.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/peak_profile.hpp>
#include <SeQuant/core/eval/scope_schedule.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

#include <catch2/catch_test_macros.hpp>

#include <string_view>
#include <utility>
#include <vector>

namespace {

using sequant::BatchModeType;
using sequant::deserialize;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::Index;
using sequant::eval::build_scope_schedule;
using sequant::eval::compute_dag_boulevard;
using sequant::eval::RichSchedule;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::SizeRegime;

// Build an EvalExpr from a single-tensor spec (e.g. "R{i_1;a_5}"); its
// canon_indices are exactly the tensor's bra+ket slots. Mirrors the helper in
// test_peak_profile.cpp / test_lifetime_mask.cpp. Named distinctly
// (scope_*) from that file's identical eval_tensor/leaf/inode helpers: under
// a Unity build both files land in the same translation unit, and same-named
// functions in file-level anonymous namespaces collide there (an unnamed
// namespace is scoped per-TU, not per-file).
EvalExpr scope_eval_tensor(std::string_view tensor) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return EvalExpr{expr->as<sequant::Tensor>()};
}

// A leaf eval node carrying the given tensor's slots.
EvalNode<EvalExpr> scope_leaf(std::string_view tensor) {
  return EvalNode<EvalExpr>{scope_eval_tensor(tensor)};
}

// An internal eval node whose OWN result slots are the given tensor's slots,
// with the two supplied child subtrees.
EvalNode<EvalExpr> scope_inode(std::string_view result, EvalNode<EvalExpr> l,
                               EvalNode<EvalExpr> r) {
  return EvalNode<EvalExpr>{scope_eval_tensor(result), std::move(l),
                            std::move(r)};
}

}  // namespace

TEST_CASE("build_scope_schedule places values at their home scope",
          "[scope-schedule]") {
  // Register the DF-aux space (Κ) on a scoped clone of the default context;
  // "i"/"a" are already present in the default context (as in
  // test_peak_profile.cpp).
  auto ctx0 = sequant::get_default_context().clone();
  auto isr = ctx0.mutable_index_space_registry();
  REQUIRE(isr != nullptr);
  sequant::mbpt::add_df_spaces(isr);  // Κ (DF aux)
  auto ctx_resetter = sequant::set_scoped_default_context(std::move(ctx0));

  Index const K{L"Κ_1"};

  // A two-term forest (P1, P2) sharing one K-carrying intermediate V. P1/P2
  // each realize a Κ batch loop (Contracted) ABOVE V: V carries Κ_1 on its
  // own slot, so the cross-occurrence meet homes it AT the Κ node; P's own
  // result slots ({i_1,a_3}) do not carry Κ (it is summed away by the loop),
  // so P -- the K-INDEPENDENT value -- homes at the root. P1 and P2 are
  // structurally identical, so they (and their identical V/W children) fold
  // by hash into one ValueCell each, exactly as in test_peak_profile.cpp's
  // CSE forest. V must be an INTERNAL node (not a leaf): stamp_lifetime_masks
  // (invoked by compute_dag_boulevard) skips leaves entirely ("leaves are not
  // stamped -- they carry no meet"), so only an internal node's sliced_modes
  // (and hence home_modes) is ever populated.
  auto make_V = [] {
    return scope_inode("V{Κ_1;i_2}", scope_leaf("V1{Κ_1;x_1}"),
                       scope_leaf("V2{x_1;i_2}"));
  };
  auto P1 = scope_inode("P{i_1;a_3}", make_V(), scope_leaf("W{a_1;a_3}"));
  auto P2 = scope_inode("P{i_1;a_3}", make_V(), scope_leaf("W{a_1;a_3}"));
  P1->set_batched_here({{K, BatchModeType::Contracted}});
  P2->set_batched_here({{K, BatchModeType::Contracted}});

  SizeRegime r;
  r.space_extent = {{L"i", 5}, {L"a", 10}, {L"Κ", 8}};
  CostModel const cm{r};
  auto const block_of = [](Index const&) -> std::size_t { return 1; };

  std::vector<EvalNode<EvalExpr>> forest{P1, P2};
  RichSchedule const rich = compute_dag_boulevard(forest, cm, block_of);
  REQUIRE(rich.cells.size() == 5);  // V1, V2, V, W, P fold pairwise

  auto const sched = build_scope_schedule(rich, /*mode_order=*/{L"Κ"});
  REQUIRE(sched.root.children.size() == 1);
  CHECK(sched.root.children[0].mode.space().base_key() == L"Κ");
  CHECK(sched.root.children[0].kind == BatchModeType::Contracted);
  // the shared K-carrying intermediate homes in the K node; a K-independent
  // value homes at root
  CHECK(!sched.root.children[0].homed_values.empty());
  CHECK(!sched.root.homed_values.empty());
  CHECK(sched.num_values == rich.cells.size());
}
