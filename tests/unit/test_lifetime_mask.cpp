#include <catch2/catch_test_macros.hpp>

#include <SeQuant/core/attr.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/lifetime_mask.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/io/shorthands.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <set>
#include <vector>

namespace {

using sequant::BatchModeType;
using sequant::EvalExpr;
using sequant::EvalNode;
using sequant::ExprPtr;
using sequant::home_scope;
using sequant::Index;
using sequant::stamp_lifetime_masks;
using sequant::stamp_seed_residency;

// A canonical eval-tree head from a two-factor product string. Two independent
// binarizations of the SAME string are structurally-equal canonical nodes (same
// hash_value + TreeNodeEqualityComparator), i.e. two occurrences of one
// canonical node -- exactly what the meet groups over.
EvalNode<EvalExpr> head(std::string_view product) {
  auto expr = sequant::deserialize<ExprPtr>(product);
  REQUIRE(static_cast<bool>(expr));
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize(expr);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  REQUIRE(!node.leaf());  // meet only stamps internal nodes
  return node;
}

std::set<Index> as_set(sequant::container::svector<Index> const& v) {
  return std::set<Index>(v.begin(), v.end());
}

std::set<Index> index_set(std::initializer_list<Index> l) {
  return std::set<Index>(l);
}

// Stamp a single External occ mode at a node.
void stamp_ext(EvalNode<EvalExpr>& n, Index ix) {
  n->set_batched_here({{std::move(ix), BatchModeType::External}});
}

// Stamp an External occ pair (i,j) at a node.
void stamp_ext_pair(EvalNode<EvalExpr>& n, Index i, Index j) {
  n->set_batched_here({{std::move(i), BatchModeType::External},
                       {std::move(j), BatchModeType::External}});
}

// Stamp a single Contracted mode at a node (mirrors stamp_ext, but tagged
// BatchModeType::Contracted instead of External).
void stamp_con(EvalNode<EvalExpr>& n, Index ix) {
  n->set_batched_here({{std::move(ix), BatchModeType::Contracted}});
}

// Stamp a Contracted occ pair (i,j) at a node (mirrors stamp_ext_pair, but
// tagged BatchModeType::Contracted instead of External).
void stamp_con_pair(EvalNode<EvalExpr>& n, Index i, Index j) {
  n->set_batched_here({{std::move(i), BatchModeType::Contracted},
                       {std::move(j), BatchModeType::Contracted}});
}

// Build an EvalExpr from a single-tensor spec (e.g. "R{i_1;a_5}"); its
// canon_indices are exactly the tensor's bra+ket slots.
EvalExpr eval_tensor(std::string_view tensor) {
  auto expr = sequant::deserialize<ExprPtr>(std::string(tensor));
  REQUIRE(static_cast<bool>(expr));
  return EvalExpr{expr->as<sequant::Tensor>()};
}

// A leaf eval node carrying the given tensor's slots.
EvalNode<EvalExpr> leaf(std::string_view tensor) {
  return EvalNode<EvalExpr>{eval_tensor(tensor)};
}

// An internal eval node whose OWN result slots are the given tensor's slots,
// with the two supplied child subtrees. The op field is immaterial to the mask
// pass (which reads only leaf-ness, batched_here, and canon_indices), so a
// tensor-derived EvalExpr with the intended result slots is a faithful stand-in
// for a real contraction result.
EvalNode<EvalExpr> inode(std::string_view result, EvalNode<EvalExpr> l,
                         EvalNode<EvalExpr> r) {
  return EvalNode<EvalExpr>{eval_tensor(result), std::move(l), std::move(r)};
}

// DFS-collect every internal (non-leaf) node reachable from a forest (the
// roots plus all their descendants), so a pass's per-node output can be
// snapshotted/compared without deep-copying the trees. Order is irrelevant
// to callers below; they only need a stable correspondence between two
// collections of the SAME (unmutated-in-shape) forest.
void collect_internal(EvalNode<EvalExpr> const& n,
                      std::vector<EvalNode<EvalExpr> const*>& out) {
  if (n.leaf()) return;
  out.push_back(&n);
  collect_internal(n.left(), out);
  collect_internal(n.right(), out);
}

std::vector<EvalNode<EvalExpr> const*> collect_forest_internal(
    std::vector<EvalNode<EvalExpr>> const& forest) {
  std::vector<EvalNode<EvalExpr> const*> out;
  for (auto const& t : forest) collect_internal(t, out);
  return out;
}

// A MIXED forest built fresh on every call (all-new deserialized ExprPtr
// instances, so two calls yield independent, unstamped trees): some
// canonical nodes carry External stamps only (one demoted across
// occurrences), some carry a Contracted pair kept in every occurrence, one
// carries BOTH an External and a Contracted mode on its own slots in every
// occurrence, and a scatter tree exercises the loop-invariant-sibling
// per-slot filter under an External ancestor.
std::vector<EvalNode<EvalExpr>> make_mixed_forest() {
  Index const i{L"i_1"}, j{L"i_2"}, k{L"i_4"};

  // External-only, demoted across occurrences (A slices, B is full).
  auto ext_A = head("F1{i_1;a_1} * F2{a_1;i_3}");
  auto ext_B = head("F1{i_1;a_1} * F2{a_1;i_3}");
  stamp_ext(ext_A, i);

  // Contracted-only pair, kept in every occurrence.
  auto con_A = head("U1{i_1;a_1} * U2{a_1;i_2}");
  auto con_B = head("U1{i_1;a_1} * U2{a_1;i_2}");
  stamp_con_pair(con_A, i, j);
  stamp_con_pair(con_B, i, j);

  // Mixed node: carries BOTH an External mode and a Contracted mode on its
  // own result slots, in every occurrence.
  auto mixed_A = head("W1{i_1;a_1} * W2{a_1;i_2}");
  auto mixed_B = head("W1{i_1;a_1} * W2{a_1;i_2}");
  auto stamp_mixed = [&](EvalNode<EvalExpr>& n) {
    n->set_batched_here(
        {{i, BatchModeType::External}, {j, BatchModeType::Contracted}});
  };
  stamp_mixed(mixed_A);
  stamp_mixed(mixed_B);

  // Scatter tree: root sliced by k (External); left child carries k on its
  // own result, right child is invariant to it.
  auto scatter_left =
      inode("A{i_4;a_3}", leaf("A1{i_4;a_1}"), leaf("A2{a_1;a_3}"));
  auto scatter_right =
      inode("B{a_3;a_5}", leaf("B1{a_3;a_4}"), leaf("B2{a_4;a_5}"));
  auto scatter_root =
      inode("R{i_4;a_5}", std::move(scatter_left), std::move(scatter_right));
  stamp_ext(scatter_root, k);

  return {ext_A, ext_B, con_A, con_B, mixed_A, mixed_B, scatter_root};
}

}  // namespace

TEST_CASE("lifetime mask cross-occurrence meet", "[lifetime_mask]") {
  Index const i{L"i_1"}, j{L"i_2"};

  // node_full: sliced by i in occurrence A, left full in occurrence B => the
  // meet demotes it to all-full (sliced /\ full = full).
  auto full_A = head("F1{i_1;a_1} * F2{a_1;i_3}");
  auto full_B = head("F1{i_1;a_1} * F2{a_1;i_3}");
  stamp_ext(full_A, i);  // A slices i; B carries no stamp (full)

  // node_pair: sliced by the same occ pair (i,j) in EVERY occurrence => the
  // meet keeps exactly {i, j}.
  auto pair_A = head("P1{i_1;a_1} * P2{a_1;i_2}");
  auto pair_B = head("P1{i_1;a_1} * P2{a_1;i_2}");
  stamp_ext_pair(pair_A, i, j);
  stamp_ext_pair(pair_B, i, j);

  std::vector<EvalNode<EvalExpr>> forest{full_A, full_B, pair_A, pair_B};
  stamp_lifetime_masks(forest);

  // A node sliced in one occurrence but full in another meets to all-full.
  CHECK(forest[0]->mask_all_full());  // full_A occurrence
  CHECK(forest[1]->mask_all_full());  // full_B occurrence

  // A node sliced by (i,j) in every occurrence meets to exactly {i, j}.
  CHECK_FALSE(forest[2]->mask_all_full());
  CHECK(as_set(forest[2]->sliced_modes()) == index_set({i, j}));
  CHECK(as_set(forest[3]->sliced_modes()) == index_set({i, j}));
}

TEST_CASE("lifetime mask OFF path is all-full (no External stamps)",
          "[lifetime_mask]") {
  // No batched_here stamps anywhere => every occurrence set empty => every mask
  // empty => all-full (the byte-identical OFF path).
  auto a = head("A1{i_1;a_1} * A2{a_1;i_2}");
  auto b = head("A1{i_1;a_1} * A2{a_1;i_2}");
  std::vector<EvalNode<EvalExpr>> forest{a, b};
  stamp_lifetime_masks(forest);
  CHECK(forest[0]->mask_all_full());
  CHECK(forest[1]->mask_all_full());
  CHECK(forest[0]->sliced_modes().empty());
}

TEST_CASE("lifetime mask expands a batched composite index proto-aware",
          "[lifetime_mask][proto]") {
  Index const i{L"i_1"}, j{L"i_2"};
  Index const a_pno{L"a_5",
                    {i, j}};  // PNO composite tied to the occ pair (i,j)

  // Both occurrences slice the SAME composite index a<i_1,i_2>; the meet
  // contributes its proto pair {i_1, i_2}, not the composite label itself.
  auto n_A = head("L1{a_3<i_1,i_2>;i_3} * M1{i_3;a_4}");
  auto n_B = head("L1{a_3<i_1,i_2>;i_3} * M1{i_3;a_4}");
  n_A->set_batched_here({{a_pno, BatchModeType::External}});
  n_B->set_batched_here({{a_pno, BatchModeType::External}});

  std::vector<EvalNode<EvalExpr>> forest{n_A, n_B};
  stamp_lifetime_masks(forest);

  CHECK_FALSE(forest[0]->mask_all_full());
  CHECK(as_set(forest[0]->sliced_modes()) == index_set({i, j}));
  // the composite index itself is NOT a sliced mode.
  CHECK(as_set(forest[0]->sliced_modes()).count(a_pno) == 0);
}

// Synthetic regression encoding the two ground-truth survey nodes
// (.superpowers/sdd/lifetime-mask-ir-survey.md), modeled directly rather than
// reproducing the heavy [.][dryrun-occ-veto] C60 forest (whose exact canonical
// hashes were captured at an ancestor commit and are fragile to reproduce).
TEST_CASE("lifetime mask survey ground-truth semantics",
          "[lifetime_mask][survey]") {
  Index const i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};

  // (1) s*C-class node (survey hash 1989507463377952644): all-full. Across its
  // occurrences it is sliced by DIFFERENT concrete occ pairs and left full in
  // others; both facts force the meet to empty (block-agnostic top cache).
  auto sC_1 = head("S1{i_1;a_1} * C1{a_1;i_1}");  // sliced by i_1
  auto sC_2 = head("S1{i_1;a_1} * C1{a_1;i_1}");  // sliced by i_2 (different)
  auto sC_3 = head("S1{i_1;a_1} * C1{a_1;i_1}");  // not sliced at all (full)
  stamp_ext(sC_1, i1);
  stamp_ext(sC_2, i2);
  // sC_3 unstamped

  // (2) g*C-class DF intermediate (survey hash 15545560759149115397): sliced by
  // the occ pair (i_1,i_2) in EVERY occurrence, carrying a free occ and a PNO
  // proto slot that differ across the two occurrence signatures. Meet = exactly
  // {i_1, i_2}; the PNO slot adds no extra sliced modes (proto-aware).
  auto gC_1 = head("G1{i_3;a_3<i_1,i_2>} * H1{a_3<i_1,i_2>;i_3}");
  auto gC_2 = head("G1{i_3;a_3<i_1,i_2>} * H1{a_3<i_1,i_2>;i_3}");
  auto gC_3 = head("G1{i_3;a_3<i_1,i_2>} * H1{a_3<i_1,i_2>;i_3}");
  stamp_ext_pair(gC_1, i1, i2);
  stamp_ext_pair(gC_2, i1, i2);
  stamp_ext_pair(gC_3, i1, i2);

  std::vector<EvalNode<EvalExpr>> forest{sC_1, sC_2, sC_3, gC_1, gC_2, gC_3};
  stamp_lifetime_masks(forest);

  // s*C is all-full under the meet.
  CHECK(forest[0]->mask_all_full());
  CHECK(forest[1]->mask_all_full());
  CHECK(forest[2]->mask_all_full());

  // g*C is sliced by exactly the external occ pair {i_1, i_2}.
  CHECK_FALSE(forest[3]->mask_all_full());
  CHECK(as_set(forest[3]->sliced_modes()) == index_set({i1, i2}));
  CHECK(as_set(forest[3]->sliced_modes()).count(i3) == 0);
}

TEST_CASE(
    "lifetime mask filters to a node's own slots (loop-invariant sibling stays "
    "all-full)",
    "[lifetime_mask]") {
  Index const i1{L"i_1"};

  // A scatter-like tree  R = A * B  under an External i_1 loop stamped at the
  // root. Only A carries i_1 on its own result; B is invariant to the i_1 loop
  // (its result slots {a_3, a_5} do not include i_1). The per-slot filter must
  // therefore leave B all-full (eligible for loop-invariant reuse) while A
  // keeps i_1 in its mask -- even though both live under the same batched
  // ancestor.
  auto A = inode("A{i_1;a_3}", leaf("A1{i_1;a_1}"), leaf("A2{a_1;a_3}"));
  auto B = inode("B{a_3;a_5}", leaf("B1{a_3;a_4}"), leaf("B2{a_4;a_5}"));
  auto R = inode("R{i_1;a_5}", std::move(A), std::move(B));
  stamp_ext(R, i1);  // External i_1 loop stamped on the root (ancestor of A, B)

  std::vector<EvalNode<EvalExpr>> forest{R};
  stamp_lifetime_masks(forest);

  auto const& root = forest[0];
  auto const& a = root.left();   // A: carries i_1
  auto const& b = root.right();  // B: invariant to i_1

  // Root carries i_1 on its own slot -> in its mask.
  CHECK_FALSE(root->mask_all_full());
  CHECK(as_set(root->sliced_modes()) == index_set({i1}));

  // Sibling A carries i_1 on its own result slot -> in its mask.
  CHECK_FALSE(a->mask_all_full());
  CHECK(as_set(a->sliced_modes()) == index_set({i1}));

  // Invariant sibling B does NOT carry i_1 -> all-full despite the batched
  // ancestor. This is the bug the per-slot filter fixes: raw ancestor+self
  // accumulation would wrongly stamp i_1 here and veto B from loop-invariant
  // run-scope residence.
  CHECK(b->mask_all_full());
  CHECK(b->sliced_modes().empty());
}

TEST_CASE("lifetime mask gates the run-scope veto", "[lifetime_mask][veto]") {
  Index const i{L"i_1"};

  auto is_volatile = [](EvalNode<EvalExpr> const&) { return false; };

  // A node whose mask is non-empty (sliced by i in EVERY occurrence) must be
  // refused run-scope residence by veto part (b) -- the mask, not a per-node
  // scalar placement level, drives that refusal. Two occurrences so
  // the CSE repeat count (>= min_repeats == 2) would otherwise admit it. The
  // node genuinely CARRIES i_1 on its own result slot (i_1 is a free bra index,
  // not contracted away), so the per-slot filter keeps it in the mask.
  auto sliced_A = head("V1{i_1;a_1} * W1{a_1;a_2}");
  auto sliced_B = head("V1{i_1;a_1} * W1{a_1;a_2}");
  stamp_ext(sliced_A, i);
  stamp_ext(sliced_B, i);

  // An all-full repeated node (no External stamp anywhere) is NOT vetoed and
  // is admitted at run scope.
  auto full_A = head("X1{i_2;a_2} * Y1{a_2;i_2}");
  auto full_B = head("X1{i_2;a_2} * Y1{a_2;i_2}");

  // No manual stamp_lifetime_masks() call here: the veto-engaged
  // sequant::cache_manager() overload below stamps the forest itself before
  // its DAG walk / veto (cache_manager.hpp), so relying on that internal
  // stamp is exactly what production callers (which never stamp themselves)
  // depend on.
  std::vector<EvalNode<EvalExpr>> forest{sliced_A, sliced_B, full_A, full_B};

  auto cache = sequant::cache_manager(forest, is_volatile);
  CHECK_FALSE(cache.exists(sliced_A));
  CHECK(cache.exists(full_A));
  // The builder's internal stamp is externally observable on the forest
  // elements it actually walked (EvalNode's copy ctor deep-copies, so
  // sliced_A/full_A above are independent of forest[0]/forest[2]).
  CHECK_FALSE(forest[0]->mask_all_full());
  CHECK(forest[2]->mask_all_full());
}

TEST_CASE(
    "stamp_seed_residency coincides with stamp_lifetime_masks when only "
    "External stamps are present",
    "[lifetime_mask][seed]") {
  Index const i{L"i_1"}, j{L"i_2"}, k{L"i_4"};

  // node_full: sliced by i in occurrence A, left full in occurrence B => the
  // meet demotes it to all-full under EITHER selector (with only External
  // stamps present, the External-only and all-batched-modes selectors
  // coincide).
  auto full_A = head("F1{i_1;a_1} * F2{a_1;i_3}");
  auto full_B = head("F1{i_1;a_1} * F2{a_1;i_3}");
  stamp_ext(full_A, i);

  // node_pair: sliced by the same occ pair (i,j) in EVERY occurrence => the
  // meet keeps exactly {i, j} under either selector.
  auto pair_A = head("P1{i_1;a_1} * P2{a_1;i_2}");
  auto pair_B = head("P1{i_1;a_1} * P2{a_1;i_2}");
  stamp_ext_pair(pair_A, i, j);
  stamp_ext_pair(pair_B, i, j);

  // Scatter tree: root sliced by k; child A carries k on its own result,
  // sibling B is invariant to it (per-slot filter exercise, same as the
  // dedicated per-slot-filter test above, but here checked under BOTH
  // selectors at once).
  auto A = inode("A{i_4;a_3}", leaf("A1{i_4;a_1}"), leaf("A2{a_1;a_3}"));
  auto B = inode("B{a_3;a_5}", leaf("B1{a_3;a_4}"), leaf("B2{a_4;a_5}"));
  auto R = inode("R{i_4;a_5}", std::move(A), std::move(B));
  stamp_ext(R, k);

  std::vector<EvalNode<EvalExpr>> forest{full_A, full_B, pair_A, pair_B, R};

  stamp_lifetime_masks(forest);
  stamp_seed_residency(forest);

  // Node-for-node equivalence: with only External stamps present, the
  // all-batched-modes selector (stamp_seed_residency) and the External-only
  // selector (stamp_lifetime_masks) walk the identical meet, so
  // seed_residency() == sliced_modes() everywhere the walk visits, including
  // descendants not directly listed in `forest`.
  auto check_same = [](EvalNode<EvalExpr> const& n) {
    CHECK(as_set(n->seed_residency()) == as_set(n->sliced_modes()));
  };
  for (auto const& n : forest) {
    check_same(n);
    if (!n.leaf()) {
      check_same(n.left());
      check_same(n.right());
    }
  }

  // Sliced-everywhere node: seed_residency is non-empty and matches the
  // known meet.
  CHECK_FALSE(forest[2]->seed_residency().empty());
  CHECK(as_set(forest[2]->seed_residency()) == index_set({i, j}));
  CHECK(as_set(forest[3]->seed_residency()) == index_set({i, j}));
}

TEST_CASE(
    "stamp_seed_residency genuinely drops the External-only filter "
    "(Contracted-mode discriminator)",
    "[lifetime_mask][seed]") {
  Index const i{L"i_1"}, j{L"i_2"};

  // node_con: sliced by a Contracted (NOT External) mode i on its own result
  // slot in EVERY occurrence. If stamp_seed_residency had accidentally kept
  // the External-only filter, its seed_residency() would be empty here (no
  // External stamps exist anywhere in this forest) -- so this case fails iff
  // that filter is present, unlike the all-External forest above where the
  // two selectors coincide regardless.
  auto con_A = head("Q1{i_1;a_1} * Q2{a_1;i_3}");
  auto con_B = head("Q1{i_1;a_1} * Q2{a_1;i_3}");
  stamp_con(con_A, i);
  stamp_con(con_B, i);

  // node_con_pair: same discriminator, but a Contracted occ PAIR kept in
  // every occurrence, exercising the multi-mode meet under the all-modes
  // selector.
  auto conpair_A = head("U1{i_1;a_1} * U2{a_1;i_2}");
  auto conpair_B = head("U1{i_1;a_1} * U2{a_1;i_2}");
  stamp_con_pair(conpair_A, i, j);
  stamp_con_pair(conpair_B, i, j);

  std::vector<EvalNode<EvalExpr>> forest{con_A, con_B, conpair_A, conpair_B};

  stamp_lifetime_masks(forest);
  stamp_seed_residency(forest);

  // stamp_seed_residency (all-batched-modes selector) keeps the Contracted
  // mode(s): the meet over every occurrence is non-empty and matches the
  // stamped set.
  CHECK(as_set(forest[0]->seed_residency()) == index_set({i}));
  CHECK(as_set(forest[1]->seed_residency()) == index_set({i}));
  CHECK(as_set(forest[2]->seed_residency()) == index_set({i, j}));
  CHECK(as_set(forest[3]->seed_residency()) == index_set({i, j}));

  // stamp_lifetime_masks (the External-only selector) sees no External
  // stamps at all here, so it must drop the Contracted mode(s) entirely --
  // every node meets to all-full under it.
  CHECK(forest[0]->mask_all_full());
  CHECK(forest[1]->mask_all_full());
  CHECK(forest[2]->mask_all_full());
  CHECK(forest[3]->mask_all_full());
  CHECK(forest[0]->sliced_modes().empty());
  CHECK(forest[2]->sliced_modes().empty());
}

TEST_CASE(
    "lifetime mask veto is byte-identical on the OFF path (no External "
    "stamps anywhere)",
    "[lifetime_mask][veto]") {
  // No batched_here stamps anywhere => every mask is empty (all-full), the
  // SAME condition the retired per-node scalar placement level satisfied for
  // every un-annotated node. A repeated canonical node is therefore
  // registered at run scope exactly as it was before this change.
  auto a = head("A1{i_1;a_1} * A2{a_1;i_2}");
  auto b = head("A1{i_1;a_1} * A2{a_1;i_2}");
  // No manual stamp_lifetime_masks() call: the veto-engaged
  // sequant::cache_manager() overload stamps the forest itself before its
  // veto reads the mask, matching production callers (which never stamp).
  std::vector<EvalNode<EvalExpr>> forest{a, b};

  auto is_volatile = [](EvalNode<EvalExpr> const&) { return false; };
  auto cache = sequant::cache_manager(forest, is_volatile);
  CHECK(cache.exists(a));
  // The builder's internal stamp leaves the OFF-path mask all-full, so the
  // veto admits exactly what it did before this change (byte-identical).
  CHECK(forest[0]->mask_all_full());
}

TEST_CASE("home_scope is an identity accessor over seed_residency",
          "[lifetime_mask][seed]") {
  // home_scope is a thin accessor: after stamp_seed_residency, home_scope(n)
  // must return exactly n->seed_residency() for every stamped node -- this
  // pins that identity so a future change to the accessor (or to what it
  // forwards to) is caught here first.
  Index const i{L"i_1"}, j{L"i_2"};

  auto pair_A = head("P1{i_1;a_1} * P2{a_1;i_2}");
  auto pair_B = head("P1{i_1;a_1} * P2{a_1;i_2}");
  stamp_con_pair(pair_A, i, j);
  stamp_con_pair(pair_B, i, j);

  std::vector<EvalNode<EvalExpr>> forest{pair_A, pair_B};
  stamp_seed_residency(forest);

  CHECK_FALSE(home_scope(forest[0]).empty());
  for (auto const& n : forest) {
    CHECK(as_set(home_scope(n)) == as_set(n->seed_residency()));
  }
  CHECK(as_set(home_scope(forest[0])) == index_set({i, j}));
}

// ---------------------------------------------------------------------
// T3: byte-unchanged-runtime guardrail.
// ---------------------------------------------------------------------
//
// stamp_seed_residency writes ONLY EvalExpr::seed_residency_; it must never
// perturb sliced_modes_ (the field the RUNTIME path -- stamp_lifetime_masks /
// place_at_this_level -- reads). These two tests turn that "purely additive"
// claim (lifetime_mask.hpp's doc comment) into an executable guarantee on a
// forest that mixes External-only, Contracted-only, and mixed-kind nodes.

TEST_CASE(
    "stamp_seed_residency never perturbs stamp_lifetime_masks' sliced_modes "
    "(byte-unchanged runtime guardrail)",
    "[lifetime_mask][seed]") {
  auto forest = make_mixed_forest();

  stamp_lifetime_masks(forest);

  // Snapshot every visited node's sliced_modes (all roots plus every
  // internal descendant), by value, before the seed pass runs.
  auto nodes = collect_forest_internal(forest);
  REQUIRE(!nodes.empty());
  std::vector<sequant::container::svector<Index>> before;
  before.reserve(nodes.size());
  for (auto const* n : nodes) before.push_back((*n)->sliced_modes());

  stamp_seed_residency(forest);

  // sliced_modes() must be byte-identical (same size, same elements, same
  // order) to its pre-seed-pass snapshot: the seed pass touches a disjoint
  // field.
  for (std::size_t k = 0; k < nodes.size(); ++k)
    CHECK((*nodes[k])->sliced_modes() == before[k]);
}

TEST_CASE(
    "stamp_lifetime_masks and stamp_seed_residency are order-independent "
    "(neither pass clobbers the other's field)",
    "[lifetime_mask][seed]") {
  // Control: stamp_lifetime_masks alone, on its own fresh forest.
  auto control = make_mixed_forest();
  stamp_lifetime_masks(control);
  auto control_nodes = collect_forest_internal(control);
  std::vector<sequant::container::svector<Index>> control_masks;
  control_masks.reserve(control_nodes.size());
  for (auto const* n : control_nodes)
    control_masks.push_back((*n)->sliced_modes());

  // Reversed order: stamp_seed_residency FIRST, then stamp_lifetime_masks,
  // on an independently constructed (but structurally identical) forest.
  auto reversed = make_mixed_forest();
  stamp_seed_residency(reversed);
  stamp_lifetime_masks(reversed);
  auto reversed_nodes = collect_forest_internal(reversed);

  REQUIRE(reversed_nodes.size() == control_nodes.size());
  for (std::size_t k = 0; k < control_nodes.size(); ++k)
    CHECK((*reversed_nodes[k])->sliced_modes() == control_masks[k]);
}

// ---------------------------------------------------------------------
// T3: against-definition tests for stamp_seed_residency, richer than T1's
// single equivalence case -- each drives the all-batched-modes meet on a
// hand-built forest and checks seed_residency() against the hand-computed
// DEEPEST meet on the node's own result slots.
// ---------------------------------------------------------------------

TEST_CASE(
    "stamp_seed_residency keeps a mode sliced (any batch kind) in EVERY "
    "occurrence",
    "[lifetime_mask][seed]") {
  Index const i{L"i_1"};

  // Occurrence A slices i as External; occurrence B slices the SAME index i
  // as Contracted. The kind differs, but the mode is batched (some kind) on
  // this node's own slot in EVERY occurrence => the all-batched-modes meet
  // keeps it regardless of kind.
  auto n_A = head("K1{i_1;a_1} * K2{a_1;i_3}");
  auto n_B = head("K1{i_1;a_1} * K2{a_1;i_3}");
  stamp_ext(n_A, i);
  stamp_con(n_B, i);

  std::vector<EvalNode<EvalExpr>> forest{n_A, n_B};
  stamp_seed_residency(forest);

  CHECK(as_set(forest[0]->seed_residency()) == index_set({i}));
  CHECK(as_set(forest[1]->seed_residency()) == index_set({i}));
}

TEST_CASE(
    "stamp_seed_residency drops a mode carried full (unbatched) in ONE "
    "occurrence",
    "[lifetime_mask][seed]") {
  Index const i{L"i_1"};

  // Occurrence A is Contracted-sliced by i on its own slot; occurrence B
  // carries i on the SAME slot but is NOT batched there at all (full) =>
  // the cross-occurrence intersection empties: the mode does not survive in
  // every occurrence, so it is dropped from the meet entirely.
  auto n_A = head("D1{i_1;a_1} * D2{a_1;i_3}");
  auto n_B = head("D1{i_1;a_1} * D2{a_1;i_3}");
  stamp_con(n_A, i);
  // n_B left unstamped (full).

  std::vector<EvalNode<EvalExpr>> forest{n_A, n_B};
  stamp_seed_residency(forest);

  CHECK(forest[0]->seed_residency().empty());
  CHECK(forest[1]->seed_residency().empty());
}

TEST_CASE(
    "stamp_seed_residency leaves a loop-invariant sibling unstamped (empty)",
    "[lifetime_mask][seed]") {
  Index const i1{L"i_1"};

  // Same scatter shape as the dedicated per-slot-filter test above, but
  // driven through stamp_seed_residency, with a Contracted (not External)
  // ancestor stamp, to also exercise the all-modes selector on this shape:
  // the left child carries i1 on its own result, the right child does not
  // (invariant to the i1 loop) and must stay unstamped even though both
  // live under the same batched ancestor.
  auto A = inode("A{i_1;a_3}", leaf("A1{i_1;a_1}"), leaf("A2{a_1;a_3}"));
  auto B = inode("B{a_3;a_5}", leaf("B1{a_3;a_4}"), leaf("B2{a_4;a_5}"));
  auto R = inode("R{i_1;a_5}", std::move(A), std::move(B));
  stamp_con(R, i1);

  std::vector<EvalNode<EvalExpr>> forest{R};
  stamp_seed_residency(forest);

  auto const& root = forest[0];
  auto const& a = root.left();
  auto const& b = root.right();

  CHECK(as_set(root->seed_residency()) == index_set({i1}));
  CHECK(as_set(a->seed_residency()) == index_set({i1}));
  CHECK(b->seed_residency().empty());
}

TEST_CASE(
    "stamp_seed_residency contains BOTH an External and a Contracted mode "
    "on a mixed node",
    "[lifetime_mask][seed]") {
  Index const i{L"i_1"}, j{L"i_2"};

  // Every occurrence stamps i as External AND j as Contracted on the SAME
  // node's own result slots => the meet keeps both, proving the union of
  // kinds (not just one selector's contribution winning).
  auto n_A = head("W1{i_1;a_1} * W2{a_1;i_2}");
  auto n_B = head("W1{i_1;a_1} * W2{a_1;i_2}");
  auto stamp_mixed = [&](EvalNode<EvalExpr>& n) {
    n->set_batched_here(
        {{i, BatchModeType::External}, {j, BatchModeType::Contracted}});
  };
  stamp_mixed(n_A);
  stamp_mixed(n_B);

  std::vector<EvalNode<EvalExpr>> forest{n_A, n_B};
  stamp_seed_residency(forest);

  CHECK(as_set(forest[0]->seed_residency()) == index_set({i, j}));
  CHECK(as_set(forest[1]->seed_residency()) == index_set({i, j}));
}

TEST_CASE("stamp_seed_residency expands a batched composite index proto-aware",
          "[lifetime_mask][seed][proto]") {
  Index const i{L"i_1"}, j{L"i_2"};
  Index const a_pno{L"a_5", {i, j}};  // PNO composite tied to the occ pair

  // Mirrors the dedicated [lifetime_mask][proto] case above, but driven
  // through stamp_seed_residency and tagged Contracted (not External) to
  // also discriminate against the External-only selector. Both occurrences
  // slice the SAME composite index a<i_1,i_2>; the meet contributes its
  // proto pair {i_1, i_2}, not the composite label itself.
  auto n_A = head("L1{a_3<i_1,i_2>;i_3} * M1{i_3;a_4}");
  auto n_B = head("L1{a_3<i_1,i_2>;i_3} * M1{i_3;a_4}");
  n_A->set_batched_here({{a_pno, BatchModeType::Contracted}});
  n_B->set_batched_here({{a_pno, BatchModeType::Contracted}});

  std::vector<EvalNode<EvalExpr>> forest{n_A, n_B};
  stamp_seed_residency(forest);

  CHECK_FALSE(forest[0]->seed_residency().empty());
  CHECK(as_set(forest[0]->seed_residency()) == index_set({i, j}));
  // the composite index itself is NOT a sliced mode.
  CHECK(as_set(forest[0]->seed_residency()).count(a_pno) == 0);
}
