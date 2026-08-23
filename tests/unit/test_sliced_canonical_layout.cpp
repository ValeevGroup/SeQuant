// Probe fixture for the sliced-value canonical layout / loop-coloring design
// (doc/dev/specs/2026-08-23-sliced-value-canonical-layout-loop-coloring-design.md,
// Summary + sec.2, mpqc4 repo). This is Task 1 of that plan: it does NOT
// implement the fix -- it pins the CURRENT (pre-fix) behavior that motivates
// it, using the substrate the design says already exists (\c occurrence_key /
// \c canonicalize_slots, occurrence_key.hpp; \c index_position,
// slicing_signature.hpp).
//
// THE BUG, MINIMIZED: a symmetric intermediate \c I(i_1,i_2) (Symmetry::Symm
// bra, so \c I(i_1,i_2) == I(i_2,i_1) as a genuine tensor-network automorphism)
// is shared -- CSE-correctly -- as the child of two different use-site
// contractions. In one use-site the surviving ("batched"/loop) mode is \c i_2
// (I's slot 1); in the other it is \c i_1 (I's slot 0) -- exactly the water-8
// pattern where the same cached cell is sliced at pos0 in some fetches and
// pos1 in others (see the design Summary, "measured on water-8 occ+aux: 29 vs
// 13").
//
// Today's router occurrence key (\c occurrence_key) colors named (batched)
// indices by SPACE, not by label, so it correctly folds these two
// occurrences -- that part is right and stays right. What it does NOT do
// today is record, as a per-VALUE fact, which physical mode the loop actually
// occupies: that "mode" is a per-*occurrence* quantity (\c index_position of
// the occurrence's own loop label), and it takes DIFFERENT values (0 vs 1) at
// the two occurrences even though they are the same value under the same
// router key. A per-cell map keyed only by the value's identity cannot store
// "pos0 here, pos1 there" -- this is the exact ambiguity sec.2 of the design
// resolves by making the loop identity a color on the canonical layout itself.
// Tasks 2 and 5 turn this probe into fold-AND-distinguish assertions once that
// coloring exists; here it is a baseline of today's reality, and it PASSES.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

namespace {

using sequant::AssertStrictBraKetSymmetry;
using sequant::bra;
using sequant::Context;
using sequant::ex;
using sequant::ExprPtr;
using sequant::ExprPtrList;
using sequant::get_default_context;
using sequant::Index;
using sequant::index_position;
using sequant::IndexSet;
using sequant::ket;
using sequant::Product;
using sequant::set_scoped_default_context;
using sequant::Symmetry;
using sequant::Tensor;
using sequant::container::svector;
using sequant::eval::occurrence_key;
using sequant::eval::RouterKeyEqual;
using sequant::eval::RouterKeyHash;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;

/// Binarizes \p expr (DryRun backend), declaring \p external as additional
/// uncontracted (batched/surviving) indices -- mirrors leaf_node/leaf helpers
/// in test_occurrence_key.cpp / test_slicing_signature.cpp.
EvalNodeDryRun probe_binarize(ExprPtr const& expr, IndexSet const& external) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(expr, external);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  return node;
}

}  // namespace

TEST_CASE(
    "sliced-canonical-layout baseline: a symmetric shared intermediate folds "
    "under occurrence_key while its loop's physical slot differs per "
    "use-site",
    "[sliced-layout]") {
  // occurrence_key's canonicalize_slots anonymizes any occ index NOT in this
  // use-site's own named/batched set (here: I's other occ index, i.e. the one
  // NOT nominated as this occurrence's loop). For a rank-2 Symm-bra tensor
  // with only one of its two same-space indices named,
  // TensorNetworkV3::create_graph's strict dummy-braket-arity check fires on
  // the anonymized one even though the network is well-formed -- a known,
  // separately-tracked pre-existing gap in occurrence_key (project
  // "occurrence_key_strict_braket" in the dev memory: "LEGITIMATELY
  // over-strict... reproduces at rank-2 CSV tensors"). Scoping the assertion
  // off is the same workaround used at test_canonicalize.cpp:270/292/331;
  // fixing occurrence_key's own strictness is out of scope for this task.
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};

  // I: the symmetric CSE intermediate. Symmetry::Symm on the (unpaired,
  // ket-less) bra makes the i1<->i2 exchange a genuine automorphism of I's own
  // tensor network -- I(i_1,i_2) and I(i_2,i_1) are literally the same value,
  // not merely two values that happen to hash the same.
  auto const I = [&] {
    return ex<Tensor>(L"I", bra(svector<Index>{i1, i2}), ket{}, Symmetry::Symm);
  };

  // Use-site A: I contracted against Y over i_1. i_1 is consumed (contracted
  // away); the mode that SURVIVES as this use-site's batched/loop index is
  // i_2 -- I's slot 1.
  auto const Y =
      ex<Tensor>(L"Y", bra(svector<Index>{i1}), ket{}, Symmetry::Nonsymm);
  auto const prodA = ex<Product>(ExprPtrList{I(), Y});
  auto const nodeA = probe_binarize(prodA, IndexSet{i2});

  // Use-site B: the SAME I, contracted against Z over i_2 instead. Now i_1 --
  // I's slot 0 -- is the one that survives as the loop index.
  auto const Z =
      ex<Tensor>(L"Z", bra(svector<Index>{i2}), ket{}, Symmetry::Nonsymm);
  auto const prodB = ex<Product>(ExprPtrList{I(), Z});
  auto const nodeB = probe_binarize(prodB, IndexSet{i1});

  // The shared operand at each use-site: I's own leaf node, as consumed by
  // that use-site's contraction (left factor, since I is listed first in both
  // products above).
  REQUIRE_FALSE(nodeA.leaf());
  REQUIRE_FALSE(nodeB.leaf());
  EvalNodeDryRun const& iA = nodeA.left();
  EvalNodeDryRun const& iB = nodeB.left();

  // --- (1) THE FOLD (correct, unchanged by this plan): the two occurrences
  // of I are the SAME value. At the leaf-node level this is literal identity
  // (both use-sites embed I(i_1,i_2) verbatim) --
  CHECK(iA->hash_value() == iB->hash_value());

  // -- and it also folds at the router-key level with each occurrence's OWN
  // loop nominated as the batched mode (occurrence_key colors named indices
  // by space, not label, so nominating i_2 in A vs i_1 in B does not break
  // the fold -- this is today's correct, existing behavior, the substrate the
  // loop-coloring design builds on):
  auto const kA = occurrence_key(iA, svector<Index>{i2});
  auto const kB = occurrence_key(iB, svector<Index>{i1});
  CHECK(RouterKeyEqual{}(kA, kB));
  CHECK(RouterKeyHash{}(kA) == RouterKeyHash{}(kB));

  // --- (2) THE AMBIGUITY (what the loop-coloring design fixes): there is no
  // per-VALUE fact today that says which physical slot the loop occupies.
  // index_position gives that answer only per-OCCURRENCE, using each
  // occurrence's own loop label -- and the two answers DIFFER even though (1)
  // just established this is the SAME value (same hash, same router key).
  auto const posA = index_position(iA, i2);  // A's loop, i_2, at I's slot 1
  auto const posB = index_position(iB, i1);  // B's loop, i_1, at I's slot 0
  REQUIRE(posA.has_value());
  REQUIRE(posB.has_value());
  CHECK(*posA == 1);
  CHECK(*posB == 0);
  // The crux: one shared value (hash_value and router key both agree it is
  // ONE thing), yet "the slot the loop occupies" is 1 at one occurrence and 0
  // at the other. A per-cell (per-value) `mode_to_level`-style map has
  // exactly one slot to record this fact and cannot hold both answers --
  // hence every attempt to resolve "which mode does the loop slice" via such
  // a map failed on precisely this case. The design's fix (sec.2) makes the
  // loop identity itself a color fed into canonicalize_slots, so the VALUE's
  // canonical storage layout pins one designated slot for the loop and each
  // occurrence carries its OWN permutation onto that layout instead of
  // requiring a single per-value slot number.
  CHECK(*posA != *posB);

  // Sanity: also confirm the two use-sites' own leaf nodes agree that BOTH
  // i_1 and i_2 are present (i.e. this is genuinely a 2-mode symmetric value,
  // not an artifact of one occurrence missing a mode).
  CHECK(index_position(iA, i1).has_value());
  CHECK(index_position(iB, i2).has_value());
}
