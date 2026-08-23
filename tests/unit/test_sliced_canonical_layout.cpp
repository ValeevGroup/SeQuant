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
#include <SeQuant/core/eval/backends/dryrun/cost_model_object.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/backends/dryrun/result.hpp>
#include <SeQuant/core/eval/backends/dryrun/size_regime.hpp>
#include <SeQuant/core/eval/cache_manager.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/eval/ordered_schedule.hpp>
#include <SeQuant/core/eval/result.hpp>
#include <SeQuant/core/eval/slicing_signature.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/tensor_network.hpp>
#include <SeQuant/core/tensor_network/typedefs.hpp>
#include <SeQuant/core/utility/macros.hpp>
#include <SeQuant/external/bliss/graph.hh>

#include <catch2/catch_test_macros.hpp>

#include <any>
#include <array>
#include <memory>

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
using sequant::TensorCanonicalizer;
using sequant::TensorNetwork;
using sequant::container::svector;
using LoopColorMap = sequant::tensor_network::NamedIndexColorMap;
using sequant::eval::in_scope_batched_on_node;
using sequant::eval::loop_colored_id;
using sequant::eval::loop_colored_layout;
using sequant::eval::LoopId;
using sequant::eval::occurrence_key;
using sequant::eval::RouterKeyEqual;
using sequant::eval::RouterKeyHash;
using sequant::eval::SlicedModeAssignment;
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

namespace {

/// Canonicalize a single-tensor network's slots with the given named (batched)
/// indices and (optional) per-index loop colors. Mirrors occurrence_key's
/// canonicalize_slots call, but exposes the loop-color knob directly so the
/// fold/distinguish/reduce probes need no eval scaffolding.
TensorNetwork::SlotCanonicalizationMetadata canon(
    ExprPtr const& tensor, IndexSet const& named,
    LoopColorMap const* colors = nullptr) {
  TensorNetwork tn{ExprPtrList{tensor}};
  return tn.canonicalize_slots(TensorCanonicalizer::cardinal_tensor_labels(),
                               &named, {}, colors);
}

bool graphs_equal(TensorNetwork::SlotCanonicalizationMetadata const& a,
                  TensorNetwork::SlotCanonicalizationMetadata const& b) {
  return bliss::ConstGraphCmp::cmp(*a.graph, *b.graph) == 0;
}

}  // namespace

// (2a) FOLD: a symmetric intermediate I whose loop is bound to i_1 in one
// occurrence and to i_2 in the other still canonicalizes to ONE form once the
// loop is colored. The two physical layouts are related by I's own i_1<->i_2
// automorphism, so coloring "the loop" (whichever physical slot it occupies)
// the same in both must NOT break the fold.
TEST_CASE(
    "sliced-canonical-layout fold: symmetric I with same loop-binding folds "
    "under loop coloring",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};
  constexpr std::size_t loopA = 1;

  auto const I = [&](Index const& a, Index const& b) {
    return ex<Tensor>(L"I", bra(svector<Index>{a, b}), ket{}, Symmetry::Symm);
  };

  // Occurrence A: physical layout I(i_1,i_2); its loop is i_1.
  LoopColorMap colorsA;
  colorsA[i1] = loopA;
  auto const kA = canon(I(i1, i2), IndexSet{i1}, &colorsA);

  // Occurrence B: physical layout I(i_2,i_1); its loop is i_2.
  LoopColorMap colorsB;
  colorsB[i2] = loopA;
  auto const kB = canon(I(i2, i1), IndexSet{i2}, &colorsB);

  // Same loop-binding-structure (loop=A at "a symm slot", free at the other)
  // -> one canonical form. The fold survives coloring.
  CHECK(graphs_equal(kA, kB));
}

// (2b) DISTINGUISH: a NON-symmetric I with both occ modes sliced, colored by
// two DIFFERENT loops, must NOT fold: binding loop A to slot 0 vs slot 1 is a
// genuinely different sliced result. To make the two slots structurally
// distinguishable (so the loop<->slot binding is observable) the two occ modes
// sit in the bra and the ket of a non-braket-symmetric I -- the natural
// non-symmetric rank-2 layout of the spec's `for i1: for i2: I(i1,i2)` case.
// Without loop coloring i_1 and i_2 are same-space and interchangeable, so the
// two orderings fold (space-only) -- it is precisely the loop color that must
// pull them apart.
TEST_CASE(
    "sliced-canonical-layout distinguish: non-symmetric I with swapped "
    "loop-bindings does not fold",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};
  constexpr std::size_t loopA = 1;
  constexpr std::size_t loopB = 2;

  // I(i_1; i_2): bra slot vs ket slot -- structurally distinct slots.
  auto const I = [&] {
    return ex<Tensor>(L"I", bra(svector<Index>{i1}), ket(svector<Index>{i2}),
                      Symmetry::Nonsymm);
  };

  // Layout 1: bra (i_1) sliced by loop A, ket (i_2) by loop B.
  LoopColorMap colors1;
  colors1[i1] = loopA;
  colors1[i2] = loopB;
  auto const k1 = canon(I(), IndexSet{i1, i2}, &colors1);

  // Layout 2: bra (i_1) sliced by loop B, ket (i_2) by loop A -- the
  // loop<->slot binding is swapped.
  LoopColorMap colors2;
  colors2[i1] = loopB;
  colors2[i2] = loopA;
  auto const k2 = canon(I(), IndexSet{i1, i2}, &colors2);

  CHECK_FALSE(graphs_equal(k1, k2));
}

// (2c) REDUCE (the plan's #1 non-regression anchor): with NO loop colors, the
// loop-colored path must be byte-identical to today's canonicalization. An
// empty color map and a null color map must both reproduce the no-color result
// exactly (same graph, same hash).
TEST_CASE(
    "sliced-canonical-layout reduce: empty loop-color map is byte-identical to "
    "today's canonicalization",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};

  // A control tensor with both indices named but no loop colors supplied.
  auto const T = [&] {
    return ex<Tensor>(L"I", bra(svector<Index>{i1, i2}), ket{}, Symmetry::Symm);
  };
  IndexSet named{i1, i2};

  auto const k_null = canon(T(), named, nullptr);  // today's path
  LoopColorMap empty;
  auto const k_empty = canon(T(), named, &empty);  // reduction path

  CHECK(graphs_equal(k_null, k_empty));
  CHECK(k_null.hash_value() == k_empty.hash_value());
}

// ===========================================================================
// Task 4: store-canonical / serve-permuted, at the CACHE granularity.
//
// The above probes (Task 1-3) establish the IDENTITY/layout side: a
// loop-colored canonicalize_slots gives one value a single canonical form.
// This section proves the other load-bearing assumption the design leans on
// (spec "Load-bearing assumptions" #1): that `CacheManager` can hold that one
// canonical-layout value ONCE and that the runtime can serve each of several
// occurrences its OWN slot order from that single resident buffer, via
// `Result::permute`. This is the exact API Task 7's executor will use:
//   - store:  `CacheManager::ensure_home_slot(key)` (or the bounded-use-count
//             overload) + `CacheManager::store(key, data)` -- `data` is built
//             once, in the value's canonical annot order (`node->annot()`).
//   - serve permuted: `CacheManager::access(key)` (or `access_at`/
//             `access_at_hops`) fetches the SAME resident buffer for every
//             occurrence; each occurrence then calls
//             `Result::permute(std::array<std::any, 2>{canonical_annot,
//             occurrence_annot})` to obtain ITS OWN distinct, correctly
//             reordered buffer. This mirrors exactly what the top-level
//             `sequant::evaluate(node, layout, leaf, cache)` overload already
//             does today (eval.hpp, the `result.pre->permute({node->annot(),
//             layout})` call) -- Task 7 only needs to supply a
//             PER-OCCURRENCE `layout` instead of one root-result layout.
// ===========================================================================
namespace {

using sequant::CacheManager;
using sequant::ResultPtr;
using sequant::eval::dryrun::CostModel;
using sequant::eval::dryrun::ExtentOverrides;
using sequant::eval::dryrun::ResultDryRun;
using sequant::eval::dryrun::SizeRegime;

/// A minimal one-space regime: the probe's assertions are about SLOT ORDER
/// and slice-width tracking through `permute()`, not about distinguishing
/// modes by extent, so a single named space suffices.
SizeRegime cache_probe_regime() {
  SizeRegime r;
  r.space_extent = {{L"i", 10}};
  return r;
}

}  // namespace

TEST_CASE(
    "sliced-canonical-layout cache: a value stored ONCE is served to two "
    "occurrences in two DIFFERENT permuted slot orders from the SAME "
    "resident buffer",
    "[sliced-layout]") {
  Index i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};

  // The value's own node identity (a leaf -- the network-canonicalization
  // machinery probed above is orthogonal to this cache-granularity probe;
  // any stable TreeNode key works). Nonsymm so its own bra order is not
  // reshuffled by TensorCanonicalizer -- the raw i_1,i_2,i_3 order is
  // preserved as this node's `annot()`.
  auto const I = ex<Tensor>(L"I", bra(svector<Index>{i1, i2, i3}), ket{},
                            Symmetry::Nonsymm);
  auto const key = probe_binarize(I, IndexSet{});
  auto const canon_order = key->annot();  // this value's OWN canonical layout
  REQUIRE(canon_order.size() == 3);

  auto const cm = std::make_shared<CostModel const>(cache_probe_regime());

  // The canonical-layout value, ALREADY narrowed on its first canonical mode
  // (as it would be if it were a batch-sliced resident value, e.g. one block
  // of an enclosing loop) -- proves the narrowed width is a per-LABEL fact
  // that must travel with its index through permute(), not a bare position.
  ExtentOverrides narrowed;
  narrowed[0] = 3;  // canon_order[0]'s width narrowed from 10 to 3
  ResultPtr const canonical =
      sequant::eval_result<ResultDryRun>(canon_order, cm, narrowed);
  auto const full_bytes = cm->memsize(canon_order);  // unsliced: 10^3
  auto const narrowed_bytes = cm->memsize(canon_order, narrowed);  // 3*10*10
  REQUIRE(narrowed_bytes < full_bytes);
  REQUIRE(canonical->size_in_bytes() == narrowed_bytes);

  // --- STORE: the value is produced and homed ONCE, resident (unbounded,
  // never drained on access -- exactly `ensure_home_slot`'s contract, used by
  // the ordered executor to home a root-scope value read by every consumer).
  auto cache = CacheManager<EvalNodeDryRun>::empty();
  cache.ensure_home_slot(key);
  auto const stored = cache.store(key, canonical);
  REQUIRE(stored);
  CHECK(stored.get() == canonical.get());  // store() is an implicit access

  // --- SERVE (canonical fetch): two INDEPENDENT occurrences each fetch the
  // resident value. Because the slot is resident (unbounded life), BOTH
  // fetches return the exact SAME buffer -- the cache stores once.
  auto const fetchA = cache.access(key);
  auto const fetchB = cache.access(key);
  REQUIRE(fetchA);
  REQUIRE(fetchB);
  CHECK(fetchA.get() == canonical.get());
  CHECK(fetchB.get() == canonical.get());
  CHECK(fetchA.get() == fetchB.get());  // same resident buffer, both reads

  // --- SERVE PERMUTED: occurrence A wants canon_order with its first two
  // slots swapped (canon_order[0] -- the NARROWED mode -- moves to
  // position 1); occurrence B wants a DIFFERENT permutation, its last two
  // slots swapped (the narrowed mode stays at position 0 this time).
  svector<Index> const layoutA{canon_order[1], canon_order[0], canon_order[2]};
  svector<Index> const layoutB{canon_order[0], canon_order[2], canon_order[1]};
  REQUIRE_FALSE(layoutA == canon_order);
  REQUIRE_FALSE(layoutB == canon_order);
  REQUIRE_FALSE(layoutA == layoutB);

  ResultPtr const readA = fetchA->permute(
      std::array<std::any, 2>{std::any{canon_order}, std::any{layoutA}});
  ResultPtr const readB = fetchB->permute(
      std::array<std::any, 2>{std::any{canon_order}, std::any{layoutB}});
  REQUIRE(readA);
  REQUIRE(readB);

  // Each permuted read is a DISTINCT buffer -- from the canonical resident
  // value AND from each other's read (this is the property `entry::holds`'s
  // doc comment at cache_manager.hpp asserts: "a sliced/permuted/phase-
  // shifted read of this entry is a DISTINCT buffer").
  CHECK(readA.get() != canonical.get());
  CHECK(readB.get() != canonical.get());
  CHECK(readA.get() != readB.get());

  // Each read's own slot order is EXACTLY its own requested layout (label
  // identity, not just a permutation count) --
  CHECK(readA->as<ResultDryRun>().indices() == layoutA);
  CHECK(readB->as<ResultDryRun>().indices() == layoutB);

  // -- and the narrowed mode's width followed ITS OWN label to its new
  // position in each occurrence's layout: canon_order[0] (width 3) is at
  // position 1 in A's layout, position 0 in B's layout. A read that got this
  // wrong (e.g. leaving the override at position 0 regardless of the actual
  // permutation) would report the WRONG occurrence's byte count here.
  auto const& ovA = readA->as<ResultDryRun>().overrides();
  auto const& ovB = readB->as<ResultDryRun>().overrides();
  REQUIRE(ovA.count(1) == 1);
  CHECK(ovA.at(1) == 3);
  CHECK(ovA.count(0) == 0);
  REQUIRE(ovB.count(0) == 1);
  CHECK(ovB.at(0) == 3);
  CHECK(ovB.count(1) == 0);

  // Both reads still report the SAME total (narrowed) size -- permutation
  // reorders modes, it does not change the value's realized volume.
  CHECK(readA->size_in_bytes() == narrowed_bytes);
  CHECK(readB->size_in_bytes() == narrowed_bytes);

  // The resident entry itself is untouched by either occurrence's permuted
  // read: a THIRD fetch still returns the one canonical buffer, unsliced-
  // layout, unpermuted.
  auto const fetchC = cache.access(key);
  REQUIRE(fetchC);
  CHECK(fetchC.get() == canonical.get());
  CHECK(fetchC->as<ResultDryRun>().indices() == canon_order);
}

// ===========================================================================
// Task 5: the loop_colored_id PRODUCER on REAL eval nodes.
//
// The canon()-based probes above exercise the coloring KNOB
// (canonicalize_slots' NamedIndexColorMap) on bare tensors. These cases
// exercise the producer `loop_colored_id(node, ctx_modes, value_id,
// assignment)`, which (a) draws the color-map keys from EXACTLY the indices
// occurrence_key nominates as named (in_scope_batched_on_node -- constraint 1),
// (b) colors each by the LoopId Task 3's assignment records for that (value,
// Index), and (c) dispatches through occurrence_key -- reducing to today's key
// byte-for-byte when nothing is colored (constraint 2). The assignment is
// hand-built plain data (loop_of only reads by_value; a real
// build_ordered_schedule is unnecessary to test the producer).
// ===========================================================================
namespace {

/// Builds the two symmetric-CSE use-site leaf nodes from the baseline probe:
/// I(i_1,i_2) symmetric, contracted at use-site A over i_1 (loop survivor
/// i_2, I's slot 1) and at use-site B over i_2 (loop survivor i_1, I's slot
/// 0). Returns the two whole product nodes; each I-leaf is `.left()`.
struct SymmProbe {
  EvalNodeDryRun nodeA, nodeB;
};
SymmProbe make_symm_probe(Index const& i1, Index const& i2) {
  auto const I = [&] {
    return ex<Tensor>(L"I", bra(svector<Index>{i1, i2}), ket{}, Symmetry::Symm);
  };
  auto const Y =
      ex<Tensor>(L"Y", bra(svector<Index>{i1}), ket{}, Symmetry::Nonsymm);
  auto const Z =
      ex<Tensor>(L"Z", bra(svector<Index>{i2}), ket{}, Symmetry::Nonsymm);
  return SymmProbe{
      probe_binarize(ex<Product>(ExprPtrList{I(), Y}), IndexSet{i2}),
      probe_binarize(ex<Product>(ExprPtrList{I(), Z}), IndexSet{i1})};
}

/// A SlicedModeAssignment with a single (value_id -> {(mode, loop)}) fact.
/// `levels` is left empty: loop_colored_id consults only `loop_of` (by_value).
SlicedModeAssignment one_fact(std::size_t vid, Index const& mode, LoopId loop) {
  SlicedModeAssignment a;
  a.by_value[vid].push_back({mode, loop});
  return a;
}

}  // namespace

// (5-FOLD) The two symmetric-value occurrences bind the loop to DIFFERENT
// physical slots (i_2 at A, i_1 at B) yet, once each is colored by the SAME
// LoopId, loop_colored_id yields the SAME id -- one batched-cache entry, no
// duplication (design sec.2). This is the exact w8 pattern (pos1 in some
// fetches, pos0 in others) folded by loop-binding-structure.
TEST_CASE(
    "loop_colored_id fold: symmetric value, loop bound to different physical "
    "slots, same LoopId -> same id",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};
  auto const probe = make_symm_probe(i1, i2);
  EvalNodeDryRun const& iA = probe.nodeA.left();  // I as consumed at A
  EvalNodeDryRun const& iB = probe.nodeB.left();  // same I, consumed at B
  constexpr LoopId loop = 7;

  // Value A's own loop is i_2; value B's own loop is i_1 -- both -> LoopId 7.
  auto const asgA = one_fact(/*vid=*/0, i2, loop);
  auto const asgB = one_fact(/*vid=*/1, i1, loop);

  auto const kA = loop_colored_id(iA, svector<Index>{i2}, /*vid=*/0, asgA);
  auto const kB = loop_colored_id(iB, svector<Index>{i1}, /*vid=*/1, asgB);

  CHECK(graphs_equal(kA, kB));
  CHECK(RouterKeyHash{}(kA) == RouterKeyHash{}(kB));

  // Constraint 1 (color actually attaches): the color-map key i_2 is exactly
  // one of the indices occurrence_key nominates as named for iA, so coloring
  // it MOVES the id off the uncolored one (had the key not matched the named
  // Index object, the color would silently no-op and the two would be equal).
  auto const uncoloredA = occurrence_key(iA, svector<Index>{i2});
  CHECK(in_scope_batched_on_node(iA, svector<Index>{i2}).count(i2) == 1);
  CHECK_FALSE(graphs_equal(kA, uncoloredA));
}

// (5-DISTINGUISH) A NON-symmetric value with both occ modes sliced by two
// DIFFERENT loops must NOT fold when the loop<->slot binding is swapped:
// binding loop 1 to the bra vs to the ket is a genuinely different sliced
// result and must be a DISTINCT cache entry (design sec.2, the `for i1: for i2:
// I(i1,i2)` case). The bra/ket-distinct slots make the binding observable; the
// loop color is what pulls the two orderings apart (space-only they would
// fold).
TEST_CASE(
    "loop_colored_id distinguish: non-symmetric value, swapped loop-bindings "
    "-> different ids",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};
  constexpr LoopId loop1 = 1, loop2 = 2;

  auto const J = ex<Tensor>(L"J", bra(svector<Index>{i1}),
                            ket(svector<Index>{i2}), Symmetry::Nonsymm);
  auto const nodeJ = probe_binarize(J, IndexSet{i1, i2});  // leaf, both named

  // Both modes are this one value's own sliced modes; only the loop binding
  // differs between the two assignments.
  SlicedModeAssignment asg1, asg2;
  asg1.by_value[0] = {{i1, loop1}, {i2, loop2}};
  asg2.by_value[0] = {{i1, loop2}, {i2, loop1}};  // swapped

  auto const k1 = loop_colored_id(nodeJ, svector<Index>{i1, i2}, 0, asg1);
  auto const k2 = loop_colored_id(nodeJ, svector<Index>{i1, i2}, 0, asg2);

  CHECK_FALSE(graphs_equal(k1, k2));

  // The canonical layout is a real per-value artifact (both slots are named).
  CHECK(loop_colored_layout(k1).size() == 2);
}

// (5-REDUCE) The plan's #1 anchor: an UNSLICED value (no assignment entry) must
// produce EXACTLY today's occurrence_key -- byte-for-byte, graph and hash. An
// empty color map dispatches through occurrence_key's null path; there is no
// sliced/unsliced branch in the value-id rule.
TEST_CASE(
    "loop_colored_id reduce: an unsliced value's id equals today's "
    "occurrence_key byte-for-byte",
    "[sliced-layout]") {
  Context ctx = get_default_context();
  ctx.set(AssertStrictBraKetSymmetry::No);
  auto const ctx_resetter = set_scoped_default_context(ctx);

  Index i1{L"i_1"}, i2{L"i_2"};
  auto const probe = make_symm_probe(i1, i2);
  EvalNodeDryRun const& iA = probe.nodeA.left();

  SlicedModeAssignment const empty;  // no by_value entries: nothing sliced
  auto const reduced =
      loop_colored_id(iA, svector<Index>{i2}, /*vid=*/0, empty);
  auto const today = occurrence_key(iA, svector<Index>{i2});

  CHECK(graphs_equal(reduced, today));
  CHECK(reduced.hash_value() == today.hash_value());

  // The extracted layout is identical to today's too (the reduction is total:
  // graph, hash, AND canonical slot order all match the pre-coloring key).
  CHECK(loop_colored_layout(reduced) == loop_colored_layout(today));
}
