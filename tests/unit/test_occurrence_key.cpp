// Phase 2 T1: tests for the router's batching-aware occurrence key
// (SeQuant/core/eval/occurrence_key.hpp). These tests pin the ONE genuinely
// new algorithmic component of the placement-router seam: canonicalizing a
// node's own sub-expression as a colored TensorNetwork with the in-scope
// batched indices as named_indices, so that:
//   - non-batched index labels never matter (flat genericization),
//   - symmetric placements of batched indices collapse (proto-symmetry),
//   - distinct batched axes never accidentally collapse,
//   - the key is identical whether computed at a read site or a store site.
// No router type and no runtime wiring are introduced here; see the Phase 2
// plan for the tasks that build on this.

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval/backends/dryrun/eval_expr.hpp>
#include <SeQuant/core/eval/eval_expr.hpp>
#include <SeQuant/core/eval/fwd.hpp>
#include <SeQuant/core/eval/occurrence_key.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>

namespace {
namespace container = sequant::container;

using sequant::aux;
using sequant::BatchModeType;
using sequant::bra;
using sequant::ColumnSymmetry;
using sequant::ex;
using sequant::ExprPtr;
using sequant::Index;
using sequant::ket;
using sequant::Symmetry;
using sequant::Tensor;
using sequant::eval::in_scope_batched_on_node;
using sequant::eval::occurrence_key;
using sequant::eval::RouterKeyEqual;
using sequant::eval::RouterKeyHash;
using sequant::eval::dryrun::EvalExprDryRun;
using sequant::eval::dryrun::EvalNodeDryRun;

/// Builds a trivial (single-leaf) DryRun eval node wrapping \p t.
EvalNodeDryRun leaf_node(ExprPtr const& t) {
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_BEGIN
  auto node = sequant::binarize<EvalExprDryRun>(t);
  SEQUANT_PRAGMA_IGNORE_DEPRECATED_END
  return node;
}

}  // namespace

TEST_CASE("occurrence_key: flat genericization renames non-batched labels",
          "[occurrence_key]") {
  // A{i1,i2,i3} and A{i1,i2,i4}, batched {i1,i2}: i3/i4 sit in the same
  // (anonymous) slot position in both tensors, so they are recolored
  // identically by space alone -- the keys must collapse. Nonsymm bra AND
  // Nonsymm columns keep every slot position-fixed, so the collapse below is
  // due ONLY to anonymous recoloring, not to an incidental permutation
  // symmetry (see the distinct-axes test below, which needs this too).
  Index i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"}, i4{L"i_4"};

  auto t1 =
      ex<Tensor>(L"A", bra(container::svector<Index>{i1, i2, i3}), ket{},
                 Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);
  auto t2 =
      ex<Tensor>(L"A", bra(container::svector<Index>{i1, i2, i4}), ket{},
                 Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);

  auto n1 = leaf_node(t1);
  auto n2 = leaf_node(t2);

  container::svector<Index> ctx{i1, i2};
  auto k1 = occurrence_key(n1, ctx);
  auto k2 = occurrence_key(n2, ctx);

  CHECK(RouterKeyEqual{}(k1, k2));
  CHECK(RouterKeyHash{}(k1) == RouterKeyHash{}(k2));
}

TEST_CASE(
    "occurrence_key: gC proto-symmetry collapse under a shared batched occ "
    "index",
    "[occurrence_key]") {
  // gC{i1,i2,mu1,K1,a1<i1,i2>} vs gC{i2,i3,mu1,K1,a1<i2,i3>}, batched {i2}:
  // the occ pair is a genuine tensor symmetry (declared here via bra
  // Symmetry::Symm; the ctor's default ColumnSymmetry::Symm on the unpaired
  // bra slots would equally suffice, as verified separately -- either axis
  // makes i1<->i2 a real graph automorphism) and a1's protoindices are (by
  // default) symmetric too, so the i1<->i2 exchange consistently relocates
  // the shared batched index i2 to a canonical slot in both instances -- the
  // keys must collapse. mu1/K1 are spectator aux legs in distinct spaces
  // ("o"/"g"): identical in both instances and never named, so they must not
  // perturb the result. (Declaring BOTH axes Nonsymm makes the two keys
  // differ, confirming the collapse is not a construction artifact.)
  Index i1{L"i_1"}, i2{L"i_2"}, i3{L"i_3"};
  Index mu1{L"o_1"}, K1{L"g_1"};
  Index a1_12(L"a_1", {i1, i2});
  Index a1_23(L"a_1", {i2, i3});

  auto gc1 = ex<Tensor>(L"gC", bra(container::svector<Index>{i1, i2}), ket{},
                        aux(container::svector<Index>{mu1, K1, a1_12}),
                        Symmetry::Symm);
  auto gc2 = ex<Tensor>(L"gC", bra(container::svector<Index>{i2, i3}), ket{},
                        aux(container::svector<Index>{mu1, K1, a1_23}),
                        Symmetry::Symm);

  auto n1 = leaf_node(gc1);
  auto n2 = leaf_node(gc2);

  container::svector<Index> ctx{i2};
  auto k1 = occurrence_key(n1, ctx);
  auto k2 = occurrence_key(n2, ctx);

  CHECK(RouterKeyEqual{}(k1, k2));
  CHECK(RouterKeyHash{}(k1) == RouterKeyHash{}(k2));
}

TEST_CASE(
    "occurrence_key: same batched slot, different batched label -> same key",
    "[occurrence_key]") {
  // The DAG-globality question (design spec 2026-08-07 sec.0): two structurally
  // identical values whose batched second (occ) slot is bound to i_3 vs i_4 --
  // the g.C legs. The occurrence key must be label-AGNOSTIC on the batched slot
  // so one DAG-global overlay serves both occurrences; the specific i_3/i_4
  // binding is resolved later, in home-scope computation.
  Index i1{L"i_1"}, i3{L"i_3"}, i4{L"i_4"};
  auto make = [](Index const& a, Index const& b) {
    return ex<Tensor>(L"B", bra(container::svector<Index>{a, b}), ket{},
                      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);
  };
  auto nA = leaf_node(make(i1, i3));
  auto nB = leaf_node(make(i1, i4));
  auto kA =
      occurrence_key(nA, container::svector<Index>{i3});  // batched slot 1
  auto kB =
      occurrence_key(nB, container::svector<Index>{i4});  // batched slot 1
  CHECK(RouterKeyEqual{}(kA, kB));
  CHECK(RouterKeyHash{}(kA) == RouterKeyHash{}(kB));
}

TEST_CASE("occurrence_key: distinct batched axes stay distinct",
          "[occurrence_key]") {
  // Same tensor, batched {i1} vs batched {i2}: Nonsymm bra AND Nonsymm
  // columns make slot position meaningful (the tensor ctor's default
  // ColumnSymmetry::Symm would otherwise let the two unpaired bra slots
  // swap, accidentally collapsing this case), so the two keys must NOT
  // collapse.
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t = ex<Tensor>(L"B", bra(container::svector<Index>{i1, i2}), ket{},
                      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);

  auto n = leaf_node(t);

  container::svector<Index> ctx1{i1};
  container::svector<Index> ctx2{i2};
  auto k1 = occurrence_key(n, ctx1);
  auto k2 = occurrence_key(n, ctx2);

  CHECK_FALSE(RouterKeyEqual{}(k1, k2));
}

TEST_CASE("occurrence_key: read-site and store-site keys agree",
          "[occurrence_key]") {
  // The correctness hazard the whole seam rests on: a consumer's (read-site)
  // batch context can carry an extra outer-loop mode (i5) that is not among
  // the node's own slots at all -- in_scope_batched_on_node must filter it
  // out so the read-site key equals the (shallower) store-site key.
  Index i1{L"i_1"}, i2{L"i_2"}, i5{L"i_5"};
  auto t = ex<Tensor>(L"A", bra(container::svector<Index>{i1, i2}), ket{},
                      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);

  auto n = leaf_node(t);

  container::svector<Index> store_ctx{i1};
  container::svector<Index> read_ctx{i1, i5};  // i5 absent from n's slots

  auto k_store = occurrence_key(n, store_ctx);
  auto k_read = occurrence_key(n, read_ctx);

  CHECK(RouterKeyEqual{}(k_store, k_read));
  CHECK(RouterKeyHash{}(k_store) == RouterKeyHash{}(k_read));
}

TEST_CASE(
    "in_scope_batched_on_node is kind-agnostic: a Contracted-origin mode on "
    "the node's own slot is in scope",
    "[occurrence_key]") {
  // in_scope_batched_on_node never inspects BatchModeType -- it filters
  // whatever Index list the caller passes as ctx_modes to the subset that
  // lives on node's own slots. This pins that a Contracted-tagged mode (not
  // just an External one) passes through when it sits on the node's own
  // slot, locking the reconciled kind-agnostic semantics documented in the
  // (fixed) header comment: this filters ALL in-scope batched modes to a
  // node's own slots, matching the unified all_batched_modes_of selector /
  // slot_modes_of in lifetime_mask.hpp -- not the External-only ext_modes_of.
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t = ex<Tensor>(L"A", bra(container::svector<Index>{i1, i2}), ket{},
                      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);

  auto n = leaf_node(t);
  // Stamp a Contracted (not External) mode on i1, which sits on n's own slot.
  n->set_node_slice_mask(container::svector<std::pair<Index, BatchModeType>>{
      {i1, BatchModeType::Contracted}});

  // Build ctx_modes the way an all-batched-modes caller would (any
  // BatchModeType, mirroring stamp_lifetime_masks's selector).
  container::svector<Index> ctx;
  for (auto const& [ix, kind] : n->node_slice_mask()) ctx.push_back(ix);

  auto named = in_scope_batched_on_node(n, ctx);
  CHECK(named.find(i1) != named.end());
}

// ---------------------------------------------------------------------------
// Pillar 1 (value identity): the loop-COLORED occurrence_key (3-arg
// NamedIndexColorMap) is the value-id substrate. These pin that the depth
// coloring DISTINGUISHES which slot a loop slices for a non-symmetric tensor,
// FOLDS it for a symmetric one, and that an empty color map is byte-identical
// to the 2-arg (space-only) key.
// ---------------------------------------------------------------------------

TEST_CASE(
    "value-id: depth coloring distinguishes two sliced slots on a "
    "non-symmetric tensor",
    "[occurrence_key][value-id]") {
  // B(i1,i2), both occ slots sliced but under loops at DIFFERENT depths. The
  // two depth assignments (i1@d0,i2@d1) vs (i1@d1,i2@d0) are genuinely
  // different sliced values for a positionally-rigid (Nonsymm) B and MUST get
  // distinct colored keys -- the space-only 2-arg key, blind to depth, cannot
  // express the difference.
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t = ex<sequant::Tensor>(
      L"B", bra(sequant::container::svector<Index>{i1, i2}), ket{},
      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);
  auto n = leaf_node(t);
  sequant::container::svector<Index> ctx{i1, i2};

  sequant::tensor_network::NamedIndexColorMap c1;
  c1.emplace(i1, 0);
  c1.emplace(i2, 1);
  sequant::tensor_network::NamedIndexColorMap c2;
  c2.emplace(i1, 1);
  c2.emplace(i2, 0);

  auto k1 = occurrence_key(n, ctx, &c1);
  auto k2 = occurrence_key(n, ctx, &c2);
  CHECK_FALSE(RouterKeyEqual{}(k1, k2));  // depth distinguishes the slicing
}

TEST_CASE("value-id: symmetric tensor folds the two depth assignments",
          "[occurrence_key][value-id]") {
  // Same as above but B is bra-symmetric: swapping the two slots is a symmetry,
  // so (i1@d0,i2@d1) and (i1@d1,i2@d0) are the SAME sliced value -- the colored
  // canonicalization must FOLD them to one key (no duplication).
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t = ex<sequant::Tensor>(
      L"B", bra(sequant::container::svector<Index>{i1, i2}), ket{},
      Symmetry::Symm, std::nullopt, ColumnSymmetry::Symm);
  auto n = leaf_node(t);
  sequant::container::svector<Index> ctx{i1, i2};

  sequant::tensor_network::NamedIndexColorMap c1;
  c1.emplace(i1, 0);
  c1.emplace(i2, 1);
  sequant::tensor_network::NamedIndexColorMap c2;
  c2.emplace(i1, 1);
  c2.emplace(i2, 0);

  auto k1 = occurrence_key(n, ctx, &c1);
  auto k2 = occurrence_key(n, ctx, &c2);
  CHECK(RouterKeyEqual{}(k1, k2));  // symmetry folds
}

TEST_CASE("value-id: empty color map is byte-identical to the 2-arg key",
          "[occurrence_key][value-id]") {
  // The #1 non-regression anchor: an EMPTY (but non-null) color map must
  // canonicalize identically to the space-only 2-arg key -- so an unsliced
  // value's id is unchanged.
  Index i1{L"i_1"}, i2{L"i_2"};
  auto t = ex<sequant::Tensor>(
      L"B", bra(sequant::container::svector<Index>{i1, i2}), ket{},
      Symmetry::Nonsymm, std::nullopt, ColumnSymmetry::Nonsymm);
  auto n = leaf_node(t);
  sequant::container::svector<Index> ctx{i1};

  sequant::tensor_network::NamedIndexColorMap empty;
  auto k_colored_empty = occurrence_key(n, ctx, &empty);
  auto k_2arg = occurrence_key(n, ctx);
  CHECK(RouterKeyEqual{}(k_colored_empty, k_2arg));
  CHECK(RouterKeyHash{}(k_colored_empty) == RouterKeyHash{}(k_2arg));
}
