#ifndef SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP

// Cross-term setup tables for the genetic optimizer: for every term and every
// subset S of its tensors, the face F(S) (the open indices, computed with the
// prototype's eff/carrier closure so composite proto-indices participate), a
// GLOBAL canonical-subnet key id (canonicalize_slots colors named indices by
// kind only, so the key identifies the array irrespective of index labels),
// and per-subset index bitmasks for nanosecond merge costing (cost.hpp).
// Two clusters -- same term or not -- share a key id iff they evaluate to
// the same array.

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/optimize/ga/genome.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <cstdint>
#include <functional>
#include <string>

namespace sequant::opt::ga {

/// Set of per-term index bits (bit k = TermTable::index_list[k]); 64 bits
/// covers the widest CCSD R2 terms (slots plus proto-indices) with headroom.
using IndexMask = std::uint64_t;

using FaceSet = TensorNetwork::NamedIndexSet;

inline constexpr std::size_t no_key = std::numeric_limits<std::size_t>::max();

/// Everything the optimizer needs to know about one summand (term).
struct TermTable {
  /// Target this term contributes to.
  std::size_t target = 0;
  /// Scalar prefactor of the input summand.
  Constant::scalar_type scalar = 1;
  /// The data-node tensors, in input factor order; bit k of a NodeMask is
  /// tensors[k].
  container::svector<ExprPtr> tensors;
  /// External (result) indices of the term.
  FaceSet ext;
  /// Per subset: global canonical-subnet key id (no_key for the empty set).
  container::vector<std::size_t> key;
  /// Per subset: the canonical ORDER of the face F(S), as index bits --
  /// `index_list[canon_face_bits[S][i]]` is slot i. Two key-equal subsets
  /// correspond axis-wise by zipping these.
  ///
  /// **Slot position IS the array axis** (`named_leaf` in emit.cpp gives the
  /// emitted intermediate all-auxiliary nonsymmetric slots precisely so that
  /// nothing downstream may legally reorder them), so this list is an ordered
  /// sequence, never a set: permuting it would relabel every axis of every
  /// shared intermediate while leaving every flop count identical.
  ///
  /// Interning to bits (T-A7) replaced a `container::svector<Index>` per
  /// subset. At C4H10/cc-pVDZ that is 1.5 M `uint8_t` instead of 1.5 M 168-byte
  /// `Index` -- and, because `svector`'s inline storage is per element type,
  /// 6 MB of vector headers instead of 220 MB. `canon_face_indices` rebuilds
  /// the `Index` form for the two consumers that need it (emission, and
  /// `Fitness::correspondences`' public form).
  container::vector<container::svector<std::uint8_t, 16>> canon_face_bits;
  /// Per subset: bits of the non-composite face indices (incl. proto-indices
  /// of composite face members -- mirrors tot_indices.outer).
  container::vector<IndexMask> outer_mask;
  /// Per subset: bits of the composite (CSV/PNO) face indices.
  container::vector<IndexMask> inner_mask;
  /// Bit v set iff tensors[v] is a volatile (amplitude-dependent) leaf, per the
  /// CostParams::is_volatile_leaf predicate handed to build_key_table. A
  /// cluster S is volatile iff `S & volatile_mask` is non-zero, which is why no
  /// propagation pass is needed: containing an amplitude leaf is inherited by
  /// every superset. Zero when no predicate is supplied, which makes replay
  /// weighting inert (see CostModel::volatile_weight).
  NodeMask volatile_mask = 0;
  /// Distinct indices of this term (slots and their proto-indices); defines
  /// the bit positions of the masks above.
  container::svector<Index> index_list;
  /// Per index bit: extent.
  container::svector<double> extent;
  /// Per index bit: kind id, dense over the whole KeyTable; two indices may
  /// correspond in an L2 face bijection only if their kinds match.
  container::svector<int> kind;

  // --- derived bit tables (T-A6) --------------------------------------------
  // The L2 pass (Fitness::find_beta and the ambient bookkeeping) used to do all
  // of its work on real `Index` objects -- grouping faces by kind, sorting
  // them, and keying a `std::map` on them -- although every question it asks is
  // a bit question about THIS term. These three tables are what let it be one.

  /// Per index bit: the position of `index_list[bit]` in the `std::sort`-by-
  /// `Index` order of `index_list`.
  ///
  /// **This is a load-bearing order, not a convenience.** `find_beta`
  /// enumerates kind-respecting bijections with `std::next_permutation` over
  /// each kind's face indices sorted by `Index`, and returns the FIRST one that
  /// satisfies (Sigma1-3); which bijection that is decides how the second
  /// operand of an extraction is renamed at emission. Bit order is the order
  /// indices were first seen while walking the term's tensors and has nothing
  /// to do with `Index::operator<`, so ordering by bit would keep every cost
  /// identical while silently changing the emitted schedule. Ordering by
  /// `index_rank` reproduces `std::sort` on `Index` exactly; that is asserted
  /// here at build time where SEQUANT_ASSERT is live and -- because the release
  /// build compiles asserts out -- pinned unconditionally by the "ga index bit
  /// tables" [ga] test case, which recomputes the sort independently.
  container::svector<std::uint8_t> index_rank;
  /// Per index bit: bits of that index's proto-indices (0 if not composite).
  /// Every proto-index of a slot is itself in `index_list` (build_term
  /// registers protos first), so (Sigma1)'s proto relation is a mask compare.
  container::svector<IndexMask> proto_mask;
  /// Per kind id (padded to the KeyTable's final kind count): bits of this
  /// term's indices of that kind. `face_mask(S) & kind_mask[k]` is the face's
  /// k-block.
  container::svector<IndexMask> kind_mask;

  std::size_t n() const { return tensors.size(); }
  NodeMask full() const { return (NodeMask{1} << n()) - 1; }
  int bit_of(Index const& ix) const;
  /// Bits of F(S). Since T-A7 this is the ONLY stored form of the face: the
  /// two masks are disjoint by construction (composite vs not) and together
  /// they name exactly the indices the dropped `container::vector<FaceSet>
  /// face` used to hold, which is what `face_set` reconstructs.
  IndexMask face_mask(NodeMask S) const {
    return outer_mask[S] | inner_mask[S];
  }
  double face_size(NodeMask S) const;

  /// F(S) as a real `Index` set, materialized on demand. Setup and emission
  /// only -- the evaluation path asks bit questions of `face_mask`. Ascending
  /// bit order is irrelevant to the result: `FaceSet` is an ordered set.
  FaceSet face_set(NodeMask S) const;
  /// The canonical face of S as `Index` objects, **in slot order** (slot i is
  /// `index_list[canon_face_bits[S][i]]`). See `canon_face_bits`.
  container::svector<Index> canon_face_indices(NodeMask S) const;
};

/// One optimization target (one ResultExpr): its name and its terms.
struct TargetBlock {
  std::wstring label;
  container::svector<std::size_t> terms;  ///< indices into KeyTable::terms
};

struct KeyTable {
  container::vector<TermTable> terms;
  container::svector<TargetBlock> targets;
  /// Whether build_key_table was given a volatile-leaf predicate. False leaves
  /// every replay weight inert: L1 through the empty volatile_masks, L2 through
  /// this flag (L2 work is replayed whenever volatility is modelled at all --
  /// see Fitness::cost_of_value -- so it cannot key off a mask).
  bool volatility_aware = false;
  std::size_t n_keys = 0;
  container::map<std::pair<IndexSpace, std::size_t>, int> kind_ids;
  /// The extent provider the table was built with, kept for callers that need
  /// to size something the table does not already carry. The optimizer itself
  /// no longer calls it: every extent it needs is in `TermTable::extent`, and
  /// since T-B1 seeding costs merges through those masks rather than through
  /// `Index` sets.
  std::function<std::size_t(Index const&)> idx_to_extent;

  /// Kind id of an index (registered during build).
  int kind_of(Index const& ix) const;
};

/// One target's input: its label and its summands, each a Product of Tensors
/// (with an optional scalar prefactor) plus that summand's external (result)
/// indices. Summands of a ResultExpr all share the result's indices, but
/// hand-assembled universes may name them per term.
struct TargetInput {
  std::wstring label;
  container::svector<ExprPtr> summands;
  container::svector<FaceSet> ext;
};

/// \param is_volatile_leaf marks a leaf tensor as amplitude-dependent (see
///        CostParams::is_volatile_leaf); fills TermTable::volatile_mask so the
///        cost model can replay-weight volatile work. Empty (default) => no
///        tensor is volatile and the cost stays volatility-blind.
KeyTable build_key_table(
    container::svector<TargetInput> const& targets,
    std::function<std::size_t(Index const&)> const& idx_to_extent,
    std::function<bool(Tensor const&)> const& is_volatile_leaf = {});

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP
