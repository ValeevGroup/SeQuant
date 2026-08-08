#ifndef SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP
#define SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP

// Cross-term setup tables for the genetic optimizer: for every term and every
// subset S of its tensors, the face F(S) (the open indices, under the
// prototype's eff/carrier closure so composite proto-indices participate), a
// GLOBAL canonical-subnet key id, and per-subset index bitmasks for merge
// costing (cost.hpp). Named indices are colored by kind only, so two clusters
// -- same term or not -- share a key id iff they evaluate to the same array.

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/expressions/tensor.hpp>
#include <SeQuant/core/optimize/ga/genome.hpp>
#include <SeQuant/core/tensor_network.hpp>

#include <cstdint>
#include <functional>
#include <string>

namespace sequant::opt::ga {

/// Set of per-term index bits (bit k = TermTable::index_list[k]).
using IndexMask = std::uint64_t;

using FaceSet = TensorNetwork::NamedIndexSet;

inline constexpr std::size_t no_key = std::numeric_limits<std::size_t>::max();

/// Everything the optimizer needs to know about one summand (term).
struct TermTable {
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
  /// nothing downstream may legally reorder them), so this is an ordered
  /// sequence, never a set: permuting it would relabel every axis of every
  /// shared intermediate while leaving every flop count identical.
  container::vector<container::svector<std::uint8_t, 16>> canon_face_bits;
  /// Per subset: bits of the non-composite face indices (incl. proto-indices
  /// of composite face members -- mirrors tot_indices.outer).
  container::vector<IndexMask> outer_mask;
  /// Per subset: bits of the composite (CSV/PNO) face indices.
  container::vector<IndexMask> inner_mask;
  /// Bit v set iff tensors[v] is a volatile (amplitude-dependent) leaf, per the
  /// CostParams::is_volatile_leaf predicate. A cluster S is volatile iff
  /// `S & volatile_mask` is non-zero -- volatility is inherited by every
  /// superset, so no propagation pass is needed. Zero (hence replay weighting
  /// inert) when no predicate is supplied.
  NodeMask volatile_mask = 0;
  /// Distinct indices of this term (slots and their proto-indices); defines
  /// the bit positions of the masks above.
  container::svector<Index> index_list;
  /// Per index bit: extent.
  container::svector<double> extent;
  /// Per index bit: kind id, dense over the whole KeyTable; two indices may
  /// correspond in an L2 face bijection only if their kinds match.
  container::svector<int> kind;

  // --- derived bit tables ---------------------------------------------------
  // Every question the L2 pass asks about a face is a bit question about THIS
  // term; these three tables are what let it be asked that way.

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
  /// `index_rank` reproduces `std::sort` on `Index` exactly; because the
  /// release build compiles asserts out, that is pinned unconditionally by the
  /// "ga index bit tables" [ga] test, which recomputes the sort independently.
  container::svector<std::uint8_t> index_rank;
  /// Per index bit: bits of that index's proto-indices (0 if not composite).
  /// Every proto-index of a slot is itself in `index_list`, so (Sigma1)'s
  /// proto relation is a mask compare.
  container::svector<IndexMask> proto_mask;
  /// Per kind id (padded to the KeyTable's final kind count): bits of this
  /// term's indices of that kind.
  container::svector<IndexMask> kind_mask;

  std::size_t n() const { return tensors.size(); }
  NodeMask full() const { return (NodeMask{1} << n()) - 1; }
  int bit_of(Index const& ix) const;
  /// Bits of F(S) -- the only stored form of the face. The two masks are
  /// disjoint by construction (composite vs not).
  IndexMask face_mask(NodeMask S) const {
    return outer_mask[S] | inner_mask[S];
  }
  double face_size(NodeMask S) const;

  /// F(S) as a real `Index` set, materialized on demand. Setup and emission
  /// only -- the evaluation path asks bit questions of `face_mask`.
  FaceSet face_set(NodeMask S) const;
  /// The canonical face of S as `Index` objects, **in slot order** (slot i is
  /// `index_list[canon_face_bits[S][i]]`). See `canon_face_bits`.
  container::svector<Index> canon_face_indices(NodeMask S) const;
};

/// One optimization target (one ResultExpr): its terms.
struct TargetBlock {
  container::svector<std::size_t> terms;  ///< indices into KeyTable::terms
};

struct KeyTable {
  container::vector<TermTable> terms;
  container::svector<TargetBlock> targets;
  /// Whether build_key_table was given a volatile-leaf predicate. False leaves
  /// every replay weight inert: L1 through the empty volatile_masks, L2 through
  /// this flag (L2 work is replayed whenever volatility is modelled at all, so
  /// it cannot key off a mask).
  bool volatility_aware = false;
  std::size_t n_keys = 0;
  container::map<std::pair<IndexSpace, std::size_t>, int> kind_ids;

  /// Kind id of an index (registered during build).
  int kind_of(Index const& ix) const;
};

/// One target's input: its label and its summands, each a Product of Tensors
/// (with an optional scalar prefactor) plus that summand's external (result)
/// indices, which hand-assembled universes may name per term.
struct TargetInput {
  std::wstring label;
  container::svector<ExprPtr> summands;
  container::svector<FaceSet> ext;
};

/// \param is_volatile_leaf marks a leaf tensor as amplitude-dependent; fills
///        TermTable::volatile_mask. Empty (default) => volatility-blind cost.
KeyTable build_key_table(
    container::svector<TargetInput> const& targets,
    std::function<std::size_t(Index const&)> const& idx_to_extent,
    std::function<bool(Tensor const&)> const& is_volatile_leaf = {});

}  // namespace sequant::opt::ga

#endif  // SEQUANT_CORE_OPTIMIZE_GA_KEY_TABLE_HPP
